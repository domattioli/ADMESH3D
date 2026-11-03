"""3D DistMesh (Persson and Strang, 2004) with sliver removal.

The functions here work on raw arrays.  Input validation, warnings and the
``Mesh`` wrapper live in ``admesh3d.api``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import Delaunay

from .domain import Domain
from .quality import tet_quality, tet_signed_volumes

Q_MIN = 0.05          # slivers below this radius ratio are removed
FSCALE = 1.1          # bar force scale factor (Persson's 3D value)
DELTAT = 0.1          # pseudo time step at the start of the run
# The time step is cooled during the last part of the run: it starts to halve
# every ANNEAL_HALF_LIFE * max_iter iterations once ANNEAL_START * max_iter
# iterations are done.  Without this, retriangulation keeps kicking nodes and
# the 3D iteration rarely settles below the dptol movement threshold.
ANNEAL_START = 0.4
ANNEAL_HALF_LIFE = 0.1
_TET_EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
_TET_FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])


@dataclass
class MeshResult:
    """Raw output of :func:`distmesh3d`."""
    nodes: np.ndarray
    tets: np.ndarray
    iterations: int = 0
    converged: bool = False
    guard_dropped: bool = False
    slivers_kept: int = 0
    notes: list = field(default_factory=list)
    move_history: list = field(default_factory=list)


def _size(size_fn, p, h_min, h_max):
    """Target size at points ``p`` (clipped to [h_min, h_max])."""
    if size_fn is None:
        return np.full(len(p), h_max)
    return np.clip(np.asarray(size_fn(p), dtype=float), h_min, h_max)


def seed_points(domain: Domain, h0: float, h_min: float, h_max: float,
                size_fn, seed: int) -> np.ndarray:
    """Jittered cubic lattice inside the domain, plus the fixed points."""
    rng = np.random.default_rng(seed)
    lo, hi = (np.array(b, dtype=float) for b in domain.bbox)
    axes = [np.arange(lo[k], hi[k] + 0.5 * h0, h0) for k in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid + rng.uniform(-0.01 * h0, 0.01 * h0, size=grid.shape)
    geps = 0.1 * h0
    p = grid[domain.sdf(grid) < geps]
    if size_fn is not None and len(p):
        h = _size(size_fn, p, h_min, h_max)
        keep = rng.uniform(size=len(p)) < (h_min / h) ** 3
        p = p[keep]
    fixed = np.asarray(domain.fixed_points(), dtype=float).reshape(-1, 3)
    if len(fixed) and len(p):
        dist = np.min(np.linalg.norm(p[:, None, :] - fixed[None, :, :], axis=2), axis=1)
        p = p[dist >= 0.5 * h0]
    return np.vstack([p, fixed]) if len(fixed) else p


def _orient(p, tets):
    """Reorder tets so every signed volume is positive; return tets, volumes."""
    vol = tet_signed_volumes(p, tets)
    tets = tets.copy()
    neg = vol < 0
    tets[neg] = tets[neg][:, [0, 1, 3, 2]]
    return tets, np.abs(vol)


def triangulate(p, domain, h0, geps):
    """Delaunay tets of ``p`` with outside centroids and flat tets removed."""
    if len(p) < 4:
        return np.empty((0, 4), dtype=np.int64)
    try:
        tets = Delaunay(p).simplices.astype(np.int64)
    except Exception:  # qhull fails on fully degenerate input
        return np.empty((0, 4), dtype=np.int64)
    cen = p[tets].mean(axis=1)
    tets = tets[domain.sdf(cen) < -geps]
    tets, vol = _orient(p, tets)
    return tets[vol > 1e-12 * h0 ** 3]


def _unique_edges(tets, n):
    e = tets[:, _TET_EDGES].reshape(-1, 2)
    e = np.sort(e, axis=1)
    key = np.unique(e[:, 0] * n + e[:, 1])
    return np.stack([key // n, key % n], axis=1)


def boundary_faces(tets):
    """Faces that belong to exactly one tet (sorted vertex triples)."""
    f = np.sort(tets[:, _TET_FACES].reshape(-1, 3), axis=1)
    uniq, inv, cnt = np.unique(f, axis=0, return_inverse=True, return_counts=True)
    return uniq[cnt == 1]


def _closed_shell(tets):
    """True when every boundary edge is shared by exactly 2 boundary faces."""
    bf = boundary_faces(tets)
    if len(bf) == 0:
        return False
    e = np.sort(bf[:, [[0, 1], [1, 2], [0, 2]]].reshape(-1, 2), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return bool(np.all(cnt == 2))


def _project(p, domain, h0):
    """Move points outside the domain back to the surface (one Newton step)."""
    d = domain.sdf(p)
    out = d > 0
    if not np.any(out):
        return p
    deps = np.sqrt(np.finfo(float).eps) * h0
    q = p[out]
    dq = d[out]
    grad = np.empty_like(q)
    for k in range(3):
        step = np.zeros(3)
        step[k] = deps
        grad[:, k] = (domain.sdf(q + step) - dq) / deps
    g2 = np.maximum(np.sum(grad ** 2, axis=1), 1e-30)
    p = p.copy()
    p[out] = q - (dq / g2)[:, None] * grad
    return p


def _score(p, tets, h0):
    """Snapshot score (lower is better) or None when it has a flat tet."""
    if len(tets) == 0:
        return None
    q = tet_quality(p, tets)
    k = max(1, int(np.ceil(0.05 * len(q))))
    worst5 = np.sort(q)[:k].mean()
    return (float(np.mean(q < Q_MIN)), -float(worst5))


def _slivers(p, tets):
    return np.flatnonzero(tet_quality(p, tets) < Q_MIN)


def _remove_boundary_slivers(p, tets, h0, order_by_quality=True):
    """Remove slivers that touch the boundary while the shell stays closed."""
    changed = True
    while changed:
        changed = False
        q = tet_quality(p, tets)
        bad = np.flatnonzero(q < Q_MIN)
        if len(bad) == 0:
            break
        bad = bad[np.argsort(q[bad])]
        bf = boundary_faces(tets)
        bkeys = {tuple(f) for f in bf.tolist()}
        for t in bad:
            faces = np.sort(tets[t][_TET_FACES], axis=1)
            if not any(tuple(f) in bkeys for f in faces.tolist()):
                continue
            trial = np.delete(tets, t, axis=0)
            if len(trial) and _closed_shell(trial):
                tets = trial
                changed = True
                break
    return tets


def _cleanup(p, tets, domain, h0, geps, fixed_mask, rng):
    """Sliver removal per research R2.  Returns (p, tets, guard_dropped, kept)."""
    tets = _remove_boundary_slivers(p, tets, h0)
    # Interior slivers: perturb their nodes and keep the better triangulation.
    for _ in range(5):
        bad = _slivers(p, tets)
        if len(bad) == 0:
            break
        nodes = np.unique(tets[bad])
        nodes = nodes[~fixed_mask[nodes]]
        trial_p = p.copy()
        trial_p[nodes] += rng.uniform(-0.1 * h0, 0.1 * h0, size=(len(nodes), 3))
        trial_p = _project(trial_p, domain, h0)
        trial = triangulate(trial_p, domain, h0, geps)
        if len(trial) == 0:
            continue
        old_q = tet_quality(p, tets)
        new_q = tet_quality(trial_p, trial)
        old_key = (np.mean(old_q < Q_MIN), -old_q.min())
        new_key = (np.mean(new_q < Q_MIN), -new_q.min())
        if new_key < old_key:
            p, tets = trial_p, trial
            tets = _remove_boundary_slivers(p, tets, h0)
    # Delete the deepest node of each remaining sliver and retriangulate.
    for _ in range(5):
        bad = _slivers(p, tets)
        if len(bad) == 0:
            break
        d = domain.sdf(p)
        drop = set()
        for t in bad:
            cand = [v for v in tets[t] if not fixed_mask[v] and d[v] < -geps]
            if cand:
                drop.add(max(cand, key=lambda v: abs(d[v])))
        if not drop:
            break
        keep = np.ones(len(p), dtype=bool)
        keep[list(drop)] = False
        remap = np.cumsum(keep) - 1
        trial_p = p[keep]
        trial = triangulate(trial_p, domain, h0, geps)
        if len(trial) == 0:
            break
        p, tets, fixed_mask = trial_p, trial, fixed_mask[keep]
        tets = _remove_boundary_slivers(p, tets, h0)
        del remap
    kept = len(_slivers(p, tets))
    # Final guard: orient and drop non-positive volumes.
    tets, vol = _orient(p, tets)
    good = vol > 0
    guard_dropped = bool(np.any(~good))
    tets = tets[good]
    return p, tets, guard_dropped, kept


def _compact(p, tets):
    """Drop unreferenced nodes and reindex."""
    used = np.unique(tets)
    remap = -np.ones(len(p), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return p[used], remap[tets]


def distmesh3d(domain: Domain, h0: float, h_min: float, h_max: float, size_fn=None,
               max_iter: int = 1000, ttol: float = 0.1, dptol: float = 0.001,
               seed: int = 0) -> MeshResult:
    """Run 3D DistMesh and return cleaned, positively oriented tets."""
    geps = 0.1 * h0
    p = seed_points(domain, h0, h_min, h_max, size_fn, seed)
    n_fixed = len(np.asarray(domain.fixed_points()).reshape(-1, 3))
    fixed_mask = np.zeros(len(p), dtype=bool)
    if n_fixed:
        fixed_mask[-n_fixed:] = True
    rng = np.random.default_rng(seed + 1)

    best = None            # (score, p, tets)
    p_old = np.full_like(p, np.inf)
    tets = None
    edges = None
    converged = False
    last_move = np.inf
    it = 0
    history = []
    for it in range(1, max_iter + 1):
        if np.max(np.linalg.norm(p - p_old, axis=1)) > ttol * h0:
            tets = triangulate(p, domain, h0, geps)
            p_old = p.copy()
            if len(tets) == 0:
                break
            edges = _unique_edges(tets, len(p))
            score = _score(p, tets, h0)
            vol = tet_signed_volumes(p, tets)
            if score is not None and np.all(np.abs(vol) >= 1e-12 * h0 ** 3):
                if best is None or score < best[0]:
                    best = (score, p.copy(), tets.copy())
        vec = p[edges[:, 0]] - p[edges[:, 1]]
        length = np.maximum(np.linalg.norm(vec, axis=1), 1e-30)
        hbars = _size(size_fn, 0.5 * (p[edges[:, 0]] + p[edges[:, 1]]), h_min, h_max)
        l0 = hbars * FSCALE * (np.sum(length ** 3) / np.sum(hbars ** 3)) ** (1.0 / 3.0)
        force = np.maximum(l0 - length, 0.0)
        fvec = (force / length)[:, None] * vec
        total = np.zeros_like(p)
        np.add.at(total, edges[:, 0], fvec)
        np.add.at(total, edges[:, 1], -fvec)
        total[fixed_mask] = 0.0
        cooling = 0.5 ** (max(it - ANNEAL_START * max_iter, 0.0) / (ANNEAL_HALF_LIFE * max_iter))
        step = DELTAT * cooling * total
        p = p + step
        p = _project(p, domain, h0)
        interior = domain.sdf(p) < -geps
        last_move = (np.max(np.linalg.norm(step[interior], axis=1)) / h0
                     if np.any(interior) else 0.0)
        history.append(last_move)
        if last_move < dptol:
            converged = True
            break

    if converged:
        tets = triangulate(p, domain, h0, geps)
        pts = p
    elif best is not None:
        pts, tets = best[1], best[2]
    else:
        tets = triangulate(p, domain, h0, geps)
        pts = p
    res = MeshResult(nodes=pts, tets=tets, iterations=it, converged=converged,
                     move_history=history)
    if len(tets) == 0:
        return res
    pts, tets, dropped, kept = _cleanup(pts, tets, domain, h0, geps, 
                                        _fixed_mask(pts, domain), rng)
    pts, tets = _compact(pts, tets)
    res.nodes, res.tets = pts, tets
    res.guard_dropped = dropped
    res.slivers_kept = kept
    return res


def _fixed_mask(p, domain):
    """Mask of nodes that coincide with the domain's fixed points."""
    fixed = np.asarray(domain.fixed_points(), dtype=float).reshape(-1, 3)
    mask = np.zeros(len(p), dtype=bool)
    if len(fixed) == 0:
        return mask
    d = np.min(np.linalg.norm(p[:, None, :] - fixed[None, :, :], axis=2), axis=1)
    return d < 1e-12
