"""3D force-based tetrahedral mesh generator.

Implements the truss-relaxation idea of Persson and Strang (SIAM Review, 2004)
in three dimensions: nodes are pushed apart along Delaunay edges toward a
target edge length, the node set is re-triangulated whenever it has moved far
enough, and nodes that leave the domain are projected back onto the boundary
with the signed distance function.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import Delaunay

from .domain import Domain
from .quality import tet_quality, tet_signed_volumes

_SLIVER_Q = 0.05
_EDGE_PAIRS = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
_FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
_PUSH = 1.1          # targets are inflated so that edges stay in compression
_BCC_CELL = 1.1     # BCC cell side over target edge length
_JITTER = 0.0        # seed jitter, in units of h0
_SNAP = 0.25         # sliver nodes this close (in h0) to the surface are snapped
_DT = 0.1            # base pseudo time step
# Under re-triangulation the 3D iteration tends to keep drifting and to grow
# slivers.  After _DECAY_START steps the time step is halved every _HALF_LIFE
# steps, which freezes the node set once its shape has settled.
_DECAY_START = 50
_HALF_LIFE = 20.0


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


def _target_size(size_fn, h_min, h_max):
    if size_fn is None:
        return lambda q: np.full(len(q), float(h_max))
    return lambda q: np.clip(np.asarray(size_fn(q), dtype=float).reshape(-1), h_min, h_max)


def _lattice(lo, hi, step):
    """Body-centred cubic lattice: a cubic grid plus the cell centres."""
    axes = [np.arange(lo[k], hi[k] + 0.5 * step, step) for k in range(3)]
    corner = np.stack([a.ravel() for a in np.meshgrid(*axes, indexing="ij")], axis=1)
    return np.vstack([corner, corner + 0.5 * step])


def _initial_nodes(domain, h0, sizes, fixed, rng, geps):
    lo = np.asarray(domain.bbox[0], float)
    hi = np.asarray(domain.bbox[1], float)
    pts = _lattice(lo, hi, _BCC_CELL * h0)
    # a small jitter breaks lattice symmetry (cospherical Delaunay ties)
    pts = pts + rng.uniform(-_JITTER, _JITTER, pts.shape) * h0
    pts = pts[domain.sdf(pts) < geps]
    if len(pts):
        h = sizes(pts)
        keep_prob = (np.min(h) / h) ** 3
        pts = pts[rng.random(len(pts)) < keep_prob]
    if len(fixed) and len(pts):
        dist = np.min(np.linalg.norm(pts[:, None, :] - fixed[None], axis=2), axis=1)
        pts = pts[dist > 0.5 * h0]
    return np.vstack([fixed, pts]) if len(fixed) else pts


def _gradient(domain, p, d, step):
    g = np.empty_like(p)
    for k in range(3):
        q = p.copy()
        q[:, k] += step
        g[:, k] = (domain.sdf(q) - d) / step
    return g


def _project(domain, p, movable, step):
    d = domain.sdf(p)
    out = (d > 0) & movable
    if np.any(out):
        q = p[out]
        for _ in range(3):
            dq = domain.sdf(q)
            g = _gradient(domain, q, dq, step)
            gg = np.maximum(np.sum(g * g, axis=1), 1e-24)
            q = q - (dq / gg)[:, None] * g
        p[out] = q
    return p


def _triangulate(domain, p, geps):
    tri = Delaunay(p).simplices
    cen = p[tri].mean(axis=1)
    tri = tri[domain.sdf(cen) < -geps]
    vol = tet_signed_volumes(p, tri)
    neg = vol < 0
    tri[neg] = tri[neg][:, [1, 0, 2, 3]]
    return tri


def _unique_edges(tri):
    e = np.sort(tri[:, _EDGE_PAIRS].reshape(-1, 2), axis=1)
    return np.unique(e, axis=0)


def _boundary_mask(tri):
    f = np.sort(tri[:, _FACES].reshape(-1, 3), axis=1)
    _, inv, cnt = np.unique(f, axis=0, return_inverse=True, return_counts=True)
    return (cnt[inv.reshape(-1)] == 1).reshape(-1, 4).any(axis=1)


def _strip_boundary_slivers(p, tri):
    while len(tri):
        bad = (tet_quality(p, tri) < _SLIVER_Q) & _boundary_mask(tri)
        if not np.any(bad):
            break
        tri = tri[~bad]
    return tri


def _final_cleanup(domain, p, movable, h0, geps, gstep, rounds=5):
    """Triangulate, strip boundary slivers, and break up interior slivers.

    A sliver that survives the boundary strip is broken through its free node
    closest to the boundary (largest sdf): snapped onto the surface when it is
    near it, deleted otherwise; then the nodes are re-triangulated.
    """
    # free nodes in a thin shell under the surface would only make flat tets
    d = domain.sdf(p)
    shell = movable & (d > -_SNAP * h0) & (d < 0)
    if np.any(shell):
        g = _gradient(domain, p[shell], d[shell], gstep)
        gg = np.maximum(np.sum(g * g, axis=1), 1e-24)
        p = p.copy()
        p[shell] -= (d[shell] / gg)[:, None] * g
        p = _project(domain, p, movable, gstep)
    dropped = int(np.sum(shell))
    for _ in range(rounds + 1):
        tri = _strip_boundary_slivers(p, _triangulate(domain, p, geps))
        bad = tri[tet_quality(p, tri) < _SLIVER_Q]
        if len(bad) == 0 or _ == rounds:
            break
        d = domain.sdf(p)
        victims = set()
        for tet in bad:
            cand = [i for i in tet if movable[i] and d[i] < -geps]
            if cand:
                victims.add(max(cand, key=lambda i: d[i]))
        if not victims:
            break
        victims = np.array(sorted(victims))
        near = d[victims] > -_SNAP * h0
        snap = victims[near]
        if len(snap):
            # near-surface node: put it on the surface, the flat tet is stripped
            g = _gradient(domain, p[snap], d[snap], gstep)
            gg = np.maximum(np.sum(g * g, axis=1), 1e-24)
            p = p.copy()
            p[snap] -= (d[snap] / gg)[:, None] * g
            p = _project(domain, p, np.ones(len(p), bool), gstep)
        deep = victims[~near]
        if len(deep):
            keep = np.ones(len(p), bool)
            keep[deep] = False
            p, movable = p[keep], movable[keep]
        dropped += len(victims)
    return tri, p, movable, dropped


def distmesh3d(domain: Domain, h0: float, h_min: float, h_max: float, size_fn=None,
               max_iter: int = 1000, ttol: float = 0.1, dptol: float = 0.001,
               seed: int = 0, callback=None, callback_every: int = 10) -> MeshResult:
    rng = np.random.default_rng(seed)
    sizes = _target_size(size_fn, h_min, h_max)
    geps = 1e-3 * h0
    gstep = np.sqrt(np.finfo(float).eps) * h0
    fixed = np.asarray(domain.fixed_points(), dtype=float).reshape(-1, 3)
    nfix = len(fixed)
    notes: list = []
    history: list = []

    p = _initial_nodes(domain, h0, sizes, fixed, rng, geps)
    movable = np.ones(len(p), bool)
    movable[:nfix] = False
    if len(p) < 4:
        notes.append("too few seed nodes")
        return MeshResult(p, np.empty((0, 4), int), notes=notes)

    if callback is not None:
        callback(0, p.copy())

    anchor = np.full_like(p, np.inf)
    tri = None
    edges = None
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        if np.max(np.linalg.norm(p - anchor, axis=1)) > ttol * h0:
            anchor = p.copy()
            tri = _triangulate(domain, p, geps)
            edges = _unique_edges(tri)

        vec = p[edges[:, 0]] - p[edges[:, 1]]
        length = np.linalg.norm(vec, axis=1)
        hbar = sizes(0.5 * (p[edges[:, 0]] + p[edges[:, 1]]))
        scale = (np.sum(length ** 3) / np.sum(hbar ** 3)) ** (1.0 / 3.0)
        rest = _PUSH * scale * hbar
        f = np.maximum(rest - length, 0.0)
        fvec = (f / np.maximum(length, 1e-300))[:, None] * vec
        force = np.zeros_like(p)
        for k in range(3):
            force[:, k] += np.bincount(edges[:, 0], fvec[:, k], len(p))
            force[:, k] -= np.bincount(edges[:, 1], fvec[:, k], len(p))
        force[~movable] = 0.0

        dt = _DT if it < _DECAY_START else _DT * 0.5 ** ((it - _DECAY_START) / _HALF_LIFE)
        step = dt * force
        p = p + step
        p = _project(domain, p, movable, gstep)

        inner = domain.sdf(p) < -geps
        move = float(np.max(np.linalg.norm(step[inner], axis=1))) if np.any(inner) else 0.0
        history.append(move)
        if callback is not None and it % callback_every == 0:
            callback(it, p.copy())
        if move < dptol * h0:
            converged = True
            break

    tri, p, movable, dropped = _final_cleanup(domain, p, movable, h0, geps, gstep)
    if dropped:
        notes.append(f"moved or removed {dropped} node(s) trapped in slivers")
    vol = tet_signed_volumes(p, tri)
    ok = vol > 1e-12 * h0 ** 3
    guard = bool(np.any(~ok))
    if guard:
        notes.append(f"guard dropped {int(np.sum(~ok))} non-positive tets")
    tri = tri[ok]

    used = np.unique(tri)
    remap = np.full(len(p), -1, int)
    remap[used] = np.arange(len(used))
    nodes = p[used]
    tets = remap[tri].astype(int)
    slivers = int(np.sum(tet_quality(nodes, tets) < _SLIVER_Q)) if len(tets) else 0
    return MeshResult(nodes, tets, iterations=it, converged=converged,
                      guard_dropped=guard, slivers_kept=slivers, notes=notes,
                      move_history=history)
