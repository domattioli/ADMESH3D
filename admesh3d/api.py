"""Public API: Mesh, tetrahedralize() and extrude()."""
from __future__ import annotations

import numbers
import warnings
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .distmesh3d import distmesh3d
from .domain import Domain
from .errors import ConvergenceWarning
from .extrude import sigma_extrude
from .quality import (QualityReport, make_report, tet_quality, tet_signed_volumes,
                      wedge_quality, wedge_signed_volumes)
from . import vtu

_K = {"tetra": 4, "wedge": 6}


@dataclass
class Mesh:
    """Volume mesh of tetrahedra or wedges with strictly positive volumes."""
    nodes: np.ndarray
    elements: np.ndarray
    element_type: str
    _volumes: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if self.element_type not in _K:
            raise ValueError(f"element_type must be 'tetra' or 'wedge', got {self.element_type!r}")
        self.nodes = np.ascontiguousarray(self.nodes, dtype=np.float64)
        self.elements = np.ascontiguousarray(self.elements, dtype=np.int64)
        if self.nodes.ndim != 2 or self.nodes.shape[1] != 3:
            raise ValueError("nodes must have shape (N, 3)")
        k = _K[self.element_type]
        if self.elements.ndim != 2 or self.elements.shape[1] != k:
            raise ValueError(f"elements must have shape (E, {k}) for {self.element_type}")
        if self.elements.size and (self.elements.min() < 0
                                   or self.elements.max() >= len(self.nodes)):
            raise ValueError("element indices out of range")
        self._volumes = self._compute_volumes()
        if self._volumes.size and not np.all(self._volumes > 0):
            bad = int(np.flatnonzero(~(self._volumes > 0))[0])
            raise ValueError(f"element {bad} has non-positive volume {self._volumes[bad]:g}")

    def _compute_volumes(self):
        if self.element_type == "tetra":
            return tet_signed_volumes(self.nodes, self.elements)
        return wedge_signed_volumes(self.nodes, self.elements)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_elements(self) -> int:
        return len(self.elements)

    @property
    def volumes(self) -> np.ndarray:
        return self._volumes

    @property
    def quality(self) -> np.ndarray:
        """Per-element quality in [0, 1]."""
        if self.element_type == "tetra":
            return tet_quality(self.nodes, self.elements)
        return wedge_quality(self.nodes, self.elements)

    def quality_report(self) -> QualityReport:
        return make_report(self.nodes, self.elements, self.element_type)

    def to_vtu(self, path) -> None:
        vtu.write_vtu(path, self.nodes, self.elements, self.element_type, self.quality)

    @classmethod
    def from_vtu(cls, path) -> "Mesh":
        nodes, elements, etype, _ = vtu.read_vtu(path)
        return cls(nodes, elements, etype)


def read_vtu(path) -> Mesh:
    """Read an ASCII .vtu file written by :meth:`Mesh.to_vtu`."""
    return Mesh.from_vtu(path)


# ---------------------------------------------------------------- validation

def _finite_pos(value, name):
    if (isinstance(value, bool) or not isinstance(value, numbers.Real)
            or not np.isfinite(value) or value <= 0):
        raise ValueError(f"{name} must be finite and > 0, got {value!r}")
    return float(value)


def _sample_grid(domain: Domain):
    lo, hi = (np.array(b, dtype=float) for b in domain.bbox)
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        raise ValueError("bbox is not finite")
    if not np.all(hi > lo):
        raise ValueError("bbox is degenerate (needs hi > lo in every coordinate)")
    axes = [np.linspace(lo[k], hi[k], 33) for k in range(3)]
    pts = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    pts = np.vstack([pts, 0.5 * (lo + hi)])
    d = np.asarray(domain.sdf(pts), dtype=float)
    if d.shape != (len(pts),):
        raise ValueError("sdf returned an array of the wrong shape")
    if not np.all(np.isfinite(d)):
        raise ValueError("sdf returned non-finite values")
    if not np.any(d < 0):
        raise ValueError("domain has no interior (sdf is never negative on the sample grid)")
    return lo, hi, pts, d


def _checked_size_fn(size_fn):
    def wrapped(p):
        h = np.asarray(size_fn(p), dtype=float)
        if h.shape != (len(p),) or not np.all(np.isfinite(h)) or np.any(h <= 0):
            raise ValueError("size_fn returned invalid values "
                             "(need a finite positive array of shape (n,))")
        return h
    return wrapped


def estimate_nodes(domain, pts, d, h_min, h_max, size_fn) -> float:
    """Estimated node count (research R9)."""
    lo, hi = (np.array(b, dtype=float) for b in domain.bbox)
    inside = d < 0
    vol = domain.volume
    if vol is None:
        vol = float(np.prod(hi - lo)) * float(np.mean(inside))
    if size_fn is None:
        tets = vol / (h_max ** 3 / (6 * np.sqrt(2)))
    else:
        h = np.clip(size_fn(pts[inside]), h_min, h_max)
        tets = float(np.mean(1.0 / h ** 3)) * vol * 6 * np.sqrt(2)
    return tets / 5.5


def tetrahedralize(domain: Domain, *, h_max: float, h_min: float | None = None,
                   size_fn: Callable[[np.ndarray], np.ndarray] | None = None,
                   max_nodes: int = 1_000_000, max_iter: int = 1000,
                   ttol: float = 0.1, dptol: float = 0.001, seed: int = 0) -> Mesh:
    """Tetrahedralize an SDF domain with 3D DistMesh.

    Without ``size_fn`` the target edge length is ``h_max``.  With ``size_fn``
    the sizes are clipped to ``[h_min, h_max]``.  A ``FunctionDomain`` thinner
    than bbox/32 everywhere may be falsely rejected as having no interior.
    Never returns an element with non-positive volume.
    """
    if not isinstance(domain, Domain):
        raise TypeError("domain must be an admesh3d Domain")
    h_max = _finite_pos(h_max, "h_max")
    h_min = h_max if h_min is None else _finite_pos(h_min, "h_min")
    if h_min > h_max:
        raise ValueError(f"h_min ({h_min}) > h_max ({h_max})")
    if isinstance(max_iter, bool) or not isinstance(max_iter, numbers.Integral) or max_iter < 1:
        raise ValueError("max_iter must be >= 1")
    for name, v in (("ttol", ttol), ("dptol", dptol)):
        if not (isinstance(v, numbers.Real) and np.isfinite(v) and v > 0):
            raise ValueError("ttol/dptol must be > 0")
    if not isinstance(max_nodes, numbers.Integral) or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    fn = _checked_size_fn(size_fn) if size_fn is not None else None

    lo, hi, pts, d = _sample_grid(domain)
    est = estimate_nodes(domain, pts, d, h_min, h_max, fn)
    if est > max_nodes:
        raise ValueError(f"estimated {est:.0f} nodes exceeds max_nodes={max_nodes}; "
                         "increase h_max or max_nodes")

    h_base = h_min if fn is not None else h_max
    cap = float(np.min(hi - lo)) / 4.0
    h0 = h_base
    if h0 > cap:
        h0 = cap
        warnings.warn(f"element size {h_base} exceeds a quarter of the smallest bbox "
                      f"side; seed spacing capped at {cap:g}", UserWarning, stacklevel=2)

    res = distmesh3d(domain, h0, h_min, h_max, size_fn=fn, max_iter=int(max_iter),
                     ttol=float(ttol), dptol=float(dptol), seed=int(seed))
    if len(res.tets) == 0:
        raise ValueError("meshing produced no elements; h_max is probably too "
                         "large for this domain")
    if not res.converged:
        warnings.warn(f"max_iter={max_iter} reached before node movement fell below "
                      "dptol*h0", ConvergenceWarning, stacklevel=2)
    if res.guard_dropped:
        warnings.warn("final guard dropped interior elements with non-positive volume",
                      ConvergenceWarning, stacklevel=2)
    if res.slivers_kept:
        warnings.warn(f"{res.slivers_kept} sliver element(s) below the quality floor "
                      "were kept", ConvergenceWarning, stacklevel=2)
    return Mesh(res.nodes, res.tets, "tetra")


def extrude(mesh2d=None, *, nodes=None, triangles=None, depth=None, layers: int) -> Mesh:
    """Extrude a 2D triangle mesh into sigma-layer wedges (z from -depth to 0)."""
    if isinstance(layers, bool) or not isinstance(layers, numbers.Integral) or layers < 1:
        raise ValueError("layers must be an integer >= 1")
    layers = int(layers)
    if mesh2d is not None:
        if nodes is not None or triangles is not None:
            raise ValueError("pass either mesh2d or nodes/triangles, not both")
        try:
            nodes, triangles = mesh2d.nodes, mesh2d.elements
        except AttributeError:
            raise ValueError("mesh2d needs .nodes and .elements attributes")
        if depth is None:
            if getattr(mesh2d, "depth", None) is not None:
                depth = mesh2d.depth
            elif getattr(mesh2d, "bathymetry", None) is not None:
                depth = -np.asarray(mesh2d.bathymetry, dtype=float)
    if nodes is None or triangles is None:
        raise ValueError("nodes and triangles are required")
    if depth is None:
        raise ValueError("depth is required (argument, mesh2d.depth or mesh2d.bathymetry)")
    nodes = np.asarray(nodes, dtype=float)
    triangles = np.asarray(triangles)
    if nodes.ndim != 2 or nodes.shape[1] not in (2, 3):
        raise ValueError("nodes must have shape (N, 2) or (N, 3)")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or len(triangles) == 0:
        raise ValueError("triangles must have shape (T, 3) with T >= 1")
    if not np.issubdtype(triangles.dtype, np.integer):
        raise ValueError("triangles must be integer indices")
    triangles = triangles.astype(np.int64)
    n = len(nodes)
    depth = np.asarray(depth, dtype=float)
    if depth.ndim == 0:
        depth = np.full(n, float(depth))
    if depth.shape != (n,):
        raise ValueError(f"shape mismatch: depth has shape {depth.shape}, expected ({n},)")
    if triangles.min() < 0 or triangles.max() >= n:
        raise ValueError("triangle indices out of range")
    bad = np.flatnonzero(~(depth > 0) | ~np.isfinite(depth))
    if len(bad):
        i = int(bad[0])
        raise ValueError(f"node {i} has depth {depth[i]:g} <= 0")
    xy = nodes[:, :2]
    a, b, c = (xy[triangles[:, j]] for j in range(3))
    area2 = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    if np.any(area2 == 0) or not np.all(np.isfinite(area2)):
        raise ValueError(f"triangle {int(np.flatnonzero(area2 == 0)[0])} is degenerate")
    cw = area2 < 0
    triangles = triangles.copy()
    triangles[cw] = triangles[cw][:, [0, 2, 1]]
    n3, el = sigma_extrude(xy, triangles, depth, layers)
    return Mesh(n3, el, "wedge")
