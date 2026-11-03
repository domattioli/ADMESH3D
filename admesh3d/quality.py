"""Element volumes, quality metrics and the quality report.

Orientation convention (same as VTK): a tetrahedron is positive when
det(p1-p0, p2-p0, p3-p0) > 0.  A wedge is positive when the right-hand normal
of its triangle (0,1,2) points toward triangle (3,4,5).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_WEDGE_TETS = np.array([[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]])


def tet_signed_volumes(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Signed volume of each tetrahedron (positive under the VTK convention)."""
    p = nodes[elements]
    a, b, c = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]
    return np.einsum("ij,ij->i", np.cross(a, b), c) / 6.0


def wedge_signed_volumes(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Signed volume of each wedge: sum of three consistently oriented tets."""
    total = np.zeros(len(elements))
    for t in _WEDGE_TETS:
        total += tet_signed_volumes(nodes, elements[:, t])
    return total


def tet_quality(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Normalized radius ratio 3*r_in/R_circ (1 for a regular tetrahedron)."""
    p = nodes[elements]
    a, b, c = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]
    vol = np.abs(np.einsum("ij,ij->i", np.cross(a, b), c)) / 6.0
    area = 0.5 * (np.linalg.norm(np.cross(a, b), axis=1)
                  + np.linalg.norm(np.cross(b, c), axis=1)
                  + np.linalg.norm(np.cross(c, a), axis=1)
                  + np.linalg.norm(np.cross(b - a, c - a), axis=1))
    num = (np.sum(a * a, 1)[:, None] * np.cross(b, c)
           + np.sum(b * b, 1)[:, None] * np.cross(c, a)
           + np.sum(c * c, 1)[:, None] * np.cross(a, b))
    with np.errstate(divide="ignore", invalid="ignore"):
        r_in = 3.0 * vol / area
        r_circ = np.linalg.norm(num, axis=1) / (12.0 * vol)
        q = 3.0 * r_in / r_circ
    return np.where(np.isfinite(q) & (vol > 0), q, 0.0)


def _tri_quality(p0, p1, p2):
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    e2 = (np.sum((p1 - p0) ** 2, 1) + np.sum((p2 - p1) ** 2, 1)
          + np.sum((p0 - p2) ** 2, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(e2 > 0, 4.0 * np.sqrt(3.0) * area / e2, 0.0)


def _wedge_parts(nodes: np.ndarray, elements: np.ndarray):
    """Return (q_tri_bottom_top_min, t, ell) for each wedge."""
    p = nodes[elements]
    q_bot = _tri_quality(p[:, 0], p[:, 1], p[:, 2])
    q_top = _tri_quality(p[:, 3], p[:, 4], p[:, 5])
    t = np.mean([np.linalg.norm(p[:, i + 3] - p[:, i], axis=1) for i in range(3)], axis=0)
    edges = []
    for base in (0, 3):
        for i, j in ((0, 1), (1, 2), (2, 0)):
            edges.append(np.linalg.norm(p[:, base + i] - p[:, base + j], axis=1))
    ell = np.mean(edges, axis=0)
    return np.minimum(q_bot, q_top), t, ell


def wedge_aspect(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Mean vertical edge length over mean triangle edge length."""
    _, t, ell = _wedge_parts(nodes, elements)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ell > 0, t / ell, 0.0)


def wedge_quality(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """min(q_tri_top, q_tri_bot) * min(t/l, l/t); 0 for non-positive volume."""
    q_tri, t, ell = _wedge_parts(nodes, elements)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where((t > 0) & (ell > 0), np.minimum(t / ell, ell / t), 0.0)
    q = q_tri * ratio
    q = np.where(wedge_signed_volumes(nodes, elements) > 0, q, 0.0)
    return np.clip(np.nan_to_num(q), 0.0, 1.0)


@dataclass(frozen=True)
class QualityReport:
    """Summary of element quality (all quality values lie in [0, 1])."""
    n_nodes: int
    n_elements: int
    element_type: str
    q_min: float
    q_mean: float
    q_worst5: float
    total_volume: float

    def __str__(self) -> str:
        rows = [
            ("element type", self.element_type),
            ("nodes", f"{self.n_nodes}"),
            ("elements", f"{self.n_elements}"),
            ("q_min", f"{self.q_min:.4f}"),
            ("q_mean", f"{self.q_mean:.4f}"),
            ("q_worst5", f"{self.q_worst5:.4f}"),
            ("total volume", f"{self.total_volume:.6g}"),
        ]
        width = max(len(k) for k, _ in rows)
        return "\n".join(f"{k:<{width}}  {v}" for k, v in rows)


def make_report(nodes, elements, element_type) -> QualityReport:
    """Build a QualityReport for a tetra or wedge mesh."""
    if element_type == "tetra":
        q = tet_quality(nodes, elements)
        vol = tet_signed_volumes(nodes, elements)
    else:
        q = wedge_quality(nodes, elements)
        vol = wedge_signed_volumes(nodes, elements)
    k = max(1, int(np.ceil(0.05 * len(q))))
    worst = np.sort(q)[:k]
    return QualityReport(
        n_nodes=len(nodes), n_elements=len(elements), element_type=element_type,
        q_min=float(q.min()), q_mean=float(q.mean()), q_worst5=float(worst.mean()),
        total_volume=float(vol.sum()))
