"""Sigma-layer prism extrusion of a 2D triangle mesh."""
from __future__ import annotations

import numpy as np


def sigma_extrude(nodes2d: np.ndarray, triangles: np.ndarray, depth: np.ndarray,
                  layers: int):
    """Extrude triangles into wedges.

    ``nodes2d`` is (N, 2), ``triangles`` (T, 3) counter-clockwise seen from +z,
    ``depth`` (N,) positive.  Node ``k*N + i`` sits at z = -depth_i * k / layers.
    Wedge (t, k) is stored at row ``k*T + t`` as
    [lower triangle (layer k+1), upper triangle (layer k)], which has a
    positive signed volume under the internal convention.
    """
    n = len(nodes2d)
    t = len(triangles)
    frac = np.arange(layers + 1, dtype=float)[:, None] / layers        # (L+1, 1)
    z = -frac * depth[None, :]                                          # (L+1, N)
    xy = np.broadcast_to(nodes2d[None, :, :], (layers + 1, n, 2))
    nodes = np.concatenate([xy, z[:, :, None]], axis=2).reshape(-1, 3)
    k = np.arange(layers, dtype=np.int64)[:, None, None]
    tri = triangles[None, :, :]
    upper = tri + k * n
    lower = tri + (k + 1) * n
    elements = np.concatenate([lower, upper], axis=2).reshape(layers * t, 6)
    return nodes, elements
