"""Signed-distance domains for 3D meshing.

A domain has a signed distance function (negative inside), a bounding box and
an optional analytic volume.  Domains can be combined with ``|`` (union),
``&`` (intersection) and ``-`` (difference).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable

import numpy as np


def _vec3(value, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have 3 components, got shape {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _as_bbox(lo: np.ndarray, hi: np.ndarray):
    return (tuple(float(v) for v in lo), tuple(float(v) for v in hi))


class Domain(ABC):
    """Abstract 3D domain described by a signed distance function."""

    @abstractmethod
    def sdf(self, p: np.ndarray) -> np.ndarray:
        """Signed distance (negative inside) for points ``p`` of shape (n, 3)."""

    @property
    @abstractmethod
    def bbox(self):
        """Bounding box ``((xmin, ymin, zmin), (xmax, ymax, zmax))``."""

    @property
    def volume(self):
        """Analytic volume if known, else ``None``."""
        return None

    def fixed_points(self) -> np.ndarray:
        """Points that must be nodes of the mesh, shape (m, 3)."""
        return np.empty((0, 3))

    def __or__(self, other: "Domain") -> "Domain":
        return Union(self, other)

    def __and__(self, other: "Domain") -> "Domain":
        return Intersection(self, other)

    def __sub__(self, other: "Domain") -> "Domain":
        return Difference(self, other)


class Box(Domain):
    """Axis-aligned box from ``lo`` to ``hi``."""

    def __init__(self, lo, hi):
        self.lo = _vec3(lo, "lo")
        self.hi = _vec3(hi, "hi")
        if not np.all(self.hi > self.lo):
            raise ValueError("Box requires hi > lo in every coordinate")

    def sdf(self, p):
        p = np.asarray(p, dtype=float)
        center = 0.5 * (self.lo + self.hi)
        half = 0.5 * (self.hi - self.lo)
        q = np.abs(p - center) - half
        outside = np.linalg.norm(np.maximum(q, 0.0), axis=1)
        inside = np.minimum(np.max(q, axis=1), 0.0)
        return outside + inside

    @property
    def bbox(self):
        return _as_bbox(self.lo, self.hi)

    @property
    def volume(self):
        return float(np.prod(self.hi - self.lo))

    def fixed_points(self):
        corners = [[x, y, z] for x in (self.lo[0], self.hi[0])
                   for y in (self.lo[1], self.hi[1])
                   for z in (self.lo[2], self.hi[2])]
        return np.array(corners, dtype=float)


class Sphere(Domain):
    """Ball with the given center and radius."""

    def __init__(self, center, radius):
        self.center = _vec3(center, "center")
        if not (np.isfinite(radius) and radius > 0):
            raise ValueError("radius must be finite and > 0")
        self.radius = float(radius)

    def sdf(self, p):
        return np.linalg.norm(np.asarray(p, dtype=float) - self.center, axis=1) - self.radius

    @property
    def bbox(self):
        return _as_bbox(self.center - self.radius, self.center + self.radius)

    @property
    def volume(self):
        return 4.0 / 3.0 * np.pi * self.radius ** 3


class Cylinder(Domain):
    """Solid cylinder centered at ``center`` with its axis along x, y or z."""

    def __init__(self, center, radius, height, axis: str = "z"):
        self.center = _vec3(center, "center")
        if not (np.isfinite(radius) and radius > 0):
            raise ValueError("radius must be finite and > 0")
        if not (np.isfinite(height) and height > 0):
            raise ValueError("height must be finite and > 0")
        if axis not in ("x", "y", "z"):
            raise ValueError("axis must be 'x', 'y' or 'z'")
        self.radius = float(radius)
        self.height = float(height)
        self.axis = axis
        self._k = "xyz".index(axis)

    def sdf(self, p):
        q = np.asarray(p, dtype=float) - self.center
        along = q[:, self._k]
        radial = np.linalg.norm(np.delete(q, self._k, axis=1), axis=1)
        d = np.stack([radial - self.radius, np.abs(along) - 0.5 * self.height], axis=1)
        outside = np.linalg.norm(np.maximum(d, 0.0), axis=1)
        inside = np.minimum(np.max(d, axis=1), 0.0)
        return outside + inside

    @property
    def bbox(self):
        half = np.full(3, self.radius)
        half[self._k] = 0.5 * self.height
        return _as_bbox(self.center - half, self.center + half)

    @property
    def volume(self):
        return float(np.pi * self.radius ** 2 * self.height)


class _Binary(Domain):
    def __init__(self, a: Domain, b: Domain):
        self.a, self.b = a, b


class Union(_Binary):
    def sdf(self, p):
        return np.minimum(self.a.sdf(p), self.b.sdf(p))

    @property
    def bbox(self):
        (alo, ahi), (blo, bhi) = self.a.bbox, self.b.bbox
        return (tuple(np.minimum(alo, blo).tolist()), tuple(np.maximum(ahi, bhi).tolist()))


class Intersection(_Binary):
    def sdf(self, p):
        return np.maximum(self.a.sdf(p), self.b.sdf(p))

    @property
    def bbox(self):
        (alo, ahi), (blo, bhi) = self.a.bbox, self.b.bbox
        return (tuple(np.maximum(alo, blo).tolist()), tuple(np.minimum(ahi, bhi).tolist()))


class Difference(_Binary):
    def sdf(self, p):
        return np.maximum(self.a.sdf(p), -self.b.sdf(p))

    @property
    def bbox(self):
        return self.a.bbox


class FunctionDomain(Domain):
    """Domain defined by a user supplied signed distance function and bbox."""

    def __init__(self, sdf: Callable[[np.ndarray], np.ndarray], bbox):
        try:
            lo, hi = bbox
            lo = np.asarray(lo, dtype=float)
            hi = np.asarray(hi, dtype=float)
        except (TypeError, ValueError):
            raise ValueError("bbox must be ((xmin, ymin, zmin), (xmax, ymax, zmax))")
        if lo.shape != (3,) or hi.shape != (3,) or not (
                np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise ValueError("bbox must be finite with 3 components per corner")
        self._sdf = sdf
        self._bbox = _as_bbox(lo, hi)

    def sdf(self, p):
        return np.asarray(self._sdf(np.asarray(p, dtype=float)), dtype=float)

    @property
    def bbox(self):
        return self._bbox
