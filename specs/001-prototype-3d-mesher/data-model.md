<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Data Model

## Domain (domain.py)
- `Domain` (ABC): `sdf(p: (n,3)) -> (n,)` negative inside; `bbox -> ((xmin,ymin,zmin),(xmax,ymax,zmax))`; `volume -> float | None` (analytic, tests).
- `Box(lo, hi)`, `Sphere(center, radius)`, `Cylinder(center, radius, height, axis="z")`.
- `Union`=min, `Intersection`=max, `Difference`=max(a,-b); operators `|`, `&`, `-`. bbox: union hull, intersection overlap, difference = a.
- `FunctionDomain(sdf, bbox)`.
- `fixed_points() -> (m,3)`: Box returns its 8 corners; all others return an empty array.
- Validation: radius>0, hi>lo, bbox finite.

## Mesh (api.py)
- dataclass `Mesh`: `nodes (N,3) float64`, `elements (E,k) int64`, `element_type: "tetra"|"wedge"` (k=4|6).
- Invariants (`__post_init__`): indices in range; all signed volumes >0 (internal convention = VTK convention, research R5).
- Properties: `n_nodes`, `n_elements`, `volumes`, `quality`.
- Methods: `quality_report()`, `to_vtu(path)`, classmethod `from_vtu(path)`.

## QualityReport (quality.py)
- frozen dataclass: `n_nodes, n_elements, element_type, q_min, q_mean, q_worst5, total_volume`; `__str__` table.

## DistMesh iteration state
seed → [Delaunay → filter → forces → move → project]* → sliver cleanup → orient → Mesh. Non-convergence: best snapshot (highest q_min, positive volumes) + `ConvergenceWarning` (defined in `admesh3d/errors.py`).
