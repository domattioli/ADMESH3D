<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Public API Contract (admesh3d)

```python
from admesh3d import (Box, Sphere, Cylinder, FunctionDomain, Domain,
                      tetrahedralize, extrude, Mesh, QualityReport,
                      read_vtu, __version__)
from admesh3d.errors import ConvergenceWarning   # also re-exported from admesh3d

tetrahedralize(domain: Domain, *, h_max: float, h_min: float | None = None,
               size_fn: Callable[[ndarray], ndarray] | None = None,
               max_nodes: int = 1_000_000, max_iter: int = 1000,
               ttol: float = 0.1, dptol: float = 0.001,
               seed: int = 0) -> Mesh                  # element_type="tetra"
```
- h_min defaults to h_max. No size_fn: h_max is the uniform target; h_min only validates and clips. With size_fn: size clipped to [h_min, h_max]. h0 = min(h_max, min bbox side/4) (warning if capped).
- ValueError (message names the input): h_max not finite or <=0 ("h_max must be finite and > 0"); h_min<=0; max_iter<1 ("max_iter must be >= 1"); ttol<=0 or dptol<=0 ("ttol/dptol must be > 0"); size_fn output wrong shape, non-finite or <=0 ("size_fn returned invalid values"); h_min>h_max ("h_min (x) > h_max (y)"); domain has no interior; non-finite or degenerate bbox; non-finite SDF on sample grid; estimated nodes > max_nodes (message includes estimate).
- h_max > domain extent → coarse valid mesh.
- ConvergenceWarning when max_iter reached with max interior move > dptol*h0, when the final guard drops interior tets, or when slivers below q_min are kept (research R1/R2). UserWarning when h0 is capped at bbox side/4. Converged → final mesh; non-converged → best snapshot; never returns non-positive volume. Same seed gives identical mesh.

```python
extrude(mesh2d=None, *, nodes=None, triangles=None, depth=None,
        layers: int) -> Mesh                            # element_type="wedge"
```
- `mesh2d`: any object with `.nodes` (N,2) or (N,3; x,y used), `.elements` (T,3). Depth: `depth=` argument, else `.depth` used as-is, else `.bathymetry` (admesh2D elevation, positive-up) with depth = -bathymetry. None or missing depth -> ValueError. admesh2D not imported.
- ValueError: layers not an int or <1; shape mismatch; depth<=0 ("node i has depth d <= 0", first bad node).
- Output: N·(L+1) nodes, T·L wedges, z in [-depth, 0]. Triangles reoriented CCW.

```python
Mesh.quality_report() -> QualityReport
Mesh.to_vtu(path) -> None   # ASCII, CellData "quality"; wedges written in internal order (same as VTK_WEDGE)
read_vtu(path) -> Mesh      # == Mesh.from_vtu
```
Quality: tetra 3·r_in/R_circ; wedge min(q_tri_top,q_tri_bot)·min(t/ℓ,ℓ/t). Both 0..1.
