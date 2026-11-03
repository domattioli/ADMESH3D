<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Research: Prototype 3D Mesher

## R1 Tet generator
- **Decision**: 3D DistMesh (Persson & Strang 2004). Seed: cubic lattice, spacing h0 = (h_min if size_fn is given else h_max), capped at min bbox side/4 (cap emits `UserWarning`), jittered by +-0.01*h0 via `np.random.default_rng(seed)`; keep d<geps; with size_fn, keep each seed point with probability (h_min/h(x))³ (h clipped to [h_min,h_max]); drop seeds within 0.5·h0 of any fixed point, then append fixed points (`Domain.fixed_points()`: Box corners; others none). Loop: scipy.spatial.Delaunay; drop tets with centroid d>-geps; unique edges; bar forces F=max(L0-L,0), L0=h(mid)*Fscale*(ΣL³/Σh³)^(1/3), Fscale=1.2; p += deltat*F (deltat=0.1); project d>0 points back via numerical gradient (deps=sqrt(eps)*h0). Defaults ttol=0.1, dptol=0.001 (both parameters of `tetrahedralize`). Retriangulate when max move >ttol*h0. Stop when interior max move <dptol*h0 or niter (default 1000). Return rule: a converged run returns the final cleaned mesh. A non-converged run returns the best snapshot. Snapshot score (lower is better) = (fraction of tets with q < q_min, then -q_worst5), computed after orienting all tets positive and dropping centroid-outside tets; a snapshot qualifies only if no tet has |volume| < 1e-12·h0³; scored only at retriangulation, vectorized; if no snapshot qualifies, the final guard is applied to the latest mesh. `ConvergenceWarning` (from `admesh3d/errors.py`) is emitted when: max_iter is reached with max interior move > dptol*h0; the final guard drops interior tets; or R2 keeps slivers below q_min. The h0 cap emits `UserWarning`, not `ConvergenceWarning`.
- **Rationale**: same family as admesh2D distmesh.py; numpy+scipy only.
- **Alternatives**: TetGen/gmsh bindings (deps, rejected); advancing front (too large).

## R2 Slivers
- **Decision**: after convergence compute radius ratio q. Remove tets with q<q_min (0.05) that have at least one boundary face, provided the boundary after removal is still a closed manifold shell; all other slivers take the interior path. Interior slivers: perturb sliver nodes ±0.1h, keep better result, retriangulate, ≤5 passes. Then, for each remaining interior sliver, delete its non-fixed node with the largest |d| (deepest interior), retriangulate and recompute; at most 5 deletion rounds. Slivers still <q_min after that are removed one at a time, only if the tet has at least one boundary face and the boundary after removal is still a closed manifold shell (every boundary edge shared by exactly 2 boundary faces); otherwise keep it and emit `ConvergenceWarning`. Final guard: drop tets with signed volume <=0, orient all positive; if it drops interior tets, emit `ConvergenceWarning`. Finally drop unreferenced nodes and reindex.
- **Rationale**: DistMesh 3D leaves slivers; needed for SC-003 (worst >0.05) and FR-004.
- **Alternatives**: sliver exudation (no weighted Delaunay in scipy).

## R3 Tet quality
- **Decision**: q = 3·r_in/R_circ (1 = regular). r_in = 3V/Σface areas; R_circ = |a²(b×c)+b²(c×a)+c²(a×b)|/(12V), edge vectors from vertex 0.

## R4 Prism quality
- **Decision**: q = min(q_tri_top, q_tri_bot) × min(t/ℓ, ℓ/t). q_tri = 4√3·A/Σedge². t = mean vertical edge length, ℓ = mean triangle edge. Non-positive volume → q=0. Report `aspect` separately (sigma layers are anisotropic by design).
- **Alternatives**: Jacobian scaled quality (harder to verify).

## R5 Prism volume
- **Decision**: split into 3 tets (0,1,2,3),(1,2,3,4),(2,3,4,5) with consistent orientation; sum signed volumes. Internal convention: tet positive iff det(p1-p0, p2-p0, p3-p0) > 0 (same as VTK_TETRA). Wedge positive iff the right-hand normal of triangle (0,1,2) points toward triangle (3,4,5). VTK_WEDGE uses the same rule (vtkWedge.h: "(0,2,1) is the base ... normal points away from (3,4,5)"; parametric coords 0=(0,0,0),1=(1,0,0),2=(0,1,0),3=(0,0,1), so J = det[p1-p0,p2-p0,p3-p0] > 0), so vtu.py writes and reads wedges in internal order with no permutation. T007 tests the internal sign; T023 tests that written wedges have det[p1-p0,p2-p0,p3-p0] > 0.

## R6 Sigma extrusion
- **Decision**: z_k(i) = -depth_i·k/L, k=0..L. Node index k·N+i. Prism (t,k) = [lower layer tri (k+1), upper layer tri (k)] = [tri+(k+1)·N, tri+k·N], with triangles CCW seen from +z; this gives positive volume under R5. Vectorized.

## R7 .vtu
- **Decision**: ASCII XML via `xml.etree.ElementTree`: Points Float64, Cells connectivity/offsets/types (VTK_TETRA=10, VTK_WEDGE=13), optional CellData "quality". Reader parses same subset.
- **Alternatives**: meshio (dep, rejected).

## R8 Domain validation
- **Decision**: SDF domains are closed by construction. Sample grid = 33×33×33 points spanning the bbox (corners included) plus the bbox center. Limitation: a FunctionDomain thinner than bbox/32 everywhere may be falsely rejected; documented in the docstring. Check bbox finite with positive extent, SDF finite (no NaN/inf) on the sample grid, and SDF negative somewhere on it; else ValueError naming the input (e.g. "domain has no interior", "bbox is not finite", "sdf returned non-finite values").

## R9 Node cap
- **Decision**: V = domain.volume if known, else bbox_volume × interior fraction from the R8 sample grid. h_eff = h_max if size_fn is None else h_min. Without size_fn: estimate = V / (h_eff³/(6√2)) / 5.5 (about 5.5 tets per node). With size_fn: tets = ∫ 1/(h³/(6√2)) dV, estimated as the mean of 1/h(x)³ over sample-grid points with d(x) < 0 (h clipped to [h_min,h_max]) × V × 6√2; nodes = tets/5.5. Exceeds max_nodes (1e6) → ValueError including estimate. T009 asserts estimate within 0.5-2x actual node count for the unit cube.

## Implementation deviations (2026-11-02)
- Fscale = 1.1 (Persson 3D value), not 1.2: at 1.2 the cube did not converge and q_mean was about 0.71.
- Time-step cooling: after 0.4·max_iter iterations, deltat halves every 0.1·max_iter (`ANNEAL_START`, `ANNEAL_HALF_LIFE` in `distmesh3d.py`). Convergence is measured by node displacement, so a converged run does not guarantee zero residual force. Open for review.
- Measured (seed 0): cube h=0.1, 659 iterations, 2.1 s, q_min 0.054, q_mean 0.887. Sphere R=1 h=0.15, 611 iterations, 2.4 s, q_min 0.066, q_mean 0.900, volume -1.3%. The q_min margin over 0.05 is thin by construction (cleanup removes q < 0.05).
