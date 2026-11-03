<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Tasks: Prototype 3D ADMESH Mesher

**Input**: Design documents from `/specs/001-prototype-3d-mesher/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/api.md, quickstart.md
**Tests**: Included (pytest, under `tests/`).
**Organization**: Phases follow user story priority: US1 tets (P1), US4 extrusion (P2), US2 quality (P3), US3 vtu (P4).

## Format: `- [ ] [ID] [P?] [Story] Description with file path`

All paths are relative to the repo root `/Users/domattioli/Projects/MADMESHing/admesh3D/`.

## Phase 1: Setup

- [X] T001 Update `pyproject.toml`: version 0.1.0 (below 1.0), `dependencies = ["numpy", "scipy"]`, optional extras `test = ["pytest"]` and `admesh2d = ["admesh2D"]`; register pytest marker `slow`
- [X] T002 [P] Add experimental-prototype banner to `README.md` (FR-009)
- [X] T003 [P] Create package skeleton `admesh3d/__init__.py` (with `__version__ = "0.1.0"`), `admesh3d/errors.py` (defines `ConvergenceWarning`), and `tests/__init__.py`

## Phase 2: Foundational (blocks all stories)

- [X] T004 [P] Write domain tests in `tests/test_domain.py` (sdf sign, bbox, analytic volume, set ops `|` `&` `-`, `Box.fixed_points()` returns 8 corners and others none, validation of radius, hi>lo, finite bbox)
- [X] T005 Implement `Domain` ABC (incl. `fixed_points()`), `Box`, `Sphere`, `Cylinder`, `Union`, `Intersection`, `Difference`, `FunctionDomain` in `admesh3d/domain.py`
- [X] T006 Implement `Mesh` dataclass (nodes, elements, element_type, `__post_init__` invariants, n_nodes, n_elements, volumes) in `admesh3d/api.py`; depends on T007
- [X] T007 [P] Write tests in `tests/test_quality.py` that a hand-built unit tet and unit wedge (triangle (0,1,2) normal toward (3,4,5)) have signed volume > 0 and their mirror images < 0 (R5 internal convention), then implement tet and wedge signed-volume helpers in `admesh3d/quality.py`

## Phase 3: User Story 1 - Mesh a simple 3D domain (P1) MVP

**Goal**: `tetrahedralize(domain, *, h_max, h_min=None, ...)` returns a valid tet mesh.
**Independent test**: unit cube at h=0.1, all volumes positive, total volume within 1% of 1.0; sphere R=1, h=0.15: q_mean >= 0.6, q_min > 0.05, volume within 3%, boundary nodes within 1% of h of the surface.

- [X] T008 [P] [US1] Write cube and sphere tests in `tests/test_distmesh3d.py`: cube volume within 1% and all volumes > 0 (SC-002, SC-004); sphere R=1 h=0.15 q_mean >= 0.6, q_min > 0.05 (SC-003 checks run under `warnings.simplefilter("error", ConvergenceWarning)`), volume within 3%, boundary nodes within 1% of h; (h_min,h_max)=(0.05,0.1) median edge ~0.1 within 25%; size_fn test (graded size, clipped to [h_min,h_max]; node density in the fine region exceeds the coarse region); same seed reproduces identical mesh; `max_iter=1` warns `ConvergenceWarning`, returns a mesh with all volumes > 0 and no unreferenced nodes (falls back to guarded latest mesh when no valid snapshot); SC-001 cube h=0.1 under 60 s marked `@pytest.mark.slow`
- [X] T009 [P] [US1] Write validation tests in `tests/test_api.py`: h_max<=0 or non-finite, h_min<=0, h_min>h_max, max_iter<1, ttol/dptol<=0, size_fn returning NaN/negative/wrong shape, non-int layers, empty domain, always-positive `FunctionDomain` (no interior), non-finite or degenerate bbox, NaN/inf SDF (errors name the input); node cap with estimate in message; node estimate within 0.5-2x actual for unit cube; h_max > cube extent gives n_nodes >= 8, volumes > 0, volume within 10% of 1; `tetrahedralize(Box(...), 0.1)` and `extrude(nodes, tris, depth, 5)` raise TypeError (keyword-only, FR-008)
- [X] T010 [US1] Implement 3D DistMesh core in `admesh3d/distmesh3d.py`: seed (h0 = h_min if size_fn else h_max, capped at min bbox side/4 with `UserWarning`; jittered cubic lattice +-0.01*h0 via `np.random.default_rng(seed)`; with size_fn keep with probability (h_min/h)³; drop seeds within 0.5·h0 of fixed points, then append `Domain.fixed_points()`), Delaunay, force relaxation (ttol=0.1, dptol=0.001), boundary projection, best-snapshot tracking (score per research R1, at retriangulation only, vectorized). Per-iteration budget at 1.5k nodes: Delaunay <= 25 ms, forces+projection <= 10 ms
- [X] T011 [US1] Add sliver removal (per research R2: boundary slivers removed if >=1 boundary face and the boundary stays a closed manifold shell; else 5 perturbation passes, then up to 5 deepest-node deletion rounds, then keep + ConvergenceWarning; drop unreferenced nodes and reindex), orientation fix, and guarantee of no non-positive volume in `admesh3d/distmesh3d.py`
- [X] T012 [US1] Implement `tetrahedralize()` in `admesh3d/api.py` with input validation (R8 incl. NaN/inf check), node-count estimate vs `max_nodes` (R9), `size_fn` clipping, `ttol`/`dptol` parameters, and `ConvergenceWarning` rule (max_iter reached with move > dptol*h0, guard drops interior tets, or slivers below q_min are kept, R1/R2)
- [X] T013 [US1] Export `Box, Sphere, Cylinder, FunctionDomain, Domain, tetrahedralize, Mesh, ConvergenceWarning` from `admesh3d/__init__.py`

## Phase 4: User Story 4 - Prism extrusion (P2)

**Goal**: `extrude()` turns a 2D triangle mesh with per-node depth into sigma-layer prisms.
**Independent test**: 2-triangle square, depth 10, 5 layers gives 10 prisms, 24 nodes, volume = area x 10.

- [X] T014 [P] [US4] Write extrusion tests in `tests/test_extrude.py`: counts N(L+1) and TL, positive volumes, total volume, CCW reorientation, depth<=0 error naming first bad node, layers<1, shape mismatch, duck-typed object with `.nodes` (N,2) or (N,3) / `.elements` / `.depth` (used as-is) and `.bathymetry` with negative elevations (depth = -bathymetry), None/missing depth raises ValueError, SC-006 10k triangles x10 layers under 5 s with volume within 0.1%
- [X] T015 [US4] Implement sigma-layer extrusion in `admesh3d/extrude.py` (wedge = [lower tri (k+1), upper tri (k)], R6)
- [X] T016 [US4] Implement `extrude()` wrapper with validation and duck-typed `mesh2d` input (no admesh2D import) in `admesh3d/api.py`
- [X] T017 [US4] Export `extrude` from `admesh3d/__init__.py`
- [X] T018 [P] [US4] Add `tests/test_admesh2d_smoke.py`: `pytest.importorskip("admesh")`, build a small admesh2D `Mesh` and pass it to `extrude()`; check counts and positive volumes

## Phase 5: User Story 2 - Quality report (P3)

**Goal**: `Mesh.quality_report()` works for tet and wedge meshes.
**Independent test**: cube mesh report has all fields, mean quality in [0,1].

- [X] T019 [P] [US2] Extend `tests/test_quality.py`: regular tet = 1.0, regular prism = 1.0, report fields (min, mean, worst-5% = mean of lowest 5%), cube and sphere (R=1, h=0.15) q_mean >= 0.6 and q_min > 0.05 per SC-003 under `warnings.simplefilter("error", ConvergenceWarning)`, both element types, SC-007 (2-triangle square, depth 10, 5 layers: prism q_min >= 0.1)
- [X] T020 [US2] Implement tet radius-ratio and wedge quality metrics in `admesh3d/quality.py`
- [X] T021 [US2] Implement frozen `QualityReport` with `__str__` table in `admesh3d/quality.py`
- [X] T022 [US2] Add `Mesh.quality`, `Mesh.quality_report()` in `admesh3d/api.py`; export `QualityReport` from `admesh3d/__init__.py`

## Phase 6: User Story 3 - Save as .vtu (P4)

**Goal**: write and read ASCII VTK unstructured grid.
**Independent test**: cube mesh round-trips with identical node and element counts.

- [X] T023 [P] [US3] Write round-trip tests for tetra and wedge meshes (counts, coordinates, element type, CellData "quality"); written wedge (0,1,2) right-hand normal points toward (3,4,5) and det[p1-p0,p2-p0,p3-p0] > 0 (VTK parametric Jacobian), in `tests/test_vtu.py`
- [X] T024 [US3] Implement ASCII `.vtu` writer and reader using `xml.etree` in `admesh3d/vtu.py`; wedges in internal order, no permutation (R5)
- [X] T025 [US3] Add `Mesh.to_vtu(path)`, `Mesh.from_vtu(path)` in `admesh3d/api.py`; export `read_vtu` from `admesh3d/__init__.py`

## Phase 7: Polish

- [X] T026 [P] Add README quickstart block (cube tets, extrusion, quality, vtu) matching `specs/001-prototype-3d-mesher/quickstart.md` in `README.md`
- [X] T027 Add `tests/test_readme.py`: extract the first ```python block of README.md and exec it unmodified with `monkeypatch.chdir(tmp_path)` (SC-005); mark `slow`; depends on T026
- [X] T028 [P] Add `tests/test_version.py`: `admesh3d.__version__` < 1.0 and equals the `pyproject.toml` version (FR-009)
- [X] T029 Run full `pytest tests/` (including `-m slow`) and confirm SC-001..SC-007 timings and thresholds pass; record cube iteration count and wall time; if > 60 s, report before changing defaults
- [ ] T030 Move specs to DomI: run `bash ~/Projects/DomI/skills/speckit-pipeline/scripts/move_specs_to_domi.sh --apply` (REQUIRES OPERATOR CONFIRMATION: pushes to DomI; depends on all other tasks; nothing may run after it)

## Dependencies

- Phase 1 then Phase 2 then stories. US1 first (MVP).
- T006 depends on T007 (Mesh uses the volume helpers); T006 is not parallel.
- US4 depends on Phase 2 only (independent of US1). US2 depends on Phase 2; its cube-mean test needs US1. US3 depends on Phase 2; its tet round-trip test needs US1, wedge test needs US4.
- T010 -> T011 -> T012 -> T013. T015 -> T016 -> T017 -> T018. T020 -> T021 -> T022. T024 -> T025. T026 -> T027.
- T006 and T012/T016/T022/T025 all edit `admesh3d/api.py`, so run those sequentially. T007, T019-T022 all edit `quality.py`/`tests/test_quality.py`: sequential.
- T030 last.

## Parallel Examples

- Phase 2: T004 with T007 (then T006 after T007).
- US1: T008 with T009. After Phase 2, US4 tests (T014) can run beside US1.
- US2/US3 tests: T019 with T023.
- Polish: T026 with T028.

## Implementation Strategy

MVP = Phase 1 + 2 + US1 (tets). Then US4, US2, US3 in order, validating each against its independent test before moving on. Polish and the DomI move last.
