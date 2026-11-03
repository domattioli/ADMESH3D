<!-- provenance: author=domattioli model=claude-sonnet-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Implementation Plan: Prototype 3D ADMESH Mesher

**Branch**: `001-prototype-3d-mesher` | **Date**: 2026-11-02 | **Spec**: [spec.md](./spec.md)

## Summary
Turn admesh3d from placeholder into prototype: (1) tet meshing of SDF domains via 3D DistMesh (Persson-Strang) on scipy Delaunay with sliver removal; (2) sigma-layer prism extrusion from a 2D triangle mesh; (3) quality report (normalized radius ratio for tets, wedge metric for prisms); (4) ASCII .vtu write/read. API mirrors admesh2D keyword names: `tetrahedralize(domain, *, h_max, h_min=None, ...)`.

## Technical Context
**Language/Version**: Python >=3.10
**Primary Dependencies**: numpy, scipy (spatial.Delaunay). Optional: admesh2D (extrusion input only, duck-typed)
**Storage**: files (.vtu ASCII XML, stdlib xml.etree)
**Testing**: pytest
**Target Platform**: any CPython (laptop)
**Project Type**: single library
**Performance Goals**: cube h=0.1 <60 s (SC-001); extrude 10k tris x10 layers <5 s (SC-006)
**Constraints**: never return inverted elements; node safety cap 1e6 (estimated before meshing); version <1.0
**Scale/Scope**: ~7 modules, prototype

## Constitution Check
`.specify/memory/constitution.md` is an unfilled template; no gates. **N/A — PASS** (re-checked after Phase 1).

## Project Structure
### Documentation
```
specs/001-prototype-3d-mesher/
├── plan.md  research.md  data-model.md  quickstart.md
├── contracts/api.md
└── tasks.md
```
### Source
```
admesh3d/
├── __init__.py      # re-exports public API, __version__ = "0.1.0"
├── domain.py        # Box, Sphere, Cylinder, Union/Intersection/Difference, FunctionDomain
├── distmesh3d.py    # 3D DistMesh: seed, Delaunay, forces, projection, sliver removal
├── extrude.py       # sigma-layer prism extrusion
├── quality.py       # tet radius ratio, wedge quality, volumes, QualityReport
├── vtu.py           # ASCII .vtu write/read
├── errors.py        # ConvergenceWarning
└── api.py           # Mesh dataclass, tetrahedralize(), extrude(), validation
tests/
├── test_domain.py
├── test_distmesh3d.py   # cube, sphere (US1, SC-002..004)
├── test_extrude.py      # US4
├── test_quality.py      # US2; regular tet/prism = 1.0
├── test_vtu.py          # US3 roundtrip
├── test_api.py          # validation errors (FR-007)
├── test_admesh2d_smoke.py  # pytest.importorskip("admesh") duck-typing smoke
├── test_readme.py       # SC-005: README block executes
└── test_version.py      # __version__ < 1.0 and equals pyproject
```
**Structure Decision**: flat package `admesh3d/` at repo root (existing layout). pyproject: `dependencies = ["numpy", "scipy"]`, optional `test = ["pytest"]` (register pytest marker `slow`), `admesh2d = ["admesh2D"]`.

## Phases
**Phase 0 – Research**: research.md.
**Phase 1 – Design**: data-model.md, contracts/api.md, quickstart.md.
**Phase 2 – Implementation order**:
1. Setup: pyproject deps, version 0.1.0, README experimental banner (FR-009).
2. Foundation: domain.py, quality.py, `Mesh` in api.py.
3. US1 (P1): distmesh3d.py + `tetrahedralize()` + cube/sphere tests.
4. US4 (P2): extrude.py + `extrude()`; duck-typed admesh2D input (`.nodes`, `.elements`, `.depth`/`.bathymetry`).
5. US2 (P3): `Mesh.quality_report()` for tet and wedge.
6. US3 (P4): vtu.py, `Mesh.to_vtu` / `read_vtu`.
7. Polish: README quickstart, perf checks SC-001/SC-006.

## Complexity Tracking
None.
