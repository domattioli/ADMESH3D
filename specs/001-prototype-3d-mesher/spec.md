<!-- provenance: author=domattioli model=claude-opus-5-5 effort=medium date=2026-11-02 skill=speckit-pipeline repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Feature Specification: Prototype 3D ADMESH Mesher

**Feature Branch**: `001-prototype-3d-mesher`
**Created**: 2026-11-02
**Status**: Draft
**Input**: User description: "a prototype 3d version of admesh. admesh3D is a PyPI name-reservation placeholder for a 3D-element extension of ADMESH (2D mesh generator, sibling repo ../ADMESH, published as admesh2D). Prototype should make admesh3d produce 3D meshes."

## Clarifications

### Session 2026-11-02

- Q: Element type? → A: Both. Tetrahedra filling a closed volume is the core path (P1). Prism layers extruded from a 2D triangle mesh is a second path (P2).
- Q: Domain input? → A: Built-in primitives (box, sphere, cylinder), set operations on them, and user-supplied signed-distance functions. No surface files.
- Q: Output format? → A: VTK unstructured grid (.vtu), readable by ParaView.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Mesh a simple 3D domain (Priority: P1)

A researcher gives the library a closed 3D domain and a maximum edge length (the uniform target), optionally with a minimum edge length and a size function. The library returns a valid 3D mesh of that domain: node coordinates and element connectivity.

**Why this priority**: This is the core claim of the package. Without it, admesh3D stays a placeholder.

**Independent Test**: Mesh a unit cube with a uniform target size. Check that every element has positive volume and that the element volumes sum to the domain volume within tolerance.

**Acceptance Scenarios**:

1. **Given** a unit cube and a uniform target edge length of 0.1, **When** the user requests a mesh, **Then** the result contains nodes and tetrahedral elements, all with positive volume, and the total volume equals 1.0 within 1%.
2. **Given** a sphere of radius 1 and h = 0.15, **When** the user requests a mesh, **Then** every boundary node lies on the sphere surface within 1% of h, q_mean >= 0.6, q_min > 0.05, and total volume is within 3% of 4/3 pi.

---

### User Story 4 - Extrude a 2D mesh into prism layers (Priority: P2)

A coastal modeler gives a 2D triangle mesh (nodes, triangles, and a depth per node) and a layer count. The library returns a 3D mesh of triangular prisms between the free surface (z = 0) and the bed (z = -depth), with layers spaced evenly in the vertical (sigma layers).

**Why this priority**: 3D shallow-water models use layered prism meshes. This links admesh3D to admesh2D output.

**Independent Test**: Extrude a 2-triangle square with constant depth 10 and 5 layers. Check 10 prisms, 24 nodes, all with positive volume, total volume = area × 10.

**Acceptance Scenarios**:

1. **Given** a 2D mesh with N nodes and T triangles and L layers, **When** the user extrudes it, **Then** the result has N×(L+1) nodes and T×L prisms, all with positive volume.
2. **Given** a node with zero or negative depth, **When** the user extrudes, **Then** the library rejects the input with an error that names the node.
3. **Given** an admesh2D mesh object, **When** admesh2D is installed, **Then** the user can pass it directly; admesh2D is not required otherwise.

---

### User Story 2 - Inspect mesh quality (Priority: P3)

The researcher asks for a quality report on the returned mesh: element count, node count, minimum, mean and worst-5% (mean of the lowest 5%) element quality on a 0 to 1 scale (1 = regular element).

**Why this priority**: ADMESH's value is element quality. The prototype must show quality numbers to be judged against admesh2D.

**Independent Test**: Run the report on the cube mesh. Check that it returns all listed fields and that mean quality is in range 0 to 1.

**Acceptance Scenarios**:

1. **Given** a returned mesh, **When** the user requests a quality report, **Then** the report lists node count, element count, and minimum, mean and worst-5% quality.

---

### User Story 3 - Save the mesh for other tools (Priority: P4)

The researcher writes the mesh to a common 3D mesh file format that standard viewers open.

**Why this priority**: Users need to view and reuse results; not needed to prove meshing works.

**Independent Test**: Save the cube mesh, read the file back, and check that node and element counts match.

**Acceptance Scenarios**:

1. **Given** a returned mesh, **When** the user saves it, **Then** a file is written that reloads with identical node and element counts.

---

### Edge Cases

- Maximum edge length larger than the domain: return a coarse valid mesh, not an error.
- Minimum edge length larger than maximum: reject with a clear error message.
- Domain with no interior, a non-finite or degenerate bounding box, or a signed-distance function returning non-finite values: reject with an error naming the input.
- SDF with no interior, non-finite or degenerate bounding box, or non-finite SDF values on the sample grid: reject with an error naming the offending input.
- Target size so small that node count would exceed a safety limit (default 1,000,000 nodes): reject before meshing, with the estimated count.
- Iteration does not converge (max_iter reached with interior movement above dptol*h0, or final guard drops interior tets): return the best mesh found with a `ConvergenceWarning`, never a mesh with inverted elements.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The library MUST accept a 3D domain given as a built-in primitive (box, sphere, cylinder), a union, intersection or difference of domains, or a user-supplied signed-distance function with a bounding box. Surface files (such as STL) are out of scope.
- **FR-002**: The library MUST accept a maximum target edge length h_max and an optional minimum h_min. Without a size function, h_max is the uniform target and h_min only validates (h_min <= h_max) and clips; the median edge length MUST be within 25% of h_max. With a size function, local size is clipped to [h_min, h_max].
- **FR-003**: The library MUST return node coordinates and element connectivity in one mesh object.
- **FR-004**: The library MUST NOT return elements with zero or negative volume.
- **FR-005**: The library MUST provide a quality report per User Story 2.
- **FR-006**: The library MUST save a mesh (tetrahedral or prism) as a VTK unstructured grid (.vtu) and read it back.
- **FR-011**: The library MUST extrude a 2D triangle mesh with per-node depth into L evenly spaced sigma layers of triangular prisms (User Story 4).
- **FR-012**: The quality report MUST work for both tetrahedral and prism meshes.
- **FR-007**: The library MUST reject invalid inputs (edge cases above) with error messages that name the bad input.
- **FR-008**: The public entry point MUST mirror admesh2D's style: one call taking a domain and edge-length bounds, with keyword names matching admesh2D (`h_min`, `h_max`, `seed`, `max_iter`, `ttol`, `dptol`). Signature: `tetrahedralize(domain, *, h_max, h_min=None, size_fn=None, max_nodes=1_000_000, max_iter=1000, ttol=0.1, dptol=0.001, seed=0)`; h_min defaults to h_max.
- **FR-009**: The package MUST state in its README that it is an experimental prototype, and the version MUST stay below 1.0.
- **FR-010**: Spatially varying size fields (physics-based sizing, bathymetry) are OUT of scope for this prototype; a uniform size or a user-supplied size function is enough.

### Key Entities

- **Domain**: the closed 3D region to mesh, with an inside/outside test.
- **Mesh**: nodes (x, y, z) plus elements (node index tuples) plus element type.
- **Quality report**: counts and quality statistics for one mesh.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Unit cube at target size 0.1 meshes in under 60 seconds on a laptop (test marked `slow`).
- **SC-002**: For cube and sphere test domains, 100% of elements have positive volume.
- **SC-003**: On cube and sphere test domains, q_min > 0.05 and q_mean >= 0.6.
- **SC-004**: Mesh volume matches the analytic domain volume within 1% (cube) and 3% (sphere).
- **SC-006**: Extruding a 10,000-triangle 2D mesh into 10 layers takes under 5 seconds, and its volume matches the depth-integrated area within 0.1%.
- **SC-005**: The README quickstart block runs unmodified (verified by a test that extracts and executes it).
- **SC-007**: Extruding the 2-triangle unit square (depth 10, 5 layers) gives prism q_min >= 0.1 (right isosceles triangles q_tri = 0.866, t = 2, mean edge ~1.14, so q ~ 0.49).

## Assumptions

- Python 3.10+ library, same language family as admesh2D; package name `admesh3d` already reserved.
- Prototype targets correctness and clarity, not speed or production scale.
- Reuse of admesh2D ideas (DistMesh-style relaxation) is allowed; no change to the ADMESH repo is in scope.
- No GUI; library calls only.
- Release to PyPI is out of scope for this feature.
