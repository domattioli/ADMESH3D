<!-- provenance: author=domattioli model=claude-haiku-5-5 effort=low date=2026-11-02 skill=speckit-checklist repo=admesh3D session=781c7800-85d3-48ff-859a-99097abe495c -->
# Mesh Quality Checklist: Prototype 3D ADMESH Mesher

**Purpose**: Unit tests for requirements writing. Checks whether the correctness and robustness requirements in spec.md are complete, clear, consistent and measurable before planning and implementation proceed.
**Created**: 2026-11-02
**Feature**: [spec.md](../spec.md)
**Domain**: mesh-quality (correctness and robustness of generated 3D meshes)
**Settings**: Depth Standard · Audience Reviewer (PR) · Focus: element validity and quality, and failure or edge-case behaviour. Defaults applied; no interactive questions asked.

## Requirement Completeness

- [ ] CHK001 - Are the element quality metrics (tetrahedral and prism) defined in the spec itself, rather than only in planning research, so FR-005 and FR-012 can be assessed from the spec alone? [Gap, Spec §FR-005, Spec §FR-012]
- [ ] CHK002 - Is a minimum acceptable quality threshold specified for prism meshes, given SC-003 sets thresholds only for tetrahedral cube and sphere domains? [Gap, Spec §SC-003]
- [ ] CHK003 - Are requirements defined for a domain that is smaller than the minimum edge length h_min, as distinct from the case where h_max exceeds the domain? [Edge Case, Gap, Spec §Edge Cases]
- [ ] CHK004 - Is rejection of an empty domain (no interior) specified in the spec, or only implied by the domain definition? [Gap, Spec §FR-001, Spec §FR-007]
- [ ] CHK005 - Are the attributes an admesh2D mesh object must expose (node coordinates, connectivity, per-node depth) specified as a requirement, not left to the implementation? [Dependency, Gap, Spec §US4]
- [ ] CHK006 - Is the user-supplied size function in scope? FR-010 says it is enough, but FR-008 defines the entry point as a domain plus two edge-length bounds only. [Conflict, Spec §FR-008, Spec §FR-010]
- [ ] CHK007 - Are determinism requirements stated (same domain and bounds must yield the same mesh across runs), or is run-to-run variation left unspecified? [Gap, Coverage]
- [ ] CHK008 - Are robustness requirements defined for degenerate sampled geometry (duplicate, coplanar or cospherical points) that the meshing process may produce? [Gap, Non-Functional]

## Requirement Clarity

- [ ] CHK009 - Is the zero-volume boundary in FR-004 quantified (a numeric tolerance for "zero"), or does "zero or negative volume" leave the threshold open? [Clarity, Spec §FR-004]
- [ ] CHK010 - Is h_min a hard lower bound on every edge, or a sizing parameter whose effect is measured only through the median edge length in FR-002? [Ambiguity, Spec §FR-002]
- [ ] CHK011 - Is the "uniform size request" in FR-002 defined when both h_min and h_max are supplied, including which target the 25% median tolerance is measured against? [Clarity, Spec §FR-002]
- [ ] CHK012 - Is "coarse valid mesh" for h_max larger than the domain quantified (for example, a minimum element count or a bound on element size)? [Clarity, Spec §Edge Cases]
- [ ] CHK013 - Are "best mesh found" and "converged" defined in the non-convergence edge case, including how candidate meshes are ranked and what counts as converged? [Ambiguity, Spec §Edge Cases]
- [ ] CHK014 - Is "worst-5% quality" defined (mean of the worst 5% of elements, or the 5th percentile value)? [Clarity, Spec §US2]
- [ ] CHK015 - Is "worst element quality" in SC-003 the same statistic as "worst-5%" in US2, or the minimum over all elements? [Consistency, Spec §SC-003, Spec §US2]
- [ ] CHK016 - Is "within 1% of the target edge length" in the sphere acceptance scenario tied to a specific size (h_min, h_max or the target)? [Ambiguity, Spec §US1]
- [ ] CHK017 - Is "depth-integrated area" in SC-006 defined, including how depth varies within a triangle? [Clarity, Spec §SC-006]
- [ ] CHK018 - Are the sigma-layer spacing rules specified for columns with differing depth, or only for constant depth? [Clarity, Spec §FR-011]
- [ ] CHK019 - Is the node-cap estimate defined with an accuracy requirement, and is the 1,000,000-node default specified as overridable or fixed? [Clarity, Spec §Edge Cases]
- [ ] CHK020 - When several nodes have zero or negative depth, does the error name the first bad node, every bad node, or either? [Clarity, Spec §US4]
- [ ] CHK021 - Is the reference hardware for the "under 60 seconds on a laptop" target (SC-001) specified, and is the timed scope defined (meshing only, or including quality report and file write)? [Measurability, Spec §SC-001]

## Requirement Consistency

- [ ] CHK022 - Are the requirement IDs ordered and continuous across the spec (FR-011 and FR-012 appear before FR-007; SC-006 appears before SC-005)? [Traceability, Spec §Requirements, Spec §Success Criteria]
- [ ] CHK023 - Does the edge case "open or self-intersecting domain surface" conflict with FR-001 and the Clarifications, which state that domains are primitives or signed-distance functions and that surface files are out of scope? [Conflict, Spec §Edge Cases, Spec §FR-001]
- [ ] CHK024 - Is the .vtu encoding (ASCII or binary) stated as a requirement, or is it left open in FR-006? [Ambiguity, Spec §FR-006]
- [ ] CHK025 - Is the sliver-removal threshold a requirement, or only a planning choice? SC-003 requires worst quality above 0.05 while planning removes elements below 0.05, so the boundary (strict or non-strict) is unstated. [Ambiguity, Spec §SC-003]
- [ ] CHK026 - Are the rules for a returned best-effort mesh (non-convergence path, which may have poor quality) consistent with FR-004 and the "never inverted" statement? [Consistency, Spec §FR-004, Spec §Edge Cases]

## Acceptance Criteria Quality

- [ ] CHK027 - Are quality thresholds (mean at least 0.6, worst above 0.05) specified for domain types beyond cube and sphere, such as cylinder and boolean combinations? [Coverage, Gap, Spec §SC-003]
- [ ] CHK028 - Is the .vtu round-trip requirement (FR-006, US3) defined in terms of coordinate precision or tolerance, or only element and node counts? [Measurability, Spec §US3]
- [ ] CHK029 - Are acceptance criteria for prism meshes (SC-006) defined for properties beyond volume, such as positive element volume and quality? [Coverage, Spec §SC-006]
- [ ] CHK030 - Is an identifier scheme for acceptance scenarios established, or can scenarios be referenced only by user story and ordinal? [Traceability, Gap]

## Scenario Coverage

- [ ] CHK031 - Are requirements defined for invalid domain parameters (zero or negative radius, zero-extent box, upper bound below lower bound) beyond the generic FR-007 statement? [Gap, Spec §FR-007]
- [ ] CHK032 - Are requirements defined for non-finite inputs (NaN or infinite bounds, edge lengths, depths or coordinates)? [Gap, Edge Case]
- [ ] CHK033 - Are recovery requirements defined when a .vtu write fails partway, including whether a partial file may remain? [Recovery, Gap, Spec §FR-006]
- [ ] CHK034 - Are requirements defined for extrusion inputs with clockwise-oriented or zero-area triangles? [Edge Case, Gap, Spec §US4]
- [ ] CHK035 - Are requirements defined for extrusion layer counts below 1 or non-integer layer counts? [Edge Case, Gap, Spec §US4]
- [ ] CHK036 - Are requirements defined when a user-supplied bounding box does not enclose the domain interior? [Edge Case, Gap, Spec §FR-001]
- [ ] CHK037 - Are requirements defined for nodes in the 2D input that are not referenced by any triangle? [Edge Case, Gap, Spec §US4]

## Non-Functional and Dependencies & Assumptions

- [ ] CHK038 - Are memory requirements specified for meshes near the node-count cap, or is only the node count limited? [Non-Functional, Gap, Spec §Edge Cases]
- [ ] CHK039 - Are the Python version and package dependency constraints stated as requirements, or only listed as assumptions? [Assumption, Spec §Assumptions]
- [ ] CHK040 - Are the convergence assumptions behind the reused DistMesh-style relaxation documented as a dependency, or left implicit in the Assumptions section? [Dependency, Spec §Assumptions]

## Notes

- Items 1 to 40 cover correctness and robustness requirements quality. Items reference spec sections or carry a quality marker ([Gap], [Ambiguity], [Conflict], [Assumption], [Traceability]).
- Check items off as resolved in spec.md: `[x]`. Add comments or findings inline.
- Items are numbered sequentially for easy reference. Each item tests the written requirement, not the implementation.
