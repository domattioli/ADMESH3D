# Changelog

All notable changes to this project will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-11-07

### Added
- `tetrahedralize`: tetrahedral meshes of signed-distance-function domains (3D DistMesh) with `Box`, `Sphere`, `Cylinder` and union, intersection and difference.
- `extrude`: sigma-layer prism extrusion of a 2D triangle mesh.
- `Mesh.quality_report` and ASCII `.vtu` output (`Mesh.to_vtu`, `read_vtu`).
- DistMesh progress callback and denser size-function seeding.

### Changed
- License changed from Apache-2.0 to PolyForm Noncommercial 1.0.0.
