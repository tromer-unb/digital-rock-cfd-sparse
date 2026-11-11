# Changelog

All notable changes are listed here. Versions follow semantic versioning, and every release is a git tag.

## Licence – 2026-11-11

The code of release 1.0.0 is licensed under the MIT licence (`LICENSE`, `pyproject.toml`, `CITATION.cff`, README).
No source file, default or number changed; the version stays 1.0.0 and the reference results of 1.0.0 remain valid.

## [1.0.0] – 2026-11-06

First release of PercoFlow.

### Workflow
- Image route:
  - pore-mask rules (green, threshold, non-zero, invert) and small-object removal;
  - inlet/outlet band cropping;
  - 8-neighbour percolation analysis with a cluster report, connected fraction and critical throat width.
- Geometry and meshing:
  - open reservoirs (b = 0.1 L_s);
  - marching-squares contours with Douglas–Peucker simplification (0.25 px default) and even–odd nesting;
  - Gmsh Frontal–Delaunay mesh with a sigmoid wall-distance size field and five level presets.
- Mesh route: Gmsh meshes with INLET/OUTLET/WALLS groups, with a contract check.
- Stokes solver:
  - Taylor–Hood and Mini elements, pressure- or velocity-driven inlet;
  - sparse LU, or MINRES with a block-diagonal preconditioner and a mandatory true-residual acceptance test;
  - mass balance and the discrete energy identity stored with every solution.
- Darcy permeability with flux-weighted probe-line pressures (images) or boundary pressures (meshes); alternative
  definitions reported.
- Descriptors: hydraulic partition, A50/A90, hydraulic diameter, domain porosity, velocity amplification, Reynolds bound.

### Independent cross-check
- OpenFOAM cross-check on by default:
  - a generic extruded polyMesh writer (pixels and triangles) and a simpleFoam case;
  - flux-weighted probe-line post-processing and sub-windows for large images;
  - OpenFOAM detection and a user-space installer.

### Verification, traceability, documentation
- Benchmarks B1–B5 with criteria declared in advance.
- pytest suite: unit and verification tests, OpenFOAM tests that are never skipped, and regression tests on the example
  inputs.
- Run manifests with SHA-256 hashes, `percoflow verify`, and `SHA256SUMS.txt` for the repository.
- Documentation: theory, inputs/outputs, OpenFOAM, benchmarks, reproducibility, audit rules, glossary, code map, tutorial.
- `requirements-lock.txt` pins every Python package of the reference environment.
