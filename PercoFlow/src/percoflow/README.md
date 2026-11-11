# `percoflow` package: code map

One module per step of the workflow. Every module starts with a docstring that states what it computes, the equation and
the reason. This page shows where each step and each equation lives, so that the code can be read against
[docs/THEORY.md](../../docs/THEORY.md).

## Data flow of an image run (`pipeline.run_image`)

```
pipeline.run_image(image, out, cfg)
 ├─ provenance.Manifest(...)                         opens manifest.json (written at the end, also on failure)
 ├─ imaging.load_mask      image -> boolean pore mask (rule: green | threshold | nonzero; invert)
 ├─ imaging.clean_mask     remove pore objects < min_object_size px
 ├─ imaging.crop_window    remove the inlet/outlet bands -> window
 ├─ percolation.analyse    label clusters, keep percolating ones -> selected mask + report (χ_c, counts, areas)
 ├─ percolation.critical_throat_width                -> w_c
 ├─ core.solve_window(selected, ...)
 │    ├─ geometry.pad_reservoirs   open strips b = round(0.1 L_s) before/after
 │    ├─ meshing.mesh_from_mask    geometry.extract_contours + contour_hierarchy -> Gmsh mesh (INLET/OUTLET/WALLS)
 │    ├─ meshing.load_mesh         read back and check the physical groups
 │    ├─ stokes.solve              assemble, solve, boundary integrals, energy identity
 │    └─ permeability.darcy_from_reservoirs   k_D with flux-weighted probe-line pressures
 ├─ descriptors.compute     partition, A50/A90, D_h, porosity, velocity amplification
 ├─ pipeline._save_fields / _figure                  fields.npz, fields.vtu, fields.png
 └─ benchmarks.openfoam.pixel.crosscheck_image       OpenFOAM on the face-connected domain (or sub-windows)
```

A mesh run (`pipeline.run_mesh`) skips the image steps. It goes from `meshing.load_mesh` to `stokes.solve`,
`permeability.darcy_from_boundaries`, `descriptors.compute` and `benchmarks.openfoam.triangles.crosscheck_mesh`.

## Modules

| module | responsibility | key functions | equation / rule (THEORY.md) |
|---|---|---|---|
| `__init__.py` | version, overview | `__version__` | — |
| `config.py` | every option, defaults, validation, YAML/JSON I/O | `RunConfig`, `MESH_LEVELS` | presets of h(d), §4 |
| `imaging.py` | image → mask → window; coordinate convention | `load_mask`, `clean_mask`, `crop_window` | §2, coordinates |
| `percolation.py` | cluster classification, χ_c, critical throat | `analyse`, `critical_throat_width` | §2 |
| `geometry.py` | reservoirs, contours, simplification, nesting | `pad_reservoirs`, `extract_contours`, `contour_hierarchy` | §3, §4 |
| `meshing.py` | Gmsh mesh with wall-distance size field; mesh contract | `mesh_from_mask`, `load_mesh` | h(d), §4 |
| `stokes.py` | Taylor–Hood / Mini assembly, LU / MINRES, true residual, boundary integrals | `solve`, `boundary_integrals` | weak form §5, solvers §6 |
| `permeability.py` | k_D definitions | `darcy_from_reservoirs`, `darcy_from_boundaries`, `probe_pressures` | k_D, §3 |
| `descriptors.py` | pore-scale descriptors, Reynolds bound, SI conversion | `compute`, `reynolds_bound`, `to_si` | §7 |
| `core.py` | the finite-element chain shared by runs and benchmarks | `solve_window` | — |
| `pipeline.py` | run directories, figures, fields, cross-check | `run_image`, `run_mesh` | — |
| `provenance.py` | manifest, hashes, environment, verification | `Manifest`, `verify_manifest`, `sha256` | — |
| `cli.py` | command line `percoflow …` | `main` | — |

## `benchmarks/`: verification, validation and the OpenFOAM cross-check

| module | content |
|---|---|
| `__init__.py` | `run_all` (B1–B5 and the JSON report), `synthetic_porous`, `openfoam_channel` (B3), `openfoam_porous` (B5) |
| `poiseuille.py` | B1: plane Poiseuille flow through the whole image pipeline, exact value of the definition |
| `cylinder_array.py` | B2: periodic cylinder array (FE) against the Hasimoto / Sangani–Acrivos series |
| `openfoam/environment.py` | locate OpenFOAM, run its utilities, record the version |
| `openfoam/polymesh.py` | write a one-layer polyMesh from 2D cells (pixels or triangles) |
| `openfoam/case.py` | simpleFoam case files, convergence assessment, reading fields, probe-line pressures |
| `openfoam/pixel.py` | cross-check of an image domain (one hexahedron per pore pixel) |
| `openfoam/triangles.py` | cross-check of a user mesh (one prism per triangle) |
| `openfoam/cylinder_cell.py` | B4: structured O-grid cylinder cell, FV against FE |

## Conventions used everywhere

* **Coordinates.** The pixel in row `i`, column `j` of an array with `N_y` rows is centred at `(x, y) = (j, N_y − 1 − i)`.
  Pixel edges lie on half-integers.
* **Directions.** `"lr"`: inlet = first column, flow along +x. `"tb"`: inlet = first row, flow along −y.
* **Units.** Lengths are in input units (px for images). With μ = 1 and p_in = 1, k_D is in px².
* **Booleans.** `True` = fluid (pore) in every mask.
