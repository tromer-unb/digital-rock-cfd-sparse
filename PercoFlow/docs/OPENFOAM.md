# OpenFOAM cross-check (always on)

## Why

Two finite-element solvers built on the same library can only verify an implementation. OpenFOAM shares neither the
discretization (cell-centred finite volumes, SIMPLEC) nor any software with `percoflow`, so agreement with it validates the
flow solution. For this reason the cross-check runs with **every** `percoflow image` / `percoflow mesh` run and is part of
`percoflow benchmark`. It can be switched off only explicitly (`--no-openfoam` or `openfoam.enabled: false`); the run
manifest then records `OpenFOAM cross-check DISABLED by the user` and the benchmark report is marked `INCOMPLETE`.

## Install

| situation | command |
|---|---|
| Ubuntu/Debian with root | `sudo apt install openfoam` (Ubuntu 24.04 ships v1912, the reference version of this release) |
| no root (Linux x86_64) | `bash scripts/install_openfoam_userspace.sh ~/opt/openfoam-v1912` |
| other | https://www.openfoam.com/download (any recent ESI version provides `simpleFoam`, `checkMesh`, `gmshToFoam`, `createPatch`) |

Then `percoflow openfoam-check`. If it is not found automatically:
`export PERCOFLOW_OPENFOAM_ENV=/path/to/etc/bashrc` (or `--openfoam-env` on the command line). Search order:
explicit option → `PERCOFLOW_OPENFOAM_ENV` → `simpleFoam` on `PATH` → `/usr/lib/openfoam/openfoam*/etc/bashrc`,
`/opt/openfoam*/etc/bashrc`, `/usr/share/openfoam/etc/bashrc`, `~/OpenFOAM/OpenFOAM-*/etc/bashrc`,
`~/opt/openfoam-sysroot/env_openfoam.sh`.

## What is run

| input | FV mesh | boundary conditions | compared with |
|---|---|---|---|
| image | one hexahedron per pore pixel of the **face-connected** percolating domain (+ the same reservoirs), one layer, `empty` front/back | `p = p_in` / `p = 0` on the reservoir faces, `noSlip` walls | FE (Taylor–Hood, LU) on the identical face-connected mask, same `k_D` definition (flux-weighted pressures on the lines 1 px inside the reservoirs) |
| mesh | the user's triangles extruded to prisms, one layer | `p = p_in` on INLET, `p = 0` on OUTLET, `noSlip` on WALLS | FE `k` from boundary pressures |
| benchmark | structured O-grid of a periodic cylinder cell, cyclic patches, body force | — | FE value and Hasimoto series |

`simpleFoam`: steady, laminar, `bounded Gauss linear` convection (negligible: ν = 10⁴, Re ≪ 1), SIMPLEC
(`consistent yes`), GAMG for p, residual targets 1e-9 (p) and 1e-10 (U). `checkMesh` must report `Mesh OK`.

**Convergence test.** simpleFoam stops when all *initial* residuals fall below their targets. A velocity component that is
identically zero in the exact solution (u_y of a straight channel) keeps a normalized residual of order one made of round-off,
so such a run ends at `endTime`. A run is therefore accepted when simpleFoam reports convergence **or** the final pressure
residual is below 1e-9; both residuals are stored (`converged_by`, `p_initial_residual_last`, `Ux/Uy_initial_residual_last`).

**Large images.** A domain with more than `openfoam.max_cells` (250 000) pore pixels is cross-checked on square sub-windows
of `openfoam.subwindow` (512) px, closest to the centre first, up to `openfoam.max_subwindows` (2); each sub-window is made
face-connected and percolating before both solvers see it. Use `openfoam.nprocs` for MPI runs.

## How to read the result

`result.json → openfoam_crosscheck` lists, per window, `k_fe_contour`, `k_fv_pixel` and `rel_diff_fv_vs_fe`.

* Images: the two numbers are **not** expected to coincide. FE follows a smooth contour through the pixel edges, FV keeps the
  staircase; their difference measures the representation of the walls. It is largest where throats are only a few pixels
  wide. Above 35 % a window is flagged `investigate`.
* Meshes: the FV scheme on triangles is less accurate than Taylor–Hood on the same mesh; a few per cent on coarse meshes is
  FV discretization error and must shrink under refinement.
* The face-connected restriction can change the domain slightly (corner contacts are removed), so `k_fe_contour` of the
  cross-check may differ from the main `k_D` of the run; both are reported.

## Verified behaviour (this repository, OpenFOAM v1912)

* channel (exact cell-centred value): relative error 4e-13, i.e. round-off (`benchmark B3`);
* synthetic rock image (examples/): FV − FE = −1.4 %; obstacle mesh (examples/): FV − FE = +1.8 %;
* cylinder cell, FV versus FE: −0.111 % (c = 0.1) and −0.125 % (c = 0.2) (`percoflow benchmark`, B4).
