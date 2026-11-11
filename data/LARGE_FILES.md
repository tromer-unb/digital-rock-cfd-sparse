# Large artefacts not stored in git

They are regenerable from the images and scripts in this repository (deterministic Gmsh meshing, direct solver).

| artefact | size | where it was produced | how to regenerate |
|---|---|---|---|
| Gmsh meshes of FR-E1 (levels 1–5, 3 configurations) | 205 MB | first workstation, `QC-20261019-v1/fix-round/FR-E1_convergence_tol0.25/mesh/` | `reproduce/01_recompute_permeability.py` (PercoFlow writes `mesh/mesh.msh` when `output.save_mesh` is true) |
| Gmsh meshes of FR-E2, FR-A, FR-C | 160 / 65 / 144 MB | same machine, sibling folders | the corresponding `reproduce/scripts_as_run/studies/run_fr_*.sh` |
| OpenFOAM cases of the cylinder cell (partial B5) | 221 MB | `fix-round/FR-G_openfoam_cell/` | `percoflow benchmark --full` |
| OpenFOAM meshes and fields of B5 final, B6, B7, B8 (`polyMesh/{points,faces,owner,neighbour}`, time directories, `processor*`) and the Gmsh `.msh` files | 4.9 GB | second workstation, `fix-round/FR-G_openfoam_cell_v2_alienware`, `FR-H_openfoam_window*`, `FR-I_subwindow_pairs_tile512`, `FR-J_pixel_convergence_t1_5`; archived by the authors | `reproduce/scripts_as_run/studies/b5_*.py`, `b6_*.py`, `b7_*.py`, `b8_*.py` with OpenFOAM v1912; the case dictionaries stored here are those of the runs |

Stored here because they are small enough (and needed by `reproduce/02_figures_from_results.py`):
`fix_round/FR-H_fields_P1/{R3H,R2V,R3V}_L5_fields.npz` (14, 12, 13 MB; SHA-256 `38a9ce8f…`, `1d4f9735…`, `f32ecd1f…`)
and the per-triangle arrays `*_tri.npz`.
