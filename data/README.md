# Data

Input images and every run output behind the paper. Files are stored **as produced** by the scripts that wrote them and are
never edited by hand; their SHA-256 are in `../SHA256SUMS.txt`. The folder and file names are those written by the scripts
(see "Naming" below), so that each file can be matched with the script that produced it.

## images/

Three segmented 4096 × 4096 px RGBA rasters. Pore pixels are pure green (rule: G ≥ 80, G ≥ R + 30, G ≥ B + 30); everything
else is solid. Lengths are in pixels and permeabilities in px².

| file | SHA-256 | sample | configurations analysed in the paper |
|---|---|---|---|
| `patch_y3800_x3800_c0_mask.png` | `852930e4…5f042b` | R1 | does not percolate in either direction (rejected by the percolation filter) |
| `patch_y7600_x19000_c0_mask.png` | `1de0a105…273490` | R2 | R2-V (flow top → bottom) |
| `patch_y7600_x53200_c0_mask.png` | `319c26fe…87395b1f` | R3 | R3-H (flow left → right), R3-V (flow top → bottom) |

These are the images of the first release of this repository (`structures.tar.gz`, July 2026), byte-identical.

## results/

### Reference results (finite elements)

| folder | content | written by |
|---|---|---|
| `fix_round/FR-E1_convergence_tol0.25/` | **reference runs**: R3-H, R2-V, R3-V on mesh levels 1–5, Taylor–Hood (`_p2`) and Mini (`_mini`), contour tolerance 0.25 px, band 8 %, reservoirs b = 0.1 L_s — Table 1, Figure 4, Figure 6, all k_D quoted in the text | `scripts_as_run/studies/run_fr_E.sh` → `fr_reservoir.py` |
| `fix_round/FR-E2_window_scan_tol0.25/` | cropping-band scan 2–16 %, level 2 — Figure 5a, A_k range | `run_fr_E.sh` |
| `fix_round/FR-A_reservoir/` | reservoir width b/L_s = 0.025–0.4, level 2, contour tolerance **1.5 px** — Figure 2b | `run_fr_all.sh`, `run_fr_A.sh` |
| `fix_round/FR-B_convergence/` | mesh convergence of the reference definition at contour tolerance 1.5 px | `run_fr_all.sh` |
| `fix_round/FR-C_geometry/` | contour tolerance 0.1–2 px — Figure 2c | `run_fr_all.sh`, `run_fr_D.sh` |
| `fix_round/FR-F_lateral_bc/` | no-slip vs free-slip lateral walls over the window scan | `run_bc_study.py` (with `drc_audit/`) |
| `fix_round/FR-H_fields_P1/` | velocity and pressure fields on level 5 (`*_L5_fields.npz`), per-triangle arrays (`*_tri.npz`) and flow descriptors (`p1_descriptors.json`) — Figures 8–11 | `fr_fields.py`, `p1_descriptors.py` |
| `fix_round/fixround_summary.json` | summary table of FR-A, FR-C, FR-E2 read by the figure scripts | `make_figures_fr.py` |

Each run JSON records its options (`level`, `tol`, `band`, `pad_frac`, `element`, `bc`), the mesh size (`n_tri`, `ndof`),
the LU residual, the flux balance and the permeability under several pressure-drop definitions. **The permeability reported
in the paper is `kD_iface_flux`** (flux-weighted pressures one pixel inside each reservoir, full width W, sample length L_s).

### Independent finite-volume code (OpenFOAM v1912)

| folder | content | written by |
|---|---|---|
| `fix_round/FR-G_openfoam_cell_v2_alienware/` | periodic cylinder cell, c = 0.05–0.4, 1/h = 25–200: case dictionaries (`system/`, `constant/`, `0/`, `polyMesh/boundary`), `checkMesh` and solver logs — Figure 3c | `scripts_as_run/studies/b5_openfoam_cell.py` |
| `fix_round/FR-H_openfoam_window*/` | the R2-V window with pixel coarsening 2–16, residual and relaxation variants | `b6_openfoam_window.py` |
| `fix_round/FR-I_subwindow_pairs_tile512/` | 11 percolating 512 px sub-windows of R2-V: contour finite-element geometry (`geometry.geo_unrolled`) and pixel finite-volume cases — Figure 5b, Figure 7, Supplementary Table 7 | `b7_subwindow_pair.py`, `b7_analyse.py` |
| `fix_round/FR-J_pixel_convergence_t1_5/` | pixel finite-volume refinement (supersampling 1–3) of sub-window t1_5 | `b8_pixel_convergence.py` |

Solver logs larger than 1 MB are stored compressed (`log.simpleFoam.xz`; read with `xzcat`).

### Benchmarks and supporting studies

| file in `benchmarks/` | content | written by |
|---|---|---|
| `b1_poiseuille.json`, `b9_poiseuille_p1.json` | plane Poiseuille channels through the whole image pipeline (first-release definition; reference definition) — Figure 3a, Suppl. Table 4 | `b1_poiseuille.py`, `b9_poiseuille_p1.py` |
| `b2_periodic_cell.json` | periodic cylinder array vs the Hasimoto/Sangani–Acrivos series — Figure 3b, Suppl. Table 6 | `b2_periodic_cell.py` |
| `b3_array_chain.json`, `b4_legacy_on_array.json` | synthetic cylinder array through the image chain; first-release pipeline on the same array | `b3_*.py`, `b4_*.py` |
| `b5_openfoam_cell.json` (`_PARTIAL`: an earlier incomplete sweep) | OpenFOAM vs Taylor–Hood on the cylinder cell — Figure 3c, Suppl. Table 6 | `b5_openfoam_cell.py` |
| `b6_openfoam_window_R2-V_coarsen{4,8,16}.json` | OpenFOAM on the R2-V window | `b6_openfoam_window.py` |
| `b7_subwindow_pair_R2-V_tile512{,_analysis}.json`, `b8_pixel_convergence_R2-V_t1_5.json` | pixel FV vs contour FE on matched sub-windows; refinement study — Figures 5b, 7 | `b7_*.py`, `b8_*.py` |
| `c2_gmsh_size_field.json` | Gmsh threshold size field — Suppl. Fig. 1 | `c2_gmsh_size_field.py` |
| `c3_wellposedness.json` | cluster census of the three images and null-space test of the discrete operator — Suppl. Table 2 | `c3_wellposedness.py` |
| `c4_minres_p1.json` | MINRES vs LU: iterations, memory, true residuals — Figures 6, 7, Suppl. Table 3 | `c4_minres_p1.py` |
| `regression_R3H_L1.json` | regression value used while developing the reference scripts | — |
| `logs/` | console logs of the benchmark runs | — |

| folder | content | written by |
|---|---|---|
| `audit_evidence/indep_sens/` | uniform-velocity vs imposed-pressure inlet at three contour tolerances, solved with an independent re-implementation on the meshes of the first-release pipeline — Figure 2a | `independent_check.py` (one run per file; it prints the JSON) |
| `audit_evidence/poiseuille/` | independent Poiseuille check of the first-release pipeline — Suppl. Fig. 2 | `independent_check.py` |
| `audit_evidence/critical_throat.json` | critical throat widths 26.9 / 60.5 / 792 px | `critical_throat.py` |
| `audit_evidence/indep_summary.json`, `claims_recomputed.json` | summary of the independent re-solutions; values of the first-release results recomputed from its raw outputs | summaries of `independent_check.py` runs |
| `legacy_convergence/` | `global_summary.csv` of the five mesh levels run with the first-release pipeline (July 2026) and its level-5 percolating domains (`.npz`, compressed; SHA-256 of each original `.npy` in the `.json`, checked on reconstruction) — Figure 1, connected fractions | `legacy/pipeline/` |
| `legacy_poiseuille_rerun/` | first-release pipeline on four Poiseuille channels — Suppl. Fig. 2, Suppl. Table 5 | `scripts_as_run/validation/run_poiseuille_validation.sh` (writes `pipe_val.py`) |
| `figure_sources/` | the numbers plotted in Figures 4 (file `figure2`) and 8–11 | `reproduce/02_figures_from_results.py` |

All scripts named above are in `../reproduce/scripts_as_run/studies/` unless a path is given.

## Naming

The names come from the scripts that wrote the files and are kept so that file and script always match:

* `FR-A` … `FR-J` are the study series of the revision of the paper (A reservoir width, B and E1 mesh convergence, C contour
  tolerance, E2 cropping window, F lateral boundary condition, G OpenFOAM cylinder cell, H fields and OpenFOAM window,
  I sub-window pairs, J pixel convergence); `fix_round` is the folder that holds them.
* `P1` (as in `FR-H_fields_P1`, `p1_descriptors.json`, `c4_minres_p1.json`) denotes the reference definition used for every
  result of the paper: pressure-driven flow between open reservoirs of width 0.1 L_s, contour tolerance 0.25 px, band 8 %,
  Taylor–Hood elements, k_D = μ Q L_s / (W ΔP) with flux-weighted pressures. It is the default of PercoFlow.
* `tol`, `band`, `pad` in file names are the contour tolerance (px), the cropping band (fraction of the side) and the
  reservoir width (fraction of L_s); `R3H`, `R2V`, `R3V` are the configurations.
* `alienware` and `audit_evidence` name the workstation and the folder where the files were first written.

## What is not stored

Meshes, OpenFOAM meshes and OpenFOAM fields (about 5.7 GB) are not stored; they can be regenerated with the scripts.
The list, with sizes and the command that regenerates each one, is in [LARGE_FILES.md](LARGE_FILES.md).
