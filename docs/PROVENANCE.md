# Provenance map: number → result file → script

This page links every display item and every number quoted in the paper to the file that holds it and to the script that
wrote it. Paths are relative to `data/results/` unless they start with `reproduce/`, `paper/` or `data/`; scripts without a
folder are in `reproduce/scripts_as_run/studies/` (runs) or `reproduce/scripts_as_run/figures/` (figures).

Figure numbers are those of the article. The figure files in `paper/figures/` are named in the order in which they were made,
so the file name is given as well (it differs for Figures 2–5). Status **R** means reproducible from this repository: the
data and the script are here, and one of the three reproduction scripts regenerates the item and compares it with the
stored one:
* `reproduce/01_recompute_permeability.py` — the permeabilities (Table 1, Figure 4 values, window scan);
* `reproduce/02_figures_from_results.py` — files figure2 and figure8–11;
* `reproduce/03_figures_and_tables_as_run.py` — files figure1, figure3–7, figureS1–S2 and every `paper/si_tables/*.tex`.

Values in the paper are taken from the result files as written; none was typed by hand.

## Display items

| display item | file | content | data | script | status |
|---|---|---|---|---|---|
| Figure 1 | figure1.pdf | samples and percolating domains | `data/images/`; level-5 percolating domains and χ_c of the first-release runs (`legacy_convergence/`) | `scripts_as_run/figures/make_fig1_fig2_split.py` | **R** |
| Figure 2 | figure3.pdf | sensitivity to inlet condition, reservoir width, contour tolerance | (a) `audit_evidence/indep_sens/`; (b) `fixround_summary.json` ← `FR-A_reservoir` (run at contour tolerance **1.5 px**); (c) same ← `FR-C_geometry` | `scripts_as_run/figures/make_fig3_fig4_sens.py` | **R** |
| Figure 3 | figure5.pdf | Poiseuille (reference definition), cylinder array, OpenFOAM vs FE | `benchmarks/b9_poiseuille_p1.json`, `b2_periodic_cell.json`, `b5_openfoam_cell.json` (cases: `fix_round/FR-G_openfoam_cell_v2_alienware/`) | `scripts_as_run/figures/make_fig5_vv.py`; runs: `scripts_as_run/studies/b9_*.py`, `b2_*.py`, `b5_*.py` | **R** |
| Figure 4 | figure2.pdf | mesh convergence (TH and Mini) | `fix_round/FR-E1_convergence_tol0.25/*.json` | `reproduce/02_figures_from_results.py` (`fig2`) | **R** |
| Table 1 | in main.tex | N_tri, N_dof, k_D by level | same | `reproduce/01_recompute_permeability.py` | **R** |
| Figure 5 | figure4.pdf | (a) crop-band scan, (b) pixel FV vs contour FE on 11 sub-windows of R2-V | (a) `fixround_summary.json` ← `FR-E2_window_scan_tol0.25`; (b) `benchmarks/b7_subwindow_pair_R2-V_tile512.json`, `b8_pixel_convergence_R2-V_t1_5.json` (cases: `FR-I_…`, `FR-J_…`) | `scripts_as_run/figures/make_fig3_fig4_sens.py`; runs: `scripts_as_run/studies/b7_*.py`, `b8_*.py` | **R** |
| Figure 6 | figure6.pdf | N_dof, solve time and memory of LU and MINRES | `FR-E1_convergence_tol0.25/*.json`, `benchmarks/c4_minres_p1.json` | `scripts_as_run/figures/make_fig6_cost.py`; run: `scripts_as_run/studies/c4_minres_p1.py` | **R** |
| Figure 7 | figure7.pdf | MINRES and Mini accuracy; matched sub-window timing | `benchmarks/c4_minres_p1.json`, `b7_subwindow_pair_R2-V_tile512.json` | `scripts_as_run/figures/make_fig6_cost.py` | **R** |
| Figure 8 | figure8.pdf | flow fields | `fix_round/FR-H_fields_P1/*_L5_fields.npz` | `reproduce/02_figures_from_results.py` (`fig8`) | **R** |
| Figure 9 | figure9.pdf | hydraulic partition | `FR-H_fields_P1/p1_descriptors.json` | `fig9` | **R** |
| Figure 10 | figure10.pdf | heterogeneity, A50/A90 | `FR-H_fields_P1/*_tri.npz`, `p1_descriptors.json` | `fig10` | **R** |
| Figure 11 | figure11.pdf | integrated response | `p1_descriptors.json` | `fig11` | **R** |
| Suppl. Fig. 1 | figureS1.pdf | Gmsh threshold size field | `benchmarks/c2_gmsh_size_field.json` | `scripts_as_run/figures/make_figures_npj.py S2`; run: `scripts_as_run/studies/c2_gmsh_size_field.py` | **R** |
| Suppl. Fig. 2 | figureS2.pdf | first-release pipeline vs Taylor–Hood on Poiseuille channels | `legacy_poiseuille_rerun/`, `audit_evidence/poiseuille/` | `scripts_as_run/studies/qc_scripts/make_figures_qc.py --fig E3` | **R** |
| Suppl. Tables 2, 3, 4, 6, 7 | `si_tables/*.tex` | clusters, MINRES, Poiseuille (reference definition), cylinders, sub-windows | `benchmarks/c3_*.json`, `c4_*.json`, `b2_*.json`, `b5_*.json`, `b9_*.json`, `b7_*.json` | `scripts_as_run/figures/make_si_tables.py` | **R** |

## Numbers quoted in the text

| quantity | value | source |
|---|---|---|
| k_D R3-H / R2-V / R3-V (TH, level 5) | 18.281 / 60.904 / 28 781 px² | `FR-E1…/{R3H,R2V,R3V}_L5_…_p2.json` → `kD_iface_flux` |
| change between levels 4 and 5 | −0.024 / −0.023 / −0.005 % | same files, levels 4 and 5 |
| Mini vs TH, level 5 | −2.5 / −1.5 / −0.35 % | `…_L5_…_mini.json` vs `_p2.json` |
| A_k (band 8 %, level 5) | 1574 | ratio of the two L5 values above |
| A_k range over bands 2–16 % (level 2) | 878–1705 | `FR-E2_window_scan_tol0.25/*.json` + `FR-E1 …_L2_…` |
| factors 1.06 / 1.83 / 1.92, jump +79 % | window scan | `FR-E2…` |
| free-slip lateral walls +0.3 / +0.3 / −0.03 % | | `FR-F_lateral_bc/*.json` |
| reservoir-width sensitivity (< 0.16 %, +1.2 %, −10.7 %) | | `FR-A_reservoir/*.json` — run at contour tolerance **1.5 px** |
| tolerance sensitivity (2.8 / 3.9 / 11.0 %, 30 % at 2 px, < 0.003 % at 0.1 px) | | `FR-C_geometry/*.json` |
| χ_c 98.6 / 53.8 / 97.8 %, porosity 59.5 / 38.4 / 52.1 % | | χ_c: `legacy_convergence/convergence_mesh_05_extra_fine/global_summary.csv` (`selected_fraction_of_cropped_pores` 0.9857 / 0.5377 / 0.9780; identical classification in the reference runs and in PercoFlow); porosity: `p1_descriptors.json`, `figure_sources/figure11_source.json` |
| critical throat 26.9 / 60.5 / 792 px | | `audit_evidence/critical_throat.json` (also PercoFlow `critical_throat_width_px`) |
| D_h 1009 / 236 / 1005 px; amplification 210 / 96 / 5.7; partition; A50/A90; p* outside [0,1] | | `p1_descriptors.json`, `figure_sources/figure10_source.json` |
| U = Q/W, reservoir share 5.7 % (R3-V) / 0.07 % (R3-H) | | `FR-H_fields_P1/*_L5_fields.json` (`U_darcy`, `dP_iface_flux`) |
| k_D·D_h and the Reynolds bound | | `p1_descriptors.json` (`kD_times_Dh_px3`) |
| Poiseuille 4.9e-7 | | `benchmarks/b9_poiseuille_p1.json` → `si_tables/poiseuille_p1.tex`; reproduced by PercoFlow benchmark B1: 4.947e-7 |
| Hasimoto 0.34 % for c ≤ 0.15 | | `benchmarks/b2_periodic_cell.json`; PercoFlow B2: 0.339 % |
| OpenFOAM vs FE 0.087–0.111 %, FV convergence 0.041 % | | `benchmarks/b5_openfoam_cell.json` (cases, logs, checkMesh: `fix_round/FR-G_openfoam_cell_v2_alienware/`) |
| 11 sub-windows, median 1.8 %, −32 %, factor 225, Richardson order 2.8, 230× wall time | | `benchmarks/b7_subwindow_pair_R2-V_tile512{,_analysis}.json`, `b8_pixel_convergence_R2-V_t1_5.json` → `si_tables/subwindows.tex` |
| MINRES iterations, memory, residuals, 14 of 15 runs | | `benchmarks/c4_minres_p1.json` → `si_tables/minres.tex`, `minres_summary.tex` |
| cluster counts of R2-V (281 / 1 / 17 / 263, 37.4 %), kernel test 1.7e-18 | | `benchmarks/c3_wellposedness.json` → `si_tables/clusters.tex` |

## Results of the first release

The values obtained with the first release of this repository (uniform inlet velocity, apparent permeability normalized by
the inlet opening, contour tolerance 1.5 px; for example k = 110.9 / 270.3 / 3.80×10⁴ px² and A_k = 342.5) come from the
pipeline in `legacy/pipeline/` and the figure scripts in `legacy/original_figure_scripts/`. They were recomputed from the raw
outputs of those runs (`legacy_convergence/`, `audit_evidence/claims_recomputed.json`). Figure 2 of the paper compares the
two definitions.

## Where the files were produced

The runs were made on two Linux workstations of the authors between July and November 2026 (the folder suffix
`_alienware` names the second one). Files were copied byte for byte; identity was checked with SHA-256 when they were
brought together. Large artefacts that are not stored (meshes, OpenFOAM meshes and fields) are listed in
`data/LARGE_FILES.md` with the command that regenerates each of them.
