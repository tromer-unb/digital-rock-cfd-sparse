# Scripts as run

These are the scripts that produced the files in `data/results/` and the figures and tables in `paper/`, **exactly as they
were run**. They are kept unmodified so that every stored output can be matched with the code that wrote it. Two
consequences:

* they refer to the directory layout and absolute paths of the machines on which they ran (for example
  `.../CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1/` and `.../CFD/PROJETO/convergencia/`), and several import the first-release
  pipeline from `legacy/pipeline/`. `../03_figures_and_tables_as_run.py` shows how to run the figure scripts from this
  repository; the study scripts need the same treatment or a re-pointing of their path constants;
* their headers and comments use the labels of the revision (see "Labels" below) and record the tool that generated them.

To recompute results, use PercoFlow and the three scripts in `../`; to see how a stored result was obtained, read the
script listed here.

## figures/ — figure and table scripts of the paper

| script | writes |
|---|---|
| `make_fig1_fig2_split.py` | `figure1.pdf` (Figure 1) |
| `make_fig3_fig4_sens.py` | `figure3.pdf` (Figure 2), `figure4.pdf` (Figure 5) |
| `make_fig5_vv.py` | `figure5.pdf` (Figure 3) |
| `make_fig6_cost.py` | `figure6.pdf`, `figure7.pdf` (Figures 6, 7) |
| `make_figures_npj.py S2` | `figureS1.pdf` (Supplementary Figure 1) |
| `make_si_tables.py` | `si_tables/*.tex` |
| `figstd.py` | common style (imported) |

Figures 4 and 8–11 are drawn by `../02_figures_from_results.py`, the repository version of `studies/make_figures_p1.py`
(identical plotting code; paths point to `data/`).

## studies/ — studies, benchmarks and their drivers

| script | writes (under `data/results/`) |
|---|---|
| `fr_reservoir.py` (driven by `run_fr_all.sh`, `run_fr_A.sh`, `run_fr_D.sh`, `run_fr_E.sh`) | `fix_round/FR-A…FR-E2/*.json`: one JSON per run of the reference definition |
| `run_bc_study.py` (with the helper package `drc_audit/`) | `fix_round/FR-F_lateral_bc/` |
| `fr_fields.py`, `p1_descriptors.py` | `fix_round/FR-H_fields_P1/` |
| `make_figures_fr.py` | `fix_round/fixround_summary.json` |
| `b1_poiseuille.py`, `b9_poiseuille_p1.py` | `benchmarks/b1_poiseuille.json`, `b9_poiseuille_p1.json` |
| `b2_periodic_cell.py`, `b3_array_chain.py`, `b4_legacy_on_array.py` | `benchmarks/b2…b4_*.json` |
| `b5_openfoam_cell.py` | `benchmarks/b5_openfoam_cell.json`, `fix_round/FR-G_openfoam_cell_v2_alienware/` |
| `b6_openfoam_window.py` | `benchmarks/b6_*.json`, `fix_round/FR-H_openfoam_window*/` |
| `b7_subwindow_pair.py`, `b7_analyse.py` | `benchmarks/b7_*.json`, `fix_round/FR-I_subwindow_pairs_tile512/` |
| `b8_pixel_convergence.py` | `benchmarks/b8_*.json`, `fix_round/FR-J_pixel_convergence_t1_5/` |
| `c2_gmsh_size_field.py`, `c3_wellposedness.py`, `c4_minres_p1.py` | `benchmarks/c2…c4_*.json` |
| `independent_check.py` | `audit_evidence/indep_sens/`, `audit_evidence/poiseuille/` (an independent Mini/Taylor–Hood solver run on the meshes of the first-release pipeline) |
| `critical_throat.py` | `audit_evidence/critical_throat.json` |
| `pipe_val.py` | `legacy_poiseuille_rerun/` (first-release pipeline with the one-line change of `validation/pipeline_patch.diff`) |
| `compare_solvers.py`, `s1_synthetic_poc.py` | exploratory checks (solver comparison; a synthetic proof of concept, labelled synthetic in its outputs) |
| `qc_scripts/make_figures_qc.py --fig E3`, `qc_scripts/qcdata.py`, `qc_scripts/templates/house_style.py` | `figureS2.pdf` (Supplementary Figure 2) and the style used by the figure scripts |
| `drc_audit/` | helper package: reservoir padding, reference Stokes solver, periodic cell, metrics, synthetic masks |

## validation/

`run_poiseuille_validation.sh` builds the four Poiseuille channel images (h = 20, 40, 80, 160 px), applies the single change
of `pipeline_patch.diff` to the first-release pipeline (minimum contour points 8 → 4, because a rectangle has 5 vertices) and
runs the five mesh levels; its output is `data/results/legacy_poiseuille_rerun/` (Supplementary Table 5, Supplementary
Figure 2).

## Labels used in the headers

| label | meaning |
|---|---|
| P1, "protocol P1" | the reference definition used for every result of the paper (see `../../data/README.md`, "Naming") |
| FR-A … FR-J, "fix round" | the study series of the revision (see `../../data/README.md`) |
| B1–B9, C2–C4 | benchmarks and checks; the acceptance criterion of each is written in its header |
| F-001 | the inlet condition (uniform velocity vs imposed pressure) and its effect on k |
| F-003 | the dependence of k on the cropping window |
| F-008 | the comparison with an independent code (OpenFOAM) |
| F-009 | the metadata of the samples (physical pixel size); results are therefore in px² |
| F-015 | a Gmsh run that did not terminate on one case, hence the per-job time limit of `run_fr_all.sh` |
