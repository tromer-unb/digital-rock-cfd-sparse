# Reproduce

Three scripts recompute the results of the paper from the data in this repository and compare them automatically with the
stored values. Run them from the repository root, in the environment made by `tools/setup_environment.sh`.

| step | script | reproduces | time (1 core) | acceptance (declared in the script before running) |
|---|---|---|---|---|
| 1 | `01_recompute_permeability.py [--quick] [--windows]` | k_D of R3-H, R2-V, R3-V on mesh levels 1–5 (Table 1, Figure 4, all k_D and A_k in the text); with `--windows` the cropping-band scan (Figure 5a, A_k range) | 4 min (`--quick`, 17 GB); full levels and windows longer, ≥ 32 GB | \|k_D / stored − 1\| < 10⁻⁶ for every run |
| 2 | `02_figures_from_results.py [outdir]` | Figures 4, 8, 9, 10, 11 from the stored results and fields | 30 s | every plotted number identical to `data/results/figure_sources/` |
| 3 | `03_figures_and_tables_as_run.py [--keep DIR]` | Figures 1, 2, 3, 5, 6, 7, Supplementary Figures 1–2 and the Supplementary Tables, by running the figure scripts exactly as run | 15 s | tables byte-identical; figures identical on the page rendered at 200 dpi |
| all | `run_all.sh [--quick]` | integrity check, PercoFlow tests and benchmarks, then steps 1–3 | 16 min with `--quick` | every step PASS |

Outputs go to `runs/` (not tracked). Each PercoFlow run of step 1 is a complete run directory with `manifest.json`
(configuration, SHA-256 of input and outputs, versions, git commit); check one with `percoflow verify <run directory>`.

## Reference output of this release

`reference_output/` holds the output of the steps run on the released commit (recorded in `reference_output/PROVENANCE.txt`
and in every manifest): `01_permeability/comparison.json` and the run manifests, `02_figures/*_source.json` and the console
logs of each step. A new run should give the same numbers (to the tolerances in `../docs/REPRODUCIBILITY.md`); timings,
dates and output hashes legitimately differ.

## How step 3 works

The figure scripts of the paper were written for the directory layout of the project in which they ran, and refer to it by
absolute paths. They are kept **verbatim** in `scripts_as_run/figures/` and `scripts_as_run/studies/qc_scripts/`, so that
what is in this repository is what produced the paper. Step 3 rebuilds that layout in a temporary directory with links into
`data/` (the level-5 domains are unpacked from their compressed archive and each is checked against the SHA-256 of the
original file), copies the scripts there, rewrites only the absolute path prefix, runs them and compares the outputs with
`paper/`.

## The scripts that produced the data

[`scripts_as_run/`](scripts_as_run/README.md) holds every script that wrote a file of `data/results/`, as run, with a map from
script to output.
