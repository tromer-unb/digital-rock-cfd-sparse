# Reproducibility

What is needed to obtain the results of the paper again, what the numbers must be, and what may legitimately differ.

## Reference environment

| item | value |
|---|---|
| software | PercoFlow 1.0.0 (`PercoFlow/`, installed by `tools/setup_environment.sh`) |
| Python | 3.12.3 |
| Python packages | `PercoFlow/requirements-lock.txt`: every package pinned, transitive dependencies included; main ones: numpy 2.5.4, scipy 1.18.1, scikit-fem 12.0.2, gmsh 4.15.2, scikit-image 0.26.0, meshio 5.3.5, matplotlib 3.11.2 |
| OpenFOAM | v1912 (Ubuntu 24.04 package `openfoam 1912.200626-2build3`); needed for the cross-check, the benchmarks B3–B5 and the finite-volume studies |
| LaTeX | TeX Live 2023 with `times`, `mathptmx`, `amsmath` (the figures use Times through `text.usetex`) |
| poppler | `pdftoppm` 24.02 (step 3 compares figures on the rendered page) |
| operating system | Linux x86_64 (Linux Mint 22.3, Ubuntu 24.04 base) |
| threads | `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` |

A different Gmsh version can produce a different (equally valid) mesh and therefore slightly different numbers; install the
pinned versions to compare with the stored results.

## Commands and expected results

| command | expected result |
|---|---|
| `sha256sum -c SHA256SUMS.txt --quiet` | no output (every file OK) |
| `cd PercoFlow && pytest -q` | `19 passed` |
| `percoflow benchmark -o runs/benchmark` | `STATUS: PASS`; B1 Taylor–Hood error ≤ 5e-7, B2 ≤ 0.34 %, B3 ~4e-13, B4 FV − FE = −0.11 / −0.13 %, B5 +1.65 % (`PercoFlow/docs/reference_results/benchmark_report.json`) |
| `python reproduce/01_recompute_permeability.py --quick` | `PASS max |rel diff| = 0.0` (12 runs) |
| `python reproduce/01_recompute_permeability.py --windows` | `PASS` (30 runs at levels 1–5 + 27 window runs) |
| `python reproduce/02_figures_from_results.py` | `PASS`, five `*_source.json` identical |
| `python reproduce/03_figures_and_tables_as_run.py` | `PASS`, 8 figures `IDENTICAL`/`SAME RASTER`, 8 tables `IDENTICAL` |

Key values (Taylor–Hood, level 5, band 8 %; px²): k_D = **18.281** (R3-H), **60.904** (R2-V), **28 781** (R3-V);
A_k = 1574. The outputs of these commands for this release are in `reproduce/reference_output/`.

## What is deterministic and what is not

| quantity | identical on re-run? | why |
|---|---|---|
| pore mask, percolation report, critical throat width | yes, bit for bit | integer image operations |
| contours and Gmsh mesh | yes, with the same Gmsh version | Gmsh is deterministic for a given version and options |
| k_D, Q, fields (finite elements, LU) | yes, to ~1e-12 relative | floating-point summation order can differ between CPUs or BLAS builds |
| OpenFOAM permeabilities | yes, to the convergence tolerance | iterative solver |
| figure PDFs | no, bytes differ | the PDF writer stores the creation date; step 3 compares the rendered page instead |
| timings, host names, dates and output hashes in `manifest.json` | no | they describe the run, not the result |

To compare a re-run with a stored result, compare numbers within these tolerances, not file hashes.

## Resources

| task | time on one core | memory |
|---|---|---|
| one image at level 1 without OpenFOAM | 25 s | 2 GB |
| `01_recompute_permeability.py --quick` | about 5 min | 4 GB |
| `01_recompute_permeability.py` (levels 1–5) | 30–40 min | 15 GB at level 5 |
| benchmark suite B1–B5 | about 15 min | 4 GB |
| steps 2 and 3 | about 3 min each | 4 GB |

Studies such as the OpenFOAM sub-window pairs (B7, B8) or the full window scans took hours on several cores; their stored
outputs are reproduced by the figure scripts, and the scripts that made them are in `reproduce/scripts_as_run/studies/`.

## When the numbers differ

1. Compare `environment.libraries` in the two `manifest.json` files; another Gmsh or scikit-fem version is the usual cause.
2. Compare `config`; every option that can change a number is recorded there.
3. Check `environment.git.dirty`: a run made from a modified working tree says so.
