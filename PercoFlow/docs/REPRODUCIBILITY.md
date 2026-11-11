# Reproducibility

This page lists what is needed to obtain the reference numbers of PercoFlow 1.0.0 again, on this or another machine.

## Reference environment

| item | value |
|---|---|
| PercoFlow | 1.0.0 (git tag `v1.0.0`) |
| Python | 3.12.3 |
| Python packages | `requirements-lock.txt`, every package pinned, transitive dependencies included |
| OpenFOAM | v1912 (Ubuntu 24.04 package `openfoam 1912.200626-2build3`) |
| operating system | Linux x86_64 (Linux Mint 22.3, kernel 7.0) |
| threads | `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` for the reference runs (the tests use 2) |

## Steps

```bash
git clone https://github.com/ribeirojr-gh/PercoFlow.git && cd PercoFlow && git checkout v1.0.0
sha256sum -c SHA256SUMS.txt --quiet                      # repository files intact
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements-lock.txt && pip install --no-deps -e .
percoflow openfoam-check                                  # OpenFOAM v1912 found
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
pytest -q                                                 # 19 passed
percoflow benchmark -o runs/bench                         # STATUS: PASS
percoflow image examples/synthetic_rock.png -o runs/rock --level 3
```

Then compare with [reference_results/](reference_results/README.md):
* the numbers of `runs/bench/benchmark_report.json` with `reference_results/benchmark_report.json`;
* `runs/rock/result.json` with `reference_results/examples/rock/result.json`.

## What is deterministic and what is not

| quantity | identical on re-run? | why |
|---|---|---|
| pore mask, percolation report, critical throat width | yes, bit for bit | integer image operations and a fixed bisection |
| contours and Gmsh mesh | yes, with the same Gmsh version | Gmsh is deterministic for a given version and options; another version can produce a different (equally valid) mesh |
| k_D, Q, fields | yes to ~1e-12 relative | sparse LU in floating point; the order of operations can differ between CPUs or BLAS builds |
| OpenFOAM k | yes to the convergence tolerance (~1e-9) | iterative solver |
| timings, host, dates in manifest.json | no | they describe the run, not the result |
| output hashes in manifest.json | no for files containing floating-point text or timestamps | `percoflow verify` checks a run against its own manifest, not against another run |

To check a re-run against a reference, compare the numbers in `result.json` within the tolerances above. Hashes are not the
right comparison between two runs. `tests/test_regression.py` follows this rule.

## When the numbers differ

1. Check `environment.libraries` in both manifests. A different Gmsh or scikit-fem version is the usual cause.
2. Check `config` in both manifests. Every option that can change a number is listed there.
3. Check `environment.git.dirty`: a run made from a modified working tree says so.
