# Reference results of PercoFlow 1.0.0

These files were produced with the released code, from a clean checkout (no uncommitted changes), in the environment of
`requirements-lock.txt` with OpenFOAM v1912 (see [../REPRODUCIBILITY.md](../REPRODUCIBILITY.md)). They are what a user should
obtain when repeating the commands below. The git commit of the code is recorded in each `manifest.json`
(`environment.git.commit`) and in `PROVENANCE.txt`.

| file | produced by |
|---|---|
| `benchmark_report.json` | `percoflow benchmark -o <dir>` (quick suite B1–B5, STATUS: PASS) |
| `pytest.txt` | `pytest -q` (all tests, OpenFOAM included) |
| `examples/rock_lr/` | `percoflow image examples/synthetic_rock.png -o <dir> --level 3 --config docs/reference_results/reference_run.yaml` |
| `examples/rock_tb/` | the same with `--direction tb` |
| `examples/obstacles/` | `percoflow mesh examples/obstacles.msh -o <dir> --config docs/reference_results/reference_run.yaml` |
| `reference_run.yaml` | defaults, except that the large mesh and field files are not written (the numbers do not depend on it) |
| `PROVENANCE.txt` | commit, date, machine and environment of these runs |

Each example directory is a complete run directory: run `percoflow verify docs/reference_results/examples/rock_lr` to check
it. Its `manifest.json` lists the configuration, the SHA-256 of the input and of every output, and the versions.

Main numbers (px², lengths in px):

| run | k_D | OpenFOAM vs PercoFlow | connected fraction χ_c | critical throat |
|---|---|---|---|---|
| rock, left → right, level 3 | 0.305888 | median \|FV/FE − 1\| = 1.395 % | 0.964 | 2.8 px |
| rock, top → bottom, level 3 | 6.24304 | 0.460 % | 0.963 | 15.6 px |
| obstacle mesh | 0.0111712 | FV/FE − 1 = +1.841 % | — | — |
