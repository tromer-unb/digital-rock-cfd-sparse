# Audit and traceability

PercoFlow is built so that any number it produces can be traced to its inputs, options and software, and checked again by
someone else. This page lists the rules the code follows, shows how to audit one result, and gives the release checklist.

## Rules the code follows

1. **No silent defaults.** Every option that can change a number is a field of `config.RunConfig`. A mistyped key raises an
   error. The complete configuration actually used is written to `manifest.json`.
2. **Hashes in, hashes out.** The SHA-256 of every input file and of every output file is recorded
   (`provenance.Manifest`). `percoflow verify` recomputes them.
3. **Software identity.** Each run records:
   - the package version and git commit, and whether the working tree had uncommitted changes;
   - the Python and library versions;
   - the OpenFOAM version and its environment script;
   - the host name, thread settings and UTC start and end times.
4. **Checks travel with the result.** Each result carries:
   - the linear residual and the mass imbalance;
   - the discrete energy identity `Φ = p_in Q`;
   - the MINRES iterations and true residual, with rejection above `solver.accept_residual`;
   - the OpenFOAM `checkMesh` status and the convergence criterion actually met.
5. **Acceptance criteria before results.** Each benchmark states its pass criterion in its module docstring and in the
   report, before any number is looked at ([BENCHMARKS.md](BENCHMARKS.md)).
6. **Independent code by default.** The OpenFOAM cross-check runs unless it is explicitly disabled, and disabling it is
   recorded in the manifest. The benchmark report is then marked `INCOMPLETE`.
7. **Failures are results too.** A run that fails still writes `manifest.json`, with `status: failed` and the error.
8. **Regression is tested.** `tests/test_regression.py` recomputes the reference numbers of the example inputs. Any change to
   the code that changes a computed number makes it fail.

## How to audit one result in five steps

1. `percoflow verify runs/X`: every output file still has the hash recorded at the end of the run.
2. Open `runs/X/manifest.json` and read `config` (options), `inputs` (hash of the image or mesh), `environment` (versions
   and commit) and `notes`.
3. Check the hash of the input file against the one in `inputs`, for instance with `sha256sum image.png`.
4. In `runs/X/result.json`, check the solver diagnostics:
   - `solver.residual` should be about 1e-12 for LU;
   - `solver.mass_imbalance` should be about 1e-13;
   - `solver.energy_identity_relerr` should be about 1e-14;
   - `openfoam_crosscheck.summary` should show both FV runs converged and `Mesh OK`.
5. Re-run the same command with the same `config` in an environment installed from `requirements-lock.txt` at the recorded
   commit, and compare `result.json`. Identical versions give identical digits ([REPRODUCIBILITY.md](REPRODUCIBILITY.md)).

## Integrity of the repository itself

`SHA256SUMS.txt` lists the SHA-256 of every tracked file. Check it with `sha256sum -c SHA256SUMS.txt`, and regenerate it with
`bash scripts/make_checksums.sh` before each release.

## Release checklist (maintainers)

* `pytest -q` passes on a clean install from `requirements-lock.txt`, with every test included (OpenFOAM and regression).
* `percoflow benchmark -o runs/bench` reports `STATUS: PASS`. Copy the report into `docs/reference_results/`.
* `CHANGELOG.md` is updated and the version is bumped in `pyproject.toml`, `src/percoflow/__init__.py` and `CITATION.cff`.
* If a computed number changed on purpose, update `tests/test_regression.py` and explain the change in the changelog.
* `bash scripts/make_checksums.sh`, then commit and tag `vX.Y.Z`.
