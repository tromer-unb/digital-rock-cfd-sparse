# How to audit this work

This page explains the rules that make every result of the paper checkable, and gives a checklist for an independent auditor.

## Rules followed

1. **Outputs are stored as produced.** Every file in `data/results/` was written by a script and has not been edited. The
   numbers in the paper are read from these files; none was typed by hand.
2. **Scripts are kept as run.** The code that wrote each output is in `reproduce/scripts_as_run/`, unmodified, even where it
   refers to the original directories. The reproduction scripts adapt paths at run time instead of editing the code.
3. **Acceptance criteria come first.** Every benchmark, study and reproduction script states its criterion in its header,
   written before the run (for example `|k_D / stored − 1| < 1e-6` in `reproduce/01_recompute_permeability.py`).
4. **Every PercoFlow run is self-describing.** It writes `manifest.json` with the full configuration, the SHA-256 of the
   input and of every output, the package version, the git commit and whether the working tree was modified, the library
   versions, the OpenFOAM version and the times. A disabled cross-check is recorded as such.
   Runs made from the unpacked Supplementary Data archive (not a git checkout) record `"commit": null`; the archive name
   carries the commit it was built from (`digital-rock-cfd-sparse-<commit>.zip`).
5. **No silent defaults.** Every option that can change a number is a field of PercoFlow's `config.RunConfig`; a mistyped
   option is an error.
6. **Integrity is checkable.** `SHA256SUMS.txt` lists the SHA-256 of every file of the repository; `PercoFlow/SHA256SUMS.txt`
   does the same for the software.
7. **Independent code.** The finite-element results are checked against OpenFOAM, a finite-volume code that shares no
   discretization and no software with PercoFlow.

## Checklist for an auditor

| # | check | how | expected |
|---|---|---|---|
| 1 | files are intact | `sha256sum -c SHA256SUMS.txt --quiet`; `cd PercoFlow && sha256sum -c SHA256SUMS.txt --quiet` | no output |
| 2 | the software is correct on exact solutions | `cd PercoFlow && pytest -q`; `percoflow benchmark -o runs/benchmark` | 19 passed; STATUS: PASS |
| 3 | the stored permeabilities come from the stated code and options | `python reproduce/01_recompute_permeability.py --quick` (or full) | PASS, \|Δk/k\| < 1e-6 |
| 4 | the figures show the stored data | `python reproduce/02_figures_from_results.py`; `python reproduce/03_figures_and_tables_as_run.py` | PASS |
| 5 | a number of the text is where the paper says | pick it in `docs/PROVENANCE.md`, open the file, read the key | same value |
| 6 | a run is what its manifest says | `percoflow verify runs/01_permeability/<run>` | all hashes OK |
| 7 | the options of a stored run | read `level`, `tol`, `band`, `pad_frac`, `element`, `bc` in its JSON | as stated in the paper |
| 8 | how a stored output was made | find the script in `reproduce/scripts_as_run/README.md` and read it | — |

## Worked example: the permeability of R3-H

1. The paper gives k_D = 18.281 px² for R3-H (Taylor–Hood, level 5).
2. `docs/PROVENANCE.md`, row "k_D R3-H / R2-V / R3-V", points to
   `data/results/fix_round/FR-E1_convergence_tol0.25/R3H_L5_pad0.1_tol0.25_band0.08_p2.json`, key `kD_iface_flux`.
3. Read it:
   ```bash
   python -c "import json; d=json.load(open('data/results/fix_round/FR-E1_convergence_tol0.25/R3H_L5_pad0.1_tol0.25_band0.08_p2.json')); print(d['kD_iface_flux'], d['level'], d['tol'], d['band'], d['element'])"
   ```
4. The file name and its keys give the options (level 5, tolerance 0.25 px, band 0.08, reservoir 0.1 L_s, element p2).
   The script that wrote it is `reproduce/scripts_as_run/studies/fr_reservoir.py`, driven by `run_fr_E.sh`.
5. Recompute it with PercoFlow: `python reproduce/01_recompute_permeability.py` compares the new value with this file (the
   finest levels need at least 32 GB of memory).

## What is not in the repository

Meshes, OpenFOAM meshes and OpenFOAM fields (about 5.7 GB) are not stored; `data/LARGE_FILES.md` lists them with the command
that regenerates each one.
