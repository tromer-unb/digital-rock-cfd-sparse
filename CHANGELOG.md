# Changelog

## 2026-11-11 — reproduction package of the revised paper

The repository was restructured so that it holds everything needed to use the method, reproduce the paper and audit it.
The history of the first release is preserved (no commit was rewritten).

* Added **PercoFlow 1.0.0** (`PercoFlow/`), the software that implements the workflow of the paper, under the MIT licence:
  a verbatim copy of PercoFlow commit `3d56815` (release 1.0.0 with the MIT licence), with its tests, benchmarks,
  documentation and reference results.
* Added the input images (`data/images/`, byte-identical to the former `structures.tar.gz`) and every run output behind the
  paper (`data/results/`).
* Added the reproduction scripts (`reproduce/01_…`, `02_…`, `03_…`, `run_all.sh`), the scripts as run
  (`reproduce/scripts_as_run/`), the display items of the paper (`paper/`) and their reference outputs.
* Added documentation: tutorial, the paper explained, provenance map, reproducibility, audit guide; READMEs in every folder.
* Added `SHA256SUMS.txt`, `CITATION.cff`, `tools/` (environment set-up, checksums, Supplementary Data archive).
* Moved the files of the first release to `legacy/original_release_2026-07/` (`git mv`, history kept); removed
  `structures.tar.gz` (its images are in `data/images/`).

## 2026-07-17 — first release

`digital_rock_cfd_sparse.py`, `visualize_cfd_outputs.py`, parameter guides, README and the three study structures.
