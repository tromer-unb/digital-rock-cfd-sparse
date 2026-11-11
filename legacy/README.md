# Legacy: the first release of this repository (July 2026)

This folder keeps, unchanged, the code that produced the results of the first version of the paper. It is here so that
those results remain traceable. **For new work use PercoFlow** ([`../PercoFlow/`](../PercoFlow/README.md)).

| path | content |
|---|---|
| `original_release_2026-07/` | the files of the first public release of this repository, moved here with their git history (`git log --follow legacy/original_release_2026-07/digital_rock_cfd_sparse.py`): `digital_rock_cfd_sparse.py` (pipeline), `visualize_cfd_outputs.py` (post-processing), the two parameter guides `run_CFD_parameters.md` and `visualize_cfd_parameters.md`, the original README (`README_original_release.md`) and the original requirements file (`requirenments.txt`, name as released) |
| `pipeline/digital_rock_cfd_pipeline_sparse.py` | the same pipeline as it was run for the first version of the paper (comments in Portuguese; md5 `57afdb1f525741bffefae221febd7586`). The release in `original_release_2026-07/` is this file with the comments translated |
| `original_figure_scripts/` | the July 2026 figure scripts and command files (`comandos`) of the first version of the paper |

## What the first release computes and how it differs from the paper

| | first release (`digital_rock_cfd_sparse.py`) | paper and PercoFlow |
|---|---|---|
| inlet condition | uniform inlet velocity | pressure difference between open reservoirs (b = 0.1 L_s) |
| permeability | apparent k from the inlet opening | k_D = μ Q L_s / (W ΔP), full width W, flux-weighted pressures |
| elements and solver | Mini, MINRES | Taylor–Hood (reference) or Mini; sparse LU or MINRES with a true-residual check |
| contour tolerance | 1.5 px | 0.25 px |

Figure 2 of the paper shows why these choices matter: the inlet condition alone changes k by up to a factor 2.2, and a
1.5 px tolerance widens the throats (+3 to +11 %). The level-5 runs of this pipeline are summarized in
`../data/results/legacy_convergence/`, and its behaviour on exact Poiseuille channels is in
`../data/results/legacy_poiseuille_rerun/` (Supplementary Figure 2).

## Running it

The images of the paper were distributed with the first release as `structures.tar.gz`; they are now in
`../data/images/` (byte-identical). To run the first release as documented in its README:

```bash
cd legacy/original_release_2026-07
mkdir -p structures && cp ../../data/images/*.png structures/
python3 digital_rock_cfd_sparse.py --help
```
