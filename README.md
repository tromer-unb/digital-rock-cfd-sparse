# digital-rock-cfd-sparse

**Code, data and complete reproduction package of the article**

> *A Percolation-Aware Sparse Finite-Element Workflow for Rapid Hydraulic Characterization of Digital Rocks*
> Raphael M. Tromer, Alysson Martins Almeida Silva, Paola Ferreira Barbosa, Roberto Ventura Santos and
> Luiz Antonio Ribeiro Junior (corresponding author, luiz.a.ribeiro@ntnu.no)

This repository holds everything needed to **use** the method of the paper, to **reproduce** every figure, table and
number of the paper, and to **audit** how each of them was obtained:

* **PercoFlow 1.0.0**, the software that implements the workflow (folder [`PercoFlow/`](PercoFlow/README.md)): percolation
  filter, contour meshing, Taylor–Hood/Mini Stokes solver, Darcy permeability, flow descriptors and an independent OpenFOAM
  cross-check, with its own tests, benchmarks and tutorial;
* the **three segmented images** analysed in the paper ([`data/images/`](data/README.md));
* **every run output** behind the paper (one JSON file per run, OpenFOAM case set-ups and logs, flow fields);
* the **scripts** that produced them and the scripts that draw every figure and Supplementary Table, exactly as run;
* **reproduction scripts** that recompute the results and compare them automatically with the stored values;
* a **provenance map** from each number of the paper to its file and script, and SHA-256 checksums of every file.

The same content is distributed as the Supplementary Data archive of the paper (`digital-rock-cfd-sparse-<version>.zip`);
the archive is made from this repository with `tools/make_supplementary_zip.sh`, so a reader can use either one.

| | |
|---|---|
| software | PercoFlow 1.0.0 (Python ≥ 3.10, Linux; Gmsh, scikit-fem, OpenFOAM v1912 for the cross-check) |
| licence | MIT ([LICENSE](LICENSE)) |
| first steps | [docs/TUTORIAL.md](docs/TUTORIAL.md) |
| documentation index | [docs/README.md](docs/README.md) |

---

## 1. The paper in one paragraph

The permeability of a rock computed from a segmented image depends not only on the rock but also on choices made before any
equation is solved. The workflow keeps only the pore clusters that connect inlet and outlet (isolated clusters make the
discrete Stokes problem singular), turns the pixel boundary into smooth polygons, meshes them with triangles refined near
walls and in throats, solves the Stokes equations with stable mixed finite elements and a sparse solver, and defines the
permeability from the measured flux and a flux-weighted pressure drop between two open reservoirs. Verification
(Poiseuille flow to 5×10⁻⁷, Hasimoto series to 0.34 %) and an independent finite-volume code (OpenFOAM, 0.09–0.11 %) show
that the discretization error is below 0.03 %, while the definition of the problem (inlet condition, contour tolerance,
cropping window, wall representation) changes the permeability by factors up to about two. A didactic summary with the
equations is in [docs/PAPER_EXPLAINED.md](docs/PAPER_EXPLAINED.md).

## 2. Quick start (10 minutes)

```bash
git clone https://github.com/tromer-unb/digital-rock-cfd-sparse.git && cd digital-rock-cfd-sparse
sha256sum -c SHA256SUMS.txt --quiet          # every file intact (no output = OK)
bash tools/setup_environment.sh              # .venv with the exact package versions + PercoFlow
source .venv/bin/activate

# permeability of one sample of the paper (R2, flow top -> bottom, coarse mesh, ~1 min)
percoflow image data/images/patch_y7600_x19000_c0_mask.png -o runs/R2-V --direction tb --level 1 --no-openfoam

# recompute the permeabilities of the paper on mesh levels 1-2 and compare with the stored results (~5 min)
python reproduce/01_recompute_permeability.py --quick           # -> PASS, max |rel diff| = 0.0
```

Inside this repository PercoFlow is installed from the folder `PercoFlow/` by `tools/setup_environment.sh`; the install
instructions in `PercoFlow/README.md` refer to the stand-alone PercoFlow repository and give the same environment.

OpenFOAM v1912 is needed only for the finite-volume cross-check and benchmarks (`--no-openfoam` switches it off and records
that in the run manifest). Installation without root: `bash PercoFlow/scripts/install_openfoam_userspace.sh`.

## 3. Reproduce the paper

| what | command | time (1 core) | acceptance |
|---|---|---|---|
| repository integrity | `sha256sum -c SHA256SUMS.txt` | seconds | all OK |
| software tests (incl. OpenFOAM) | `cd PercoFlow && pytest -q` | 1 min | 19 passed |
| verification benchmarks B1–B5 | `percoflow benchmark -o runs/bench` | 15 min | STATUS: PASS |
| Table 1, Figure 4 (k_D, levels 1–5) | `python reproduce/01_recompute_permeability.py` | 30–40 min, 15 GB | \|Δk/k\| < 10⁻⁶ |
| Figure 5a, A_k range (window scan) | `python reproduce/01_recompute_permeability.py --windows` | +10 min | \|Δk/k\| < 10⁻⁶ |
| Figures 4, 8, 9, 10, 11 | `python reproduce/02_figures_from_results.py` | 3 min | plotted numbers identical |
| Figures 1, 2, 3, 5, 6, 7, Suppl. Figs 1–2, Suppl. Tables | `python reproduce/03_figures_and_tables_as_run.py` | 3 min | raster/byte identical |
| everything above, in order | `bash reproduce/run_all.sh` (`--quick` for levels 1–2) | ~1 h | every step PASS |

The acceptance criteria are written in each script before any result. Step-by-step explanations are in
[docs/TUTORIAL.md](docs/TUTORIAL.md); the environment and the expected numbers are in
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). The reference outputs of these commands for this release are in
[`reproduce/reference_output/`](reproduce/README.md).

## 4. Audit a number

Every number of the paper can be followed from the text to the file that holds it and the script that wrote it:

```
paper (e.g. k_D of R3-H = 18.281 px²)
  -> docs/PROVENANCE.md         row "k_D R3-H / R2-V / R3-V"
  -> data/results/fix_round/FR-E1_convergence_tol0.25/R3H_L5_pad0.1_tol0.25_band0.08_p2.json   key kD_iface_flux
  -> reproduce/scripts_as_run/studies/run_fr_E.sh  (script that wrote it)
  -> reproduce/01_recompute_permeability.py         (recomputes it with PercoFlow and compares)
```

[docs/AUDIT.md](docs/AUDIT.md) describes the rules that make this possible (stored outputs are never edited, scripts are kept
as run, acceptance criteria are declared in advance, every PercoFlow run writes a manifest with hashes and versions) and gives
a checklist for an independent auditor.

## 5. Repository map

```
README.md                  this page
PercoFlow/                 the software (PercoFlow 1.0.0, MIT): source, tests, benchmarks, tutorial, its own docs
data/
  images/                  the three segmented 4096 x 4096 px images of the paper
  results/                 every run output behind the paper (JSON, OpenFOAM set-ups and logs, fields)
reproduce/
  01_recompute_permeability.py     step 1: recompute k_D with PercoFlow and compare with data/results
  02_figures_from_results.py       step 2: Figures 4, 8-11 from the stored results and fields
  03_figures_and_tables_as_run.py  step 3: the other figures and the Supplementary Tables, scripts run verbatim
  run_all.sh                       all steps in order
  reference_output/                outputs of the steps for this release
  scripts_as_run/                  the scripts exactly as run for the paper (studies, benchmarks, figures)
paper/
  figures/, si_tables/     the figure files and Supplementary Tables of the paper (targets of the comparisons)
docs/                      tutorial, the paper explained, provenance map, reproducibility, audit guide
legacy/                    the first public release of this repository (July 2026) and the original pipeline
tools/                     environment set-up, checksum and archive generators
SHA256SUMS.txt             SHA-256 of every file
CITATION.cff, CHANGELOG.md, LICENSE
```

## 6. Relation to the first release of this repository

The first version of this repository (July 2026, `digital_rock_cfd_sparse.py`) was a velocity-driven pipeline with Mini
elements and MINRES. The paper now uses the reference definition implemented in PercoFlow (pressure-driven flow between open
reservoirs, Taylor–Hood elements, contour tolerance 0.25 px), for the reasons shown in Figure 2 of the paper. The first
release is kept unchanged in [`legacy/`](legacy/README.md), with its git history, so that its results remain traceable; use
PercoFlow for new work.

## 7. Scope

Two-dimensional, single-phase, incompressible, steady creeping flow of a Newtonian fluid. Lengths are in pixels and
permeabilities in px²; multiply by the square of the pixel size to obtain m². A 2D section can miss connections that exist in
3D, and the permeability belongs to the analysed window and to its definition, which must be reported with it.

## 8. Citation and licence

Please cite the article and the software; the metadata are in [CITATION.cff](CITATION.cff). Code, data and documentation in
this repository are released under the MIT licence ([LICENSE](LICENSE)).
Contact: Luiz Antonio Ribeiro Junior (luiz.a.ribeiro@ntnu.no) and Raphael M. Tromer.
