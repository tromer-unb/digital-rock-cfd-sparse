# Tutorial

This tutorial goes from an empty Linux machine to (1) a permeability computed on an image of the paper, (2) the reproduction
of every result of the paper and (3) an analysis of your own images. Each step shows the command, what it prints and how to
read it. Commands are run from the repository root.

Contents
1. What you need
2. Get the repository and check it
3. Install the environment
4. Your first run: sample R2, flow top → bottom
5. Read the outputs
6. Add the independent OpenFOAM cross-check
7. Reproduce the paper
8. Follow one number from the paper to its source
9. Use PercoFlow on your own images
10. Troubleshooting

---

## 1. What you need

* Linux x86_64 with Python ≥ 3.10 (the reference is 3.12.3) and about 5 GB of free disk space.
* Memory: about 17 GB for any run on the 4096 × 4096 px images of the paper (the image stage dominates), more for the
  finest mesh levels (use a machine with at least 32 GB for levels 3–5).
* For the figures: LaTeX with the Times fonts and `pdftoppm` (Debian/Ubuntu:
  `sudo apt install texlive-latex-extra texlive-fonts-recommended dvipng cm-super poppler-utils`).
* For the finite-volume cross-check and benchmarks: OpenFOAM v1912 (step 3.2).

## 2. Get the repository and check it

From GitHub:

```bash
git clone https://github.com/tromer-unb/digital-rock-cfd-sparse.git
cd digital-rock-cfd-sparse
```

or from the Supplementary Data archive of the paper:

```bash
unzip digital-rock-cfd-sparse-<version>.zip
cd digital-rock-cfd-sparse
```

Both have the same content. Check that every file is intact:

```bash
sha256sum -c SHA256SUMS.txt --quiet      # prints nothing when every file is OK
```

## 3. Install the environment

### 3.1 Python packages and PercoFlow

```bash
bash tools/setup_environment.sh          # creates .venv with the exact versions used for the paper
source .venv/bin/activate
percoflow --version                      # percoflow 1.0.0
```

The script installs `PercoFlow/requirements-lock.txt` (every package pinned) and then PercoFlow itself from `PercoFlow/`.
Activate the environment (`source .venv/bin/activate`) in every new terminal.

For runs that should match the stored results exactly, use one thread for the linear algebra:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
```

### 3.2 OpenFOAM (optional)

| situation | command |
|---|---|
| Ubuntu/Debian with root | `sudo apt install openfoam` (Ubuntu 24.04 ships v1912) |
| no root | `bash PercoFlow/scripts/install_openfoam_userspace.sh` (extracts the packages under `~/opt`) |

If OpenFOAM is not on the `PATH`, tell PercoFlow where its environment script is and check it:

```bash
export PERCOFLOW_OPENFOAM_ENV=$HOME/opt/openfoam-sysroot/env_openfoam.sh   # path printed by the installer
percoflow openfoam-check                                                    # "version": "v1912"
```

Without OpenFOAM, add `--no-openfoam` to the `percoflow image` and `percoflow mesh` commands; the manifest of the run
records that the cross-check was switched off.

## 4. Your first run: sample R2, flow top → bottom

The three images of the paper are in `data/images/` (see [data/README.md](../data/README.md)). Sample R2 percolates only from
top to bottom (configuration R2-V of the paper). Compute its permeability on the coarsest mesh level:

```bash
percoflow image data/images/patch_y7600_x19000_c0_mask.png -o runs/R2-V --direction tb --level 1 --no-openfoam
```

After about 30 seconds it prints:

```
k_D = 60.767   (Q = 70.9339, dP = 0.980358, L = 3440, W = 4096)
OpenFOAM cross-check: disabled_by_user
```

What the numbers mean (lengths in pixels, viscosity μ = 1):

| symbol | meaning | value |
|---|---|---|
| `k_D` | Darcy permeability k_D = μ Q L_s / (W ΔP), in px² | 60.767 |
| `Q` | flux through the sample (integrated, not imposed) | 70.93 |
| `dP` | difference of the flux-weighted pressures one pixel inside the inlet and outlet reservoirs | 0.980 |
| `L` | sample length L_s after removing the 8 % bands at the inlet and outlet sides | 3440 px |
| `W` | full width of the window | 4096 px |

The paper gives 60.904 px² for R2-V at level 5. The difference from 60.767 is the mesh error of level 1 (0.2 %); Table 1 and
Figure 4 of the paper show how the value converges with the mesh level (`--level 1` … `--level 5`).

This run reproduces a stored result of the paper exactly:

```bash
python -c "import json; print(json.load(open('runs/R2-V/result.json'))['permeability']['k_D'])"
python -c "import json; print(json.load(open('data/results/fix_round/FR-E1_convergence_tol0.25/R2V_L1_pad0.1_tol0.25_band0.08_p2.json'))['kD_iface_flux'])"
# both print 60.76696061428415
```

## 5. Read the outputs

```
runs/R2-V/
  manifest.json       configuration, SHA-256 of the input and of every output, versions, git commit, times
  result.json         every number of the run
  mask_selected.png   the pore space that was actually solved (only the percolating cluster)
  mesh/mesh.msh       the Gmsh triangle mesh
  fields.npz, .vtu    velocity and pressure at the mesh vertices (open the .vtu in ParaView)
  fields.png          pressure and speed maps
```

The main blocks of `result.json`:

| block | what to look at | R2-V, level 1 |
|---|---|---|
| `percolation` | how many clusters, which ones were kept, connected fraction χ_c | 281 clusters: 1 percolating, 17 dead ends, 263 isolated; χ_c = 0.538 |
| `critical_throat_width_px` | width of the narrowest passage that the flow must cross | 60.5 px (as in the paper) |
| `permeability` | `k_D` and alternative definitions (`k_D_length_weighted`, `k_D_outer_faces`), `U_darcy`, reservoir pressure share | k_D = 60.767 |
| `solver` | element, unknowns, residual, mass balance, energy identity (dissipation = Δp Q) | 336 060 unknowns, residual 2e-14, energy identity 8e-16 |
| `descriptors` | active/stagnant partition, A50/A90, hydraulic diameter, velocity amplification, Reynolds bound | — |
| `openfoam_crosscheck` | finite-volume value and its difference from the finite-element value | `disabled_by_user` here |

χ_c = 53.8 % is the connected fraction quoted in the paper for R2: almost half of the pore area of R2 does not connect the
inlet to the outlet, and solving it would make the linear system singular (263 isolated clusters → 263 null vectors).

Check that the run directory is exactly what its manifest says:

```bash
percoflow verify runs/R2-V
```

A key-by-key description of every output is in [PercoFlow/docs/INPUTS_OUTPUTS.md](../PercoFlow/docs/INPUTS_OUTPUTS.md).

## 6. Add the independent OpenFOAM cross-check

Run the same case without `--no-openfoam`:

```bash
percoflow image data/images/patch_y7600_x19000_c0_mask.png -o runs/R2-V_of --direction tb --level 1
```

PercoFlow writes an OpenFOAM case of the same percolating pixels (one hexahedral cell per pixel, staircase walls), solves it
with `simpleFoam` and post-processes it with the same flux-weighted pressure definition. An image of the paper has millions
of pore pixels, so the cross-check is made on square sub-windows of 512 px nearest the centre (options
`openfoam.max_cells`, `openfoam.subwindow`, `openfoam.max_subwindows`); each sub-window is solved by both codes on the same
face-connected percolating mask. The result is in the block `openfoam_crosscheck` of `result.json` (per window:
`k_fe_contour`, `k_fv_pixel`, `rel_diff_fv_vs_fe`), and the OpenFOAM cases are kept in `runs/R2-V_of/openfoam/`.

On one core this run takes considerably longer than the finite-element run alone, because simpleFoam iterates to a tight
residual on two sub-windows of 512 × 512 px (170 000 and 217 000 pore cells); on the reference machine the two OpenFOAM
solutions took more than 25 minutes each.

The two values are not expected to be equal: the finite elements solve the flow in a domain with smooth contour walls, the
finite volumes in a domain with pixel-staircase walls. Figure 5b of the paper measures this difference on 11 matched
sub-windows (median 1.8 %, up to 32 %). Agreement to within this range shows that both solutions are correct for their own
domain; the paper validates the solver itself on the cylinder cell (Figure 3c, 0.09–0.11 %).

## 7. Reproduce the paper

The three reproduction scripts compare their results with the stored ones and print PASS or FAIL
([reproduce/README.md](../reproduce/README.md)).

**Software checks** (about 12 minutes):

```bash
(cd PercoFlow && pytest -q)                  # 19 passed
percoflow benchmark -o runs/benchmark        # B1-B5 ... STATUS: PASS
```

B1 runs exact Poiseuille channels through the whole image pipeline, B2 compares with the Hasimoto series, B3–B5 compare with
OpenFOAM. Each criterion is in the docstring of its module ([PercoFlow/docs/BENCHMARKS.md](../PercoFlow/docs/BENCHMARKS.md)).

**Step 1 — permeabilities** (Table 1, Figure 4, the k_D and A_k values of the text):

```bash
python reproduce/01_recompute_permeability.py --quick      # levels 1-2, about 4 minutes
```

Each line compares a recomputed value with the stored one:

```
R3-H L1 p2   k_D = ...  archived ...  rel +0.0e+00
...
PASS max |rel diff| = 0.0
```

Without `--quick` it runs levels 1–5 (longer, and more memory than level 2); `--windows` adds the cropping-window scan of
Figure 5a.

**Step 2 — Figures 4 and 8–11** from the stored results and fields:

```bash
python reproduce/02_figures_from_results.py                 # figures in runs/02_figures/, then PASS
```

**Step 3 — the other figures and the Supplementary Tables**, running the figure scripts of the paper verbatim:

```bash
python reproduce/03_figures_and_tables_as_run.py            # IDENTICAL / SAME RASTER per item, then PASS
```

**Everything in order:** `bash reproduce/run_all.sh --quick` (16 minutes on the reference machine) or `bash reproduce/run_all.sh` (all levels and the window scan).

## 8. Follow one number from the paper to its source

Take the permeability of R3-H, 18.281 px².

1. [docs/PROVENANCE.md](PROVENANCE.md), row "k_D R3-H / R2-V / R3-V", gives the file
   `data/results/fix_round/FR-E1_convergence_tol0.25/R3H_L5_pad0.1_tol0.25_band0.08_p2.json` and the key `kD_iface_flux`.
2. The file name gives the options: level 5, reservoir 0.1 L_s (`pad0.1`), tolerance 0.25 px, band 8 %, Taylor–Hood (`p2`).
   They are also stored inside the file.
3. [reproduce/scripts_as_run/README.md](../reproduce/scripts_as_run/README.md) names the script that wrote it
   (`fr_reservoir.py`, driven by `run_fr_E.sh`).
4. `reproduce/01_recompute_permeability.py` recomputes it with PercoFlow and compares.

The same path exists for every figure, table and quoted number. [docs/AUDIT.md](AUDIT.md) gives a checklist for a complete
audit.

## 9. Use PercoFlow on your own images

PercoFlow takes any segmented 2D image. The default pore rule is the one of the paper (pure green pixels); grey-level and
binary images are handled with a threshold or the non-zero rule. Options go on the command line or in a YAML file:

```bash
cp PercoFlow/examples/config_default.yaml my_config.yaml     # every option, with its default and a comment
# edit my_config.yaml, e.g.   image: {pore_rule: threshold, threshold: 128}
percoflow image my_image.png -o runs/my_image --config my_config.yaml --direction lr --level 3
```

From Python:

```python
from percoflow.config import RunConfig
from percoflow.pipeline import run_image

cfg = RunConfig.from_dict({"domain": {"direction": "lr"}, "mesh": {"level": 3}})
res = run_image("my_image.png", "runs/my_image", cfg)
print(res["permeability"]["k_D"], res["percolation"]["connected_fraction"])
```

Points to keep in mind, from the paper:

* report the permeability **with its definition**: inlet condition, reservoirs, pressure probe, contour tolerance, window and
  wall representation (all are in `manifest.json`);
* check mesh convergence by repeating the run at two or three levels;
* `k_D` is in px²; multiply by the square of the pixel size for m²;
* a Stokes solution does not fix the Reynolds number: `descriptors.reynolds_bound` gives the largest pressure gradient for
  which creeping flow holds for a given fluid and pixel size.

A tutorial on two small synthetic inputs, including a mesh input (`percoflow mesh`), is in
[PercoFlow/examples/README.md](../PercoFlow/examples/README.md).

## 10. Troubleshooting

| symptom | cause and fix |
|---|---|
| `percoflow: command not found` | the environment is not active: `source .venv/bin/activate` |
| `OpenFOAM not found` | set `PERCOFLOW_OPENFOAM_ENV` (step 3.2) or add `--no-openfoam` |
| LaTeX errors in steps 2 and 3 | install the LaTeX packages of section 1 (the figures use Times through LaTeX) |
| `pdftoppm: not found` in step 3 | install `poppler-utils` |
| a recomputed value differs at the 1e-12 level | floating-point order on another CPU; within the tolerance of [REPRODUCIBILITY.md](REPRODUCIBILITY.md) |
| a recomputed value differs at the 1e-4 level or more | another Gmsh or scikit-fem version: compare `environment.libraries` in the two manifests and reinstall with `tools/setup_environment.sh` |
| memory error | the images of the paper need about 17 GB, the finest levels more; use `--quick` or a larger machine |
