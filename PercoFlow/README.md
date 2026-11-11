# PercoFlow

**Percolation-aware flow solver for the permeability of two-dimensional porous images and meshes, with an independent
OpenFOAM cross-check in every run.**

PercoFlow takes a segmented image of a porous material (or a triangle mesh of a pore space) and returns its Darcy
permeability, the full pore-scale velocity and pressure fields, and a set of descriptors of how the flow is organized.
It is meant for rock and soil images, filters, membranes, scaffolds and micromodels, or any 2D creeping flow through a pore
space. Every number it writes can be traced to its inputs, options and software versions, and checked again.

| | |
|---|---|
| version | 1.0.0 (2026-11-06) |
| language | Python ≥ 3.10 (Linux), with Gmsh, scikit-fem and OpenFOAM |
| licence | MIT (see [Licence](#licence)) |
| documentation | [docs/README.md](docs/README.md) (reading order and index) |

---

## 1. What the name says

**Perco** stands for *percolation*. Before any mesh is made, PercoFlow labels the pore clusters of the image and keeps only
those that join the inlet and the outlet. This is not only a saving. An isolated cluster makes the discrete Stokes system
singular, and a dead end makes a velocity-driven problem unsolvable ([docs/THEORY.md](docs/THEORY.md), Section 2).

**Flow** stands for the steady Stokes flow solved on what remains. The solver uses stable mixed finite elements: Taylor–Hood
as the reference and Mini for screening. The same domain is then solved again with OpenFOAM, a finite-volume code that
shares no discretization and no software with PercoFlow.

## 2. How it works, in one picture

```
 segmented image                                                                   user mesh (.msh)
      │                                                                                   │
      ▼  imaging.py       pore mask (green / threshold / non-zero rule)                   │
      ▼  imaging.py       remove a band at the inlet and at the outlet side               │
      ▼  percolation.py   keep only clusters that join inlet and outlet; χ_c, w_c         │
      ▼  geometry.py      add open reservoirs (b = 0.1 L_s) before and after the sample   │
      ▼  geometry.py      contours: marching squares + Douglas–Peucker (0.25 px)          │
      ▼  meshing.py       Gmsh triangle mesh, small elements near walls and in throats    │
      ▼  stokes.py        Stokes flow, Taylor–Hood (P2/P1) or Mini, LU or MINRES  ◄───────┘
      ▼  permeability.py  k_D = μ Q L_s / (W ΔP), pressures weighted by the flux
      ▼  descriptors.py   active/dead partition, A50/A90, D_h, velocity amplification
      ▼  benchmarks/openfoam/   OpenFOAM solution of the same domain, FV vs FE difference
      ▼  pipeline.py, provenance.py   result.json + manifest.json (hashes, versions, options)
```

Each step is one module of [`src/percoflow/`](src/percoflow/README.md). Each module starts with a docstring that states
what it computes, the equation it implements and why it does it that way.

## 3. Install (Linux)

```bash
git clone https://github.com/ribeirojr-gh/PercoFlow.git && cd PercoFlow
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements-lock.txt      # the exact versions used to validate this release
pip install -e ".[dev]"
percoflow --version                        # -> percoflow 1.0.0
```

OpenFOAM is needed for the cross-check, which is on by default:

| situation | command |
|---|---|
| Ubuntu/Debian with root | `sudo apt install openfoam` (Ubuntu 24.04 ships OpenFOAM v1912, the reference version) |
| no root | `bash scripts/install_openfoam_userspace.sh` (extracts the Ubuntu packages under `~/opt`) |

Then run `percoflow openfoam-check`. Details are in [docs/OPENFOAM.md](docs/OPENFOAM.md).

## 4. First run (5 minutes)

```bash
python examples/make_examples.py                                   # two small synthetic inputs (deterministic)
percoflow image examples/synthetic_rock.png -o runs/rock --level 3  # image route, flow left -> right
percoflow mesh  examples/obstacles.msh      -o runs/obstacles       # mesh route
percoflow verify runs/rock                                          # recompute every output hash of the run
```

The first command prints, in about one minute on one core:

```
k_D = 0.305888   (Q = 0.302659, dP = 0.99934, L = 202, W = 200)
OpenFOAM cross-check: full domain, median |FV/FE - 1| = 1.395%
```

`k_D` is in px² (lengths are pixels). Multiply it by the square of the measured pixel size to obtain m². The tutorial in
[examples/README.md](examples/README.md) explains every line of the output and every file of the run directory.

From Python:

```python
from percoflow.config import RunConfig
from percoflow.pipeline import run_image

cfg = RunConfig.from_dict({"domain": {"direction": "lr"}, "mesh": {"level": 3}})
res = run_image("my_segmented_image.png", "runs/my_image", cfg)
print(res["permeability"]["k_D"], res["percolation"]["connected_fraction"], res["openfoam_crosscheck"]["summary"])
```

## 5. What a run writes

```
runs/rock/
  manifest.json        what was run, on what, with which software: configuration, SHA-256 of inputs and outputs,
                       package version and git commit, library versions, OpenFOAM version, host, UTC times
  result.json          permeability (and alternative definitions), percolation report, critical throat width,
                       solver diagnostics (residual, mass balance, energy identity), descriptors, OpenFOAM cross-check
  mask_selected.png    the domain that was actually solved
  mesh/mesh.msh        the finite-element mesh
  fields.npz/.vtu      velocity and pressure at the mesh vertices (open the .vtu in ParaView)
  fields.png           pressure and speed maps
  openfoam/openfoam_crosscheck.json
```

Every key is described in [docs/INPUTS_OUTPUTS.md](docs/INPUTS_OUTPUTS.md).

## 6. Why the numbers can be trusted

| check | where | result for this release |
|---|---|---|
| Poiseuille flow through the whole image pipeline (exact solution) | benchmark B1 | Taylor–Hood error ≤ 5e-7 |
| periodic cylinder array against the Hasimoto series | benchmark B2 | ≤ 0.34 % for solid fraction ≤ 0.15 |
| OpenFOAM pixel route on a channel (exact cell-centred value) | benchmark B3 | 4e-13 (round-off) |
| OpenFOAM against PercoFlow on the cylinder cell | benchmark B4 | FV − FE = −0.11 % and −0.13 % |
| OpenFOAM against PercoFlow on a synthetic porous window | benchmark B5 | both converge; difference reported |
| reference numbers of the example inputs | `tests/test_regression.py` | unchanged |
| discrete energy identity Φ = p_in Q | every run (`energy_identity_relerr`) | ~1e-14 |

Run them with `percoflow benchmark -o runs/bench` (about 15 min on one core) and `pytest -q`. The acceptance criteria are
written in each benchmark module **before** any result. The report of this release is in
[docs/reference_results/](docs/reference_results/README.md), and [docs/BENCHMARKS.md](docs/BENCHMARKS.md) explains each case.

## 7. Traceable, reproducible, auditable

* **Traceable.** Every run writes `manifest.json`: the complete configuration, SHA-256 of every input and output file,
  versions and git commit. A disabled cross-check is recorded as such.
* **Reproducible.** `requirements-lock.txt` pins every Python package. [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)
  lists the reference machine, the commands and the numbers they must give.
* **Auditable.** `percoflow verify <run>` recomputes every output hash, and `SHA256SUMS.txt` lists the hash of every file of
  the repository. Every option that can change a number is a field of `config.RunConfig`, and a mistyped option is an error,
  never a silent default. The rules are in [docs/AUDIT_AND_TRACEABILITY.md](docs/AUDIT_AND_TRACEABILITY.md).

## 8. Repository map

```
README.md                 this page
docs/                     theory, inputs/outputs, OpenFOAM, benchmarks, reproducibility, audit rules, glossary
  reference_results/      benchmark report and example results of this release
src/percoflow/            the package (one module per step; see src/percoflow/README.md)
tests/                    pytest suite: unit + verification, OpenFOAM (never skipped), regression
examples/                 synthetic inputs, default configuration, tutorial
scripts/                  OpenFOAM user-space installer, checksum generator
requirements-lock.txt     exact Python package versions of the reference environment
SHA256SUMS.txt            SHA-256 of every tracked file
CITATION.cff, AUTHORS.md, CHANGELOG.md, LICENSE
```

## 9. Scope and limitations

PercoFlow handles two-dimensional, single-phase, incompressible, steady creeping flow of a Newtonian fluid. It does not
handle inertial, transient, multiphase, non-Newtonian or reactive flow.

Three cautions apply to how the results are read:
* A 2D section can miss connections that exist in 3D.
* The permeability belongs to the analysed window and its definition, and both must be reported with it.
* A Stokes solution does not fix a Reynolds number. `descriptors.reynolds_bound` gives the largest pressure gradient for
  which the creeping-flow assumption holds for a given fluid and pixel size.

## 10. Authors and citation

Raphael M. Tromer, Alysson Martins Almeida Silva, Paola Ferreira Barbosa, Roberto Ventura Santos and Luiz Antonio Ribeiro
Junior (corresponding author, luiz.a.ribeiro@ntnu.no). Affiliations are in [AUTHORS.md](AUTHORS.md); the citation metadata
is in [CITATION.cff](CITATION.cff).

## Licence

Copyright © 2026 the PercoFlow authors. PercoFlow is released under the MIT licence; see [LICENSE](LICENSE).
