# Digital Rock CFD Sparse

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Gmsh](https://img.shields.io/badge/Mesh-Gmsh-orange.svg)](https://gmsh.info/)
[![Finite Elements](https://img.shields.io/badge/FEM-scikit--fem-green.svg)](https://github.com/kinnala/scikit-fem)
[![Research Software](https://img.shields.io/badge/status-research%20software-lightgrey.svg)](#scientific-scope-and-limitations)

**Digital Rock CFD Sparse** is an automated, percolation-aware, and
memory-efficient workflow that converts segmented two-dimensional pore images
into directional permeability estimates, resolved Stokes-flow fields, and
pore-scale hydraulic diagnostics with minimal user intervention.

The framework integrates:

- segmented-image preprocessing;
- automatic hydraulic-connectivity analysis;
- extraction of inlet-to-outlet percolating pore domains;
- adaptive unstructured triangular meshing with Gmsh;
- mixed finite-element solution of the incompressible Stokes equations;
- sparse MINRES solution of the velocity-pressure system;
- apparent-permeability calculation;
- active, intermediate, and nearly stagnant flow-zone analysis;
- VTK and CSV export;
- lightweight CFD post-processing and visualization.

The software is designed for rapid and reproducible analysis of
**two-dimensional, single-phase, incompressible, low-Reynolds-number flow**
through segmented porous structures.

> [!IMPORTANT]
> This software is not intended to replace general-purpose CFD platforms such
> as OpenFOAM. It provides a focused and rapidly deployable workflow for
> digital-rock screening, directional-permeability estimation, hydraulic
> connectivity analysis, and visualization of preferential flow pathways.

---

## Contents

- [Overview](#overview)
- [Main features](#main-features)
- [Computational workflow](#computational-workflow)
- [Governing model](#governing-model)
- [Repository structure](#repository-structure)
- [Study structures](#study-structures)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Running the three flow configurations](#running-the-three-flow-configurations)
- [CFD post-processing](#cfd-post-processing)
- [Input-image requirements](#input-image-requirements)
- [Main simulation parameters](#main-simulation-parameters)
- [Output files](#output-files)
- [Study demonstration](#study-demonstration)
- [Reproducibility](#reproducibility)
- [Scientific scope and limitations](#scientific-scope-and-limitations)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Overview

Image-based pore-scale simulation can require several independent operations,
including connectivity analysis, geometry reconstruction, mesh generation,
boundary identification, numerical-solver configuration, convergence
verification, and post-processing.

Digital Rock CFD Sparse combines these stages into a single parameterized
Python workflow.

The central idea is to solve the CFD problem only in the pore components that
establish a continuous hydraulic connection between the prescribed inlet and
outlet boundaries. Disconnected pores and non-conducting cavities are excluded
before mesh generation and therefore do not consume finite-element degrees of
freedom.

The retained pore geometry is discretized using an adaptive triangular mesh.
Smaller elements are concentrated near pore walls and narrow throats, whereas
larger elements are used in wider pore bodies. The incompressible Stokes
equations are then solved using a stable mixed finite-element formulation and
a sparse iterative solver.

The workflow produces not only an apparent-permeability value, but also the
pressure field, velocity field, preferential pathways, nearly stagnant
regions, hydraulic diameter, Reynolds number, solver diagnostics, and
standardized visualization files.

---

## Main features

### Automated hydraulic-domain extraction

Only pore components connecting the selected inlet and outlet planes are
retained. Non-percolating samples are detected automatically and reported
without attempting to construct an invalid through-flow simulation.

### Directional analysis

Three execution modes are available:

| Option | Meaning |
|---|---|
| `--direction lr` | Horizontal flow from left to right |
| `--direction tb` | Vertical flow from top to bottom |
| `--direction both` | Runs both horizontal and vertical configurations |

### Adaptive unstructured meshing

The pore-solid interface is reconstructed from the segmented image and
discretized using Gmsh. Mesh resolution is controlled by the distance from the
pore walls, allowing narrow throats to receive greater resolution than wide
pore bodies.

### Sparse mixed finite elements

The velocity-pressure Stokes system is discretized using a mixed
finite-element formulation. The global saddle-point system remains in sparse
format and is solved iteratively using MINRES with block-oriented diagonal
preconditioning.

### Memory-controlled processing

The workflow avoids dense global matrix factorization. Vertex interpolation
is performed in batches, and a configurable upper limit can be applied to the
total number of degrees of freedom.

### Integrated hydraulic characterization

The code automatically evaluates quantities including:

- volumetric flow rate;
- Darcy velocity;
- inlet and outlet mean pressure;
- pressure drop;
- apparent permeability;
- hydraulic diameter;
- Reynolds number;
- connected-pore fraction;
- active-flow area;
- intermediate-flow area;
- dead or nearly stagnant area;
- maximum velocity amplification;
- solver residual;
- number of degrees of freedom;
- number of sparse nonzero coefficients;
- estimated CSR-system storage.

### Automated visualization

A separate lightweight post-processing script reads the generated VTK files
and creates:

- normalized velocity maps;
- pressure maps;
- velocity-vector maps;
- streamline visualizations;
- active, intermediate, and dead-flow-zone maps;
- combined CFD overview figures.

---

## Computational workflow

```text
Segmented pore image
        │
        ▼
Image reading and binary-mask cleaning
        │
        ▼
Definition of fixed inlet and outlet planes
        │
        ▼
Connected-component and percolation analysis
        │
        ▼
Extraction of the complete hydraulic domain
        │
        ▼
Pore-boundary reconstruction
        │
        ▼
Adaptive triangular mesh generation with Gmsh
        │
        ▼
Mixed finite-element Stokes discretization
        │
        ▼
Sparse MINRES solution
        │
        ▼
Pressure drop, flow rate, and permeability
        │
        ▼
Hydraulic metrics, VTK files, CSV reports, and figures
```

The complete procedure is executed independently for every image and requested
flow direction.

---

## Governing model

The current implementation considers steady, incompressible creeping flow:

```math
-\mu \nabla^2 \mathbf{u} + \nabla p = \mathbf{0},
```

```math
\nabla \cdot \mathbf{u} = 0,
```

where:

- \(\mathbf{u}\) is the velocity field;
- \(p\) is pressure;
- \(\mu\) is the dynamic viscosity.

The main boundary conditions are:

- prescribed uniform velocity at the inlet;
- zero pressure at the outlet;
- no-slip condition at pore-solid walls.

Apparent permeability is calculated from Darcy's law:

```math
k_{\mathrm{app}}
=
\frac{\mu U_D L}{\Delta P},
```

where:

- \(U_D\) is the Darcy velocity;
- \(L\) is the hydraulic-domain length;
- \(\Delta P\) is the inlet-to-outlet pressure drop.

The hydraulic diameter is evaluated from the connected fluid area and wetted
perimeter, and the corresponding Reynolds number is reported as an
a posteriori verification of the creeping-flow assumption.

---

## Repository structure

```text
digital-rock-cfd-sparse/
│
├── digital_rock_cfd_sparse.py
├── visualize_cfd_outputs.py
│
├── structures.tar.gz
│   ├── patch_y3800_x3800_c0_mask.png
│   ├── patch_y7600_x19000_c0_mask.png
│   ├── patch_y7600_x53200_c0_mask.png
│
│
├── run_parameter.md
├── visualize_CFD_parameter.md
├── requirements.txt
├── CITATION.cff
├── LICENSE
├── CONTRIBUTING.md
├── CHANGELOG.md
├── .gitignore
└── README.md
```

The `figures/` directory is optional and may be used to distribute the scripts
employed to reproduce the article figures.

---

## Study structures

The repository includes the three segmented structures used in the numerical
study.

| Rock | File | Hydraulic behavior observed with the study parameters |
|---|---|---|
| R1 | `patch_y3800_x3800_c0_mask.png` | Non-percolating |
| R2 | `patch_y7600_x19000_c0_mask.png` | Vertically percolating |
| R3 | `patch_y7600_x53200_c0_mask.png` | Horizontally and vertically percolating |

The connectivity classification depends on the selected inlet and outlet
planes and on the preprocessing parameters. Users should therefore report the
complete configuration used in each analysis.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/digital-rock-cfd-sparse.git
cd digital-rock-cfd-sparse
```

### 2. Create a virtual environment

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

### 3. Install the dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Alternatively, install the principal dependencies directly:

```bash
python -m pip install \
  numpy \
  scipy \
  matplotlib \
  pillow \
  scikit-image \
  gmsh \
  meshio \
  "scikit-fem[all]"
```

The lightweight post-processing script requires:

```bash
python -m pip install numpy matplotlib meshio
```

### 4. Check the command-line interface

```bash
python digital_rock_cfd_sparse.py --help
```

```bash
python visualize_cfd_outputs.py --help
```

---

## Quick start

Place the segmented images in:

```text
structures/
```

Then run both hydraulic directions:

```bash
python digital_rock_cfd_sparse.py \
  --structures structures \
  --output results_both \
  --direction both \
  --boundary-band-fraction 0.08 \
  --closing-radius 0 \
  --h-min 4 \
  --h-max 30 \
  --distance-min 4 \
  --distance-max 60 \
  --viscosity 1.0 \
  --density 1.0 \
  --inlet-velocity 0.001 \
  --solver-tolerance 1e-8
```

The main output directory will contain:

```text
results_both/
├── lr/
├── tb/
├── global_summary.csv
└── run_configuration.json
```

Every successful sample-direction case is stored in its own subdirectory.

---

## Running the three flow configurations

Complete command descriptions are available in
[`run_parameter.md`](run_parameter.md).

### Horizontal flow

```bash
python digital_rock_cfd_sparse.py \
  --structures structures \
  --output results_horizontal \
  --direction lr \
  --boundary-band-fraction 0.08 \
  --closing-radius 0 \
  --h-min 4 \
  --h-max 30 \
  --distance-min 4 \
  --distance-max 60 \
  --viscosity 1.0 \
  --density 1.0 \
  --inlet-velocity 0.001 \
  --solver-tolerance 1e-8
```

### Vertical flow

```bash
python digital_rock_cfd_sparse.py \
  --structures structures \
  --output results_vertical \
  --direction tb \
  --boundary-band-fraction 0.08 \
  --closing-radius 0 \
  --h-min 4 \
  --h-max 30 \
  --distance-min 4 \
  --distance-max 60 \
  --viscosity 1.0 \
  --density 1.0 \
  --inlet-velocity 0.001 \
  --solver-tolerance 1e-8
```

### Horizontal and vertical flow

```bash
python digital_rock_cfd_sparse.py \
  --structures structures \
  --output results_both \
  --direction both \
  --boundary-band-fraction 0.08 \
  --closing-radius 0 \
  --h-min 4 \
  --h-max 30 \
  --distance-min 4 \
  --distance-max 60 \
  --viscosity 1.0 \
  --density 1.0 \
  --inlet-velocity 0.001 \
  --solver-tolerance 1e-8
```

---

## CFD post-processing

Complete post-processing examples are available in
[`visualize_CFD_parameter.md`](visualize_CFD_parameter.md).

### Visualize horizontal results

```bash
python visualize_cfd_outputs.py \
  --results-root results_horizontal \
  --direction lr \
  --output-dir cfd_visualizations_horizontal \
  --format png \
  --dpi 220
```

### Visualize vertical results

```bash
python visualize_cfd_outputs.py \
  --results-root results_vertical \
  --direction tb \
  --output-dir cfd_visualizations_vertical \
  --format png \
  --dpi 220
```

### Visualize both directions

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_both \
  --format png \
  --dpi 220 \
  --max-points 60000 \
  --max-triangles 120000 \
  --grid-resolution 220 \
  --seed-density 24 \
  --streamline-density 1.05 \
  --max-arrows 1600 \
  --speed-percentile 99.5
```

For each successful CFD case, the post-processing script generates:

```text
normalized_speed.png
pressure.png
velocity_vectors.png
streamlines.png
flow_zones.png
overview.png
```

A global post-processing manifest is also created:

```text
visualization_manifest.csv
```

The visualization script may reduce the mesh used for plotting to control
memory consumption and output-file size. This reduction affects only the
figures. The original finite-element solution, pressure drop, flow rate, and
permeability are not modified.

---

## Input-image requirements

Supported image extensions include:

```text
.png
.jpg
.jpeg
.tif
.tiff
.bmp
```

For binary or grayscale images:

- white or high-intensity pixels represent the pore or fluid phase;
- black or low-intensity pixels represent the solid matrix.

For the clearest and most reproducible input, use a binary lossless image such
as PNG or TIFF.

Recommended requirements:

- use one segmented sample per file;
- avoid lossy JPEG compression;
- preserve the original image dimensions;
- avoid scale bars, labels, borders, or annotations inside the mask;
- ensure that the pore phase is represented consistently;
- record the physical pixel size when permeability is required in physical
  units.

When:

```text
--pixel-size 1.0
--unit-name pixel
```

the reported permeability unit is:

```text
pixel^2
```

To report permeability in a physical unit, supply the physical length
represented by one pixel and a corresponding unit name.

---

## Main simulation parameters

| Parameter | Description | Default |
|---|---|---:|
| `--structures` | Directory containing segmented images | `structures` |
| `--output` | Main result directory | `results_complete` |
| `--direction` | `lr`, `tb`, or `both` | `lr` |
| `--boundary-band-px` | Boundary-band width in pixels | `0` |
| `--boundary-band-fraction` | Fraction removed from each hydraulic end | `0.08` |
| `--closing-radius` | Morphological-closing radius | `0` |
| `--closing-iterations` | Number of closing operations | `1` |
| `--keep-percolating` | Retain `all` or only the `largest` percolating component | `all` |
| `--min-percolating-area` | Minimum connected-component area in pixels | `0` |
| `--min-object-size` | Removes smaller pore objects | `20` |
| `--fill-holes-smaller-than` | Fills enclosed holes below the specified size | `0` |
| `--pixel-size` | Physical length represented by one pixel | `1.0` |
| `--unit-name` | Spatial-unit name | `pixel` |
| `--thickness` | Out-of-plane model thickness | `1.0` |
| `--h-min` | Minimum mesh size near walls | `4.0` |
| `--h-max` | Maximum mesh size in wider regions | `30.0` |
| `--distance-min` | Distance below which the finest size is applied | `4.0` |
| `--distance-max` | Distance above which the coarsest size is applied | `60.0` |
| `--contour-tolerance` | Polygon-simplification tolerance | `1.5` |
| `--min-contour-area` | Minimum retained contour area | `20.0` |
| `--mesh-algorithm` | Gmsh two-dimensional meshing algorithm | `6` |
| `--smoothing` | Gmsh mesh-smoothing steps | `5` |
| `--viscosity` | Dynamic viscosity | `1.0` |
| `--density` | Fluid density | `1.0` |
| `--inlet-velocity` | Prescribed uniform inlet velocity | `0.001` |
| `--active-threshold-relative` | Active-zone threshold relative to maximum speed | `0.01` |
| `--dead-threshold-relative` | Dead-zone threshold relative to maximum speed | `0.001` |
| `--solver-tolerance` | Relative MINRES tolerance | `1e-8` |
| `--solver-max-iterations` | Maximum MINRES iterations | `5000` |
| `--interpolation-batch-size` | Number of interpolation points per batch | `2048` |
| `--max-total-dofs` | Safety limit for total degrees of freedom | `2000000` |
| `--dpi` | Resolution of automatically generated previews | `220` |

Display every available parameter with:

```bash
python digital_rock_cfd_sparse.py --help
```

---

## Output files

### Global files

The main output directory contains:

```text
global_summary.csv
run_configuration.json
```

#### `global_summary.csv`

Contains one row for every image-direction combination and records whether the
case:

- completed successfully;
- was classified as non-percolating;
- failed during geometry, meshing, or numerical solution.

For successful cases, it includes permeability, pressure drop, flow metrics,
connectivity information, mesh size, solver convergence, and sparse-memory
diagnostics.

#### `run_configuration.json`

Stores the complete configuration used in the run. This file should be
preserved with published results to support reproducibility.

### Sample-level files

A successful case follows the structure:

```text
RESULTS_ROOT/
└── DIRECTION/
    └── SAMPLE_NAME/
        ├── solution.vtk
        ├── summary.csv
        ├── full_percolating_domain.png
        ├── full_percolating_domain.npy
        └── generated mesh and diagnostic files
```

#### `solution.vtk`

Contains the computational mesh and pore-scale fields required for external
visualization and post-processing, including:

- velocity;
- velocity magnitude;
- pressure.

The VTK file can be opened in software such as ParaView or processed using the
included `visualize_cfd_outputs.py` script.

#### `summary.csv`

Contains the hydraulic, geometric, numerical, and solver metrics for one
sample-direction case.

#### `full_percolating_domain.png`

Provides a directly viewable image of the pore space that connects the selected
inlet and outlet boundaries.

#### `full_percolating_domain.npy`

Stores the same percolating-domain mask as a NumPy Boolean array for scientific
reuse.

---

## Study demonstration

Using the converged meshes adopted in the associated study, the workflow
identified three hydraulically distinct configurations:

| Configuration | Apparent permeability |
|---|---:|
| R3 horizontal | \(110.9\ \mathrm{pixel}^2\) |
| R2 vertical | \(270.3\ \mathrm{pixel}^2\) |
| R3 vertical | \(\approx 3.80\times10^4\ \mathrm{pixel}^2\) |

For R3, the resulting directional permeability ratio was:

```math
A_k
=
\frac{k_{\mathrm{R3-V}}}{k_{\mathrm{R3-H}}}
\approx 342.5.
```

The resolved fields showed that:

- R3 horizontal flow was controlled by a localized hydraulic bottleneck;
- R2 vertical flow followed a tortuous pathway constrained by a central neck;
- R3 vertical flow was distributed through a broad conductive corridor.

These results demonstrate why porosity alone is insufficient to describe
hydraulic behavior. Connectivity, directional organization, bottlenecks, and
the spatial distribution of active flow strongly control apparent
permeability.

> [!NOTE]
> These values correspond to the study configuration and converged meshes.
> Different preprocessing, mesh, physical-property, or boundary parameters may
> produce different numerical values.

---

## Reproducibility

For reproducible scientific use, preserve the following items together:

1. the original segmented images;
2. the exact version or release of the code;
3. `run_configuration.json`;
4. `global_summary.csv`;
5. the sample-level `summary.csv` files;
6. the final `solution.vtk` files;
7. the Python environment or dependency lock file;
8. the mesh-convergence settings.

To record the exact Python environment:

```bash
python -m pip freeze > requirements-lock.txt
```

To reproduce a published analysis, use a tagged software release rather than
the continuously changing development branch.

Recommended release tag:

```text
v1.0.0
```

Recommended GitHub release title:

```text
Digital Rock CFD Sparse — Initial Research Release
```

---

## Scientific scope and limitations

Digital Rock CFD Sparse is intentionally specialized.

The present implementation is intended for:

- two-dimensional segmented porous structures;
- steady flow;
- single-phase flow;
- incompressible Newtonian fluids;
- low-Reynolds-number or creeping-flow conditions;
- directional apparent-permeability estimation;
- preferential-pathway and stagnant-region analysis.

The present implementation does not currently provide:

- three-dimensional tomographic simulation;
- turbulent flow;
- transient flow;
- multiphase displacement;
- non-Newtonian rheology;
- heat transfer;
- reactive transport;
- deformation or fluid-structure interaction;
- moving interfaces;
- distributed-memory parallelism.

General-purpose CFD, finite-volume, finite-difference, and lattice-Boltzmann
platforms remain more appropriate when these capabilities are required.

The advantage of this repository is its focused workflow: it substantially
reduces image-to-simulation setup for a well-defined class of pore-scale flow
problems while retaining resolved pressure and velocity fields.

Results should always be evaluated through:

- mesh-convergence analysis;
- solver-residual inspection;
- Reynolds-number verification;
- sensitivity to segmentation;
- sensitivity to boundary-plane placement;
- comparison with experiments or an independent numerical method when
  available.

This repository is research software and is provided without a guarantee that
it is suitable for safety-critical, regulatory, clinical, or commercial
decision-making.

---

## Citation

A `CITATION.cff` file is included in the repository. GitHub can use this file
to display the recommended software citation.

### Software citation

```text
YOUR_NAME. Digital Rock CFD Sparse: an automated, percolation-aware,
and memory-efficient workflow for two-dimensional pore-scale flow simulation.
Version 1.0.0, YEAR. DOI_OR_REPOSITORY_URL.
```

### Associated article

```text
YOUR_AUTHORS. ARTICLE_TITLE. JOURNAL, YEAR.
https://doi.org/ARTICLE_DOI
```

### BibTeX 

```bibtex
@software{digital_rock_cfd_sparse,
  author  = {Raphael M. Tromer and Luiz Ribeiro},
  title   = {Digital Rock CFD Sparse: An Automated, Percolation-Aware,
             and Memory-Efficient Workflow for Two-Dimensional
             Pore-Scale Flow Simulation},
  year    = {YEAR},
  version = {1.0.0},
  url     = {https://github.com/YOUR_USERNAME/digital-rock-cfd-sparse},
  doi     = {SOFTWARE_DOI}
}
```

Replace the placeholder fields after the software release and article
publication.

---

## Contributing

Contributions that improve robustness, documentation, validation, testing, or
computational performance are welcome.

Suggested contribution areas include:

- automated tests;
- additional image formats;
- improved mesh-quality checks;
- advanced saddle-point preconditioners;
- parallel sparse assembly;
- three-dimensional extensions;
- benchmark datasets;
- experimental validation;
- additional visualization tools;
- improved operating-system support.

Before submitting a contribution:

1. create a dedicated branch;
2. make the proposed changes;
3. test the solver on at least one percolating and one non-percolating case;
4. document any new parameters;
5. submit a pull request with a clear description of the change.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for additional guidance.

---

## License

This project is distributed under the MIT License.

See [`LICENSE`](LICENSE) for the complete license text.

The segmented study structures may be redistributed under the repository
license only when the authors have the legal right to release the underlying
image data. Any third-party dataset must retain its original attribution and
license information.

---

## Acknowledgments

This project was developed using open scientific-computing tools, including:

- NumPy;
- SciPy;
- Matplotlib;
- Pillow;
- scikit-image;
- Gmsh;
- meshio;
- scikit-fem.

The authors acknowledge the developers and maintainers of these projects.

---

## Contact
**Raphal M. Tromer and Luiz Ribeiro**  
  


For reproducibility questions, bug reports, or feature requests, open a GitHub
issue in this repository.
