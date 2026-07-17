# CFD Post-Processing and Visualization

This file contains example commands for post-processing the CFD results generated
by `digital_rock_cfd_sparse.py`.

The visualization script reads the `solution.vtk` and `summary.csv` files created
for each successful simulation and automatically generates velocity, pressure,
streamline, vector-field and hydraulic-zone figures.

## Requirements

Before running the post-processing commands, confirm that:

- `visualize_cfd_outputs.py` is in the repository root;
- the CFD simulations have already been completed;
- each successful case contains `solution.vtk`;
- Python and the required visualization packages are installed.

Install the post-processing dependencies with:

```bash
python -m pip install numpy matplotlib meshio
```

The recommended repository structure is:

```text
digital-rock-cfd-sparse/
├── digital_rock_cfd_sparse.py
├── visualize_cfd_outputs.py
├── structures/
├── results_horizontal/
├── results_vertical/
├── results_both/
└── visualize_CFD_parameter.md
```

---

## 1. Visualize horizontal-flow results

This command processes only the horizontal simulations stored in
`results_horizontal/`.

```bash
python visualize_cfd_outputs.py \
  --results-root results_horizontal \
  --direction lr \
  --output-dir cfd_visualizations_horizontal \
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

The generated images will be stored in:

```text
cfd_visualizations_horizontal/
```

---

## 2. Visualize vertical-flow results

This command processes only the vertical simulations stored in
`results_vertical/`.

```bash
python visualize_cfd_outputs.py \
  --results-root results_vertical \
  --direction tb \
  --output-dir cfd_visualizations_vertical \
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

The generated images will be stored in:

```text
cfd_visualizations_vertical/
```

---

## 3. Visualize both flow directions

This command processes the horizontal and vertical simulations stored in
`results_both/`.

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

The generated images will be stored in:

```text
cfd_visualizations_both/
```

---

## 4. Process horizontal and vertical result folders together

Use repeated `--results-root` arguments to process separate horizontal and
vertical result directories in a single command.

```bash
python visualize_cfd_outputs.py \
  --results-root results_horizontal \
  --results-root results_vertical \
  --direction both \
  --output-dir cfd_visualizations_all \
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

The generated images will be stored in:

```text
cfd_visualizations_all/
```

---

## 5. Automatic result discovery

The script can automatically search the current directory for result folders
whose names begin with `results_complete`.

For example:

```text
results_complete_lr/
results_complete_tb/
results_complete_both/
```

Run automatic discovery with:

```bash
python visualize_cfd_outputs.py
```

To specify the base directory used during automatic discovery:

```bash
python visualize_cfd_outputs.py \
  --base-dir . \
  --output-dir cfd_visualizations
```

> Automatic discovery only searches for directories whose names start with
> `results_complete`. For folders named `results_horizontal`,
> `results_vertical`, or `results_both`, use `--results-root` explicitly.

---

## 6. Generate publication-quality PNG figures

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_publication \
  --format png \
  --dpi 350 \
  --max-points 60000 \
  --max-triangles 120000 \
  --grid-resolution 260 \
  --seed-density 28 \
  --streamline-density 1.15 \
  --max-arrows 1800 \
  --speed-percentile 99.5
```

This configuration increases image resolution and streamline detail while
maintaining controlled plotting memory.

---

## 7. Generate SVG figures

Use SVG when editable vector text, axes and streamline elements are required.

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_svg \
  --format svg \
  --dpi 300 \
  --max-points 50000 \
  --max-triangles 100000 \
  --grid-resolution 240 \
  --seed-density 26 \
  --streamline-density 1.10 \
  --max-arrows 1400 \
  --speed-percentile 99.5
```

The scalar CFD fields are rasterized internally when necessary to avoid
excessively large SVG files.

---

## 8. Generate PDF figures

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_pdf \
  --format pdf \
  --dpi 300 \
  --max-points 50000 \
  --max-triangles 100000 \
  --grid-resolution 240 \
  --seed-density 26 \
  --streamline-density 1.10 \
  --max-arrows 1400 \
  --speed-percentile 99.5
```

---

## 9. Lightweight visualization for very large meshes

Use stronger mesh reduction when the original CFD mesh contains hundreds of
thousands or millions of elements.

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_light \
  --format png \
  --dpi 180 \
  --max-points 30000 \
  --max-triangles 60000 \
  --grid-resolution 180 \
  --seed-density 20 \
  --streamline-density 0.95 \
  --max-arrows 900 \
  --speed-percentile 99.5
```

This option affects only visualization. It does not modify the original
finite-element solution or the calculated permeability.

---

## 10. High-detail streamline visualization

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_streamlines \
  --format png \
  --dpi 350 \
  --max-points 70000 \
  --max-triangles 140000 \
  --grid-resolution 320 \
  --seed-density 36 \
  --streamline-density 1.35 \
  --max-arrows 1600 \
  --speed-percentile 99.5
```

Higher values of `--grid-resolution`, `--seed-density`, and
`--streamline-density` produce more detailed streamline figures but increase
post-processing time and memory use.

---

## 11. Visualize one specific sample

The `--sample` option must receive the complete sample-directory name.

Example:

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --sample patch_y7600_x53200_c0_mask \
  --output-dir cfd_visualizations_R3 \
  --format png \
  --dpi 300
```

To process more than one selected sample, repeat `--sample`:

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --sample patch_y7600_x19000_c0_mask \
  --sample patch_y7600_x53200_c0_mask \
  --output-dir cfd_visualizations_selected \
  --format png \
  --dpi 300
```

---

## 12. Remove triangular mesh lines

Use `--no-mesh-lines` to generate cleaner scalar-field and streamline figures
without displaying the reduced triangular mesh.

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_clean \
  --format png \
  --dpi 300 \
  --no-mesh-lines
```

---

## Generated post-processing files

For every successful CFD case, the script creates:

```text
normalized_speed.png
pressure.png
velocity_vectors.png
streamlines.png
flow_zones.png
overview.png
```

When another output format is selected, the corresponding extension changes to
`.svg` or `.pdf`.

### `normalized_speed`

Normalized velocity-magnitude field:

```text
|u| / U0
```

where `U0` is the prescribed inlet velocity.

### `pressure`

Resolved pressure field obtained from the finite-element solution.

### `velocity_vectors`

Normalized velocity field with spatially distributed vector arrows.

### `streamlines`

Interpolated velocity-magnitude field with streamline trajectories initiated
near the hydraulic inlet.

### `flow_zones`

Classification of the connected pore domain into:

- active-flow regions;
- intermediate-flow regions;
- dead or nearly stagnant regions.

### `overview`

Combined four-panel visualization containing:

- normalized velocity magnitude;
- pressure field;
- velocity vectors;
- streamlines.

---

## Output-directory structure

The visualization script preserves the hierarchy of the CFD result folders.

Example:

```text
cfd_visualizations_both/
└── results_both/
    ├── lr/
    │   └── patch_y7600_x53200_c0_mask/
    │       ├── normalized_speed.png
    │       ├── pressure.png
    │       ├── velocity_vectors.png
    │       ├── streamlines.png
    │       ├── flow_zones.png
    │       └── overview.png
    └── tb/
        ├── patch_y7600_x19000_c0_mask/
        └── patch_y7600_x53200_c0_mask/
```

A global processing report is also created:

```text
visualization_manifest.csv
```

The manifest records:

- result directory;
- flow direction;
- sample name;
- processing status;
- original mesh size;
- reduced visualization-mesh size;
- number of generated images;
- output location;
- error message, if a case fails.

---

## Main post-processing parameters

| Parameter | Description | Default |
|---|---|---:|
| `--base-dir` | Directory searched for `results_complete*` folders | `.` |
| `--results-root` | Explicit CFD result directory | none |
| `--output-dir` | Main visualization-output directory | `cfd_visualizations` |
| `--direction` | Direction filter: `auto`, `lr`, `tb`, or `both` | `auto` |
| `--sample` | Complete sample-directory name used as a filter | all samples |
| `--format` | Output format: `png`, `svg`, or `pdf` | `png` |
| `--dpi` | Figure resolution | `220` |
| `--max-points` | Maximum vertices retained for plotting | `60000` |
| `--max-triangles` | Maximum triangles retained for plotting | `120000` |
| `--grid-resolution` | Regular-grid resolution used for streamlines | `220` |
| `--seed-density` | Base number of streamline seeds | `24` |
| `--streamline-density` | Streamline density passed to Matplotlib | `1.05` |
| `--max-arrows` | Maximum number of velocity vectors | `1600` |
| `--speed-percentile` | Upper percentile used for the velocity color scale | `99.5` |
| `--no-mesh-lines` | Removes reduced-mesh lines from figures | disabled |

---

## Recommended commands

### Fast inspection

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_quick \
  --dpi 160 \
  --max-points 25000 \
  --max-triangles 50000 \
  --grid-resolution 160 \
  --max-arrows 700
```

### Standard analysis

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_standard \
  --dpi 220 \
  --max-points 60000 \
  --max-triangles 120000 \
  --grid-resolution 220 \
  --seed-density 24 \
  --streamline-density 1.05 \
  --max-arrows 1600
```

### Publication-quality analysis

```bash
python visualize_cfd_outputs.py \
  --results-root results_both \
  --direction both \
  --output-dir cfd_visualizations_publication \
  --format png \
  --dpi 350 \
  --max-points 70000 \
  --max-triangles 140000 \
  --grid-resolution 300 \
  --seed-density 32 \
  --streamline-density 1.25 \
  --max-arrows 1800 \
  --speed-percentile 99.5 \
  --no-mesh-lines
```

---

## Notes

- Run the commands from the root directory of the repository.
- The visualization script does not rerun the CFD simulation.
- The original `solution.vtk` files are never modified.
- Mesh reduction is applied only to plotting and does not change numerical
  permeability, pressure-drop or flow-rate results.
- Streamlines are generated from an interpolated regular grid.
- Velocity, pressure and permeability calculations remain based on the original
  finite-element solution.
- The backslash `\` at the end of each command line continues the command on
  Linux and macOS.
- On Windows PowerShell, replace each backslash with a backtick `` ` `` or place
  the complete command on a single line.
