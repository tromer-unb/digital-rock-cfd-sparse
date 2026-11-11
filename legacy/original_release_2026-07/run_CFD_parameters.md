# Run Parameters

This file contains example commands for running `digital_rock_cfd_sparse.py`
from the root directory of the repository.

Before running the commands, confirm that:

- `digital_rock_cfd_sparse.py` is in the repository root;
- the segmented rock images are stored in the `structures/` directory;
- all required Python dependencies are installed.

---

## 1. Horizontal flow

Runs the simulations from left to right using `--direction lr`.

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

The results will be stored in:

```text
results_horizontal/
```

---

## 2. Vertical flow

Runs the simulations from top to bottom using `--direction tb`.

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

The results will be stored in:

```text
results_vertical/
```

---

## 3. Horizontal and vertical flow

Runs both flow directions sequentially using `--direction both`.

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

The results will be stored in:

```text
results_both/
```

---

## Main parameters

| Parameter | Description | Value |
|---|---|---:|
| `--structures` | Directory containing the segmented pore images | `structures` |
| `--output` | Main output directory | case dependent |
| `--direction lr` | Horizontal flow from left to right | `lr` |
| `--direction tb` | Vertical flow from top to bottom | `tb` |
| `--direction both` | Runs horizontal and vertical flow | `both` |
| `--boundary-band-fraction` | Fraction removed from each hydraulic boundary | `0.08` |
| `--closing-radius` | Morphological closing radius; zero disables it | `0` |
| `--h-min` | Minimum mesh-element size | `4` |
| `--h-max` | Maximum mesh-element size | `30` |
| `--distance-min` | Distance below which `h-min` is applied | `4` |
| `--distance-max` | Distance above which `h-max` is applied | `60` |
| `--viscosity` | Dynamic viscosity | `1.0` |
| `--density` | Fluid density | `1.0` |
| `--inlet-velocity` | Prescribed uniform inlet velocity | `0.001` |
| `--solver-tolerance` | Relative MINRES convergence tolerance | `1e-8` |

## Notes

- White or `True` pixels represent the pore space.
- Black or `False` pixels represent the solid matrix.
- With `pixel_size = 1`, permeability is reported in `pixel²`.
- The model assumes steady, incompressible Stokes flow at low Reynolds number.
- The backslash `\` at the end of each line continues the command on Linux and macOS.
- Run these commands from the root directory of the repository.
