# Run Parameters

This file contains example commands for running the `digital_rock_cfd_sparse.py`
workflow from the root directory of the repository.

Before running the commands, confirm that:

- `digital_rock_cfd_sparse.py` is in the repository root;
- the segmented rock images are stored in `structures/`;
- all Python dependencies are installed.

## Horizontal flow

Runs the simulations from left to right (`lr`).

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
