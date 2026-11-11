#!/usr/bin/env bash
# Poiseuille verification of the author's pipeline (generated during the manuscript revision).
# Straight channels of width h (px) and length L between the inlet/outlet planes, exact k = h^2/12.
# Requires: gmsh, scikit-fem, scikit-image, meshio, numpy, scipy, pillow, matplotlib.
# The only change to the author's code is min_contour_points 8 -> 4 (rectangles have 5 vertices),
# see pipeline_patch.diff.
set -euo pipefail
PIPE=${1:-digital_rock_cfd_pipeline_sparse.py}
python3 - <<'PY'
import numpy as np
from PIL import Image
import os
os.makedirs('structs2', exist_ok=True)
def make(h, Lflow, name, pad=40):
    W = int(round(Lflow / 0.84)); H = h + 2 * pad
    img = np.zeros((H, W), np.uint8); img[pad:pad + h, :] = 255
    Image.fromarray(img, 'L').save(f'structs2/{name}.png')
make(20, 400, 'h020'); make(40, 800, 'h040'); make(80, 800, 'h080'); make(160, 1600, 'h160')
PY
sed 's/min_contour_points: int = 8/min_contour_points: int = 4/' "$PIPE" > pipe_val.py
LV=("01_coarse 12 70 12 140" "02_intermediate 8 50 8 100" "03_fine 4 30 4 60" "04_very_fine 2 15 2 30" "05_extra_fine 1.5 12 1.5 24")
for lvl in "${LV[@]}"; do set -- $lvl
  python3 pipe_val.py --structures structs2 --direction lr --output v_$1 --boundary-band-fraction 0.08 \
    --closing-radius 0 --keep-percolating all --min-object-size 20 --contour-tolerance 1.5 --min-contour-area 20 \
    --mesh-algorithm 6 --smoothing 5 --pixel-size 1.0 --unit-name pixel --viscosity 1.0 --density 1.0 \
    --inlet-velocity 0.0001 --solver-tolerance 1e-7 --h-min $2 --h-max $3 --distance-min $4 --distance-max $5 \
    --solver-max-iterations 20000 --max-total-dofs 6000000
done
