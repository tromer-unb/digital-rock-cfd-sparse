#!/usr/bin/env bash

set -euo pipefail

SCRIPT="digital_rock_cfd_pipeline_sparse.py"
STRUCTURES="structures"

COMMON_ARGS=(
  --structures "$STRUCTURES"
  --direction both
  --boundary-band-fraction 0.08
  --closing-radius 0
  --keep-percolating all
  --min-object-size 20
  --fill-holes-smaller-than 0
  --contour-tolerance 1.5
  --min-contour-area 20
  --mesh-algorithm 6
  --smoothing 5
  --pixel-size 1.0
  --unit-name pixel
  --thickness 1.0
  --viscosity 1.0
  --density 1.0
  --inlet-velocity 0.0001
  --active-threshold-relative 0.01
  --dead-threshold-relative 0.001
  --solver-tolerance 1e-7
)

python3 "$SCRIPT" \
  "${COMMON_ARGS[@]}" \
  --output convergence_mesh_01_coarse \
  --h-min 12 \
  --h-max 70 \
  --distance-min 12 \
  --distance-max 140 \
  --solver-max-iterations 10000 \
  --interpolation-batch-size 1024 \
  --max-total-dofs 2000000

python3 "$SCRIPT" \
  "${COMMON_ARGS[@]}" \
  --output convergence_mesh_02_intermediate \
  --h-min 8 \
  --h-max 50 \
  --distance-min 8 \
  --distance-max 100 \
  --solver-max-iterations 10000 \
  --interpolation-batch-size 1024 \
  --max-total-dofs 2000000

python3 "$SCRIPT" \
  "${COMMON_ARGS[@]}" \
  --output convergence_mesh_03_fine \
  --h-min 4 \
  --h-max 30 \
  --distance-min 4 \
  --distance-max 60 \
  --solver-max-iterations 10000 \
  --interpolation-batch-size 1024 \
  --max-total-dofs 2500000

python3 "$SCRIPT" \
  "${COMMON_ARGS[@]}" \
  --output convergence_mesh_04_very_fine \
  --h-min 2 \
  --h-max 15 \
  --distance-min 2 \
  --distance-max 30 \
  --solver-max-iterations 15000 \
  --interpolation-batch-size 512 \
  --max-total-dofs 4000000

python3 "$SCRIPT" \
  "${COMMON_ARGS[@]}" \
  --output convergence_mesh_05_extra_fine \
  --h-min 1.5 \
  --h-max 12 \
  --distance-min 1.5 \
  --distance-max 24 \
  --solver-max-iterations 20000 \
  --interpolation-batch-size 256 \
  --max-total-dofs 6000000

echo "Estudo de convergência concluído."
