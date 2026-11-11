#!/usr/bin/env bash
# FR-D: finer contour tolerances and a denser scan of the crop band (window placement), level 2, b=10%, P2/P1.
set -u
QC="$(cd "$(dirname "$0")/.." && pwd)"; PY=${PY:?}; FR="$QC/scripts/fr_reservoir.py"; OUT="$QC/fix-round/FR-C_geometry"
run() { timeout 1500 nice -n 10 $PY "$FR" "$@" --out "$OUT" 2>&1 | grep -v "^Info\|^Warning\|^$" >> "$OUT/log.txt"; }
for c in R3-H R2-V R3-V; do
  for tol in 0.25 0.1; do run --case $c --level 2 --pad 0.10 --tol $tol --elements p2; done
  for band in 0.02 0.04 0.06 0.10 0.14; do run --case $c --level 2 --pad 0.10 --band $band --elements p2; done
done
echo DONE > "$QC/fix-round/DONE_D.txt"
