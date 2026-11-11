#!/usr/bin/env bash
# FR-E: production configuration of the fix round (pressure-driven, reservoirs b=10%, contour tolerance 0.25 px):
#   E1 mesh-convergence series (levels 1-5, P2/P1 and Mini); E2 window (crop-band) scan at level 2.
set -u
QC="$(cd "$(dirname "$0")/.." && pwd)"; PY=${PY:?}; FR="$QC/scripts/fr_reservoir.py"
run() { local out=$1; shift; mkdir -p "$out"; timeout 1800 nice -n 10 $PY "$FR" "$@" --out "$out" 2>&1 | grep -v "^Info\|^Warning\|^$" >> "$out/log.txt"; }
for lv in 1 2 3 4 5; do for c in R3-H R2-V R3-V; do run "$QC/fix-round/FR-E1_convergence_tol0.25" --case $c --level $lv --pad 0.10 --tol 0.25 --elements p2 mini; done; done
for c in R3-H R2-V R3-V; do for band in 0.02 0.04 0.06 0.07 0.09 0.10 0.12 0.14 0.16; do
  run "$QC/fix-round/FR-E2_window_scan_tol0.25" --case $c --level 2 --pad 0.10 --tol 0.25 --band $band --elements p2; done; done
echo DONE > "$QC/fix-round/DONE_E.txt"
