#!/usr/bin/env bash
# Fix-round F-001: pressure-driven + open reservoirs. A = reservoir width; B = mesh convergence (b=10%); C = geometry sensitivity (b=10%).
# Sequential, niced, per-job timeout (the Gmsh hang of F-015 must not stall the queue).
set -u
QC="$(cd "$(dirname "$0")/.." && pwd)"; PY=${PY:?set PY to the venv python}; FR="$QC/scripts/fr_reservoir.py"
run() { local out=$1; shift; mkdir -p "$out"; timeout 1500 nice -n 10 $PY "$FR" "$@" --out "$out" 2>&1 | grep -v "^Info\|^Warning\|^$" >> "$out/log.txt"; }
for c in R3-H R2-V R3-V; do for pad in 0 0.025 0.05 0.10 0.20 0.40; do run "$QC/fix-round/FR-A_reservoir" --case $c --level 2 --pad $pad --elements p2; done; done
for lv in 1 2 3 4 5; do for c in R3-H R2-V R3-V; do run "$QC/fix-round/FR-B_convergence" --case $c --level $lv --pad 0.10 --elements p2 mini; done; done
for c in R3-H R2-V R3-V; do
  for tol in 0.5 1.0 2.0; do run "$QC/fix-round/FR-C_geometry" --case $c --level 2 --pad 0.10 --tol $tol --elements p2; done
  for band in 0.04 0.12 0.16; do run "$QC/fix-round/FR-C_geometry" --case $c --level 2 --pad 0.10 --band $band --elements p2; done
done
echo DONE > "$QC/fix-round/DONE.txt"
