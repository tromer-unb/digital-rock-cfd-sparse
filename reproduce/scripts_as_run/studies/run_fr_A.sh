#!/usr/bin/env bash
# FR-A: reservoir-thickness study (level 2, P2/P1, pressure-driven). Sequential, niced, with timeout.
set -u
QC="$(cd "$(dirname "$0")/.." && pwd)"; PY=${PY:?set PY to the venv python}
OUT="$QC/fix-round/FR-A_reservoir"; mkdir -p "$OUT"
for case in R3-H R2-V R3-V; do for pad in 0 0.025 0.05 0.10 0.20; do
  timeout 1500 nice -n 10 $PY "$QC/scripts/fr_reservoir.py" --case $case --level 2 --pad $pad --elements p2 --out "$OUT" 2>&1 | grep -v "^Info\|^Warning" >> "$OUT/log.txt"
done; done
echo DONE >> "$OUT/log.txt"
