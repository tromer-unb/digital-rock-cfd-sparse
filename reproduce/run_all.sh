#!/usr/bin/env bash
# Run the whole reproduction of the paper in order and stop at the first failure.
# Usage (from the repository root, with the environment of tools/setup_environment.sh active):
#   bash reproduce/run_all.sh            # full: mesh levels 1-5 and the window scan (about 1 h, 15 GB RAM)
#   bash reproduce/run_all.sh --quick    # mesh levels 1-2 only (about 25 min)
# Every output goes to runs/ (not tracked by git). OpenFOAM v1912 is needed for steps 2 and 3
# (set PERCOFLOW_OPENFOAM_ENV if it is not on the PATH, see PercoFlow/docs/OPENFOAM.md).
set -euo pipefail
cd "$(dirname "$0")/.."
QUICK="${1:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
step() { echo; echo "=== $1"; }

step "0. integrity of the repository"
sha256sum -c SHA256SUMS.txt --quiet && echo "SHA256SUMS: OK"

step "1. PercoFlow tests (unit, verification, OpenFOAM, regression)"
(cd PercoFlow && python -m pytest -q)

step "2. verification benchmarks B1-B5"
percoflow benchmark -o runs/benchmark | tail -3

step "3. permeabilities of the paper (Table 1, Figure 4$([ "$QUICK" = --quick ] || echo ', Figure 5a'))"
if [ "$QUICK" = --quick ]; then
  python reproduce/01_recompute_permeability.py --quick
else
  python reproduce/01_recompute_permeability.py --windows
fi

step "4. Figures 4, 8-11 from the stored results"
python reproduce/02_figures_from_results.py runs/02_figures

step "5. Figures 1-3, 5-7, Supplementary Figures 1-2 and Supplementary Tables, scripts as run"
python reproduce/03_figures_and_tables_as_run.py

echo; echo "ALL STEPS PASS"
