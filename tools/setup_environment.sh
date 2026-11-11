#!/usr/bin/env bash
# Create the Python environment of this repository in .venv:
#   1. the exact package versions used for the paper (PercoFlow/requirements-lock.txt, every package pinned);
#   2. PercoFlow itself, installed in editable mode from PercoFlow/ (no dependency resolution, so nothing else changes).
# Usage:  bash tools/setup_environment.sh [python-executable]      (default: python3; the reference is Python 3.12.3)
# Afterwards:  source .venv/bin/activate && percoflow --version      -> percoflow 1.0.0
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${1:-python3}"
"$PY" -m venv .venv
.venv/bin/python -m pip install --quiet -r PercoFlow/requirements-lock.txt
.venv/bin/python -m pip install --quiet --no-deps -e PercoFlow
.venv/bin/percoflow --version
echo "Environment ready: source .venv/bin/activate"
echo "OpenFOAM (optional, for the cross-check and benchmarks): .venv/bin/percoflow openfoam-check"
