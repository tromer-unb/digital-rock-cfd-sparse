"""Shared test settings. The suite has two parts:
    tests/test_core.py      finite-element workflow and verification (no OpenFOAM), about 1 min on one core
    tests/test_openfoam.py  OpenFOAM cross-check (never skipped: fails with installation instructions if OpenFOAM is missing)
    tests/test_regression.py  reference numbers of the example inputs (detects any change of the computed results)
"""
import os

# keep the test suite polite on shared machines
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
