"""Regression tests: the reference numbers of the two example inputs.

These values were computed with PercoFlow 1.0.0 in the environment of requirements-lock.txt (Python 3.12, Linux x86_64) and
are also stored, with their manifests, in docs/reference_results/. Any change of the code that changes a computed number
makes these tests fail. If such a change is intended, update the numbers here and explain it in CHANGELOG.md.

Tolerances: integer quantities (triangle counts) must match exactly; permeabilities to a relative 1e-8, which leaves room
only for round-off of the sparse LU factorization on another CPU. With library versions other than the locked ones the mesh
generator may produce a slightly different mesh: these tests then fail by design, and the difference should be small
(of the order of the discretization error, < 1e-3 relative at mesh level 3).

Run: pytest -q -m regression   (about 20 s, no OpenFOAM needed)
"""
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from percoflow.config import RunConfig
from percoflow.pipeline import run_image, run_mesh

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
pytestmark = pytest.mark.regression

# reference values (PercoFlow 1.0.0, requirements-lock.txt)
ROCK_LR_LEVEL3 = {"k_D": 0.3058875747211206, "Q": 0.3026592104006075, "n_triangles": 5880,
                  "connected_fraction": 0.9639289678135405, "clusters_isolated": 2, "critical_throat_width_px": 2.81640625}
OBSTACLES = {"k_D": 0.011171228780279464, "Q": 0.002792641501519497, "n_triangles": 5409}


def _no_openfoam(level=None):
    d = {"openfoam": {"enabled": False}, "output": {"save_vtk": False, "figures": False}}
    if level:
        d["mesh"] = {"level": level}
    return RunConfig.from_dict(d)


def test_example_image_is_the_published_one(tmp_path):
    """examples/synthetic_rock.png is regenerated identically by examples/make_examples.py (pixel by pixel)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("make_examples", EXAMPLES / "make_examples.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.synthetic_rock(tmp_path / "rock.png")
    a = np.asarray(Image.open(tmp_path / "rock.png"))
    b = np.asarray(Image.open(EXAMPLES / "synthetic_rock.png"))
    assert a.shape == b.shape and np.array_equal(a, b)


def test_rock_image_reference_values(tmp_path):
    r = run_image(EXAMPLES / "synthetic_rock.png", tmp_path / "rock", _no_openfoam(level=3))
    ref = ROCK_LR_LEVEL3
    assert r["status"] == "ok"
    assert r["solver"]["n_triangles"] == ref["n_triangles"]
    assert r["percolation"]["clusters_isolated"] == ref["clusters_isolated"]
    assert r["percolation"]["connected_fraction"] == pytest.approx(ref["connected_fraction"], rel=1e-12)
    assert r["critical_throat_width_px"] == pytest.approx(ref["critical_throat_width_px"], rel=1e-12)
    assert r["permeability"]["k_D"] == pytest.approx(ref["k_D"], rel=1e-8)
    assert r["permeability"]["Q"] == pytest.approx(ref["Q"], rel=1e-8)


def test_obstacle_mesh_reference_values(tmp_path):
    r = run_mesh(EXAMPLES / "obstacles.msh", tmp_path / "obst", _no_openfoam())
    assert r["solver"]["n_triangles"] == OBSTACLES["n_triangles"]
    assert r["permeability"]["k_D"] == pytest.approx(OBSTACLES["k_D"], rel=1e-8)
    assert r["permeability"]["Q"] == pytest.approx(OBSTACLES["Q"], rel=1e-8)
