"""OpenFOAM cross-check tests. They are NOT skipped when OpenFOAM is missing: they fail with installation instructions,
because the cross-check is part of the method. Run: pytest -q -m openfoam"""
import json
from pathlib import Path

import numpy as np
import pytest

from percoflow.benchmarks import openfoam_channel
from percoflow.benchmarks.openfoam import find_openfoam
from percoflow.benchmarks.openfoam.pixel import face_connected, fv_pixel
from percoflow.config import OpenFOAMOptions, RunConfig
from percoflow.pipeline import run_image, run_mesh
from percoflow.provenance import verify_manifest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
pytestmark = pytest.mark.openfoam


def test_openfoam_is_installed():
    of = find_openfoam()
    assert of.version


def test_face_connectivity_rule():
    w = np.zeros((4, 4), bool)
    w[0:2, 0:2] = True
    w[2:4, 2:4] = True               # touches the first block only at a corner
    assert not face_connected(w).any()


def test_channel_exact(tmp_path):
    r = openfoam_channel(tmp_path, OpenFOAMOptions())
    assert r["pass"], r["rel_err"]
    assert r["fv"]["mesh_ok"] and r["fv"]["converged"]


def test_image_pipeline_with_crosscheck(tmp_path):
    cfg = RunConfig.from_dict({"mesh": {"level": 3}})
    r = run_image(EXAMPLES / "synthetic_rock.png", tmp_path / "rock", cfg)
    oc = r["openfoam_crosscheck"]
    assert r["status"] == "ok" and oc["status"] == "completed"
    assert oc["summary"]["all_fv_converged"] and oc["summary"]["all_mesh_ok"]
    assert abs(oc["summary"]["median_abs_rel_diff"]) < 0.10          # representation difference, a few per cent here
    assert verify_manifest(tmp_path / "rock")["ok"]
    man = json.loads((tmp_path / "rock" / "manifest.json").read_text())
    assert man["openfoam"]["version"] and man["status"] == "completed"


def test_mesh_pipeline_with_crosscheck(tmp_path):
    r = run_mesh(EXAMPLES / "obstacles.msh", tmp_path / "obst")
    oc = r["openfoam_crosscheck"]
    assert oc["status"] == "completed" and oc["converged"] and oc["mesh_ok"]
    assert abs(oc["rel_diff_fv_vs_fe"]) < 0.05


def test_disabling_is_recorded(tmp_path):
    cfg = RunConfig.from_dict({"mesh": {"level": 1}, "openfoam": {"enabled": False}})
    r = run_image(EXAMPLES / "synthetic_rock.png", tmp_path / "off", cfg)
    man = json.loads((tmp_path / "off" / "manifest.json").read_text())
    assert r["openfoam_crosscheck"]["status"] == "disabled_by_user"
    assert any("DISABLED" in n for n in man["notes"])
