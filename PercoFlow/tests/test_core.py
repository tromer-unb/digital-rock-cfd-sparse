"""Fast unit and verification tests (no OpenFOAM). Run: pytest -q -m "not openfoam" """
import numpy as np
import pytest

from percoflow import imaging, percolation, stokes
from percoflow.benchmarks import cylinder_array, poiseuille
from percoflow.config import DomainOptions, MeshOptions, RunConfig, SolverOptions
from percoflow.core import solve_window
from percoflow.geometry import contour_hierarchy, extract_contours, pad_reservoirs


def test_config_roundtrip_and_typos():
    cfg = RunConfig.from_dict({"mesh": {"level": 3}, "solver": {"element": "mini"}})
    assert cfg.mesh.sizes() == (4.0, 30.0, 4.0, 60.0) and cfg.solver.element == "mini"
    assert RunConfig.from_dict(cfg.to_dict()).to_dict() == cfg.to_dict()
    with pytest.raises(KeyError):
        RunConfig.from_dict({"mesh": {"levle": 3}})
    with pytest.raises(ValueError):
        RunConfig.from_dict({"domain": {"reservoir_fraction": 0.0}})      # pressure-driven needs reservoirs


def test_band_and_crop():
    m = np.ones((100, 50), bool)
    w, info = imaging.crop_window(m, DomainOptions(direction="lr", band_fraction=0.1))
    assert info["band_px"] == 5 and w.shape == (100, 40)
    w, info = imaging.crop_window(m, DomainOptions(direction="tb", band_fraction=0.1))
    assert info["band_px"] == 10 and w.shape == (80, 50)


def test_percolation_classification():
    w = np.zeros((20, 30), bool)
    w[5:8, :] = True                 # percolating channel
    w[12:14, 0:10] = True            # dead end (touches the inlet only)
    w[16:19, 12:18] = True           # isolated cluster
    sel, rep = percolation.analyse(w, "lr")
    assert (rep["clusters_percolating"], rep["clusters_dead_end"], rep["clusters_isolated"]) == (1, 1, 1)
    assert sel.sum() == 3 * 30 and rep["percolates"]
    assert not percolation.analyse(w[:, :15][:, ::-1][10:, :], "lr")[1]["percolates"]


def test_critical_throat_of_a_channel():
    w = np.ones((11, 40), bool)
    wc = percolation.critical_throat_width(w, "lr")
    assert 10.0 <= wc <= 12.0            # channel of 11 px between lateral walls


def test_contours_square_with_hole():
    m = np.ones((40, 40), bool)
    m[15:25, 15:25] = False
    cs = extract_contours(m, tolerance=0.25, min_points=4, min_area=1.0)
    depth, children = contour_hierarchy(cs)
    assert sorted(depth) == [0, 1] and sum(len(c) for c in children) == 1


def test_reservoirs():
    padded, Ls, b = pad_reservoirs(np.zeros((10, 50), bool), "lr", 0.1)
    assert (Ls, b) == (50, 5) and padded.shape == (10, 60) and padded[:, :5].all() and padded[:, -5:].all()


def test_poiseuille_exact_through_whole_pipeline(tmp_path):
    r = poiseuille.run(tmp_path, widths=(20,), levels=(1, 3), aspect=20)
    assert r["pass"], r["max_abs_rel_err_taylor_hood"]
    th = [x for x in r["rows"] if x["element"] == "p2"]
    assert all(abs(x["energy_identity_relerr"]) < 1e-10 for x in th)       # discrete energy identity Phi = p_in Q
    assert all(x["rel_err"] < 0 for x in r["rows"] if x["element"] == "mini")


def test_tb_equals_rotated_lr(tmp_path):
    """Flow direction is a convention: a top-to-bottom run equals the left-to-right run of the transposed window."""
    from percoflow.benchmarks import synthetic_porous
    w = synthetic_porous(n=80, seed=7)
    k_lr = solve_window(w, "lr", tmp_path / "lr", MeshOptions(level=2), SolverOptions())["darcy"]["k_D"]
    k_tb = solve_window(np.ascontiguousarray(w.T), "tb", tmp_path / "tb", MeshOptions(level=2), SolverOptions())["darcy"]["k_D"]
    assert abs(k_tb / k_lr - 1) < 2e-3            # different but equivalent meshes


def test_minres_true_residual_gate(tmp_path):
    w = np.ones((10, 60), bool)
    ok = solve_window(w, "lr", tmp_path / "a", MeshOptions(level=1), SolverOptions(element="mini", method="minres"))
    ref = solve_window(w, "lr", tmp_path / "b", MeshOptions(level=1), SolverOptions(element="mini", method="lu"))
    assert abs(ok["darcy"]["k_D"] / ref["darcy"]["k_D"] - 1) < 1e-3
    assert ok["solution"].info["residual"] > 0 and "minres_iterations" in ok["solution"].info
    with pytest.raises(stokes.SolverNotConverged):
        solve_window(w, "lr", tmp_path / "c", MeshOptions(level=1),
                     SolverOptions(element="mini", method="minres", minres_tol=1e-2, accept_residual=1e-12))


def test_cylinder_array_vs_hasimoto(tmp_path):
    r = cylinder_array.run(tmp_path, cs=(0.1,), h=0.02)
    assert abs(r["rows"][0]["rel_diff"]) < 5e-3
