"""OpenFOAM cross-check of an image domain: one hexahedral finite-volume cell per pore pixel (the staircase of the image),
versus the contour-based finite-element solution of the SAME boolean mask.

Protocol:
1. The domain is restricted to the clusters connected through pixel FACES (4-neighbour) that join both planes: two pore
   pixels touching only at a corner have no contact area in a finite-volume mesh. The identical restricted mask is given to
   the finite-element side, so that both codes solve the same domain.
2. Both sides get the same reservoirs (b = 0.1 L_s), the same imposed pressure and the same k_D definition (flux-weighted
   pressures on lines one pixel inside the reservoirs).
3. The two numbers are NOT expected to coincide: the FE mesh follows a smooth contour through the pixel edges, the FV mesh
   keeps the staircase. The difference is largest where throats are only a few pixels wide. Both solvers are validated on
   the cylinder array (benchmarks B2 and B4); their difference measures the representation of the walls, which a report of
   image-based permeability should state. Differences above ``INVESTIGATE`` (35 %) are flagged.

Exactness check of the FV discretization (used by the benchmark suite): an open channel of N pixels has the cell-centred
solution k = (N^2 + 2)/12 for the face-to-face pressure drop; with probe lines L_s + 2 apart the definition gives
k = (N^2 + 2)/12 * L_s / (L_s + 2).
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from ...config import MeshOptions, OpenFOAMOptions, SolverOptions
from ...geometry import pad_reservoirs
from . import case as foam_case
from .environment import find_openfoam
from .polymesh import pixel_cells, write_polymesh

INVESTIGATE = 0.35
_CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8)


def to_lr(window: np.ndarray, direction: str) -> np.ndarray:
    """Rotate a 'tb' window so that the flow runs left to right (inlet row -> first column). k is unchanged."""
    return np.ascontiguousarray(window.T) if direction == "tb" else window


def face_connected(window_lr: np.ndarray) -> np.ndarray:
    """Pore pixels of the 4-connected clusters that join the first and the last column."""
    labels, _ = ndi.label(window_lr, structure=_CROSS)
    keep = (set(np.unique(labels[:, 0])) & set(np.unique(labels[:, -1]))) - {0}
    return np.isin(labels, sorted(keep)) if keep else np.zeros_like(window_lr, dtype=bool)


def fv_pixel(window_lr: np.ndarray, case_dir, of, opts: OpenFOAMOptions, reservoir_fraction: float = 0.10, p_in: float = 1.0) -> dict:
    """Finite-volume k_D of a left-to-right window (True = fluid, already face-connected)."""
    padded, Ls, b = pad_reservoirs(window_lr, "lr", reservoir_fraction)
    ny, nx = padded.shape
    pts, cells, _ = pixel_cells(padded)

    def edge_patch(mid):
        # boundary edges on x = -0.5 are the inlet, on x = nx - 0.5 the outlet, all others walls
        names = np.full(len(mid), "walls", dtype=object)
        names[np.abs(mid[:, 0] - (-0.5)) < 1e-9] = "inlet"
        names[np.abs(mid[:, 0] - (nx - 0.5)) < 1e-9] = "outlet"
        return names
    case_dir = Path(case_dir)
    if case_dir.exists():
        shutil.rmtree(case_dir)
    pm = write_polymesh(case_dir, pts, cells, edge_patch, {"inlet": "patch", "outlet": "patch", "walls": "wall"}, dz=1.0)
    foam_case.write_case(case_dir, nu=opts.nu, p_in=p_in, max_iterations=opts.max_iterations, nprocs=opts.nprocs)
    diag = foam_case.run(of, case_dir, nprocs=opts.nprocs)
    # probe lines one pixel inside each reservoir: the same definition of k_D as on the finite-element side
    lines = foam_case.line_pressures(case_dir, pm["_internal"], axis=0, coords={"in": b - 1.5, "out": b + Ls + 0.5})
    dP = lines["in"]["p_flux"] - lines["out"]["p_flux"]
    Q = 0.5 * (lines["in"]["Q"] + lines["out"]["Q"])
    W = window_lr.shape[0]
    return {"k_D": opts.nu * Q * Ls / (W * dP), "Q": Q, "dP": dP, "L_s": int(Ls), "W": int(W), "b": int(b),
            "n_cells": pm["n_cells"], "flux_mismatch_in_out": abs(lines["in"]["Q"] - lines["out"]["Q"]) / Q,
            "Re_window_scale": float(Q / opts.nu), **diag}


def _subwindows(window_lr: np.ndarray, size: int, max_n: int):
    """Square tiles of ``size`` px from the window, closest to the centre first, keeping those that percolate."""
    ny, nx = window_lr.shape
    if ny < size or nx < size:
        return [((0, 0), window_lr)]
    tiles = [(r, c) for r in range(0, ny - size + 1, size) for c in range(0, nx - size + 1, size)]
    centre = np.array([ny / 2, nx / 2])
    tiles.sort(key=lambda rc: np.hypot(rc[0] + size / 2 - centre[0], rc[1] + size / 2 - centre[1]))
    out = []
    for r, c in tiles:
        sub = face_connected(window_lr[r: r + size, c: c + size])
        if sub.any():
            out.append(((r, c), sub))
        if len(out) >= max_n:
            break
    return out


def crosscheck_image(window: np.ndarray, direction: str, workdir, mesh_opts: MeshOptions, solver_opts: SolverOptions,
                     of_opts: OpenFOAMOptions, reservoir_fraction: float = 0.10) -> dict:
    """Run FE and FV on the face-connected domain (or on sub-windows when it is too large). Returns a JSON-able report."""
    from ...core import solve_window          # local import: core imports the FE stack only
    of = find_openfoam(of_opts.env_script)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    w_lr = face_connected(to_lr(window, direction))
    if not w_lr.any():
        return {"status": "not_run", "reason": "no face-connected percolating path (only corner contacts)"}
    targets = [((0, 0), w_lr)] if w_lr.sum() <= of_opts.max_cells else _subwindows(w_lr, of_opts.subwindow, of_opts.max_subwindows)
    fe_solver = SolverOptions(**{**solver_opts.__dict__, "element": "p2", "method": "lu", "inlet": "pressure"})
    pairs = []
    for k, ((r, c), sub) in enumerate(targets):
        tag = f"window_r{r}_c{c}"
        fv = fv_pixel(sub, workdir / f"{tag}_openfoam", of, of_opts, reservoir_fraction, fe_solver.p_in)
        fe = solve_window(sub, "lr", workdir / f"{tag}_fe", mesh_opts, fe_solver, reservoir_fraction)
        k_fe = fe["darcy"]["k_D"]
        rel = fv["k_D"] / k_fe - 1.0
        pairs.append({"tag": tag, "origin_rc_in_lr_frame": [int(r), int(c)], "shape": list(sub.shape), "pore_pixels": int(sub.sum()),
                      "k_fe_contour": k_fe, "k_fv_pixel": fv["k_D"], "rel_diff_fv_vs_fe": rel, "fv": fv,
                      "fe_mesh": {k2: v for k2, v in fe["mesh"].items() if k2 != "mesh_path"},
                      "flag": "investigate" if abs(rel) > INVESTIGATE else "ok"})
        if not of_opts.keep_cases:
            shutil.rmtree(workdir / f"{tag}_openfoam", ignore_errors=True)
    rels = np.array([p["rel_diff_fv_vs_fe"] for p in pairs])
    report = {"status": "completed", "openfoam": of.describe(), "scope": "full domain" if len(targets) == 1 and targets[0][0] == (0, 0)
              and w_lr.sum() <= of_opts.max_cells else f"{len(targets)} sub-window(s) of {of_opts.subwindow} px",
              "face_connected_pixels": int(w_lr.sum()), "pairs": pairs,
              "summary": {"n": len(pairs), "median_abs_rel_diff": float(np.median(np.abs(rels))), "mean_rel_diff": float(rels.mean()),
                          "min_rel_diff": float(rels.min()), "max_rel_diff": float(rels.max()),
                          "all_fv_converged": all(p["fv"]["converged"] for p in pairs),
                          "all_mesh_ok": all(p["fv"]["mesh_ok"] for p in pairs)},
              "interpretation": "FV keeps the pixel staircase, FE a smooth contour; the difference measures the wall representation"}
    (workdir / "openfoam_crosscheck.json").write_text(json.dumps(report, indent=1, default=float))
    return report
