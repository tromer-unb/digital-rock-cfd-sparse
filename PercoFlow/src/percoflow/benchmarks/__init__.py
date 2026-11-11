"""Verification and validation suite. ``run_all`` executes every benchmark, OpenFOAM included, and writes one JSON report.

| id | what it tests | reference | criterion (declared in advance) |
|----|---------------|-----------|---------------------------------|
| B1 poiseuille         | whole image pipeline, FE | exact Poiseuille value of the k_D definition | Taylor-Hood error < 1e-5 |
| B2 cylinder_array     | FE solver, analytic geometry | Hasimoto / Sangani-Acrivos series | <= 0.5 % for c <= 0.15 |
| B3 openfoam_channel   | OpenFOAM pixel route + probe lines | exact cell-centred channel value | error < 1e-6 |
| B4 openfoam_cylinder  | independent code (FV) vs FE | FE value of B2 | FV converged < 1 %, FV vs FE < 2 % |
| B5 openfoam_porous    | FV pixel vs FE contour on a synthetic porous window | none (representation difference) | both converge; difference reported |

OpenFOAM benchmarks are always part of the suite. If OpenFOAM is missing the suite fails with installation instructions;
``skip_openfoam=True`` (CLI ``--no-openfoam``) runs the rest and marks the report "INCOMPLETE".
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..config import MeshOptions, OpenFOAMOptions, SolverOptions


def synthetic_porous(n: int = 160, porosity: float = 0.62, seed: int = 7) -> np.ndarray:
    """Reproducible synthetic porous window: overlapping discs of solid on a fluid background, kept percolating."""
    from ..percolation import analyse
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    solid = np.zeros((n, n), bool)
    while 1 - solid.mean() > porosity:
        cy, cx, r = rng.uniform(0, n), rng.uniform(0, n), rng.uniform(4, 11)
        solid |= (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r
    sel, rep = analyse(~solid, "lr")
    if not rep["percolates"]:
        raise RuntimeError("synthetic window does not percolate; change the seed")
    return sel


def openfoam_channel(workdir, opts: OpenFOAMOptions, h: int = 16, Ls: int = 80) -> dict:
    """B3: open channel of h x L_s pixels on the OpenFOAM pixel route versus its exact cell-centred value."""
    from .openfoam.environment import find_openfoam
    from .openfoam.pixel import fv_pixel
    from .poiseuille import exact_fv
    of = find_openfoam(opts.env_script)
    # u_y = 0 exactly here, so simpleFoam's own stopping test cannot trigger (see case.residual_history); 5000 SIMPLE
    # iterations are far more than needed for this 1536-cell case.
    from dataclasses import replace
    r = fv_pixel(np.ones((h, Ls), bool), Path(workdir) / "openfoam_channel", of, replace(opts, max_iterations=5000))
    ex = exact_fv(h, Ls)
    return {"benchmark": "openfoam_channel_exactness", "openfoam": of.describe(), "k_fv": r["k_D"], "exact": ex,
            "rel_err": r["k_D"] / ex - 1.0, "fv": r, "criterion": "|rel_err| < 1e-6", "pass": abs(r["k_D"] / ex - 1.0) < 1e-6}


def openfoam_porous(workdir, opts: OpenFOAMOptions) -> dict:
    """B5: synthetic porous window, OpenFOAM pixel route versus finite elements on contours (difference reported)."""
    from .openfoam.pixel import crosscheck_image
    window = synthetic_porous()
    rep = crosscheck_image(window, "lr", Path(workdir) / "openfoam_porous", MeshOptions(level=5), SolverOptions(), opts)
    s = rep["summary"]
    return {"benchmark": "openfoam_porous_window", "report": rep, "criterion": "FV converged, mesh OK; difference reported",
            "pass": bool(s["all_fv_converged"] and s["all_mesh_ok"])}


def run_all(workdir, skip_openfoam: bool = False, of_opts: OpenFOAMOptions | None = None, quick: bool = True) -> dict:
    """Run B1-B5 and write ``workdir/benchmark_report.json``. Status PASS only if every criterion passes."""
    from . import cylinder_array, poiseuille
    from .openfoam import cylinder_cell
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    of_opts = of_opts or OpenFOAMOptions()
    res = {"B1_poiseuille": poiseuille.run(workdir / "B1"),
           "B2_cylinder_array": cylinder_array.run(workdir / "B2", h=0.01, cs=(0.05, 0.1, 0.15, 0.2) if quick else (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5))}
    if skip_openfoam:
        res["B3_openfoam_channel"] = res["B4_openfoam_cylinder"] = res["B5_openfoam_porous"] = {"pass": None, "status": "SKIPPED by user"}
    else:
        res["B3_openfoam_channel"] = openfoam_channel(workdir / "B3", of_opts)
        res["B4_openfoam_cylinder"] = cylinder_cell.run(workdir / "B4", opts=of_opts,
                                                        hs=(0.04, 0.02, 0.01) if quick else (0.04, 0.02, 0.01, 0.005),
                                                        cs=(0.1, 0.2) if quick else (0.05, 0.1, 0.2, 0.4))
        res["B5_openfoam_porous"] = openfoam_porous(workdir / "B5", of_opts)
    passes = [v.get("pass") for v in res.values()]
    status = "INCOMPLETE (OpenFOAM skipped)" if skip_openfoam else ("PASS" if all(passes) else "FAIL")
    out = {"status": status, "quick": quick, "results": res}
    (workdir / "benchmark_report.json").write_text(json.dumps(out, indent=1, default=float))
    return out
