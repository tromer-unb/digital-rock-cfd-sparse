"""Plane Poiseuille flow through the WHOLE image pipeline: an exactness test.

A straight channel image (all pixels fluid, N_y = h rows) is processed like any rock image: reservoirs, contours, Gmsh mesh,
Taylor-Hood solve, flux-weighted probe pressures. The quadratic velocity and linear pressure of Poiseuille flow lie inside
the Taylor-Hood spaces and satisfy the open-boundary conditions exactly, so the computed k_D must equal the exact value of the
definition on ANY mesh, to round-off:

    walls at y = -0.5 and y = h - 0.5 (channel width h), probe lines L_s + 2 apart, reservoirs are channel too:
    k_D = (h^2 / 12) * L_s / (L_s + 2)

This single number checks assembly, boundary conditions, flux integration and pressure probing together.
Pass criterion (declared in advance): |k_D / exact - 1| < 1e-5 for Taylor-Hood. The Mini element is reported, not tested
(it cannot represent the parabola and lies below the exact value).

The same channel solved by OpenFOAM on the pixel mesh is checked against the exact cell-centred value
k = (h^2 + 2)/12 * L_s / (L_s + 2) (criterion 1e-6) in ``benchmarks/__init__.py``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..config import MeshOptions, SolverOptions
from ..core import solve_window


def exact_fe(h: int, Ls: int) -> float:
    """Exact k_D of the definition for a channel of h px (finite elements reproduce it to round-off)."""
    return h ** 2 / 12.0 * Ls / (Ls + 2.0)


def exact_fv(h: int, Ls: int) -> float:
    """Exact k_D of the cell-centred finite-volume solution on the pixel mesh of the same channel."""
    return (h ** 2 + 2.0) / 12.0 * Ls / (Ls + 2.0)


def run(workdir, widths=(20, 40), levels=(1, 3, 5), aspect: int = 20) -> dict:
    """B1 for each channel width, mesh level and element; pass if every Taylor-Hood error is below 1e-5."""
    workdir = Path(workdir)
    rows = []
    for h in widths:
        Ls = aspect * h
        mask = np.ones((h, Ls), dtype=bool)
        for level in levels:
            for element in ("p2", "mini"):
                r = solve_window(mask, "lr", workdir / f"h{h}_L{level}_{element}", MeshOptions(level=level),
                                 SolverOptions(element=element))
                k = r["darcy"]["k_D"]
                rows.append({"h": h, "L_s": Ls, "level": level, "element": element, "k_D": k, "exact": exact_fe(h, Ls),
                             "rel_err": k / exact_fe(h, Ls) - 1.0, "n_dof": r["solution"].info["n_dof"],
                             "energy_identity_relerr": r["solution"].info["energy_identity_relerr"]})
    th = [abs(r["rel_err"]) for r in rows if r["element"] == "p2"]
    return {"benchmark": "poiseuille_exactness", "rows": rows, "criterion": "Taylor-Hood |rel_err| < 1e-5",
            "max_abs_rel_err_taylor_hood": max(th), "pass": max(th) < 1e-5}
