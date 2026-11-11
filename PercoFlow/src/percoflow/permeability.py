"""Darcy permeability of the analysed sample.

Definition (reference protocol; docs/THEORY.md, Section 3):

    k_D = mu Q L_s / (W dP),     dP = Pbar_in - Pbar_out,     Pbar = int_S p u_s ds / int_S u_s ds

Q   integrated flux through the inlet (|flux|),
L_s sample length along the flow (reservoirs excluded),
W   full transverse width of the window (NOT the open length of the inlet: that would give an interstitial quantity),
S   transverse probe line one pixel inside the inlet or outlet reservoir, sampled every half pixel,
u_s velocity component along the flow.

Why flux weighting: with it, Q * dP equals the mechanical power delivered between the two probe lines (up to the small
viscous flux through them), so k_D is the permeability of the homogeneous Darcy medium that dissipates the same power at the
same flux. Measuring dP at the outer reservoir faces instead would add the access resistance of the reservoirs.

For a user mesh (no reservoirs known to the code) the same formula is evaluated with the flux-weighted pressures ON the
INLET and OUTLET curves; L_s and W default to the bounding-box extents along and across the flow.
"""
from __future__ import annotations

import numpy as np
from skfem import Basis, ElementTriP1

from .stokes import StokesSolution


def probe_lines(direction: str, b: int, Ls: int, eps: float = 1.0):
    """Coordinates of the inlet and outlet probe lines on the flow axis for a padded window (reservoirs of width b)."""
    if direction == "lr":
        return {"in": b - 0.5 - eps, "out": b + Ls - 0.5 + eps}
    top = (b + Ls + b) - 1.0       # y of the centre of the first (inlet-side) row of the padded array
    return {"in": top - (b - 0.5) + eps, "out": top - (b + Ls - 0.5) - eps}


def probe_pressures(sol: StokesSolution, direction: str, b: int, Ls: int, W: int, eps: float = 1.0, dt: float = 0.5) -> dict:
    """Length- and flux-weighted pressures on the two probe lines (all probe points lie in the open reservoirs)."""
    fu = sol.basis_u.interpolator(sol.u)
    fp = Basis(sol.mesh, ElementTriP1(), intorder=4).interpolator(sol.p)
    t = np.arange(0.0, W - 1 + 1e-9, dt)
    out = {}
    for key, c in probe_lines(direction, b, Ls, eps).items():
        pts = np.vstack([np.full_like(t, c), t]) if direction == "lr" else np.vstack([t, np.full_like(t, c)])
        uu, pp = fu(pts), fp(pts)
        un = uu[0] if direction == "lr" else -uu[1]
        out[key] = {"p_length": float(np.mean(pp)), "p_flux": float(np.sum(pp * un) / np.sum(un)), "line": float(c)}
    return out


def darcy_from_reservoirs(sol: StokesSolution, direction: str, b: int, Ls: int, W: int, mu: float = 1.0) -> dict:
    """k_D of an image-based run with reservoirs (the default route)."""
    pr = probe_pressures(sol, direction, b, Ls, W)
    Q = sol.info["Q"]
    dP = pr["in"]["p_flux"] - pr["out"]["p_flux"]                  # reported definition
    dP_len = pr["in"]["p_length"] - pr["out"]["p_length"]          # alternative: plain line averages
    dP_outer = sol.info["inlet"]["p_mean_flux"] - sol.info["outlet"]["p_mean_flux"]   # alternative: outer reservoir faces
    return {"k_D": mu * Q * Ls / (W * dP), "Q": Q, "dP": dP, "Pbar_in": pr["in"]["p_flux"], "Pbar_out": pr["out"]["p_flux"],
            "L_s": int(Ls), "W": int(W), "b": int(b), "U_darcy": Q / W,
            "k_D_length_weighted": mu * Q * Ls / (W * dP_len),
            "k_D_outer_faces": mu * Q * Ls / (W * dP_outer),
            "reservoir_pressure_share": float(1.0 - dP / dP_outer) if dP_outer else float("nan"),
            "definition": "k_D = mu Q L_s / (W dP), dP between flux-weighted pressures 1 px inside the reservoirs"}


def darcy_from_boundaries(sol: StokesSolution, direction: str, L=None, W=None, mu: float = 1.0) -> dict:
    """k_D of a user mesh: flux-weighted pressures on INLET and OUTLET; L and W default to the bounding-box extents."""
    ax = 0 if direction == "lr" else 1
    L = float(L) if L is not None else float(np.ptp(sol.mesh.p[ax]))
    W = float(W) if W is not None else float(np.ptp(sol.mesh.p[1 - ax]))
    Q = sol.info["Q"]
    dP = sol.info["inlet"]["p_mean_flux"] - sol.info["outlet"]["p_mean_flux"]
    return {"k_D": mu * Q * L / (W * dP), "Q": Q, "dP": dP, "L_s": L, "W": W, "U_darcy": Q / W,
            "definition": "k_D = mu Q L / (W dP), dP between flux-weighted pressures on the INLET and OUTLET curves"}
