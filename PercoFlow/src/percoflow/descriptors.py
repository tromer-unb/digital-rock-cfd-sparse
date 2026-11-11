"""Pore-scale descriptors derived from a converged solution (docs/THEORY.md, Section 7).

All descriptors are evaluated on the triangles of the SAMPLE (centroid between the two reservoir planes); the reservoirs are
excluded. Speeds are normalized by the mean Darcy velocity U = Q / W.

* Hydraulic partition. Each triangle takes the mean speed of its three vertices, u_T. With u_max the largest vertex speed in
  the sample: active if u_T >= alpha u_max, dead if u_T <= beta u_max, intermediate otherwise (alpha = 0.01, beta = 0.001).
  The thresholds are relative to the maximum of each solution, so the shares compare spatial organization, not mobility.
* Transport localization. Rank triangles by decreasing u_T; a_n = cumulative area fraction, tau_n = cumulative u_T |T|
  fraction. A_50 and A_90 are the smallest a_n with tau_n >= 0.5 and 0.9 (uniform transport: A_50 = 0.5).
* Hydraulic diameter D_h = 4 A_f / P_w (P_w = wall length inside the sample, lateral window edges included).
* Domain porosity phi = A_f / (L_s * transverse extent of the domain) (bounding rectangle of the fluid domain).
* Velocity amplification u_max / U.
* Reynolds number. A Stokes solution is linear in the driving pressure, so it does NOT fix Re. With a mean gradient G,
  U = k_D G / mu and Re = rho U D_h / mu = rho k_D D_h G / mu^2. ``reynolds_bound`` returns the gradient at which Re = 1
  for a given fluid and pixel size; creeping flow requires G well below it.
"""
from __future__ import annotations

import numpy as np

from .stokes import StokesSolution


def triangle_areas(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Areas of the triangles ``tris`` (3 x n vertex indices) of the vertex array ``points`` (2 x N)."""
    a, b, c = points[:, tris[0]], points[:, tris[1]], points[:, tris[2]]
    return 0.5 * np.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def sample_triangles(sol: StokesSolution, direction: str, b: float, Ls: float) -> np.ndarray:
    """Boolean selector of the triangles whose centroid lies in the sample (between the reservoir planes)."""
    ax = 0 if direction == "lr" else 1
    cen = sol.mesh.p[:, sol.mesh.t].mean(axis=1)
    lo, hi = b - 0.5, b + Ls - 0.5          # identical for lr and tb in the package's coordinate convention
    return (cen[ax] > lo) & (cen[ax] < hi)


def compute(sol: StokesSolution, direction: str, darcy: dict, b: float, Ls: float, alpha=0.01, beta=0.001) -> dict:
    """All descriptors of the module docstring for one solution (b = 0: the whole mesh is the sample)."""
    mesh = sol.mesh
    P, T = mesh.p, mesh.t
    ax = 0 if direction == "lr" else 1
    sel = sample_triangles(sol, direction, b, Ls) if b > 0 else np.ones(T.shape[1], bool)
    Ts = T[:, sel]
    area = triangle_areas(P, Ts)
    A_f = float(area.sum())
    ux, uy = sol.vertex_velocity()
    speed = np.hypot(ux, uy)
    verts = np.unique(Ts)
    u_max = float(speed[verts].max())
    U = darcy["U_darcy"]
    u_T = speed[Ts].mean(axis=0)
    active = float(area[u_T >= alpha * u_max].sum() / A_f)
    dead = float(area[u_T <= beta * u_max].sum() / A_f)
    order = np.argsort(-u_T)
    a_c = np.cumsum(area[order]) / A_f
    w = u_T[order] * area[order]
    tau = np.cumsum(w) / w.sum()
    walls = mesh.facets[:, mesh.boundaries["WALLS"]]
    fm = P[:, walls].mean(axis=1)
    inside = (fm[ax] > b - 0.5) & (fm[ax] < b + Ls - 0.5) if b > 0 else np.ones(walls.shape[1], bool)
    P_w = float(np.linalg.norm(P[:, walls[1, inside]] - P[:, walls[0, inside]], axis=0).sum())
    extent = float(np.ptp(P[1 - ax, verts]))
    D_h = 4.0 * A_f / P_w if P_w > 0 else float("nan")
    return {
        "area_fluid": A_f, "porosity_domain": A_f / (Ls * extent) if extent > 0 else float("nan"),
        "wall_length": P_w, "hydraulic_diameter": D_h,
        "speed_max": u_max, "velocity_amplification": u_max / U,
        "partition": {"active": active, "dead": dead, "intermediate": 1.0 - active - dead, "alpha": alpha, "beta": beta},
        "A50": float(a_c[np.searchsorted(tau, 0.5)]), "A90": float(a_c[np.searchsorted(tau, 0.9)]),
        "kD_times_Dh": darcy["k_D"] * D_h,
        "n_triangles_sample": int(sel.sum()),
    }


def reynolds_bound(k_D_px2: float, D_h_px: float, pixel_size_m: float, mu: float = 1.0e-3, rho: float = 998.0) -> float:
    """Mean pressure gradient [Pa/m] at which Re = rho U D_h / mu reaches 1 (water by default)."""
    return mu ** 2 / (rho * k_D_px2 * D_h_px * pixel_size_m ** 3)


def to_si(k_px2: float, pixel_size_m: float) -> dict:
    """Convert k [px^2] to m^2, darcy and millidarcy (1 darcy = 0.986923e-12 m^2). Use a MEASURED pixel size."""
    k_m2 = k_px2 * pixel_size_m ** 2
    return {"k_m2": k_m2, "k_darcy": k_m2 / 0.986923e-12, "k_mD": 1e3 * k_m2 / 0.986923e-12}
