"""Finite-element core shared by the pipeline and the benchmarks: percolating window -> reservoirs -> mesh -> Stokes -> k_D."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from . import meshing, permeability, stokes
from .config import MeshOptions, SolverOptions
from .geometry import pad_reservoirs


def solve_window(selected: np.ndarray, direction: str, outdir, mesh_opts: MeshOptions, solver_opts: SolverOptions,
                 reservoir_fraction: float = 0.10) -> dict:
    """Solve the Stokes problem on a percolating boolean window (True = fluid) and return
    {"solution", "darcy", "mesh", "b", "L_s", "W", "timings"}."""
    if not selected.any():
        raise ValueError("empty domain: nothing percolates")
    # 1. open reservoirs of width b = round(f * L_s) before the inlet and after the outlet
    padded, Ls, b = pad_reservoirs(selected, direction, reservoir_fraction)
    W = selected.shape[0] if direction == "lr" else selected.shape[1]     # full transverse width of the window
    # 2. contours -> Gmsh triangle mesh with physical groups INLET / OUTLET / WALLS
    t0 = time.time()
    info = meshing.mesh_from_mask(padded, direction, Path(outdir), mesh_opts)
    t_mesh = time.time() - t0
    mesh = meshing.load_mesh(info["mesh_path"])
    # 3. Stokes solution (velocity, pressure, boundary fluxes, residual, energy identity)
    sol = stokes.solve(mesh, direction, solver_opts)
    # 4. permeability: probe lines inside the reservoirs (reference) or, without reservoirs, the boundary curves
    if b > 0 and solver_opts.inlet == "pressure":
        darcy = permeability.darcy_from_reservoirs(sol, direction, b, Ls, W, solver_opts.mu)
    else:
        darcy = permeability.darcy_from_boundaries(sol, direction, L=Ls, W=W, mu=solver_opts.mu)
    return {"solution": sol, "darcy": darcy, "mesh": info, "b": b, "L_s": Ls, "W": W,
            "timings": {"mesh_s": t_mesh, "assembly_s": sol.info["t_assembly_s"], "solve_s": sol.info["t_solve_s"]}}
