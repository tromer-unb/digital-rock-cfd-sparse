"""OpenFOAM cross-check of a user triangle mesh: the same triangles extruded to one layer of prisms.

The finite-volume solution uses the identical geometry and boundary groups (INLET -> inlet, OUTLET -> outlet, WALLS -> walls)
with p = p_in on the inlet and p = 0 on the outlet. Both sides report k = mu Q L / (W dP) with the flux-weighted pressures on
the inlet and outlet curves (FV: the imposed values), L and W from the bounding box unless given. A cell-centred FV scheme on
triangles is less accurate than Taylor-Hood elements on the same mesh, so differences of a few per cent on coarse meshes are
discretization error of the FV side; they must shrink when the mesh is refined.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ...config import OpenFOAMOptions
from . import case as foam_case
from .environment import find_openfoam
from .polymesh import write_polymesh

GROUP_TO_PATCH = {"INLET": "inlet", "OUTLET": "outlet", "WALLS": "walls"}


def crosscheck_mesh(mesh, fe_darcy: dict, workdir, of_opts: OpenFOAMOptions, p_in: float = 1.0) -> dict:
    """Solve the extruded user mesh with simpleFoam and compare its k with the finite-element value ``fe_darcy``."""
    of = find_openfoam(of_opts.env_script)
    workdir = Path(workdir)
    case_dir = workdir / "openfoam_mesh"
    if case_dir.exists():
        shutil.rmtree(case_dir)
    mids, names = [], []
    for g, patch in GROUP_TO_PATCH.items():
        f = mesh.facets[:, mesh.boundaries[g]]
        mids.append(mesh.p[:, f].mean(axis=1).T)
        names += [patch] * f.shape[1]
    tree = cKDTree(np.vstack(mids))
    names = np.array(names, dtype=object)

    def edge_patch(mid):
        d, i = tree.query(mid)
        if np.max(d) > 1e-6:
            raise ValueError("a boundary edge of the mesh carries no INLET/OUTLET/WALLS group")
        return names[i]
    pm = write_polymesh(case_dir, mesh.p.T, mesh.t.T, edge_patch, {"inlet": "patch", "outlet": "patch", "walls": "wall"}, dz=1.0)
    foam_case.write_case(case_dir, nu=of_opts.nu, p_in=p_in, max_iterations=of_opts.max_iterations, nprocs=of_opts.nprocs)
    diag = foam_case.run(of, case_dir, nprocs=of_opts.nprocs)
    t = foam_case.latest_time(case_dir)
    q_out = float(np.sum(foam_case.read_patch_values(t / "phi", "outlet")))
    q_in = float(np.sum(foam_case.read_patch_values(t / "phi", "inlet")))
    Q = 0.5 * (abs(q_in) + abs(q_out))
    L, W = fe_darcy["L_s"], fe_darcy["W"]
    k_fv = of_opts.nu * Q * L / (W * p_in)
    rep = {"status": "completed", "openfoam": of.describe(), "n_cells": pm["n_cells"], "k_fv_prisms": k_fv, "k_fe": fe_darcy["k_D"],
           "rel_diff_fv_vs_fe": k_fv / fe_darcy["k_D"] - 1.0, "flux_mismatch_in_out": abs(abs(q_in) - abs(q_out)) / Q, **diag}
    if not of_opts.keep_cases:
        shutil.rmtree(case_dir, ignore_errors=True)
    (workdir / "openfoam_crosscheck.json").write_text(json.dumps(rep, indent=1, default=float))
    return rep
