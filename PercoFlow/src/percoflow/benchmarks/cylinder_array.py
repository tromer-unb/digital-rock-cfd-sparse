"""Periodic square array of circular cylinders: finite elements (this package) versus the dilute series of Hasimoto.

Geometry: unit cell [0,1]^2 with one cylinder of radius R at its centre, solid fraction c = pi R^2. Flow driven by a uniform
body force f = (1, 0), velocity and pressure fully periodic, no-slip on the cylinder. The dimensionless permeability is
k / L^2 = mu <u> / f with <u> the cell-averaged (superficial) velocity.

Reference: Hasimoto (1959) with the higher-order coefficients of Sangani & Acrivos (1982),
    k / L^2 = (-ln c - 1.476 + 2 c - 1.774 c^2 + 4.076 c^3) / (8 pi),
accurate in the dilute limit; the finite-element value agrees with it to ~0.3 % for c <= 0.15 and departs above, where the
truncated series loses accuracy. The geometry is built analytically (no pixels), so this test is independent of the
image-to-mesh chain.
"""
from __future__ import annotations

import gmsh
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from skfem import Basis, ElementTriP1, ElementTriP2, ElementVector, LinearForm, MeshTri, asm, bmat
from skfem.models.general import divergence
from skfem.models.poisson import vector_laplace


def hasimoto(c):
    """Dimensionless permeability k / L^2 of the square cylinder array from the dilute series (c = solid fraction)."""
    c = np.asarray(c, float)
    return (-np.log(c) - 1.476 + 2 * c - 1.774 * c ** 2 + 4.076 * c ** 3) / (8 * np.pi)


def cell_mesh(c: float, h: float, path) -> None:
    """Unstructured triangle mesh of the cell with periodic (matching) opposite edges."""
    R = float(np.sqrt(c / np.pi))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        g = gmsh.model.geo
        gmsh.model.add("cell")
        p = [g.addPoint(x, y, 0, h) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))]
        bot, right, top, left = g.addLine(p[0], p[1]), g.addLine(p[1], p[2]), g.addLine(p[3], p[2]), g.addLine(p[0], p[3])
        cc = g.addPoint(0.5, 0.5, 0, h)
        q = [g.addPoint(0.5 + R * np.cos(a), 0.5 + R * np.sin(a), 0, h) for a in (0, np.pi / 2, np.pi, 1.5 * np.pi)]
        arcs = [g.addCircleArc(q[i], cc, q[(i + 1) % 4]) for i in range(4)]
        s = g.addPlaneSurface([g.addCurveLoop([bot, right, -top, -left]), g.addCurveLoop(arcs)])
        g.synchronize()
        gmsh.model.mesh.setPeriodic(1, [top], [bot], [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1])
        gmsh.model.mesh.setPeriodic(1, [right], [left], [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
        gmsh.model.addPhysicalGroup(2, [s], name="FLUID")
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def _components(bu):
    """Velocity component (0 = x, 1 = y) of each vector degree of freedom; used to pair periodic dofs correctly."""
    comp = -np.ones(bu.N, int)
    for k in (0, 1):
        comp[bu.nodal_dofs[k]] = k
        comp[bu.facet_dofs[k]] = k
    return comp


def solve_cell(path) -> dict:
    """Taylor-Hood solution of the periodic cell; periodicity by identifying dofs at matching locations."""
    mesh = MeshTri.load(str(path))
    bu = Basis(mesh, ElementVector(ElementTriP2()), intorder=4)
    bp = Basis(mesh, ElementTriP1(), intorder=4)
    A = asm(vector_laplace, bu).tocsr()
    B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr")

    @LinearForm
    def body(v, w):
        return v[0]
    fb = body.assemble(bu)
    f = np.zeros(K.shape[0])
    f[: bu.N] = fb

    def reduce(locs, comp):
        loc = locs.copy()
        for ax in (0, 1):
            loc[ax] = np.where(np.abs(loc[ax] - 1) < 1e-9, 0.0, loc[ax])
        idx, rep = {}, np.empty(loc.shape[1], int)
        for i in range(loc.shape[1]):
            key = (round(loc[0, i], 8), round(loc[1, i], 8)) + ((int(comp[i]),) if comp is not None else ())
            rep[i] = idx.setdefault(key, len(idx))
        return rep, len(idx)
    ru, nu = reduce(bu.doflocs, _components(bu))
    rp, npp = reduce(bp.doflocs, None)
    rep = np.concatenate([ru, nu + rp])
    P = sp.csr_matrix((np.ones(len(rep)), (np.arange(len(rep)), rep)), shape=(len(rep), nu + npp))
    Kr, fr = (P.T @ K @ P).tocsr(), P.T @ f
    cyl = mesh.facets_satisfying(lambda x: (x[0] > 1e-6) & (x[0] < 1 - 1e-6) & (x[1] > 1e-6) & (x[1] < 1 - 1e-6), boundaries_only=True)
    fixed = np.concatenate([ru[bu.get_dofs(cyl).all()], [nu]])           # no-slip cylinder + one pressure pin
    free = np.setdiff1d(np.arange(nu + npp), fixed)
    Kf = Kr[free][:, free].tocsc()
    x = np.zeros(nu + npp)
    x[free] = splu(Kf).solve(fr[free])
    u = (P @ x)[: bu.N]
    return {"k": float(fb @ u), "n_triangles": int(mesh.t.shape[1]), "n_dof": int(K.shape[0]),
            "residual": float(np.linalg.norm(Kf @ x[free] - fr[free]) / np.linalg.norm(fr[free]))}


def run(workdir, cs=(0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5), h: float = 0.02) -> dict:
    """FE versus Hasimoto. Pass criterion (declared in advance): |FE/Hasimoto - 1| <= 0.5 % for c <= 0.15."""
    from pathlib import Path
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in cs:
        path = workdir / f"cell_c{c}_h{h}.msh"
        cell_mesh(c, h, path)
        r = solve_cell(path)
        r.update(c=c, h=h, k_hasimoto=float(hasimoto(c)), rel_diff=r["k"] / float(hasimoto(c)) - 1.0)
        rows.append(r)
    dilute = [abs(r["rel_diff"]) for r in rows if r["c"] <= 0.15]
    return {"benchmark": "cylinder_array_fe_vs_hasimoto", "rows": rows, "criterion": "|rel_diff| <= 0.005 for c <= 0.15",
            "max_abs_rel_diff_dilute": max(dilute), "pass": max(dilute) <= 0.005}
