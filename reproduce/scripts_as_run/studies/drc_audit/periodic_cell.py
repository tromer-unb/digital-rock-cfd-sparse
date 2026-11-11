"""Periodic unit cell of a square array of circular cylinders: Stokes flow driven by a body force, fully periodic velocity and pressure.

Used as a literature benchmark (Hasimoto 1959, dilute series): dimensionless permeability k/L^2 of the square array as a function of
the solid fraction c = pi R^2 / L^2.  k = mu <u> / f with <u> the superficial (cell-average) velocity and f the body-force density (= -grad p).
Independent of the digital-rock mesher: the geometry is built analytically (no pixels) in Gmsh with periodic edge meshes.
"""
from __future__ import annotations
import numpy as np
import gmsh
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from skfem import (Basis, ElementTriP1, ElementTriP2, ElementVector, LinearForm, MeshTri, asm, bmat, Functional)
from skfem.models.general import divergence
from skfem.models.poisson import vector_laplace


def hasimoto_k(c):
    """k/L^2 of the square array of cylinders, dilute-limit series of Hasimoto (1959), c = solid fraction (O(c^4) remainder).
    Drag per unit length F = 8 pi mu U / (-ln c - 1.476 + 2c - 1.774 c^2 + 4.076 c^3) with U the superficial velocity, as quoted by
    Kallemov et al. (arXiv:1505.07865, Eq. 25: 4 pi/(-ln sqrt(phi) - 0.738 + phi - 0.887 phi^2 + 2.038 phi^3)); k = mu U L^2 / F."""
    c = np.asarray(c, float)
    return (-np.log(c) - 1.476 + 2 * c - 1.774 * c ** 2 + 4.076 * c ** 3) / (8 * np.pi)


def make_mesh(c, h, path):
    R = float(np.sqrt(c / np.pi))
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    try:
        g = gmsh.model.geo; gmsh.model.add("cell")
        p1, p2, p3, p4 = g.addPoint(0, 0, 0, h), g.addPoint(1, 0, 0, h), g.addPoint(1, 1, 0, h), g.addPoint(0, 1, 0, h)
        bot, right, top, left = g.addLine(p1, p2), g.addLine(p2, p3), g.addLine(p4, p3), g.addLine(p1, p4)   # all oriented +x / +y
        cc = g.addPoint(0.5, 0.5, 0, h)
        pts = [g.addPoint(0.5 + R * np.cos(a), 0.5 + R * np.sin(a), 0, h) for a in (0, np.pi / 2, np.pi, 1.5 * np.pi)]
        arcs = [g.addCircleArc(pts[i], cc, pts[(i + 1) % 4]) for i in range(4)]
        outer = g.addCurveLoop([bot, right, -top, -left]); hole = g.addCurveLoop(arcs)
        s = g.addPlaneSurface([outer, hole]); g.synchronize()
        gmsh.model.mesh.setPeriodic(1, [top], [bot], [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1])
        gmsh.model.mesh.setPeriodic(1, [right], [left], [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2); gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return R


def solve_cell(mesh_path, element="p2"):
    mesh = MeshTri.load(mesh_path)
    assert element == "p2", "only Taylor-Hood P2/P1 is implemented for the periodic cell"
    ev = ElementVector(ElementTriP2())
    bu, bp = Basis(mesh, ev, intorder=4), Basis(mesh, ElementTriP1(), intorder=4)
    A = asm(vector_laplace, bu).tocsr(); B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr")
    nu, npp = bu.N, bp.N

    @LinearForm
    def body(v, w):
        return v[0]                      # f = (1, 0)
    f = np.zeros(K.shape[0]); f[:nu] = body.assemble(bu)

    # periodic identification by dof location and component
    def red_map(locs, comp):
        loc = locs.copy()
        loc[0] = np.where(np.abs(loc[0] - 1) < 1e-9, 0.0, loc[0]); loc[1] = np.where(np.abs(loc[1] - 1) < 1e-9, 0.0, loc[1])
        keys = [tuple(np.round(loc[:, i], 8)) + ((int(comp[i]),) if comp is not None else ()) for i in range(loc.shape[1])]
        idx, rep = {}, np.empty(len(keys), int)
        for i, k in enumerate(keys):
            rep[i] = idx.setdefault(k, len(idx))
        return rep, len(idx)
    comp_u = _components(bu)
    rep_u, nru = red_map(bu.doflocs, comp_u)
    rep_p, nrp = red_map(bp.doflocs, None)
    rep = np.concatenate([rep_u, nru + rep_p]); nred = nru + nrp
    P = sp.csr_matrix((np.ones(len(rep)), (np.arange(len(rep)), rep)), shape=(len(rep), nred))
    Kr = (P.T @ K @ P).tocsr(); fr = P.T @ f
    # no-slip on the cylinder: boundary facets that do not lie on the unit-square edges
    cyl = mesh.facets_satisfying(lambda x: (x[0] > 1e-6) & (x[0] < 1 - 1e-6) & (x[1] > 1e-6) & (x[1] < 1 - 1e-6), boundaries_only=True)
    d_cyl = rep[bu.get_dofs(cyl).all()]
    pin = nru                              # pin the first pressure dof (pressure level is arbitrary)
    free = np.setdiff1d(np.arange(nred), np.concatenate([d_cyl, [pin]]))
    Kf = Kr[free][:, free].tocsc()
    sol = np.zeros(nred); sol[free] = splu(Kf).solve(fr[free])
    full = P @ sol
    u = full[:nu]
    ux_int = float(body.assemble(bu) @ u)                                # integral of u_x over the cell (area 1)
    return dict(k=ux_int, n_tri=int(mesh.t.shape[1]), n_dof=int(K.shape[0]), residual=float(np.linalg.norm(Kf @ sol[free] - fr[free]) / np.linalg.norm(fr[free])))


def _components(bu):
    """Vector component (0/1) of every velocity dof, from the nodal/facet dof layout of ElementVector(P2)."""
    comp = -np.ones(bu.N, int)
    for cidx in (0, 1):
        comp[bu.nodal_dofs[cidx]] = cidx
        comp[bu.facet_dofs[cidx]] = cidx
    assert (comp >= 0).all()
    return comp
