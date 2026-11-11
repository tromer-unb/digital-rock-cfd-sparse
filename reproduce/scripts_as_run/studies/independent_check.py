#!/usr/bin/env python3
"""
TEMPLATE-ID: none (audit script, not a figure template)
DESCRIPTION: Independent re-solution of the author's Stokes problem on the author's own Gmsh meshes.
  Solvers: (i) Mini P1b/P1 + sparse direct LU (isolates MINRES algebraic error),
           (ii) Taylor-Hood P2/P1 + direct LU (isolates element error).
  No pressure pin (outlet natural condition fixes the pressure level).
  Reports k under several definitions of the pressure drop and of the Darcy flux.
INPUT: <mesh.msh> (physical groups FLUID, INLET, OUTLET, WALLS), <direction lr|tb>
USAGE: independent_check.py <mesh.msh> <lr|tb> [--element mini|p2] [--u0 1e-4]
"""
import argparse, json, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from skfem import (LinearForm, Basis, ElementTriMini, ElementTriP1, ElementTriP2, ElementVector,
                   FacetBasis, MeshTri, asm, bmat, Functional)
from skfem.helpers import dot
from skfem.models.general import divergence
from skfem.models.poisson import vector_laplace


def solve(mesh_path, direction, element, u0=1e-4, mu=1.0, bc="velocity", p_in=1.0, return_fields=False):
    mesh = MeshTri.load(mesh_path)
    ev = ElementVector(ElementTriMini() if element == "mini" else ElementTriP2())
    bu, bp = Basis(mesh, ev, intorder=4), Basis(mesh, ElementTriP1(), intorder=4)
    A = (mu * asm(vector_laplace, bu)).tocsr()
    B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr")
    n_u = bu.N
    fin = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)
    sign = 1.0 if direction == "lr" else -1.0
    comp = 0 if direction == "lr" else 1
    def vel(x):
        v = np.zeros((2,) + x.shape[1:])
        v[comp] = sign * u0
        return v
    proj = fin.project(vel)
    d_in = bu.get_dofs("INLET").all(); d_w = bu.get_dofs("WALLS").all()
    x = np.zeros(K.shape[0])
    if bc == "velocity":
        x[d_in] = proj[d_in]; x[d_w] = 0.0
        pres = np.unique(np.concatenate([d_in, d_w]))
    else:  # pressure-driven: traction -p_in n on INLET, traction-free OUTLET, no-slip walls only
        pres = np.unique(d_w)
    free = np.setdiff1d(np.arange(K.shape[0]), pres)
    rhs = -K[free][:, pres] @ x[pres]
    if bc == "pressure":
        @LinearForm
        def trac(v, w):
            return -p_in * dot(w.n, v)
        f = np.zeros(K.shape[0]); f[:n_u] = trac.assemble(fin)
        rhs = rhs + f[free]
    Kf = K[free][:, free].tocsc()
    t0 = time.time(); lu = splu(Kf); xs = lu.solve(rhs); t_lu = time.time() - t0
    x[free] = xs
    res = np.linalg.norm(Kf @ xs - rhs) / np.linalg.norm(rhs)
    u, p = x[:n_u], x[n_u:]
    # boundary integrals
    out = {}
    for name in ("INLET", "OUTLET"):
        fb_u = FacetBasis(mesh, ev, facets=mesh.boundaries[name], intorder=4)
        fb_p = FacetBasis(mesh, ElementTriP1(), facets=mesh.boundaries[name], intorder=4)
        @Functional
        def flux(w):  # u.n
            return dot(w["u"], w.n)
        @Functional
        def one(w):
            return 1.0 + 0 * w.x[0]
        @Functional
        def pint(w):
            return w["p"]
        @Functional
        def pflux(w):
            return w["p"] * dot(w["u"], w.n)
        uu = fb_u.interpolate(u)
        # p on a facet basis that shares the same quadrature: build with same element-mesh facets
        pp = fb_p.interpolate(p)
        L = one.assemble(fb_u)
        Qn = flux.assemble(fb_u, u=uu)
        out[name] = dict(L=float(L), Q=float(Qn), p_len=float(pint.assemble(fb_p, p=pp) / L))
        dofs = bp.get_dofs(name).all()
        out[name]["p_nodal"] = float(np.mean(p[dofs]))
    # flux-weighted pressure via a common facet basis (element = vector P2/Mini, pressure interpolated separately)
    for name in ("INLET", "OUTLET"):
        fb = FacetBasis(mesh, ev, facets=mesh.boundaries[name], intorder=4)
        fbp = FacetBasis(mesh, ElementTriP1(), facets=mesh.boundaries[name], intorder=4)
        @Functional
        def pf(w):
            return w["p"] * dot(w["u"], w.n)
        out[name]["p_flux"] = float(pf.assemble(fb, u=fb.interpolate(u), p=fbp.interpolate(p)) / out[name]["Q"])
    phi = float(mu * u @ (asm(vector_laplace, bu) @ u))  # total viscous dissipation integral |grad u|^2
    tot_Q = abs(out["INLET"]["Q"])
    flow_len = float(mesh.p[0 if direction == "lr" else 1].max() - mesh.p[0 if direction == "lr" else 1].min())
    W = float(mesh.p[1 if direction == "lr" else 0].max() - mesh.p[1 if direction == "lr" else 0].min())
    Lin = out["INLET"]["L"]
    k = {}
    for tag, dP in (("nodal", abs(out["INLET"]["p_nodal"] - out["OUTLET"]["p_nodal"])),
                    ("length", abs(out["INLET"]["p_len"] - out["OUTLET"]["p_len"])),
                    ("flux", abs(out["INLET"]["p_flux"] - out["OUTLET"]["p_flux"])),
                    ("energy", phi / tot_Q)):
        k[tag] = dict(dP=dP,
                      k_author_def=mu * (tot_Q / Lin) * flow_len / dP,   # Q / inlet opening
                      k_darcy_total=mu * (tot_Q / W) * flow_len / dP)     # Q / full transverse width
    res_dict = dict(bc=bc, element=element, n_tri=int(mesh.t.shape[1]), n_dof=int(K.shape[0]), residual_lu=float(res),
                t_lu_s=t_lu, flow_len=flow_len, W=W, inlet=out["INLET"], outlet=out["OUTLET"],
                Q_balance=(out["OUTLET"]["Q"] + out["INLET"]["Q"]) / tot_Q, dissipation=phi, k=k)
    if return_fields:
        res_dict['_fields'] = dict(mesh=mesh, u=u, p=p, bu=bu)
    return res_dict


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh"); ap.add_argument("direction")
    ap.add_argument("--element", default="mini", choices=["mini", "p2"])
    ap.add_argument("--u0", type=float, default=1e-4)
    ap.add_argument("--bc", default="velocity", choices=["velocity", "pressure"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    r = solve(a.mesh, a.direction, a.element, a.u0, bc=a.bc)
    s = json.dumps(r, indent=1)
    if a.out: open(a.out, "w").write(s)
    print(s)
