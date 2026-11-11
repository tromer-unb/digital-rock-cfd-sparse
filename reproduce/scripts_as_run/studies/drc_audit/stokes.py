"""Reference Stokes solver (direct sparse LU) on the Gmsh meshes produced by the author's mesher.

Differences to the author's solver (legacy pipeline):
  * direct factorisation (no MINRES), so the algebraic error is at round-off level;
  * Taylor-Hood P2/P1 (default) or Mini P1b/P1;
  * inlet: uniform velocity ("velocity", as in the legacy code) or imposed pressure ("pressure": traction -p n);
  * lateral boundaries of the window: no-slip ("noslip", as in the legacy code) or free-slip ("freeslip": u.n = 0, zero shear traction),
    which is mathematically the same as a mirror-symmetric (or periodic mirrored) cell;
  * no pressure pin: the traction-free outlet fixes the pressure level.
Units: whatever the mesh uses (pixels in this project); mu = 1 unless given.
"""
from __future__ import annotations
import time
import numpy as np
from scipy.sparse.linalg import splu
from skfem import (LinearForm, Basis, ElementTriMini, ElementTriP1, ElementTriP2, ElementVector, FacetBasis, MeshTri, asm, bmat, Functional)
from skfem.helpers import dot
from skfem.models.general import divergence
from skfem.models.poisson import vector_laplace


def _component_dofs(bu, facets, comp):
    dv = bu.get_dofs(facets)
    key = f"u^{comp + 1}"
    out = [np.asarray(dv.nodal[key], int)]
    if key in dv.facet:
        out.append(np.asarray(dv.facet[key], int))
    return np.unique(np.concatenate(out))


def lateral_facets(mesh, direction, tol=0.75):
    """Boundary facets of group WALLS lying on the two lateral edges of the window (transverse coordinate at its extremes)."""
    t = 1 if direction == "lr" else 0
    lo, hi = mesh.p[t].min(), mesh.p[t].max()
    walls = np.asarray(mesh.boundaries["WALLS"], int)
    nodes = mesh.facets[:, walls]                      # (2, nf)
    c = mesh.p[t][nodes]
    on = (np.abs(c - lo) <= tol).all(axis=0) | (np.abs(c - hi) <= tol).all(axis=0)
    return walls[on], walls[~on]


def solve(mesh_path, direction, element="p2", bc="pressure", lateral="noslip", u0=1e-4, mu=1.0, p_in=1.0, return_fields=False, mesh=None):
    assert bc in ("velocity", "pressure") and lateral in ("noslip", "freeslip") and element in ("p2", "mini")
    mesh = mesh if mesh is not None else MeshTri.load(mesh_path)
    ev = ElementVector(ElementTriMini() if element == "mini" else ElementTriP2())
    bu, bp = Basis(mesh, ev, intorder=4), Basis(mesh, ElementTriP1(), intorder=4)
    A = (mu * asm(vector_laplace, bu)).tocsr()
    B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr")
    n_u = bu.N
    fin = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)
    sign = 1.0 if direction == "lr" else -1.0
    comp_flow = 0 if direction == "lr" else 1
    comp_trans = 1 - comp_flow

    def vel(x):
        v = np.zeros((2,) + x.shape[1:]); v[comp_flow] = sign * u0; return v

    d_in = bu.get_dofs("INLET").all()
    if lateral == "freeslip":
        lat, noslip = lateral_facets(mesh, direction)
        d_w = bu.get_dofs(noslip).all() if len(noslip) else np.array([], int)
        d_slip = _component_dofs(bu, lat, comp_trans)
        d_w = np.unique(np.concatenate([d_w, d_slip]))      # corners shared with no-slip facets are already fully fixed
        corner = bu.get_dofs(noslip).all() if len(noslip) else np.array([], int)
        full_fixed = np.unique(corner)                       # both components zero
    else:
        d_w = bu.get_dofs("WALLS").all(); full_fixed = d_w
    x = np.zeros(K.shape[0])
    if bc == "velocity":
        x[d_in] = fin.project(vel)[d_in]
        x[full_fixed] = 0.0; x[d_w] = 0.0
        fixed = np.unique(np.concatenate([d_in, d_w]))
    else:
        fixed = np.unique(d_w)
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    rhs = -K[free][:, fixed] @ x[fixed]
    if bc == "pressure":
        @LinearForm
        def trac(v, w):
            return -p_in * dot(w.n, v)
        f = np.zeros(K.shape[0]); f[:n_u] = trac.assemble(fin)
        rhs = rhs + f[free]
    Kf = K[free][:, free].tocsc()
    t0 = time.time(); xs = splu(Kf).solve(rhs); t_lu = time.time() - t0
    x[free] = xs
    res = float(np.linalg.norm(Kf @ xs - rhs) / np.linalg.norm(rhs))
    u, p = x[:n_u], x[n_u:]

    out = {}
    for name in ("INLET", "OUTLET"):
        fb = FacetBasis(mesh, ev, facets=mesh.boundaries[name], intorder=4)
        fbp = FacetBasis(mesh, ElementTriP1(), facets=mesh.boundaries[name], intorder=4)

        @Functional
        def one(w): return 1.0 + 0 * w.x[0]

        @Functional
        def flux(w): return dot(w["u"], w.n)

        @Functional
        def pint(w): return w["p"]

        @Functional
        def pf(w): return w["p"] * dot(w["u"], w.n)
        L = one.assemble(fb); Qn = flux.assemble(fb, u=fb.interpolate(u))
        pp = fbp.interpolate(p)
        out[name] = dict(L=float(L), Q=float(Qn), p_len=float(pint.assemble(fbp, p=pp) / L),
                         p_flux=float(pf.assemble(fb, u=fb.interpolate(u), p=pp) / Qn),
                         p_nodal=float(np.mean(p[bp.get_dofs(name).all()])))      # unweighted nodal mean: the legacy definition
    phi = float(mu * u @ (asm(vector_laplace, bu) @ u))
    Q = abs(out["INLET"]["Q"])
    flow_axis = 0 if direction == "lr" else 1
    flow_len = float(mesh.p[flow_axis].max() - mesh.p[flow_axis].min())
    Wd = float(mesh.p[1 - flow_axis].max() - mesh.p[1 - flow_axis].min())
    res_d = dict(bc=bc, lateral=lateral, element=element, n_tri=int(mesh.t.shape[1]), n_dof=int(K.shape[0]), residual_lu=res, t_lu_s=t_lu,
                 flow_len=flow_len, W=Wd, inlet=out["INLET"], outlet=out["OUTLET"], Q=Q,
                 Q_balance=float((out["OUTLET"]["Q"] + out["INLET"]["Q"]) / Q), dissipation=phi,
                 dP_flux=abs(out["INLET"]["p_flux"] - out["OUTLET"]["p_flux"]), dP_len=abs(out["INLET"]["p_len"] - out["OUTLET"]["p_len"]),
                 dP_nodal=abs(out["INLET"]["p_nodal"] - out["OUTLET"]["p_nodal"]), dP_energy=phi / Q)
    if return_fields:
        res_d["_fields"] = dict(mesh=mesh, u=u, p=p, bu=bu)
    return res_d


def interface_pressure_drop(fields, direction, b, Ls, W, eps=1.0, dy=0.5):
    """dP between the two sample faces, flux-weighted pressure probed `eps` px inside the open reservoirs (all probe points are fluid).
    b: reservoir width [px]; Ls: sample length [px]; W: transverse width [px]. Returns dict(len=, flux=)."""
    mesh, u, p, bu = fields["mesh"], fields["u"], fields["p"], fields["bu"]
    bp = Basis(mesh, ElementTriP1(), intorder=4)
    fu, fp = bu.interpolator(u), bp.interpolator(p)
    t = np.arange(0.0, W - 1 + 1e-9, dy)
    if direction == "lr":
        planes = {"in": b - 0.5 - eps, "out": b + Ls - 0.5 + eps}
    else:   # rows run top to bottom, y runs upward: inlet reservoir at large y
        top = (b + Ls + b) - 1.0
        planes = {"in": top - (b - 0.5) + eps, "out": top - (b + Ls - 0.5) - eps}
    v = {}
    for k, c in planes.items():
        pts = np.vstack([np.full_like(t, c), t]) if direction == "lr" else np.vstack([t, np.full_like(t, c)])
        uu, pp = fu(pts), fp(pts)
        un = uu[0] if direction == "lr" else uu[1]
        v[k] = (float(np.mean(pp)), float(np.sum(pp * un) / np.sum(un)))
    return dict(len=abs(v["in"][0] - v["out"][0]), flux=abs(v["in"][1] - v["out"][1]))
