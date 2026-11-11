#!/usr/bin/env python3
"""C3 - the percolation filter as a well-posedness condition, not only a saving (npj revision, 2026-11-03).

Theory (Supplementary Note of the npj draft): let the fluid domain split into connected components Omega_i.
  (i) A component whose boundary is entirely no-slip wall (an ISOLATED cluster) admits u = 0, p = const_i as a non-trivial
      solution of the homogeneous discrete problem: the vector x = (0, 1 on the pressure dofs of Omega_i) is in the kernel
      of the reduced saddle-point matrix, because b(v, 1) = int div v = int_dOmega_i v.n = 0 for every admissible v.
      The system is singular (one kernel vector per isolated cluster), whatever the inlet condition.
  (ii) With a prescribed inlet velocity, a component that touches the inlet but not the outlet (a DEAD END) receives a
      non-zero imposed flux with nowhere to leave. The same constant-pressure vector y of that component lies in the
      kernel of the symmetric reduced matrix, and y.rhs = -(imposed flux into the component) != 0, so the discrete system
      is inconsistent and has no solution.
  (iii) After the percolation filter every component touches both planes and the pressure-driven problem has a unique
      solution.
Test geometry: 120 x 600 px mask with a straight percolating channel, an isolated pore blob and a dead end attached to
the inlet plane only.
Declared criteria (before running):
  A  pressure-driven + reservoirs, unfiltered: ||K x|| / ||K|| ||x|| < 1e-12 for the isolated-cluster vector x;
  B  velocity-driven, unfiltered: |y.rhs| / (||y|| ||rhs||) > 1e-6 for the dead-end vector y (inconsistency), and the
     MINRES relative residual stays above 1e-6 after 20,000 iterations;
  C  pressure-driven + reservoirs, filtered: LU relative residual < 1e-10 and mass imbalance < 1e-10.
Also reported: cluster statistics of the three rock images (8-connectivity, band 8%), computed with the pipeline's own
functions.
OUTPUT: results/c3_wellposedness.json
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu, minres
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import reservoir                                                    # noqa: E402
from drc_audit.legacy import load_pipeline                                        # noqa: E402
from skfem import (Basis, ElementTriMini, ElementTriP1, ElementVector, FacetBasis, MeshTri, asm, bmat, LinearForm)  # noqa: E402
from skfem.helpers import dot                                                     # noqa: E402
from skfem.models.general import divergence                                       # noqa: E402
from skfem.models.poisson import vector_laplace                                   # noqa: E402

m = load_pipeline(); tmp = Path(tempfile.mkdtemp(prefix="c3_"))
H, L = 120, 600
mask = np.zeros((H, L), bool)
mask[40:80, :] = True            # percolating channel
mask[5:25, 200:300] = True       # isolated blob
mask[92:112, 0:150] = True       # dead end, touches the inlet plane (x = 0) only


def assemble(mesh_path, bc):
    mesh = MeshTri.load(str(mesh_path))
    ev = ElementVector(ElementTriMini()); bu = Basis(mesh, ev, intorder=4); bp = Basis(mesh, ElementTriP1(), intorder=4)
    A = asm(vector_laplace, bu).tocsr(); B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr"); nu = bu.N
    fin = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)
    x = np.zeros(K.shape[0]); d_w = bu.get_dofs("WALLS").all()
    if bc == "velocity":
        d_in = bu.get_dofs("INLET").all()
        def vel(xx):
            v = np.zeros((2,) + xx.shape[1:]); v[0] = 1e-4; return v
        x[d_in] = fin.project(vel)[d_in]; x[d_w] = 0.0; fixed = np.unique(np.concatenate([d_in, d_w]))
    else:
        fixed = np.unique(d_w)
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    rhs = -K[free][:, fixed] @ x[fixed]
    if bc == "pressure":
        @LinearForm
        def trac(v, w):
            return -1.0 * dot(w.n, v)
        f = np.zeros(K.shape[0]); f[:nu] = trac.assemble(fin); rhs = rhs + f[free]
    # connected components of the mesh (vertex graph) -> pressure dofs of P1 are the vertices
    t = mesh.t; nv = mesh.p.shape[1]
    rows = np.concatenate([t[0], t[1], t[2]]); cols = np.concatenate([t[1], t[2], t[0]])
    G = sp.coo_matrix((np.ones(rows.size), (rows, cols)), shape=(nv, nv))
    ncomp, lab = connected_components(G, directed=False)
    return dict(mesh=mesh, K=K, free=free, rhs=rhs, nu=nu, ncomp=ncomp, lab=lab)


def comp_vector(s, comp):
    """constant pressure 1 on component comp, zero velocity, restricted to the free dofs"""
    full = np.zeros(s["K"].shape[0]); full[s["nu"] + np.where(s["lab"] == comp)[0]] = 1.0
    return full[s["free"]]


def classify_components(s):
    mesh = s["mesh"]; out = {}
    xs = mesh.p[0]
    for c in range(s["ncomp"]):
        v = s["lab"] == c
        out[c] = dict(n_vert=int(v.sum()), touches_in=bool(np.any(xs[v] <= xs.min() + 0.5)),
                      touches_out=bool(np.any(xs[v] >= xs.max() - 0.5)), y_mean=float(mesh.p[1][v].mean()))
    return out


res = {}
cfg = reservoir.build_cfg(m, ".", tmp, level=3, tol=0.25)
# A: pressure-driven + reservoirs, unfiltered
padded, Ls, b = reservoir.pad_reservoirs(mask, "lr", 0.10)
mi = m.create_gmsh_mesh(padded, "lr", tmp / "A", cfg) if (tmp / "A").mkdir() is None else None
sA = assemble(mi["mesh_path"], "pressure"); cA = classify_components(sA)
iso = [c for c, d in cA.items() if not d["touches_in"] and not d["touches_out"]]
KfA = sA["K"][sA["free"]][:, sA["free"]]
normK = sp.linalg.norm(KfA)
A_rows = []
for c in iso:
    xv = comp_vector(sA, c); A_rows.append(dict(component=c, y_mean=cA[c]["y_mean"], rel_Kx=float(np.linalg.norm(KfA @ xv) / (normK * np.linalg.norm(xv)))))
res["A_pressure_unfiltered"] = dict(n_components=sA["ncomp"], components=cA, isolated=A_rows,
                                    PASS=bool(len(A_rows) == 1 and all(r["rel_Kx"] < 1e-12 for r in A_rows)))
# B: velocity-driven, no reservoirs, unfiltered
(tmp / "B").mkdir(); mi = m.create_gmsh_mesh(mask, "lr", tmp / "B", cfg)
sB = assemble(mi["mesh_path"], "velocity"); cB = classify_components(sB)
dead = [c for c, d in cB.items() if d["touches_in"] and not d["touches_out"]]
KfB = sB["K"][sB["free"]][:, sB["free"]]
B_rows = []
for c in dead:
    y = comp_vector(sB, c)
    B_rows.append(dict(component=c, y_mean=cB[c]["y_mean"], rel_y_rhs=float(abs(y @ sB["rhs"]) / (np.linalg.norm(y) * np.linalg.norm(sB["rhs"]))),
                       rel_Ky=float(np.linalg.norm(KfB @ y) / (sp.linalg.norm(KfB) * np.linalg.norm(y)))))
xs, info = minres(KfB, sB["rhs"], rtol=1e-12, maxiter=20000)
rres = float(np.linalg.norm(KfB @ xs - sB["rhs"]) / np.linalg.norm(sB["rhs"]))
res["B_velocity_unfiltered"] = dict(n_components=sB["ncomp"], components=cB, dead_ends=B_rows, minres_info=int(info), minres_rel_residual=rres,
                                    PASS=bool(len(B_rows) == 1 and B_rows[0]["rel_y_rhs"] > 1e-6 and rres > 1e-6))
# C: pressure-driven + reservoirs, filtered with the pipeline's percolation selection
sel, info_sel = m.select_full_percolating_components(mask, "lr", cfg)
padded, Ls, b = reservoir.pad_reservoirs(sel, "lr", 0.10)
(tmp / "C").mkdir(); mi = m.create_gmsh_mesh(padded, "lr", tmp / "C", cfg)
sC = assemble(mi["mesh_path"], "pressure")
KfC = sC["K"][sC["free"]][:, sC["free"]].tocsc(); xc = splu(KfC).solve(sC["rhs"])
rC = float(np.linalg.norm(KfC @ xc - sC["rhs"]) / np.linalg.norm(sC["rhs"]))
r = reservoir.run_mask(sel, "lr", tmp / "C2", level=3, tol=0.25, pad=0.10, elements=("p2",), m=m)[0]
res["C_pressure_filtered"] = dict(n_components=sC["ncomp"], lu_rel_residual_mini=rC, p2_residual_lu=r["residual_lu"], p2_Q_balance=r["Q_balance"],
                                  kD_face=r["kD_face"], kD_channel_expected=40 ** 2 / 12 * (L / (L + 2.0)) * 40 / H,
                                  PASS=bool(rC < 1e-10 and r["residual_lu"] < 1e-10 and abs(r["Q_balance"]) < 1e-10))

# cluster statistics of the rock images (pipeline functions, 8-connectivity, band 8%)
stats = {}
SD = ROOT / "data" / "structures"
for name, stem in (("R1", "patch_y3800_x3800"), ("R2", "patch_y7600_x19000"), ("R3", "patch_y7600_x53200")):
    cfg8 = reservoir.build_cfg(m, SD, tmp, level=2, tol=0.25, band=0.08)
    pore = m.clean_mask(m.load_pore_mask(SD / f"{stem}_c0_mask.png", cfg8), cfg8)
    for direction in ("lr", "tb"):
        crop, _ = m.crop_between_fixed_planes(pore, direction, cfg8)
        lab, n = ndi.label(crop, structure=np.ones((3, 3), int)); cnt = np.bincount(lab.ravel())
        a_in = set(np.unique(lab[:, 0] if direction == "lr" else lab[0, :])) - {0}
        a_out = set(np.unique(lab[:, -1] if direction == "lr" else lab[-1, :])) - {0}
        perc = a_in & a_out; only_in = a_in - a_out; only_out = a_out - a_in; iso_l = set(range(1, n + 1)) - a_in - a_out
        tot = cnt[1:].sum()
        def area(s):
            return float(sum(cnt[i] for i in s) / tot) if tot else 0.0
        stats[f"{name}-{'H' if direction == 'lr' else 'V'}"] = dict(n_clusters=int(n), n_percolating=len(perc), n_dead_end_inlet=len(only_in),
                                                                    n_dead_end_outlet=len(only_out), n_isolated=len(iso_l),
                                                                    area_frac_percolating=area(perc), area_frac_dead_end=area(only_in | only_out),
                                                                    area_frac_isolated=area(iso_l), pore_px=int(tot))
res["rock_cluster_statistics"] = stats
res["PASS_all"] = bool(res["A_pressure_unfiltered"]["PASS"] and res["B_velocity_unfiltered"]["PASS"] and res["C_pressure_filtered"]["PASS"])
json.dump(res, open(ROOT / "results" / "c3_wellposedness.json", "w"), indent=1, default=float)
print(json.dumps({k: (v if k == "rock_cluster_statistics" else {kk: vv for kk, vv in v.items() if kk != "components"}) if isinstance(v, dict) else v
                  for k, v in res.items()}, indent=1, default=float))
