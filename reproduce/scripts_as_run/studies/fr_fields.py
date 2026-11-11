#!/usr/bin/env python3
"""
TEMPLATE-ID: none (audit/revision script, not a figure template)
DESCRIPTION: P1 fields for the manuscript (resolves the "[P1 PENDING: recompute from P1 fields]" flags).
  Re-solves the protocol-P1 reference (pressure-driven, open reservoirs b = 0.1 L_s, contour tolerance 0.25 px,
  band 8 %, Taylor-Hood P2/P1, sparse LU) on the SAME level-5 meshes of FR-E1 (fix-round/FR-E1_convergence_tol0.25/mesh)
  and stores vertex fields for the descriptor and figure scripts.
  Acceptance (declared before running): k_D (flux-weighted face pressures) reproduces the stored FR-E1 value
  kD_iface_flux to a relative 1e-9, residual_lu < 1e-10, |Q_balance| < 1e-10. Otherwise the run is rejected.
INPUT:  FR-E1 level-5 meshes and result JSONs.
OUTPUT: fix-round/FR-H_fields_P1/<case>_L5_fields.npz and <case>_L5_fields.json
USAGE:  fr_fields.py --case R3-H
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from independent_check import solve          # noqa: E402
from fr_reservoir import interface_pressures  # noqa: E402

QC = HERE.parent
E1 = QC / "fix-round/FR-E1_convergence_tol0.25"
OUT = QC / "fix-round/FR-H_fields_P1"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--case", required=True, choices=["R3-H", "R2-V", "R3-V"])
    a = ap.parse_args()
    tag = f"{a.case.replace('-', '')}_L5_pad0.1_tol0.25_band0.08"
    ref = json.load(open(E1 / f"{tag}_p2.json"))
    direction, b, Ls, W = ref["direction"], ref["pad_px"], ref["Ls"], ref["W"]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    r = solve(str(E1 / "mesh" / tag / "mesh.msh"), direction, "p2", bc="pressure", return_fields=True)
    f = r.pop("_fields")
    iface = interface_pressures(f, direction, b, Ls, W)
    Q = abs(r["inlet"]["Q"])
    kD = Q * Ls / (W * iface["iface_flux"])
    mesh, u, p, bu = f["mesh"], f["u"], f["p"], f["bu"]
    nd = bu.nodal_dofs                      # vector P2: rows = (u_x, u_y) at vertices
    ux, uy = u[nd[0]], u[nd[1]]
    pv = p[:mesh.p.shape[1]]
    # spot check of the vertex extraction against the interpolator
    # (probe points pulled 1e-6 of the way from a vertex towards the centroid of one of its triangles)
    rng = np.random.default_rng(0); it = rng.choice(mesh.t.shape[1], 50, replace=False)
    vid = mesh.t[0, it]; cen = mesh.p[:, mesh.t[:, it]].mean(axis=1)
    fu = bu.interpolator(u); chk = fu(mesh.p[:, vid] + 1e-6 * (cen - mesh.p[:, vid]))
    err_vertex = float(np.max(np.abs(chk[0] - ux[vid])) / max(np.abs(ux).max(), 1e-300))
    # probe-line flux-weighted pressures (same as interface_pressures, kept for p* normalisation)
    from skfem import Basis, ElementTriP1
    bp = Basis(mesh, ElementTriP1(), intorder=4); fp = bp.interpolator(p[:bp.N])
    tt = np.arange(0.0, W - 1 + 1e-9, 0.5)
    if direction == "lr":
        planes = {"in": b - 0.5 - 1.0, "out": b + Ls - 0.5 + 1.0}
    else:
        top = (b + Ls + b) - 1.0
        planes = {"in": top - (b - 0.5) + 1.0, "out": top - (b + Ls - 0.5) - 1.0}
    pbar = {}
    for k_, c in planes.items():
        pts = np.vstack([np.full_like(tt, c), tt]) if direction == "lr" else np.vstack([tt, np.full_like(tt, c)])
        uu = fu(pts); pp = fp(pts); un = uu[0] if direction == "lr" else uu[1]
        pbar[k_] = float(np.sum(pp * un) / np.sum(un))
    ok = (abs(kD / ref["kD_iface_flux"] - 1) < 1e-9 and r["residual_lu"] < 1e-10 and abs(r["Q_balance"]) < 1e-10)
    rec = dict(case=a.case, tag=tag, direction=direction, pad_px=b, Ls=Ls, W=W, n_tri=r["n_tri"], n_dof=r["n_dof"],
               residual_lu=r["residual_lu"], Q_balance=r["Q_balance"], Q=Q, U_darcy=Q / W, dP_iface_flux=iface["iface_flux"],
               kD=kD, kD_stored=ref["kD_iface_flux"], p_in_bar=pbar["in"], p_out_bar=pbar["out"],
               dissipation=r["dissipation"], vertex_extraction_relerr=err_vertex, t_total_s=time.time() - t0,
               acceptance="PASS" if ok else "FAIL")
    np.savez_compressed(OUT / f"{a.case.replace('-', '')}_L5_fields.npz", points=mesh.p, triangles=mesh.t,
                        ux=ux, uy=uy, p=pv, wall_facets=mesh.facets[:, mesh.boundaries["WALLS"]])
    json.dump(rec, open(OUT / f"{a.case.replace('-', '')}_L5_fields.json", "w"), indent=1)
    print(json.dumps(rec, indent=1), flush=True)


if __name__ == "__main__":
    main()
