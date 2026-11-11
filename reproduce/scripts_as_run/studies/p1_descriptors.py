#!/usr/bin/env python3
"""
TEMPLATE-ID: none (revision script, not a figure template)
DESCRIPTION: Pore-scale descriptors of the manuscript (Results "Resolved fields", "Hydraulic partition",
  "Integrated response"; Methods "Pore-scale descriptors"; Supplementary Note 9) evaluated on the protocol-P1 fields
  (fix-round/FR-H_fields_P1). Definitions follow the manuscript and the legacy pipeline (calculate_metrics):
    - only triangles of the sample (centroid between the two reservoir/sample planes) enter, reservoirs excluded;
    - vertex speed |u| from the Taylor-Hood vertex values; triangle speed = mean of its three vertex speeds;
    - |u|_max = max vertex speed in the sample; active: >= 0.01 |u|_max, dead: <= 0.001 |u|_max;
    - U = Q/W (mean Darcy velocity), p* = (p - Pbar_out)/dP with the flux-weighted probe-line pressures of k_D;
    - D_h = 4 A_f / P_w, P_w = wall length inside the sample (pore walls + lateral window edges);
    - phi = A_f / (L_s * transverse extent of the domain), the bounding-rectangle porosity of the legacy pipeline; A_50, A_90 from triangles ranked by decreasing speed, weight |u|_T |T|;
    - velocity amplification = |u|_max / U.
  Acceptance (declared before running): field JSON acceptance == PASS for all three configurations.
OUTPUT: fix-round/FR-H_fields_P1/p1_descriptors.json and <case>_tri.npz (per-triangle arrays for the figures)
"""
import json
from pathlib import Path
import numpy as np

QC = Path(__file__).resolve().parents[1]
D = QC / "fix-round/FR-H_fields_P1"
CASES = ["R3-H", "R2-V", "R3-V"]


def tri_area(P, T):
    a, b, c = P[:, T[0]], P[:, T[1]], P[:, T[2]]
    return 0.5 * np.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def run(case):
    cc = case.replace("-", "")
    meta = json.load(open(D / f"{cc}_L5_fields.json"))
    assert meta["acceptance"] == "PASS", case
    z = np.load(D / f"{cc}_L5_fields.npz")
    P, T, ux, uy, p, wf = z["points"], z["triangles"], z["ux"], z["uy"], z["p"], z["wall_facets"]
    b, Ls, W, d = meta["pad_px"], meta["Ls"], meta["W"], meta["direction"]
    ax = 0 if d == "lr" else 1
    lo, hi = b - 0.5, b + Ls - 0.5                 # sample planes on the flow axis (both orientations)
    cen = P[:, T].mean(axis=1)
    ins = (cen[ax] > lo) & (cen[ax] < hi)
    Ts = T[:, ins]
    area = tri_area(P, Ts)
    sp = np.hypot(ux, uy)
    vmask = np.zeros(P.shape[1], bool); vmask[Ts.ravel()] = True
    umax = float(sp[vmask].max())
    U = meta["U_darcy"]
    ut = sp[Ts].mean(axis=0)
    A = area.sum()
    act = area[ut >= 0.01 * umax].sum() / A; dead = area[ut <= 0.001 * umax].sum() / A
    # wall length inside the sample
    fm = P[:, wf].mean(axis=1)
    win = (fm[ax] > lo) & (fm[ax] < hi)
    Pw = float(np.linalg.norm(P[:, wf[1, win]] - P[:, wf[0, win]], axis=0).sum())
    Dh = 4 * A / Pw
    # transport localization
    o = np.argsort(-ut); a_c = np.cumsum(area[o]) / A; tw = ut[o] * area[o]; tau = np.cumsum(tw) / tw.sum()
    A50 = float(a_c[np.searchsorted(tau, 0.5)]); A90 = float(a_c[np.searchsorted(tau, 0.9)])
    # normalised pressure and axial coordinate (s = 0 at the inlet plane)
    dP = meta["dP_iface_flux"]; pstar_v = (p - meta["p_out_bar"]) / dP
    if d == "lr":
        s = cen[0, ins] - lo
    else:
        s = hi - cen[1, ins]
    pstar_t = pstar_v[Ts].mean(axis=0)
    vin = np.where(vmask)[0]
    frac_lo = float(np.mean(pstar_v[vin] < 0)); frac_hi = float(np.mean(pstar_v[vin] > 1))
    # location of the speed maximum in the cropped frame (flow axis shifted by b)
    iv = vin[np.argmax(sp[vin])]; xmax = P[:, iv].copy(); xmax[ax] -= b
    # domain porosity as in the legacy pipeline: A_f over the bounding rectangle of the domain (flow length L_s)
    tr = 1 - ax; ext = float(P[tr, vin].max() - P[tr, vin].min())
    phi_bbox = float(A / (Ls * ext))
    ps_area_below = float(area[pstar_t < 0].sum() / A); ps_area_above = float(area[pstar_t > 1].sum() / A)
    ps_far = float(np.mean((pstar_v[vin] < -0.05) | (pstar_v[vin] > 1.05)))
    np.savez_compressed(D / f"{cc}_tri.npz", ut_over_U=ut / U, area=area, s_star=s / Ls, pstar=pstar_t,
                        umax_over_U=umax / U)
    return dict(case=case, kD=meta["kD"], Q=meta["Q"], U=U, Ls=Ls, W=W, n_tri_sample=int(ins.sum()),
                A_f=float(A), porosity_full_window=float(A / (Ls * W)), porosity=phi_bbox, transverse_extent=ext,
                pstar_area_below0=ps_area_below, pstar_area_above1=ps_area_above, pstar_vertex_beyond_5pct=ps_far, wall_length=Pw, D_h=float(Dh),
                umax_over_U=umax / U, active=float(act), dead=float(dead), intermediate=float(1 - act - dead),
                A50=A50, A90=A90, pstar_vertex_below0=frac_lo, pstar_vertex_above1=frac_hi,
                pstar_vertex_outside=frac_lo + frac_hi, umax_position_cropped_frame=xmax.tolist(),
                kD_times_Dh_px3=float(meta["kD"] * Dh))


def main():
    out = {c: run(c) for c in CASES}
    out["A_k_R3"] = out["R3-V"]["kD"] / out["R3-H"]["kD"]
    json.dump(out, open(D / "p1_descriptors.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
