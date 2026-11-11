#!/usr/bin/env python3
"""B9 - Poiseuille channels under the exact protocol P1 (npj revision, 2026-11-03).

Why: Figure 8a of the Nature Communications draft took its "pressure-driven" channel data from evidence/poiseuille,
which re-solved the ORIGINAL pipeline meshes (contour tolerance 1.5 px, no reservoirs, nodal-mean pressures on the
boundary). Douglas-Peucker at 1.5 px chamfers the channel ends (area deficit 1.3% for h = 20), so those deviations are
a geometry defect, not a property of the P1 definition. This script repeats the channels with the P1 protocol itself.

Geometry: fully open channel of width h (W = h) and length L_s (L_s/h = 20 for h = 20, 40 and 10 for h = 80, 160, as in
the draft), lateral edges are no-slip walls, open reservoirs of width b = 0.1 L_s, contour tolerance 0.25 px.

Exact value of the P1 definition: the probe planes sit 1 px inside each reservoir, so they are L_s + 2 px apart while
k_D divides by L_s. For a straight channel the reservoirs are channel too, hence
    k_D,exact = (h^2/12) * L_s / (L_s + 2)          (face definition, flux- or length-weighted)
    k_D,outer,exact = (h^2/12) * L_s / (L_s + 2 b)   (outer reservoir faces)
Plane Poiseuille flow (quadratic velocity, linear pressure) lies in the Taylor-Hood P2/P1 space and satisfies the
do-nothing outflow condition of the Laplacian form exactly, so P2/P1 must reproduce it to round-off on any mesh.

Declared acceptance criterion (written before running): |k_D,face(P2/P1) / k_D,exact - 1| < 1e-6 at every level and
width. Mini (P1b/P1) cannot represent the parabola: its error and observed order are reported, not tested.
Also reported: the velocity-driven variant without reservoirs (uniform U0 on the inlet), k from outer-face pressures.

OUTPUT: results/b9_poiseuille_p1.json
"""
import json, sys, tempfile, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import reservoir                                                    # noqa: E402
from drc_audit.legacy import load_pipeline                                        # noqa: E402

m = load_pipeline()
tmp = Path(tempfile.mkdtemp(prefix="b9_"))
rows = []
for h, ratio in ((20, 20), (40, 20), (80, 10), (160, 10)):
    Ls = ratio * h
    mask = np.ones((h, Ls), bool)
    k0 = h * h / 12.0
    for level in (1, 2, 3, 4, 5):
        row = dict(h=h, Ls=Ls, level=level, k_poiseuille=k0, k_exact_face=k0 * Ls / (Ls + 2.0))
        for el in ("p2", "mini"):
            t0 = time.time()
            r = reservoir.run_mask(mask, "lr", tmp / f"p_{h}_{level}_{el}", level=level, tol=0.25, pad=0.10,
                                   elements=(el,), bc="pressure", lateral="noslip", m=m)[0]
            b = r["pad_px"]
            row[f"P1_{el}"] = dict(n_tri=r["n_tri"], n_dof=r["n_dof"], residual_lu=r["residual_lu"], Q_balance=r["Q_balance"],
                                   kD_face=r["kD_face"], kD_face_len=r["kD_face_len"], kD_outer=r["kD_outer"],
                                   rel_face=r["kD_face"] / row["k_exact_face"] - 1,
                                   rel_face_len=r["kD_face_len"] / row["k_exact_face"] - 1,
                                   rel_outer=r["kD_outer"] / (k0 * Ls / (Ls + 2.0 * b)) - 1,
                                   rel_vs_poiseuille=r["kD_face"] / k0 - 1, pad_px=b, wall_s=time.time() - t0)
            r = reservoir.run_mask(mask, "lr", tmp / f"v_{h}_{level}_{el}", level=level, tol=0.25, pad=0.0,
                                   elements=(el,), bc="velocity", lateral="noslip", m=m)[0]
            row[f"velocity_{el}"] = dict(n_tri=r["n_tri"], kD_outer=r["kD_outer"], rel_vs_poiseuille=r["kD_outer"] / k0 - 1,
                                         Q_balance=r["Q_balance"])
        rows.append(row)
        print(json.dumps(dict(h=h, level=level, p2=row["P1_p2"]["rel_face"], mini=row["P1_mini"]["rel_face"],
                              vel_p2=row["velocity_p2"]["rel_vs_poiseuille"]), default=float), flush=True)

worst = max(abs(r["P1_p2"]["rel_face"]) for r in rows)
criteria = dict(declared="|kD_face(P2/P1)/kD_exact - 1| < 1e-6 at every level and width", worst_p2=worst, PASS=bool(worst < 1e-6))
json.dump(dict(protocol="P1: pressure-driven, reservoirs b = 0.1 Ls, tol 0.25 px, face planes 1 px inside reservoirs",
               rows=rows, criteria=criteria), open(ROOT / "results" / "b9_poiseuille_p1.json", "w"), indent=1, default=float)
print(json.dumps(criteria, indent=1, default=float))
