#!/usr/bin/env python3
"""B3 - end-to-end check of protocol P1 on a pixelised cylinder array (image -> author's mask/mesh functions -> reservoirs -> reference solver -> k_D).
Reference: periodic unit cell solved on the analytic geometry (benchmarks/b2, independent of pixels and of the mesher), evaluated at the actual
pixel solid fraction. End effects are removed by the increment between arrays with N and N+2 cells:
    k_inc = (L_{N+2} - L_N) / ( W * ( dP_{N+2}/Q_{N+2} - dP_N/Q_N ) )     [dP: face pressure drop; mu = 1]
Pass criterion (declared before the run): |k_inc / k_cell - 1| < 2% with free-slip lateral walls.
USAGE: b3_array_chain.py [--level 3]
OUTPUT: results/b3_array_chain.json
"""
import argparse, json, sys, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import periodic_cell as pc, reservoir, synthetic
from drc_audit.legacy import load_pipeline

ap = argparse.ArgumentParser(); ap.add_argument("--level", type=int, default=3); a_ = ap.parse_args()
m = load_pipeline()
a, M = 160, 2
tmp = Path(tempfile.mkdtemp())
res = []
for c in (0.10, 0.20):
    runs = {}
    for lat in ("freeslip", "noslip"):
        for N in (4, 6):
            mask, R = synthetic.cylinder_array_mask(a, c, M, N)
            rec = reservoir.run_mask(mask, "lr", tmp / f"c{c}_{lat}_{N}", level=a_.level, tol=0.25, pad=0.10, elements=("p2",), bc="pressure", lateral=lat, m=m)[0]
            runs[(lat, N)] = dict(rec=rec, c_pix=float(1 - mask.mean()), mask_shape=mask.shape)
    c_pix = runs[("freeslip", 4)]["c_pix"]
    pc.make_mesh(c_pix, 0.01, tmp / "cell.msh"); kcell = pc.solve_cell(tmp / "cell.msh")["k"] * a ** 2          # px^2
    row = dict(c_nominal=c, c_pixel=c_pix, k_cell_px2=kcell, k_hasimoto_px2=float(pc.hasimoto_k(c_pix)) * a ** 2, level=a_.level)
    for lat in ("freeslip", "noslip"):
        r4, r6 = runs[(lat, 4)]["rec"], runs[(lat, 6)]["rec"]
        W = r4["W"]
        # resistance dP/Q from kD_face = Q Ls / (W dP)  ->  dP/Q = Ls / (W kD_face)
        R4, R6 = r4["Ls"] / (W * r4["kD_face"]), r6["Ls"] / (W * r6["kD_face"])
        kinc = (r6["Ls"] - r4["Ls"]) / (W * (R6 - R4))
        row[lat] = dict(k_D_N4=r4["kD_face"], k_D_N6=r6["kD_face"], k_inc=kinc, rel_to_cell=kinc / kcell - 1, n_tri=[r4["n_tri"], r6["n_tri"]],
                        residual=max(r4["residual_lu"], r6["residual_lu"]))
    row["PASS_freeslip_2pct"] = bool(abs(row["freeslip"]["rel_to_cell"]) < 0.02)
    res.append(row)
    print(json.dumps(row, indent=1))
json.dump(res, open(ROOT / "results" / "b3_array_chain.json", "w"), indent=1)
