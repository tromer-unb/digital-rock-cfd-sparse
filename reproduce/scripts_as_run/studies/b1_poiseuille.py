#!/usr/bin/env python3
"""B1 - Poiseuille channel (analytic): k_D = h^2/12 for a straight channel of width h (W = h) with no-slip walls.
Pass criterion (declared before the run): |k_D / (h^2/12) - 1| < 0.5% for P2/P1 with protocol P1 (pressure-driven, reservoirs 10%).
Also reported (information only): Mini P1b/P1, and the uniform-velocity inlet WITHOUT reservoirs (legacy-type inlet, k from outer-face pressure).
OUTPUT: results/b1_poiseuille.json
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import reservoir
from drc_audit.legacy import load_pipeline

m = load_pipeline(); tmp = Path(tempfile.mkdtemp()); rows = []
for h in (100, 200):
    mask = np.ones((h, 5 * h), bool)                      # fully open channel, lateral edges are the walls
    ref = h ** 2 / 12
    row = dict(h=h, k_ref=ref)
    for el in ("p2", "mini"):
        r = reservoir.run_mask(mask, "lr", tmp / f"p1_{h}_{el}", level=3, tol=0.25, pad=0.10, elements=(el,), bc="pressure", lateral="noslip", m=m)[0]
        row[f"P1_{el}"] = dict(k_face=r["kD_face"], rel=r["kD_face"] / ref - 1, k_outer=r["kD_outer"], rel_outer=r["kD_outer"] / ref - 1)
    r = reservoir.run_mask(mask, "lr", tmp / f"v_{h}", level=3, tol=0.25, pad=0.0, elements=("p2",), bc="velocity", lateral="noslip", m=m)[0]
    row["velocity_inlet_noreservoir_p2"] = dict(k_outer=r["kD_outer"], rel_outer=r["kD_outer"] / ref - 1)
    row["PASS_P1_p2_0.5pct"] = bool(abs(row["P1_p2"]["rel"]) < 0.005)
    rows.append(row); print(json.dumps(row, indent=1))
json.dump(rows, open(ROOT / "results" / "b1_poiseuille.json", "w"), indent=1)
