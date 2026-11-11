#!/usr/bin/env python3
"""B2 - literature benchmark (independent of the digital-rock mesher): periodic square array of cylinders.
Pass criteria (declared before looking at the result):
  (i)  mesh convergence: P2/P1 k changes < 1% between the two finest meshes;
  (ii) agreement with the Hasimoto (1959) series within 0.5% for solid fraction c <= 0.15;
  (iii) the remainder (series - numerics) scales as c^4 for c >= 0.2 (the next term of the series), i.e. remainder/c^4 constant within 10%.
OUTPUT: results/b2_periodic_cell.json
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drc_audit import periodic_cell as pc

out = {"cases": [], "criteria": {}}
tmp = Path(tempfile.mkdtemp()) / "cell.msh"
for c in (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5):
    ks = []
    for h in (0.08, 0.04, 0.02, 0.01):
        pc.make_mesh(c, h, tmp); r = pc.solve_cell(tmp); ks.append((h, r["n_tri"], r["k"]))
    kref = float(pc.hasimoto_k(c))
    br = lambda k: 8 * np.pi * k                       # numeric bracket
    poly = -np.log(c) - 1.476 + 2 * c - 1.774 * c ** 2 + 4.076 * c ** 3
    out["cases"].append(dict(c=c, mesh=[dict(h=h, n_tri=n, k=k) for h, n, k in ks], k_num=ks[-1][2], k_hasimoto=kref,
                             rel_diff=float(ks[-1][2] / kref - 1), conv_last_step=float(ks[-1][2] / ks[-2][2] - 1), remainder_over_c4=float((br(ks[-1][2]) - poly) / c ** 4)))
cs = out["cases"]
out["criteria"]["conv_max_last_step"] = float(max(abs(x["conv_last_step"]) for x in cs))
out["criteria"]["dilute_max_rel_diff(c<=0.15)"] = float(max(abs(x["rel_diff"]) for x in cs if x["c"] <= 0.15))
rem = [x["remainder_over_c4"] for x in cs if x["c"] >= 0.2]
out["criteria"]["remainder_over_c4_range(c>=0.2)"] = [float(min(rem)), float(max(rem))]
out["criteria"]["PASS_i"] = bool(out["criteria"]["conv_max_last_step"] < 0.01)
out["criteria"]["PASS_ii"] = bool(out["criteria"]["dilute_max_rel_diff(c<=0.15)"] < 0.005)
out["criteria"]["PASS_iii"] = bool((max(rem) - min(rem)) / abs(np.mean(rem)) < 0.10)
Path(__file__).resolve().parents[1].joinpath("results").mkdir(exist_ok=True)
json.dump(out, open(Path(__file__).resolve().parents[1] / "results" / "b2_periodic_cell.json", "w"), indent=1)
for x in cs: print(x["c"], x["k_num"], x["k_hasimoto"], round(x["rel_diff"] * 100, 2), "%", round(x["conv_last_step"] * 100, 3), "%", round(x["remainder_over_c4"], 2))
print(out["criteria"])
