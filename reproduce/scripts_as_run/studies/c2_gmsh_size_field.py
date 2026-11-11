#!/usr/bin/env python3
"""C2 - which size law does the pipeline's Gmsh Threshold field (Sigmoid = 1) impose? (npj revision, 2026-11-03)

The pipeline sets Distance -> Threshold(SizeMin = h_min, SizeMax = h_max, DistMin = d_min, DistMax = d_max, Sigmoid = 1)
as background mesh. Two candidate laws, with r = clip((d - d_min)/(d_max - d_min), 0, 1):
  sigmoid: h = h_min (1 - s) + h_max s,  s = exp(12 r - 6) / (1 + exp(12 r - 6))
  linear : h = h_min + (h_max - h_min) r
Test: mesh a 4000 x 1000 px rectangle with the field measured from the bottom edge only, d_min = 50, d_max = 450,
h_min = 4, h_max = 40, and compare the RMS edge length of the triangles in 20-px strips of d with both laws.
Declared criterion (before running): the sigmoid law must fit the strip means with a smaller RMS relative misfit than the
linear law over d_min < d < d_max; otherwise the paper describes the transition without a formula.
OUTPUT: results/c2_gmsh_size_field.json
"""
import json
from pathlib import Path
import numpy as np
import gmsh

ROOT = Path(__file__).resolve().parents[1]
HMIN, HMAX, DMIN, DMAX, LX, LY = 4.0, 40.0, 50.0, 450.0, 4000.0, 1000.0
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("c2")
p = [gmsh.model.geo.addPoint(x, y, 0, HMAX) for x, y in ((0, 0), (LX, 0), (LX, LY), (0, LY))]
c = [gmsh.model.geo.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(c)]); gmsh.model.geo.synchronize()
f = gmsh.model.mesh.field.add("Distance"); gmsh.model.mesh.field.setNumbers(f, "CurvesList", [c[0]])
gmsh.model.mesh.field.setNumber(f, "Sampling", 400)
t = gmsh.model.mesh.field.add("Threshold")
for k, v in (("InField", f), ("SizeMin", HMIN), ("SizeMax", HMAX), ("DistMin", DMIN), ("DistMax", DMAX), ("Sigmoid", 1)):
    gmsh.model.mesh.field.setNumber(t, k, v)
gmsh.model.mesh.field.setAsBackgroundMesh(t)
for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
             ("Mesh.Algorithm", 6), ("Mesh.Smoothing", 5)):
    gmsh.option.setNumber(k, v)
gmsh.model.mesh.generate(2)
tags, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3); idx = {int(g): i for i, g in enumerate(tags)}
et, _, en = gmsh.model.mesh.getElements(2); tri = np.array([idx[int(g)] for g in en[0]]).reshape(-1, 3)
gmsh.finalize()
P = xyz[tri][:, :, :2]
edge = np.sqrt(np.mean([np.sum((P[:, i] - P[:, (i + 1) % 3]) ** 2, axis=1) for i in range(3)], axis=0))
yc = P[:, :, 1].mean(axis=1); xc = P[:, :, 0].mean(axis=1)
keep = (xc > 200) & (xc < LX - 200) & (yc < LY - 100)          # away from the other three edges
def sig(d):
    r = np.clip((d - DMIN) / (DMAX - DMIN), 0, 1); s = np.exp(12 * r - 6) / (1 + np.exp(12 * r - 6)); return HMIN * (1 - s) + HMAX * s
def lin(d):
    r = np.clip((d - DMIN) / (DMAX - DMIN), 0, 1); return HMIN + (HMAX - HMIN) * r
bins = np.arange(0, 600, 20.0); rows = []
for a, b in zip(bins[:-1], bins[1:]):
    s = keep & (yc >= a) & (yc < b)
    if s.sum() < 20:
        continue
    d = 0.5 * (a + b); rows.append(dict(d=d, n=int(s.sum()), h_meas=float(np.mean(edge[s])), h_sig=float(sig(d)), h_lin=float(lin(d))))
tr = [r for r in rows if DMIN < r["d"] < DMAX]
mis = {k: float(np.sqrt(np.mean([(r["h_meas"] / r[k] - 1) ** 2 for r in tr]))) for k in ("h_sig", "h_lin")}
# a measured mean edge differs from the target size by a constant factor; fit that factor before comparing shapes
fac = {k: float(np.median([r["h_meas"] / r[k] for r in rows if r["d"] < DMIN or r["d"] > DMAX])) for k in ("h_sig", "h_lin")}
mis_scaled = {k: float(np.sqrt(np.mean([(r["h_meas"] / (fac[k] * r[k]) - 1) ** 2 for r in tr]))) for k in ("h_sig", "h_lin")}
res = dict(params=dict(h_min=HMIN, h_max=HMAX, d_min=DMIN, d_max=DMAX), n_tri=int(len(tri)), rows=rows, rms_misfit=mis,
           plateau_factor=fac, rms_misfit_after_plateau_factor=mis_scaled,
           PASS_sigmoid_better=bool(mis_scaled["h_sig"] < mis_scaled["h_lin"]))
json.dump(res, open(ROOT / "results" / "c2_gmsh_size_field.json", "w"), indent=1)
for r in rows:
    print(f"d={r['d']:6.1f} n={r['n']:5d} meas={r['h_meas']:6.2f} sig={r['h_sig']:6.2f} lin={r['h_lin']:6.2f}")
print(json.dumps({k: res[k] for k in ("rms_misfit", "plateau_factor", "rms_misfit_after_plateau_factor", "PASS_sigmoid_better")}, indent=1))
