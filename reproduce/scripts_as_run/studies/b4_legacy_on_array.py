#!/usr/bin/env python3
"""B4 - the author's pipeline (legacy/, unmodified, run as its own CLI with the manuscript's parameters) on the same pixelised cylinder array,
compared with the periodic-cell reference and with protocol P1.  Shows how the definition of k in the legacy code (nominal flux over the inlet
opening, uniform inlet velocity, no-slip lateral walls, 8% band) maps onto a case whose answer is known.
USAGE: b4_legacy_on_array.py [--level 2]
OUTPUT: results/b4_legacy_on_array.json
"""
import argparse, csv, json, subprocess, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import periodic_cell as pc, reservoir, synthetic
from drc_audit.legacy import load_pipeline, pipeline_path

ap = argparse.ArgumentParser(); ap.add_argument("--level", type=int, default=2); A = ap.parse_args()
hmin, hmax, dmin, dmax = reservoir.LEVELS[A.level]
a, M, N = 160, 6, 6
tmp = Path(tempfile.mkdtemp()); m = load_pipeline(); rows = []
for c in (0.10, 0.20):
    mask, R = synthetic.cylinder_array_mask(a, c, M, N)
    sdir = tmp / f"s{c}"; sdir.mkdir(); synthetic.save_green_png(mask, sdir / "array_c0_mask.png")
    c_pix = float(1 - mask.mean())
    pc.make_mesh(c_pix, 0.01, tmp / "cell.msh"); kcell = pc.solve_cell(tmp / "cell.msh")["k"] * a ** 2
    H, L = mask.shape
    out = dict(c_pixel=c_pix, k_cell_px2=kcell, level=A.level)
    for band in (0.08, 0.10, 0.001):   # 0.001 (= minimum band, 1 px) leaves the inlet plane fully open: tests the legacy failure mode
        od = tmp / f"o{c}_{band}"
        cmd = [sys.executable, str(pipeline_path()), "--structures", str(sdir), "--direction", "lr", "--boundary-band-fraction", str(band), "--closing-radius", "0",
               "--keep-percolating", "all", "--min-object-size", "20", "--fill-holes-smaller-than", "0", "--contour-tolerance", "1.5", "--min-contour-area", "20",
               "--mesh-algorithm", "6", "--smoothing", "5", "--pixel-size", "1.0", "--unit-name", "pixel", "--thickness", "1.0", "--viscosity", "1.0", "--density", "1.0",
               "--inlet-velocity", "0.0001", "--active-threshold-relative", "0.01", "--dead-threshold-relative", "0.001", "--solver-tolerance", "1e-7",
               "--h-min", str(hmin), "--h-max", str(hmax), "--distance-min", str(dmin), "--distance-max", str(dmax), "--solver-max-iterations", "10000",
               "--interpolation-batch-size", "1024", "--max-total-dofs", "2000000", "--output", str(od)]
        subprocess.run(cmd, check=True, capture_output=True, timeout=3000)
        r = list(csv.DictReader(open(od / "global_summary.csv")))[0]
        if "apparent_permeability" not in r:
            out[f"legacy_band{band}"] = dict(status=r.get("status"), error=r.get("error"), exit_code_of_cli=0); continue
        k_auth = float(r["apparent_permeability"])
        open_frac = float(mask[:, int(r["x_offset_px"])].mean())        # pore fraction on the cropped inlet plane = L_in / W
        out[f"legacy_band{band}"] = dict(k_author=k_auth, inlet_open_fraction=open_frac, ratio_to_cell=k_auth / kcell, k_author_times_open_fraction=k_auth * open_frac,
                                         ratio_darcy_to_cell=k_auth * open_frac / kcell, residual=float(r["solver_relative_residual"]))
    rows.append(out); print(json.dumps(out, indent=1))
json.dump(rows, open(ROOT / "results" / "b4_legacy_on_array.json", "w"), indent=1)
