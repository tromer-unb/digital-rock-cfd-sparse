#!/usr/bin/env python3
"""Lateral-boundary study on the manuscript's windows: free-slip (symmetry) versus no-slip lateral walls, window (crop band) scan, protocol P1
(level 2, P2/P1, tolerance 0.25 px, reservoirs 10%).  Free-slip on the lateral edges is the mirror-symmetric (equivalently, mirrored-periodic)
condition; a genuine periodic condition cannot be imposed on these windows because the pore openings on opposite edges do not match.
USAGE: run_bc_study.py --structures data/structures --out results/bc_study [--laterals freeslip noslip]
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from drc_audit import reservoir
from drc_audit.legacy import load_pipeline

ap = argparse.ArgumentParser()
ap.add_argument("--structures", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--laterals", nargs="+", default=["freeslip"]); ap.add_argument("--level", type=int, default=2)
ap.add_argument("--bands", nargs="+", type=float, default=[0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.16])
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
m = load_pipeline()
for lat in a.laterals:
    for case in reservoir.CASES:
        for band in a.bands:
            tag = f"{case.replace('-', '')}_{lat}_L{a.level}_band{band:g}"
            f = out / f"{tag}.json"
            if f.exists(): continue
            try:
                rec = reservoir.run_case(case, a.structures, out / "mesh" / tag, level=a.level, tol=0.25, band=band, pad=0.10, bc="pressure", lateral=lat, m=m)[0]
            except RuntimeError as e:
                rec = dict(status="not_percolating", case=case, band=band, error=str(e))
            rec["tag"] = tag
            json.dump(rec, open(f, "w"), indent=1)
            print(tag, rec.get("kD_face"), flush=True)
