#!/usr/bin/env python3
"""Run BOTH solvers on one case and print a side-by-side benchmark table.
  (1) legacy pipeline, run unchanged as its own CLI with the manuscript's parameters (MINRES, Mini, uniform inlet velocity, tol 1.5 px, band 8%);
  (2) reference solver on the SAME Gmsh mesh with the SAME boundary conditions (verification of the solver implementation: expect agreement of the
      legacy-definition permeability within the MINRES tolerance);
  (3) protocol P1 (pressure-driven, reservoirs, P2/P1, tolerance 0.25 px) on the same image (the recommended definition).
USAGE: compare_solvers.py --case R3-H --level 1 [--structures data/structures] [--out results/compare]
"""
import argparse, csv, json, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import reservoir, stokes
from drc_audit.legacy import load_pipeline, pipeline_path

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True, choices=list(reservoir.CASES)); ap.add_argument("--level", type=int, default=1)
ap.add_argument("--structures", default=str(ROOT / "data" / "structures")); ap.add_argument("--out", default=str(ROOT / "results" / "compare"))
a = ap.parse_args()
stem, direction = reservoir.CASES[a.case]
out = Path(a.out) / f"{a.case.replace('-', '')}_L{a.level}"; out.mkdir(parents=True, exist_ok=True)
one = out / "structure"; one.mkdir(exist_ok=True); shutil.copy(Path(a.structures) / f"{stem}_c0_mask.png", one)
hmin, hmax, dmin, dmax = reservoir.LEVELS[a.level]
cmd = [sys.executable, str(pipeline_path()), "--structures", str(one), "--direction", direction, "--boundary-band-fraction", "0.08", "--closing-radius", "0", "--keep-percolating", "all",
       "--min-object-size", "20", "--fill-holes-smaller-than", "0", "--contour-tolerance", "1.5", "--min-contour-area", "20", "--mesh-algorithm", "6", "--smoothing", "5",
       "--pixel-size", "1.0", "--unit-name", "pixel", "--thickness", "1.0", "--viscosity", "1.0", "--density", "1.0", "--inlet-velocity", "0.0001",
       "--active-threshold-relative", "0.01", "--dead-threshold-relative", "0.001", "--solver-tolerance", "1e-7", "--h-min", str(hmin), "--h-max", str(hmax),
       "--distance-min", str(dmin), "--distance-max", str(dmax), "--solver-max-iterations", "10000", "--interpolation-batch-size", "1024",
       "--max-total-dofs", "2000000", "--output", str(out / "legacy")]
subprocess.run(cmd, check=True, capture_output=True)
row = list(csv.DictReader(open(out / "legacy" / "global_summary.csv")))[0]
assert row["status"] == "success", row
k_leg = float(row["apparent_permeability"]); mesh = row["mesh_path"]
ref = stokes.solve(mesh, direction, "mini", bc="velocity", lateral="noslip")
k_ref_legdef = 1e-4 * ref["flow_len"] / ref["dP_nodal"]                  # mu U0 L / dP_nodal, the legacy definition (Q_nom / L_in = U0)
p1 = reservoir.run_case(a.case, a.structures, out / "p1", level=a.level, tol=0.25, band=0.08, pad=0.10, elements=("p2",), bc="pressure", lateral="noslip")[0]
table = dict(case=a.case, level=a.level,
             legacy_solver_k_legacy_definition=k_leg, legacy_solver_residual=float(row["solver_relative_residual"]),
             reference_solver_same_mesh_k_legacy_definition=k_ref_legdef, rel_diff_reference_vs_legacy=k_ref_legdef / k_leg - 1,
             reference_solver_same_mesh_Q_actual_over_nominal=ref["Q"] / (1e-4 * ref["inlet"]["L"]),
             protocol_P1_kD_face_px2=p1["kD_face"], protocol_P1_kD_outer_px2=p1["kD_outer"], ratio_P1_over_legacy=p1["kD_face"] / k_leg)
json.dump(table, open(out / "compare.json", "w"), indent=1); print(json.dumps(table, indent=1))
