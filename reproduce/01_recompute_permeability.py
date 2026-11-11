#!/usr/bin/env python3
"""Step 1 of the reproduction: recompute the permeabilities of the paper with PercoFlow and compare them with the archived
run outputs in data/results/.

What is recomputed: the Darcy permeability k_D of the three percolating configurations (R3-H, R2-V, R3-V) on mesh levels
1-5 with Taylor-Hood (p2) and Mini elements, i.e. Table 1 and Figure 4 of the paper, and optionally the cropping-window scan
of Figure 5a. Options are the reference definition of the paper (pressure-driven flow between open reservoirs b = 0.1 L_s,
contour tolerance 0.25 px, boundary band 8 %, flux-weighted pressures, full width W), which are the PercoFlow defaults.

Acceptance (declared before running): every recomputed k_D agrees with the archived value to a relative 1e-6
(same images, same options, same mesher and solver logic; differences only from floating-point summation order).

Usage (inside the PercoFlow environment, see docs/TUTORIAL.md):
    python reproduce/01_recompute_permeability.py --quick       # mesh levels 1-2, about 5 min on one core
    python reproduce/01_recompute_permeability.py               # levels 1-5 (Table 1, Figure 4), 30-40 min, ~15 GB RAM
    python reproduce/01_recompute_permeability.py --windows     # + crop-band scan (Figure 5a, A_k range), +10 min
The OpenFOAM cross-check is switched off here because these runs reproduce the finite-element numbers; the finite-volume
numbers of the paper are reproduced by `percoflow benchmark` and by the scripts in reproduce/scripts_as_run/studies/.
Writes <out>/comparison.json (default out: runs/01_permeability) and one PercoFlow run directory per case, each with its
manifest.json; exit status 0 on PASS.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGES = ROOT / "data" / "images"
ARCH = ROOT / "data" / "results" / "fix_round"
OUT = ROOT / "runs" / "01_permeability"
# configuration -> (image, flow direction). Mapping of samples R1-R3 to the files: data/README.md.
CASES = {"R3-H": ("patch_y7600_x53200_c0_mask.png", "lr"), "R2-V": ("patch_y7600_x19000_c0_mask.png", "tb"),
         "R3-V": ("patch_y7600_x53200_c0_mask.png", "tb")}


def run(case, level, element, band=0.08):
    from percoflow.config import RunConfig
    from percoflow.pipeline import run_image
    img, d = CASES[case]
    cfg = RunConfig.from_dict({"domain": {"direction": d, "band_fraction": band}, "mesh": {"level": level},
                               "solver": {"element": element}, "openfoam": {"enabled": False},
                               "output": {"save_mesh": False, "save_fields": False, "save_vtk": False, "figures": False}})
    r = run_image(IMAGES / img, OUT / f"{case}_L{level}_{element}_band{band}", cfg)
    return r


def archived(case, level, element, band=0.08):
    tag = case.replace("-", "")
    if band == 0.08:
        f = ARCH / "FR-E1_convergence_tol0.25" / f"{tag}_L{level}_pad0.1_tol0.25_band0.08_{element}.json"
    else:
        f = ARCH / "FR-E2_window_scan_tol0.25" / f"{tag}_L{level}_pad0.1_tol0.25_band{band:g}_{element}.json"
    return json.loads(f.read_text()) if f.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--windows", action="store_true")
    ap.add_argument("--out", help="output directory (default runs/01_permeability)")
    a = ap.parse_args()
    global OUT
    if a.out:
        OUT = Path(a.out).resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    levels = (1, 2) if a.quick else (1, 2, 3, 4, 5)
    rows = []
    for case in CASES:
        for level in levels:
            for el in ("p2", "mini"):
                r = run(case, level, el)
                ref = archived(case, level, el)
                k = r["permeability"]["k_D"]
                kref = ref["kD_iface_flux"] if ref else None
                rows.append({"case": case, "level": level, "element": el, "band": 0.08, "k_D": k, "k_D_archived": kref,
                             "rel_diff": (k / kref - 1) if kref else None, "n_tri": r["solver"]["n_triangles"],
                             "n_tri_archived": ref["n_tri"] if ref else None})
                print(f"{case} L{level} {el:4s} k_D = {k:.6f}  archived {kref:.6f}  rel {k / kref - 1:+.1e}", flush=True)
    if a.windows:
        for case in CASES:
            for band in (0.02, 0.04, 0.06, 0.07, 0.09, 0.10, 0.12, 0.14, 0.16):
                r = run(case, 2, "p2", band)
                ref = archived(case, 2, "p2", band)
                k = r.get("permeability", {}).get("k_D")
                kref = ref.get("kD_iface_flux") if ref else None
                rows.append({"case": case, "level": 2, "element": "p2", "band": band, "k_D": k, "k_D_archived": kref,
                             "rel_diff": (k / kref - 1) if (k and kref) else None, "status": r.get("status")})
                print(f"{case} band {band:.2f} k_D = {k}  archived {kref}", flush=True)
    diffs = [abs(x["rel_diff"]) for x in rows if x["rel_diff"] is not None]
    ok = bool(diffs) and max(diffs) < 1e-6
    rep = {"criterion": "|k_D / archived - 1| < 1e-6", "max_abs_rel_diff": max(diffs) if diffs else None, "pass": ok, "rows": rows}
    (OUT / "comparison.json").write_text(json.dumps(rep, indent=1))
    print("PASS" if ok else "FAIL", "max |rel diff| =", rep["max_abs_rel_diff"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
