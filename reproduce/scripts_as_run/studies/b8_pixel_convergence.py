#!/usr/bin/env python3
"""B8 - is the FV-vs-FEM gap of the b7 pairs a discretisation error of the pixel mesh, or a disagreement about
the physics? (finding F-008, follow-up to b7_subwindow_pair.py)

On one sub-window the two solvers differ by up to 32%. Four explanations were tested on the b7 results and all
four were refuted:
  1. wall resolution (+2/N^2, the law B6-S measured on an open channel): Pearson r of -0.11..+0.28 against three
     aperture measures, and the predicted magnitudes are ~0.05% against an observed median of -0.61%;
  2. flow paths lost to face connectivity (a diagonal pinch carries nothing in a finite-volume mesh): restricting
     both sides to face-connected components changes k_FEM by 1.3e-7 - it fixed two GAMG crashes but moved no number;
  3. solid inclusions below cfg.min_contour_area_px2 = 20 px^2, which the FEM geometry drops and the pixel mesh
     keeps: the dropped area is 0.004% to 0.036% of the pore area, orders of magnitude too small;
  4. the contour simplification opening passages the pixel mask does not have: re-solving the FEM at contour
     tolerance 0.1 and 0.0 px moves k_FEM by at most 1e-5 relative.
What remains is the discretisation of the pixel mesh itself. This script tests exactly that, by SUPERSAMPLING:
each pore pixel is split into s x s finite-volume cells. The staircase boundary stays byte-identical - the walls
still sit on the original pixel edges - so the geometry is unchanged and only the cell size shrinks, by 1/s.

Declared expectation (written before running): if the gap is pixel-mesh discretisation error, k_FV(s) must
converge monotonically towards k_FEM as s grows, at the second-order rate the scheme has, and the Richardson
extrapolation of k_FV(s) to s -> infinity must land within a per cent or so of k_FEM. If instead k_FV(s)
converges to a value well away from k_FEM, the two codes disagree about the same geometry and the discrepancy is
a finding against one of them, which would have to be chased further before F-008 can be closed on rock windows.

ENV: B8_CASE, B8_TILE, B8_TAG (e.g. t1_5), B8_SUPER (default "1,2,3"), B8_NPROC, B8_WORK, B8_LEVEL (FEM, default 5)
OUTPUT: results/b8_pixel_convergence_<case>_<tag>.json
"""
import importlib.util, json, os, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("b7", Path(__file__).resolve().parent / "b7_subwindow_pair.py")
b7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b7)
from drc_audit import legacy, reservoir                                            # noqa: E402

CASE = os.environ.get("B8_CASE", "R2-V")
TILE = int(os.environ.get("B8_TILE", "512"))
TAG = os.environ["B8_TAG"]
SUPER = [int(x) for x in os.environ.get("B8_SUPER", "1,2,3").split(",")]
NPROC = int(os.environ.get("B8_NPROC", "20"))
LEVEL = int(os.environ.get("B8_LEVEL", "5"))
BAND, PAD, TOL = 0.08, 0.10, 0.25
STRUCT = ROOT / "data/structures"


def main():
    work = Path(os.environ["B8_WORK"])
    work.mkdir(parents=True, exist_ok=True)
    m = legacy.load_pipeline()
    stem, direction = reservoir.CASES[CASE]
    cfg = reservoir.build_cfg(m, STRUCT, work / "_cfg", level=LEVEL, tol=TOL, band=BAND)
    sel, _ = reservoir.select_window(m, cfg, STRUCT / f"{stem}_c0_mask.png", direction)
    i_f, i_t = (int(x) for x in TAG.lstrip("t").split("_"))
    sl_f, sl_t = slice(i_f * TILE, (i_f + 1) * TILE), slice(i_t * TILE, (i_t + 1) * TILE)
    sub = sel[sl_f, sl_t] if direction == "tb" else sel[sl_t, sl_f]
    tsel, _ = m.select_full_percolating_components(sub.copy(), direction, cfg)
    tsel4, dropped = b7.face_connected(tsel, direction)
    print(f"[{CASE} {TAG}] mask {tsel4.shape} dir={direction} pore_px={int(tsel4.sum()):_} "
          f"porosity={tsel4.mean():.4f} dropped_by_face_conn={dropped}", flush=True)

    out = {"case": CASE, "tag": TAG, "tile_px": TILE, "direction": direction, "supersampling": SUPER,
           "fem_level": LEVEL, "contour_tol_px": TOL, "pore_px": int(tsel4.sum()), "fv": []}

    r = reservoir.run_mask(tsel4, direction, work / "fem", level=LEVEL, tol=TOL, pad=PAD, elements=("p2",),
                           bc="pressure", lateral="noslip", m=m, structures_dir=STRUCT)[0]
    out["fem"] = r
    k_fem = r["kD_face"]
    print(f"  FEM P2/P1 level {LEVEL}: k = {k_fem:.4f} px^2  (n_dof {r['n_dof']:_}, LU residual {r['residual_lu']:.1e})",
          flush=True)

    for s in SUPER:
        up = np.repeat(np.repeat(tsel4, s, axis=0), s, axis=1)      # identical staircase, cells s times smaller
        padded, Ls, b = reservoir.pad_reservoirs(up, direction, PAD)
        W = up.shape[0] if direction == "lr" else up.shape[1]
        t0 = time.time()
        rec = b7.solve_fv(padded, direction, Ls, b, W, work / f"fv_s{s}", NPROC)
        rec["supersampling"] = s
        rec["cell_px"] = 1.0 / s
        rec["kD_face_px2"] = rec["kD_face"] / s ** 2                # k is an area: cell units -> px^2
        rec["rel_vs_fem"] = rec["kD_face_px2"] / k_fem - 1
        rec["wall_s"] = time.time() - t0
        out["fv"].append(rec)
        print(f"  FV  s={s} (cell = {1 / s:.3f} px): k = {rec['kD_face_px2']:10.4f} px^2  "
              f"{100 * rec['rel_vs_fem']:+7.2f}% vs FEM   cells={rec['n_cells']:_}  it={rec['iterations']}  "
              f"ranks={rec['nproc']}  {rec['wall_s']:.0f} s  res_ok={rec['res_ok']}", flush=True)

    ks = [x["kD_face_px2"] for x in out["fv"]]
    if len(ks) >= 3:
        # Richardson on the last three, assuming the observed order p from the three values
        k1, k2, k3 = ks[-3:]
        s1, s2, s3 = (x["supersampling"] for x in out["fv"][-3:])
        denom = (k2 - k1)
        ratio = (k3 - k2) / denom if denom else float("nan")
        out["observed_order"] = float(-np.log(abs(ratio)) / np.log(s3 / s2)) if ratio and abs(ratio) > 0 else None
        out["richardson_limit"] = float(k3 + (k3 - k2) * ratio / (1 - ratio)) if ratio and abs(ratio) < 1 else None
        if out["richardson_limit"] is not None:
            out["richardson_vs_fem"] = out["richardson_limit"] / k_fem - 1
            print(f"\n  observed order {out['observed_order']:.2f}; Richardson limit "
                  f"{out['richardson_limit']:.4f} px^2 = {100 * out['richardson_vs_fem']:+.2f}% vs FEM {k_fem:.4f}")
        else:
            print(f"\n  Richardson not applicable (successive differences do not shrink: ratio {ratio})")
    p = ROOT / f"results/b8_pixel_convergence_{CASE}_{TAG}.json"
    json.dump(out, open(p, "w"), indent=1)
    print("->", p)


if __name__ == "__main__":
    main()
