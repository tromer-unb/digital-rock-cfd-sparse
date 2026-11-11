#!/usr/bin/env python3
"""B7 analysis: does the FV-vs-FEM deviation of the matched sub-window pairs follow the +2/N^2 law that the
B6-S self-test measured on an open channel?

B6-S established that the pixel finite-volume discretisation of a channel of N px returns k (1 + 2/N^2) exactly.
If the deviation seen on real pore space is that same wall-resolution effect, then for each sub-window
  rel_diff = k_FV/k_FEM - 1  ~=  2 / N_eff^2
for some effective aperture N_eff, and N_eff must be a pore-geometry quantity: it should correlate with the
sub-window's aperture statistics and NOT with the sub-window size. This script extracts N_eff from each pair and
compares it with three geometric apertures measured on the identical mask:
  r_c       percolation-threshold erosion radius (the bottleneck half-width of the best path, as scripts/critical_throat.py)
  a_mean    mean local aperture 2*d over the pore space (d = Euclidean distance to the nearest wall)
  a_harm    flow-weighted harmonic mean aperture, sum(d)/sum(1/(2d)) - the average that a series of resistances sees,
            which is the one a permeability should follow.
Printed, not thresholded: this is criterion (iii) of b7_subwindow_pair.py (interpretation), not a pass/fail test.
USAGE: b7_analyse.py results/b7_subwindow_pair_*.json
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np
import scipy.ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import legacy, reservoir                                            # noqa: E402

spec = importlib.util.spec_from_file_location("b7", Path(__file__).resolve().parent / "b7_subwindow_pair.py")
_b7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_b7)

STRUCT = ROOT / "data/structures"


def apertures(mask, direction):
    """r_c (percolation-threshold erosion radius) and the mean / flow-weighted-harmonic local apertures, in px."""
    d = ndi.distance_transform_edt(mask)

    def perc(thr):
        lab, _ = ndi.label(d > thr, structure=np.ones((3, 3)))
        a, b = (lab[:, 0], lab[:, -1]) if direction == "lr" else (lab[0, :], lab[-1, :])
        return bool((set(np.unique(a)) & set(np.unique(b))) - {0})

    lo, hi = 0.0, float(d.max())
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if perc(mid) else (lo, mid)
    dd = d[mask]
    dd = dd[dd > 0]
    return dict(r_c_px=lo, a_crit_px=2 * lo, a_mean_px=float(2 * dd.mean()),
                a_harm_px=float(dd.size / (1.0 / (2 * dd)).sum()))


def dropped_inclusions(mask, min_area):
    """Solid inclusions that the FEM geometry does not have.

    The authors' extract_contours discards any contour enclosing less than cfg.min_contour_area_px2 (20 px^2 in the
    audited configuration), so a solid grain smaller than that is simply absent from the meshed domain - while the
    pixel mesh, built from the mask itself, keeps every solid pixel. The FEM then sees more open space than the FV,
    which should push k_FEM above k_FV. Counted here: 8-connected solid components that lie strictly inside the
    window (so they are holes of the pore polygon, not part of the outer boundary) and are below the threshold.
    """
    solid = ~mask
    lab, n = ndi.label(solid, structure=np.ones((3, 3)))
    if n == 0:
        return dict(n_small_solid=0, area_small_solid_px=0.0, frac_small_solid=0.0,
                    n_inner_solid=0, area_inner_solid_px=0.0)
    border = set(np.unique(np.concatenate([lab[0, :], lab[-1, :], lab[:, 0], lab[:, -1]]))) - {0}
    areas = np.bincount(lab.ravel())
    inner = [i for i in range(1, n + 1) if i not in border]
    small = [i for i in inner if areas[i] < min_area]
    return dict(n_small_solid=len(small), area_small_solid_px=float(sum(areas[i] for i in small)),
                frac_small_solid=float(sum(areas[i] for i in small) / mask.sum()),
                n_inner_solid=len(inner), area_inner_solid_px=float(sum(areas[i] for i in inner)))


def wall_excess(mask, m, cfg):
    """How much more pore wall the pixel mesh carries than the FEM contour geometry.

    B6-S showed the pixel finite-volume discretisation is exact for an axis-aligned channel, so aperture alone
    cannot explain a deviation. What a staircase does carry is extra wall: a pore boundary running at an angle to
    the grid is represented by steps whose total length is up to sqrt(2) times the length of the smooth contour
    through the same pixels, and every bit of that wall is a no-slip surface that dissipates. The ratio below is
    therefore the natural predictor of a negative deviation, and it is a property of the wall orientation, not of
    the aperture.
      P_pixel   : number of pore/solid faces in the pixel mesh, 1 px each (lateral window edges excluded)
      P_contour : perimeter of the contours the FEM actually meshes (after approximate_polygon at cfg tolerance)
    """
    px = ((mask[:, :-1] & ~mask[:, 1:]).sum() + (mask[:, 1:] & ~mask[:, :-1]).sum()
          + (mask[:-1, :] & ~mask[1:, :]).sum() + (mask[1:, :] & ~mask[:-1, :]).sum())
    p_pixel = float(px)                                  # each such face is 1 px long
    p_contour = 0.0
    for c in m.extract_contours(mask, cfg):              # closed polygons, already simplified at cfg tolerance
        p_contour += float(np.linalg.norm(np.diff(c, axis=0), axis=1).sum())
    return dict(p_pixel_px=p_pixel, p_contour_px=p_contour,
                perimeter_ratio=p_pixel / p_contour if p_contour else float("nan"))


def main(path):
    res = json.load(open(path))
    case, direction, tile = res["case"], res["direction"], res["tile_px"]
    m = legacy.load_pipeline()
    cfg = reservoir.build_cfg(m, STRUCT, "/tmp/b7an", level=5, tol=res["contour_tol_px"], band=res["band"])
    stem, _ = reservoir.CASES[case]
    sel, _ = reservoir.select_window(m, cfg, STRUCT / f"{stem}_c0_mask.png", direction)

    rows = []
    for t in res["tiles"]:
        if t.get("status") != "ok" or not t.get("sound"):
            continue
        sl_f = slice(t["i_flow"] * tile, (t["i_flow"] + 1) * tile)
        sl_t = slice(t["i_trans"] * tile, (t["i_trans"] + 1) * tile)
        sub = sel[sl_f, sl_t] if direction == "tb" else sel[sl_t, sl_f]
        tsel, _ = m.select_full_percolating_components(sub.copy(), direction, cfg)
        tsel, _ = _b7.face_connected(tsel, direction)      # the mask both solvers actually saw
        ap = apertures(tsel, direction)
        ap.update(dropped_inclusions(tsel, res.get("min_contour_area_px2", 20.0)))
        ap.update(wall_excess(tsel, m, cfg))
        rd = t["rel_diff"]
        n_eff = float(np.sqrt(2.0 / rd)) if rd > 0 else None
        rows.append(dict(tag=t["tag"], phi=t["porosity_percolating"], k_fem=t["fem"]["5"]["kD_face"],
                         k_fv=t["fv"]["kD_face"], rel_diff=rd, n_eff_px=n_eff, **ap))

    print(f"\n{case}, sub-windows of {tile} px, {len(rows)} sound pairs")
    print(f"{'tag':>7} {'phi':>6} {'k_FEM':>11} {'k_FV':>11} {'rel':>8} "
          f"{'a_crit':>7} {'a_mean':>7} | {'frac_sm':>9} | {'P_pixel':>9} {'P_contour':>10} {'P ratio':>8}")
    for r in rows:
        print(f"{r['tag']:>7} {r['phi']:6.3f} {r['k_fem']:11.4f} {r['k_fv']:11.4f} "
              f"{100 * r['rel_diff']:+7.2f}% {r['a_crit_px']:7.1f} {r['a_mean_px']:7.1f} | "
              f"{100 * r['frac_small_solid']:8.4f}% | {r['p_pixel_px']:9.0f} {r['p_contour_px']:10.0f} "
              f"{r['perimeter_ratio']:8.4f}")
    if len(rows) >= 3:
        rd = np.array([r["rel_diff"] for r in rows])
        print("\nhypothesis 1 - wall resolution, the +2/N^2 law of B6-S:")
        for key in ("a_crit_px", "a_mean_px", "a_harm_px"):
            a = np.array([r[key] for r in rows])
            pred = 2.0 / a ** 2
            print(f"  {key:>10}: Pearson r(rel_diff, 2/a^2) = {np.corrcoef(rd, pred)[0, 1]:+.3f}; "
                  f"predicted median {100 * np.median(pred):+.3f}% vs observed median {100 * np.median(rd):+.3f}%")
        print("hypothesis 3 - excess no-slip wall carried by the staircase (wall orientation, not aperture):")
        pr = np.array([r["perimeter_ratio"] for r in rows])
        print(f"  {'perimeter_ratio':>20}: Pearson r(rel_diff, ratio) = {np.corrcoef(rd, pr)[0, 1]:+.3f}   "
              f"range {pr.min():.4f} to {pr.max():.4f}")
        fit = np.polyfit(pr, rd, 1)
        print(f"  {'':>20}  least squares: rel_diff = {fit[0]:+.4f} * ratio {fit[1]:+.4f}; "
              f"ratio = 1 would give {fit[0] + fit[1]:+.4f}")
        print("hypothesis 2 - solid inclusions below min_contour_area_px2, absent from the FEM geometry:")
        for key in ("frac_small_solid", "n_small_solid", "n_inner_solid", "area_small_solid_px"):
            a = np.array([float(r[key]) for r in rows])
            r_ = np.corrcoef(rd, a)[0, 1] if a.std() > 0 else float("nan")
            print(f"  {key:>20}: Pearson r(rel_diff, x) = {r_:+.3f}   values {np.round(a, 4).tolist()}")
    out = Path(path).with_name(Path(path).stem + "_analysis.json")
    json.dump(rows, open(out, "w"), indent=1)
    print("\n->", out)


if __name__ == "__main__":
    main(sys.argv[1])
