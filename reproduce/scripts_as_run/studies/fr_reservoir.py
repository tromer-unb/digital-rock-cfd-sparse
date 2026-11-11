#!/usr/bin/env python3
"""
TEMPLATE-ID: none (audit/fix-round script, not a figure template)
DESCRIPTION: Fix-round FR (finding F-001). Pressure-driven Stokes flow with open inlet/outlet reservoirs.
  Re-uses the author's pipeline functions UNMODIFIED (imported read-only): mask loading, cleaning, cropping,
  percolating-component selection, contour extraction, Gmsh meshing. Between the percolating-component selection and
  the meshing, the mask is padded with fully open pore columns/rows (reservoirs) of width b = pad_frac * L_s on both
  flow-direction sides. The reviewer solver (independent_check.solve, direct LU, pressure-driven) is then applied.
  Darcy permeability:  k_D = mu * Q * L_s / (W * dP),  L_s = sample length (without reservoirs), W = full transverse width.
INPUT:  author masks (PROJETO/convergencia/github/structures/*.png)
OUTPUT: <out>/<tag>.json (one file per run) + mesh in <out>/mesh/<tag>/
USAGE:  fr_reservoir.py --case R3-H --level 2 --pad 0.05 --elements p2 mini --out <dir> [--tol 1.5] [--band 0.08]
"""
import argparse, importlib.util, json, sys, time
from pathlib import Path
import numpy as np

QC = Path(__file__).resolve().parents[1]
PROJ = QC.parents[1] / "convergencia"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_check import solve  # noqa: E402

LEVELS = {1: (12, 70, 12, 140), 2: (8, 50, 8, 100), 3: (4, 30, 4, 60), 4: (2, 15, 2, 30), 5: (1.5, 12, 1.5, 24)}
CASES = {"R3-H": ("patch_y7600_x53200", "lr"), "R2-V": ("patch_y7600_x19000", "tb"), "R3-V": ("patch_y7600_x53200", "tb")}


def load_pipeline():
    spec = importlib.util.spec_from_file_location("drc", PROJ / "digital_rock_cfd_pipeline_sparse.py")
    m = importlib.util.module_from_spec(spec); sys.modules["drc"] = m; spec.loader.exec_module(m)
    return m


def build_cfg(m, level, tol, band, outdir):
    hmin, hmax, dmin, dmax = LEVELS[level]
    return m.Config(structures_dir=PROJ / "github/structures", output_dir=Path(outdir), directions=("lr", "tb"),
                    boundary_band_fraction=band, h_min=hmin, h_max=hmax, distance_min=dmin, distance_max=dmax,
                    contour_tolerance_px=tol, min_contour_area_px2=20.0, min_object_size=20, closing_radius=0,
                    keep_percolating="all", mesh_algorithm=6, smoothing_steps=5, viscosity=1.0, density=1.0)


def prepare(m, cfg, patch, direction):
    img = PROJ / "github/structures" / (patch + "_c0_mask.png")
    pore = m.clean_mask(m.load_pore_mask(img, cfg), cfg)
    cropped, info = m.crop_between_fixed_planes(pore, direction, cfg)
    sel, _ = m.select_full_percolating_components(cropped, direction, cfg)
    return sel, info


def pad_reservoirs(sel, direction, frac):
    Ls = sel.shape[1] if direction == "lr" else sel.shape[0]
    b = int(round(frac * Ls))
    if b == 0:
        return sel, Ls, 0
    if direction == "lr":
        pad = np.ones((sel.shape[0], b), bool); out = np.concatenate([pad, sel, pad], axis=1)
    else:
        pad = np.ones((b, sel.shape[1]), bool); out = np.concatenate([pad, sel, pad], axis=0)
    return out, Ls, b


def interface_pressures(f, direction, b, Ls, W, eps=1.0, dy=0.5):
    """Pressure drop between the two sample faces, probed EPS px inside the open reservoirs (all points are fluid).
    Returns flux-weighted and length-weighted plane means. Excludes the spreading drop inside the reservoirs."""
    from skfem import Basis, ElementTriP1
    mesh, u, p, bu = f["mesh"], f["u"], f["p"], f["bu"]
    bp = Basis(mesh, ElementTriP1(), intorder=4)
    n_u = bu.N
    fu = bu.interpolator(u); fp = bp.interpolator(p[:bp.N] if p.size == bp.N else p)
    t = np.arange(0.0, W - 1 + 1e-9, dy)
    res = {}
    sign = 1.0 if direction == "lr" else -1.0
    planes = {}
    if direction == "lr":
        planes["in"] = b - 0.5 - eps; planes["out"] = b + Ls - 0.5 + eps
    else:   # tb: y axis is flipped (inlet at the top of the padded array = largest y)
        top = (b + Ls + b) - 1.0   # y of last row centre
        planes["in"] = top - (b - 0.5) + eps; planes["out"] = top - (b + Ls - 0.5) - eps
    vals = {}
    for k_, c in planes.items():
        pts = np.vstack([np.full_like(t, c), t]) if direction == "lr" else np.vstack([t, np.full_like(t, c)])
        uu = fu(pts); pp = fp(pts)
        un = (uu[0] if direction == "lr" else uu[1])
        vals[k_] = (float(np.mean(pp)), float(np.sum(pp * un) / np.sum(un)), float(np.sum(un) * dy))
    res["iface_len"] = abs(vals["in"][0] - vals["out"][0])
    res["iface_flux"] = abs(vals["in"][1] - vals["out"][1])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=list(CASES)); ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--pad", type=float, default=0.05); ap.add_argument("--tol", type=float, default=1.5)
    ap.add_argument("--band", type=float, default=0.08); ap.add_argument("--elements", nargs="+", default=["p2"])
    ap.add_argument("--bc", default="pressure", choices=["pressure", "velocity"])
    ap.add_argument("--out", required=True); ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    patch, direction = CASES[a.case]
    tag = a.tag or f"{a.case.replace('-', '')}_L{a.level}_pad{a.pad:g}_tol{a.tol:g}_band{a.band:g}"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    m = load_pipeline()
    mdir = out / "mesh" / tag; mdir.mkdir(parents=True, exist_ok=True)
    cfg = build_cfg(m, a.level, a.tol, a.band, mdir)
    t0 = time.time()
    try:
        sel, info = prepare(m, cfg, patch, direction)
    except RuntimeError as e:
        json.dump(dict(tag=tag, case=a.case, status="not_percolating", error=str(e)), open(out / f"{tag}.json", "w"))
        print(tag, "not percolating"); return
    if not np.any(sel):
        json.dump(dict(tag=tag, case=a.case, status="not_percolating", band=a.band, tol=a.tol), open(out / f"{tag}.json", "w"))
        print(tag, "not percolating (empty selection)"); return
    padded, Ls, b = pad_reservoirs(sel, direction, a.pad)
    W = sel.shape[0] if direction == "lr" else sel.shape[1]
    mi = m.create_gmsh_mesh(padded, direction, mdir, cfg)
    t_mesh = time.time() - t0
    for el in a.elements:
        r = solve(str(mi["mesh_path"]), direction, el, bc=a.bc, return_fields=(b > 0))
        Q = abs(r["inlet"]["Q"])
        iface = interface_pressures(r.pop("_fields"), direction, b, Ls, W) if b > 0 else {}
        rec = dict(tag=tag, case=a.case, status="ok", level=a.level, pad_frac=a.pad, pad_px=b, tol=a.tol, band=a.band,
                   bc=a.bc, element=el, direction=direction, Ls=int(Ls), W=int(W), n_tri=r["n_tri"], n_dof=r["n_dof"],
                   residual_lu=r["residual_lu"], t_lu_s=r["t_lu_s"], t_mesh_s=t_mesh, Q=Q, Q_balance=r["Q_balance"],
                   selected_pore_px=int(sel.sum()), porosity_sel=float(sel.mean()), band_px=int(info["band_px"]))
        for key, kk in r["k"].items():
            dP = kk["dP"]
            rec[f"dP_{key}"] = dP
            rec[f"kD_{key}"] = Q * Ls / (W * dP)   # mu = 1
        for key, dP in iface.items():
            rec[f"dP_{key}"] = dP
            rec[f"kD_{key}"] = Q * Ls / (W * dP)
        json.dump(rec, open(out / f"{tag}_{el}.json", "w"), indent=1)
        print(tag, el, "n_dof", rec["n_dof"], "kD_flux %.4f" % rec["kD_flux"], "t_lu %.0fs" % rec["t_lu_s"], flush=True)


if __name__ == "__main__":
    main()
