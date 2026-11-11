#!/usr/bin/env python3
"""S1 - SYNTHETIC proof of concept for the "k per window" analysis (finding F-009 workaround; NOT real rock data).
A Boolean model of overlapping solid grains (log-normal radius) with a prescribed porosity is generated with a fixed seed; the pixel size is an
ASSUMED value (default 2 um, i.e. a fine-grained sandstone-like image at 2 mm field of view) so that physical units are plausible.
Checks: (1) the workflow runs on masks in the author's green-pore PNG format; (2) k_D per window and its dependence on the window size (REV-type
analysis) for an isotropic medium, where k(lr) ~ k(tb) is expected; (3) Kozeny-Carman indicator K = phi^3/(k S^2) >= 3; (4) k(phi) is monotone.
Every output carries synthetic=True.  USAGE: s1_synthetic_poc.py [--pixel-size-um 2] [--level 4]
OUTPUT: results/s1_synthetic_poc.json, data/synthetic/SYNTHETIC_*.png
"""
import argparse, json, sys, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from drc_audit import reservoir, synthetic, metrics
from drc_audit.legacy import load_pipeline

ap = argparse.ArgumentParser(); ap.add_argument("--pixel-size-um", type=float, default=2.0); ap.add_argument("--level", type=int, default=4); A = ap.parse_args()
m = load_pipeline(); tmp = Path(tempfile.mkdtemp()); (ROOT / "data" / "synthetic").mkdir(parents=True, exist_ok=True)
a_m = A.pixel_size_um * 1e-6
res = dict(synthetic=True, pixel_size_um_ASSUMED=A.pixel_size_um, level=A.level, porosity_series=[], rev=[])


def window_k(sub, direction, tag):
    from drc_audit import reservoir as R
    # keep only the component(s) connecting the two flow-direction faces, exactly like the author's selection
    cfg = R.build_cfg(m, ".", tmp / tag, A.level, 0.25, 0.08)
    sel, _ = m.select_full_percolating_components(sub, direction, cfg)
    if not np.any(sel): return None
    r = R.run_mask(sel, direction, tmp / tag, level=A.level, tol=0.25, pad=0.10, elements=("p2",), bc="pressure", lateral="freeslip", m=m)[0]
    K, phi, S = metrics.kozeny_carman_K(r["kD_face"], sel)
    return dict(k_px2=r["kD_face"], **metrics.to_si(r["kD_face"], a_m), phi_percolating=phi, K_KC=K, residual=r["residual_lu"], n_tri=r["n_tri"])


# (a) porosity series on 1024^2 media, both directions
for phi_t, seed in ((0.25, 1), (0.40, 2), (0.55, 3)):
    mask, phi = synthetic.boolean_disks_mask((1024, 1024), phi_t, seed=seed)
    synthetic.save_green_png(mask, ROOT / "data" / "synthetic" / f"SYNTHETIC_boolean_phi{phi_t:.2f}_seed{seed}.png")
    row = dict(phi_target=phi_t, phi=phi, seed=seed)
    for d in ("lr", "tb"):
        row[d] = window_k(mask, d, f"p{phi_t}{d}")
    res["porosity_series"].append(row); print("phi", phi_t, row["lr"] and row["lr"]["k_px2"], row["tb"] and row["tb"]["k_px2"], flush=True)

# (b) window-size (REV-type) analysis on the phi = 0.40 medium (1024^2)
mask, phi = synthetic.boolean_disks_mask((1024, 1024), 0.40, seed=2)
for size in (256, 512, 1024):
    n = 1024 // size; ks = {"lr": [], "tb": []}
    for i in range(n):
        for j in range(n):
            sub = mask[i * size:(i + 1) * size, j * size:(j + 1) * size]
            for d in ("lr", "tb"):
                r = window_k(sub, d, f"w{size}_{i}{j}{d}")
                ks[d].append(None if r is None else r["k_px2"])
    row = dict(window_px=size, n_windows=n * n)
    for d in ("lr", "tb"):
        v = np.array([x for x in ks[d] if x is not None], float)
        row[d] = dict(n_percolating=int(len(v)), k_values=v.tolist(), mean=float(v.mean()) if len(v) else None,
                      cv=float(v.std(ddof=1) / v.mean()) if len(v) > 1 else None)
    res["rev"].append(row); print("window", size, row["lr"]["mean"], row["lr"]["cv"], row["tb"]["mean"], row["tb"]["cv"], flush=True)
json.dump(res, open(ROOT / "results" / "s1_synthetic_poc.json", "w"), indent=1)
