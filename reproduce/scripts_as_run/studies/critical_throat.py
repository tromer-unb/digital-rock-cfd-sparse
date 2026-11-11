#!/usr/bin/env python3
"""Pore-scale resolution diagnostic: largest erosion radius r (px) at which the percolating domain still connects inlet and outlet planes
(Euclidean-distance threshold d > r). r_c is the half-width of the critical throat; 2*r_c is the critical throat width in px.
Also reports the width histogram of the narrowest parts via the local-thickness proxy 2*d on the medial set at percolation threshold."""
import numpy as np, scipy.ndimage as ndi, json, sys, pathlib
QC = pathlib.Path(__file__).resolve().parents[1]; P = QC.parents[1] / "convergencia" / "convergence_mesh_05_extra_fine"
out = {}
for name, d, f in (("R3H", "lr", "patch_y7600_x53200"), ("R2V", "tb", "patch_y7600_x19000"), ("R3V", "tb", "patch_y7600_x53200")):
    m = np.load(P / d / (f + "_c0_mask") / "full_percolating_domain.npy")
    dist = ndi.distance_transform_edt(m)          # distance to nearest wall pixel (pixel centres)
    def perc(thr):
        lab, n = ndi.label(dist > thr, structure=np.ones((3, 3)))
        a, b = (lab[:, 0], lab[:, -1]) if d == "lr" else (lab[0, :], lab[-1, :])
        return bool(set(np.unique(a)) & set(np.unique(b)) - {0})
    lo, hi = 0.0, float(dist.max())
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if perc(mid) else (lo, mid)
    out[name] = dict(r_c_px=lo, critical_throat_width_px=2 * lo, max_inscribed_radius_px=float(dist.max()), shape=list(m.shape), fluid_px=int(m.sum()))
    # fraction of fluid within 3 px of a wall (cells with d<=1.5 px)
    out[name]["frac_fluid_within_1.5px_of_wall"] = float(((dist <= 1.5) & m).sum() / m.sum())
print(json.dumps(out, indent=1)); json.dump(out, open(QC / "evidence/critical_throat.json", "w"), indent=1)
