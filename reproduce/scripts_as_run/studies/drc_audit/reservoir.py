"""Protocol P1 driver: author's mask pipeline + open reservoirs + reference solver + Darcy permeability k_D.

k_D = mu Q L_s / (W dP)       Q: integrated boundary flux, L_s: sample length (without reservoirs), W: full transverse width,
dP: difference of flux-weighted pressures at the two sample faces (probed 1 px inside the reservoirs) = "face"; "outer" uses the
outer reservoir faces.  mu = 1, so k_D is in px^2 when the mesh is in px.
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
from . import stokes
from .legacy import load_pipeline

# mesh-size parameters of the five levels of the manuscript (script.sh): h_min, h_max, distance_min, distance_max [px]
LEVELS = {1: (12, 70, 12, 140), 2: (8, 50, 8, 100), 3: (4, 30, 4, 60), 4: (2, 15, 2, 30), 5: (1.5, 12, 1.5, 24)}
# case -> (mask file stem, direction). Mapping of R1..R3 to files is inferred (see docs/AUDIT-TRAIL.md, finding F-009).
CASES = {"R3-H": ("patch_y7600_x53200", "lr"), "R2-V": ("patch_y7600_x19000", "tb"), "R3-V": ("patch_y7600_x53200", "tb")}


def build_cfg(m, structures_dir, outdir, level=2, tol=0.25, band=0.08):
    hmin, hmax, dmin, dmax = LEVELS[level]
    return m.Config(structures_dir=Path(structures_dir), output_dir=Path(outdir), directions=("lr", "tb"), boundary_band_fraction=band,
                    h_min=hmin, h_max=hmax, distance_min=dmin, distance_max=dmax, contour_tolerance_px=tol, min_contour_area_px2=20.0,
                    min_object_size=20, closing_radius=0, keep_percolating="all", mesh_algorithm=6, smoothing_steps=5, viscosity=1.0, density=1.0)


def select_window(m, cfg, mask_png, direction):
    """Author's functions, unchanged: load pores (green), clean, crop flow-direction bands, keep percolating components."""
    pore = m.clean_mask(m.load_pore_mask(Path(mask_png), cfg), cfg)
    cropped, info = m.crop_between_fixed_planes(pore, direction, cfg)
    sel, _ = m.select_full_percolating_components(cropped, direction, cfg)
    return sel, info


def pad_reservoirs(sel, direction, frac):
    Ls = sel.shape[1] if direction == "lr" else sel.shape[0]
    b = int(round(frac * Ls))
    if b == 0:
        return sel, Ls, 0
    if direction == "lr":
        pad = np.ones((sel.shape[0], b), bool); return np.concatenate([pad, sel, pad], axis=1), Ls, b
    pad = np.ones((b, sel.shape[1]), bool); return np.concatenate([pad, sel, pad], axis=0), Ls, b


def run_mask(sel, direction, outdir, level=2, tol=0.25, pad=0.10, elements=("p2",), bc="pressure", lateral="noslip", m=None, structures_dir="."):
    """Mesh a (percolating) boolean window with reservoirs and solve. Returns a list of records, one per element."""
    m = m or load_pipeline()
    cfg = build_cfg(m, structures_dir, outdir, level, tol)
    if not np.any(sel):
        return [dict(status="not_percolating")]
    padded, Ls, b = pad_reservoirs(sel, direction, pad)
    W = sel.shape[0] if direction == "lr" else sel.shape[1]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    t0 = time.time(); mi = m.create_gmsh_mesh(padded, direction, Path(outdir), cfg); t_mesh = time.time() - t0
    recs = []
    for el in elements:
        r = stokes.solve(str(mi["mesh_path"]), direction, el, bc=bc, lateral=lateral, return_fields=(b > 0))
        face = stokes.interface_pressure_drop(r.pop("_fields"), direction, b, Ls, W) if b > 0 else {}
        rec = dict(status="ok", element=el, bc=bc, lateral=lateral, level=level, tol=tol, pad_frac=pad, pad_px=b, Ls=int(Ls), W=int(W), n_tri=r["n_tri"],
                   n_dof=r["n_dof"], residual_lu=r["residual_lu"], t_lu_s=r["t_lu_s"], t_mesh_s=t_mesh, Q=r["Q"], Q_balance=r["Q_balance"],
                   porosity=float(sel.mean()), kD_outer=r["Q"] * Ls / (W * r["dP_flux"]))
        if face:
            rec["kD_face"] = r["Q"] * Ls / (W * face["flux"]); rec["kD_face_len"] = r["Q"] * Ls / (W * face["len"])
        recs.append(rec)
    return recs


def run_case(case, structures_dir, outdir, level=2, tol=0.25, band=0.08, pad=0.10, elements=("p2",), bc="pressure", lateral="noslip", m=None):
    m = m or load_pipeline()
    stem, direction = CASES[case]
    cfg = build_cfg(m, structures_dir, outdir, level, tol, band)
    sel, info = select_window(m, cfg, Path(structures_dir) / f"{stem}_c0_mask.png", direction)
    recs = run_mask(sel, direction, outdir, level, tol, pad, elements, bc, lateral, m, structures_dir)
    for r in recs:
        r.update(case=case, band=band, band_px=int(info["band_px"]))
    return recs
