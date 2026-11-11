#!/usr/bin/env python3
"""B7 - independent-code benchmark on real digital-rock geometry, MATCHED PAIRS (finding F-008,
docs/EXTERNAL-BENCHMARK-PLAN.md step 1 and 3).

Why pairs instead of the whole window. Running OpenFOAM on a full 4096 px window at native pixel resolution is
23-38 h of wall time (measured: SIMPLE outer iterations grow as N^1.8 and the per-iteration cost as N^2; see
results/b6_openfoam_window_R2-V_coarsen*.json and MEMORY-README Sec. 5.2), and the only way to make it cheap -
coarsening the pixels - opens the throats and inflates k by tens of per cent, so it measures nothing physical.
Instead this script cuts the SAME percolating window into native-resolution square sub-windows and, on each one,
solves the identical boolean mask twice:
  FEM  reference: the audit's scikit-fem Taylor-Hood P2/P1 on the authors' Gmsh contour mesh (contour tol 0.25 px),
                  through drc_audit.reservoir.run_mask - i.e. the exact protocol P1 machinery that produced the
                  corrected k of the fix round;
  FV   challenger: OpenFOAM v1912 simpleFoam (SIMPLEC) on a pixel hex polyMesh written directly by
                  benchmarks/b6_openfoam_window.py - a different code, a different method (finite volume vs mixed
                  finite elements) AND a different geometry discretisation (staircase pixels vs 0.25 px contour).
One declared preprocessing step is applied to BOTH sides before solving: the mask is restricted to the
face-connected (4-neighbour) components that span inlet to outlet. The authors' percolation selection labels with
8-connectivity, so a diagonal pinch - two pore pixels touching only at a corner - counts as a flow path; a
finite-volume mesh connects cells only through shared faces, where that pinch has zero contact area. Keeping it
would leave islands with no fixed-pressure boundary (singular pressure matrix, GAMG dies in
DICPreconditioner::calcReciprocalD) and would hand the two sides different flow topologies. The pixel count removed
and the resulting k difference (8-connected vs face-connected, on the FEM side) are both recorded per sub-window.

Both sides use protocol P1 unchanged: open reservoirs b = 0.10 L_s on the flow sides, pressure inlet / traction-free
outlet, no-slip on the pore walls and on the lateral window walls (decision 7), and
k_D = mu Q L_s / (W dP) with dP the difference of the flux-weighted plane pressures 1 px inside the reservoirs.

Pass criteria (declared before looking at any result):
  (i)   per sub-window soundness: FV mass imbalance < 1e-6 with converged p and flow-direction U residuals;
        FEM LU residual < 1e-10 and |Q_balance| < 1e-10. A sub-window failing (i) is excluded and reported, not silently dropped.
  (ii)  paired agreement: median |k_FV / k_FEM - 1| over the sound sub-windows < 0.10, and every sound sub-window
        < 0.25. Rationale for the size of the bound: B6-S shows the FV side reproduces the exact discrete result
        k = (N^2+2)/12 to 1e-10 on an open channel, so the whole residual deviation here is the staircase
        representation of the pore walls, which scales as 2/N^2 per throat of N px. The claim under test is that
        an independent code and an independent geometry discretisation confirm the protocol-P1 permeability to
        within that spread - the manuscript's values sit a factor 4.4 (R2-V) and 6.1 (R3-H) away from P1, so even
        the 0.25 bound separates P1 from the manuscript by more than a factor of 3.
  (iii) reported, no threshold: (a) the FEM mesh-level 4 vs 5 change on each sub-window, which shows the FEM
        reference is itself converged; (b) the sign of the deviation, which the staircase argument predicts to be
        negative (a staircase wall is longer than the smooth contour, so it dissipates more); (c) the spread of k
        across sub-windows, which is the window dependence of finding F-003 measured at a smaller scale.

Requires the OpenFOAM v1912 environment on PATH (source /usr/share/openfoam/etc/bashrc on this host).
ENV: B7_CASE (R3-H|R2-V|R3-V), B7_TILE (px, default 1024), B7_NPROC, B7_WORK, B7_LEVELS (default "5,4"), B7_MAX_TILES.
OUTPUT: results/b7_subwindow_pair_<case>_tile<TILE>.json
"""
import importlib.util, json, os, re, shutil, sys, time
from pathlib import Path
import numpy as np
import scipy.ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(QC / "scripts"))

spec = importlib.util.spec_from_file_location("b6", Path(__file__).resolve().parent / "b6_openfoam_window.py")
b6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b6)

from drc_audit import legacy, reservoir                                            # noqa: E402

CASE = os.environ.get("B7_CASE", "R2-V")
TILE = int(os.environ.get("B7_TILE", "1024"))
NPROC = int(os.environ.get("B7_NPROC", "20"))
BAND, PAD, TOL = 0.08, 0.10, 0.25
LEVELS = [int(x) for x in os.environ.get("B7_LEVELS", "5,4").split(",")]
MAX_TILES = int(os.environ.get("B7_MAX_TILES", "99"))
MIN_CELLS_PER_RANK = int(os.environ.get("B7_MIN_CELLS_PER_RANK", "25000"))
STRUCT = ROOT / "data/structures"


def face_connected(sel, direction):
    """Restrict a pore mask to the FACE-connected (4-neighbour) components that span inlet to outlet.

    The authors' `select_full_percolating_components` labels with 8-connectivity, so two pore pixels touching only
    at a corner count as connected. A finite-volume mesh connects cells only through shared faces, so such a
    diagonal pinch has zero contact area and transmits nothing. Leaving it in has two consequences, both observed:
    a pore cluster reachable only through a diagonal pinch becomes an island with no fixed-pressure boundary, which
    makes the pressure matrix singular and kills GAMG in DICPreconditioner::calcReciprocalD (signal 8); and a flow
    path that exists for the contour mesh does not exist for the pixel mesh, which biases k downwards.
    So the mask is restricted here, and - so that the pair stays an exact pair - the SAME restricted mask is given
    to the FEM side as well. The number of pixels this removes is recorded: it is the geometric difference between
    a corner-connected and a face-connected representation of the same segmented image, not a numerical detail.
    """
    cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
    lab, n = ndi.label(sel, structure=cross)
    if n == 0:
        return np.zeros_like(sel), 0
    a, bnd = (lab[:, 0], lab[:, -1]) if direction == "lr" else (lab[0, :], lab[-1, :])
    keep = (set(np.unique(a)) & set(np.unique(bnd))) - {0}
    out = np.isin(lab, list(keep)) if keep else np.zeros_like(sel)
    return out, int(sel.sum() - out.sum())


def solve_fv(mask, direction, Ls, b, W, d, nproc):
    """OpenFOAM pixel-mesh side, using the b6 primitives unchanged."""
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    geom = b6.write_poly_mesh(mask, direction, d)
    b6.write_case(d, geom["n_cells"])
    t_mesh = time.time() - t0
    b6.run(["checkMesh"], d, "log.checkMesh")
    cm = (d / "log.checkMesh").read_text()
    t1 = time.time()
    # A geometric decomposition of a sparse pore mask can hand a subdomain (almost) no cells, and GAMG's coarsest-level
    # DIC preconditioner then divides by zero (signal 8 in DICPreconditioner::calcReciprocalD). So cap the rank count by
    # the cell count and verify decomposePar's own per-processor tally before launching the solver.
    nproc_eff = max(1, min(nproc, geom["n_cells"] // MIN_CELLS_PER_RANK))
    while nproc_eff > 1:
        for stale in d.glob("processor*"):
            shutil.rmtree(stale)
        b6.write_decompose(d, nproc_eff, mask.shape)
        b6.run(["decomposePar", "-force"], d, "log.decomposePar")
        counts = [int(x) for x in re.findall(r"Number of cells = (\d+)", (d / "log.decomposePar").read_text())]
        if counts and min(counts) >= MIN_CELLS_PER_RANK // 10:
            break
        nproc_eff = nproc_eff // 2
    if nproc_eff > 1:
        b6.run(["mpirun", "-np", str(nproc_eff), "simpleFoam", "-parallel"], d, "log.simpleFoam")
        b6.run(["reconstructPar", "-latestTime"], d, "log.reconstructPar")
    else:
        for stale in d.glob("processor*"):
            shutil.rmtree(stale)
        b6.run(["simpleFoam"], d, "log.simpleFoam")
    t_solve = time.time() - t1
    log = (d / "log.simpleFoam").read_text()
    res = b6.final_residuals(log)
    u_main = "Ux" if direction == "lr" else "Uy"
    t = max((q for q in d.iterdir() if re.fullmatch(r"\d+", q.name) and q.name != "0"), key=lambda q: int(q.name))
    qi, qo = abs(b6.patch_sum(t / "phi", "inlet")), abs(b6.patch_sum(t / "phi", "outlet"))
    Q = 0.5 * (qi + qo) / b6.DZ
    pl = b6.plane_pressures(t, mask, direction, b, Ls, geom["n_internal"], geom)
    dP = abs(pl["in"]["p_flux"] - pl["out"]["p_flux"])
    return dict(n_cells=geom["n_cells"], n_faces=geom["n_faces"], mesh_ok="Mesh OK" in cm, iterations=int(t.name),
                final_residuals=res, residual_field=u_main,
                res_ok=bool(res.get("p", 1.0) < 1e-8 and res.get(u_main, 1.0) < 1e-8),
                Q=Q, mass_imbalance=abs(qi - qo) / qi, dP_face=dP, kD_face=Q * Ls / (W * dP),
                kD_outer=Q * Ls / (W * b6.DP), Re=(Q / W) * Ls / b6.NU,
                t_mesh_s=t_mesh, t_solve_s=t_solve, nproc=nproc_eff, nproc_asked=nproc,
                decomp_cells=counts if nproc_eff > 1 else None, case_dir=str(d))


def main():
    work = Path(os.environ["B7_WORK"])
    work.mkdir(parents=True, exist_ok=True)
    m = legacy.load_pipeline()
    stem, direction = reservoir.CASES[CASE]
    cfg = reservoir.build_cfg(m, STRUCT, work / "_cfg", level=LEVELS[0], tol=TOL, band=BAND)
    sel, info = reservoir.select_window(m, cfg, STRUCT / f"{stem}_c0_mask.png", direction)
    print(f"[{CASE}] full percolating window {sel.shape} dir={direction} porosity={sel.mean():.4f} "
          f"band_px={int(info['band_px'])}", flush=True)

    nflow, ntrans = (sel.shape[1], sel.shape[0]) if direction == "lr" else (sel.shape[0], sel.shape[1])
    n_f, n_t = nflow // TILE, ntrans // TILE
    print(f"[{CASE}] tiling into {n_f} x {n_t} = {n_f * n_t} sub-windows of {TILE} px "
          f"({nflow - n_f * TILE} px of flow length and {ntrans - n_t * TILE} px of width left over)", flush=True)

    out = {"case": CASE, "direction": direction, "tile_px": TILE, "band": BAND, "pad_frac": PAD, "contour_tol_px": TOL,
           "fem": "scikit-fem Taylor-Hood P2/P1 on the authors' Gmsh contour mesh (drc_audit.reservoir.run_mask)",
           "fv": "OpenFOAM v1912 simpleFoam (SIMPLEC) on a pixel hex polyMesh (benchmarks/b6_openfoam_window.py)",
           "full_window_shape": list(sel.shape), "full_window_porosity": float(sel.mean()),
           "k_p2p1_full_window_ref": json.load(open(QC / "fix-round/fixround_summary.json"))["convergence"][CASE]["k_final_iface"],
           "tiles": []}
    done = 0
    for i_f in range(n_f):
        for i_t in range(n_t):
            if done >= MAX_TILES:
                break
            sl_f = slice(i_f * TILE, (i_f + 1) * TILE)
            sl_t = slice(i_t * TILE, (i_t + 1) * TILE)
            tile = sel[sl_f, sl_t] if direction == "tb" else sel[sl_t, sl_f]
            tag = f"t{i_f}_{i_t}"
            rec = {"tag": tag, "i_flow": i_f, "i_trans": i_t, "porosity_raw": float(tile.mean())}
            # authors' own percolation selection, applied to the sub-window
            try:
                tsel, _ = m.select_full_percolating_components(tile.copy(), direction, cfg)
            except RuntimeError as e:
                rec.update(status="not_percolating", error=str(e))
                out["tiles"].append(rec); print(f"  {tag}: not percolating ({e})", flush=True); continue
            if not np.any(tsel):
                rec.update(status="not_percolating")
                out["tiles"].append(rec); print(f"  {tag}: not percolating (empty selection)", flush=True); continue
            rec["porosity_percolating"] = float(tsel.mean())
            tsel4, n_dropped = face_connected(tsel, direction)
            rec["px_dropped_by_face_connectivity"] = n_dropped
            rec["frac_dropped_by_face_connectivity"] = n_dropped / float(tsel.sum())
            if not np.any(tsel4):
                rec.update(status="not_percolating_face_connected")
                out["tiles"].append(rec)
                print(f"  {tag}: percolates only through diagonal pinches - no face-connected path", flush=True)
                continue

            # --- FEM reference: on the face-connected mask (the paired value, each mesh level) and on the authors'
            #     8-connected mask (the protocol value), so the cost of the connectivity choice is measured too
            rec["fem"], rec["fem8"] = {}, {}
            for lvl in LEVELS:
                rec["fem"][str(lvl)] = reservoir.run_mask(
                    tsel4, direction, work / tag / f"fem4_L{lvl}", level=lvl, tol=TOL, pad=PAD,
                    elements=("p2",), bc="pressure", lateral="noslip", m=m, structures_dir=STRUCT)[0]
            rec["fem8"][str(LEVELS[0])] = reservoir.run_mask(
                tsel, direction, work / tag / f"fem8_L{LEVELS[0]}", level=LEVELS[0], tol=TOL, pad=PAD,
                elements=("p2",), bc="pressure", lateral="noslip", m=m, structures_dir=STRUCT)[0]
            fem = rec["fem"][str(LEVELS[0])]
            if fem.get("status") != "ok":
                rec["status"] = "fem_failed"
                out["tiles"].append(rec); print(f"  {tag}: FEM failed ({fem})", flush=True); continue
            f8 = rec["fem8"][str(LEVELS[0])]
            rec["fem_8conn_vs_4conn"] = (f8["kD_face"] / fem["kD_face"] - 1) if f8.get("status") == "ok" else None

            # --- FV challenger on the identical (face-connected) mask
            padded, Ls, b = reservoir.pad_reservoirs(tsel4, direction, PAD)
            W = tsel4.shape[0] if direction == "lr" else tsel4.shape[1]
            try:
                rec["fv"] = solve_fv(padded, direction, Ls, b, W, work / tag / "fv", NPROC)
            except Exception as e:                                   # one bad sub-window must not kill the sweep
                rec.update(status="fv_failed", error=str(e)[:800])
                out["tiles"].append(rec)
                print(f"  {tag}: FV FAILED - {str(e)[:200]}", flush=True)
                continue

            k_fem, k_fv = fem["kD_face"], rec["fv"]["kD_face"]
            rec["rel_diff"] = k_fv / k_fem - 1
            rec["fem_level_change"] = (rec["fem"][str(LEVELS[0])]["kD_face"] / rec["fem"][str(LEVELS[1])]["kD_face"] - 1
                                       if len(LEVELS) > 1 and rec["fem"][str(LEVELS[1])].get("status") == "ok" else None)
            rec["sound"] = bool(rec["fv"]["res_ok"] and rec["fv"]["mass_imbalance"] < 1e-6 and rec["fv"]["mesh_ok"]
                                and fem["residual_lu"] < 1e-10 and abs(fem["Q_balance"]) < 1e-10)
            rec["status"] = "ok"
            out["tiles"].append(rec)
            done += 1
            lvlch = rec["fem_level_change"]
            lvlch_s = "n/a" if lvlch is None else f"{100 * lvlch:+.3f}%"
            print(f"  {tag}: phi={rec['porosity_percolating']:.3f}  k_FEM={k_fem:10.4f}  k_FV={k_fv:10.4f}  "
                  f"{100 * rec['rel_diff']:+7.2f}%  FEM L{LEVELS[0]}vsL{LEVELS[1]} {lvlch_s:>9}  "
                  f"cells={rec['fv']['n_cells']:_}  it={rec['fv']['iterations']}  sound={rec['sound']}", flush=True)

    sound = [t for t in out["tiles"] if t.get("sound")]
    devs = sorted(abs(t["rel_diff"]) for t in sound)
    out["criteria"] = {
        "n_tiles": len(out["tiles"]),
        "n_percolating_8conn": len([t for t in out["tiles"] if t.get("status") in
                                    ("ok", "fv_failed", "fem_failed", "not_percolating_face_connected")]),
        "n_percolating_face_connected": len([t for t in out["tiles"] if t.get("status") in ("ok", "fv_failed", "fem_failed")]),
        "n_solved": len([t for t in out["tiles"] if t.get("status") == "ok"]),
        "n_failed": len([t for t in out["tiles"] if t.get("status") in ("fv_failed", "fem_failed")]),
        "n_sound": len(sound),
        "median_abs_rel_diff": float(np.median(devs)) if devs else None,
        "max_abs_rel_diff": float(devs[-1]) if devs else None,
        "mean_rel_diff_signed": float(np.mean([t["rel_diff"] for t in sound])) if sound else None,
        "k_fem_spread_ratio": (max(t["fem"][str(LEVELS[0])]["kD_face"] for t in sound)
                               / min(t["fem"][str(LEVELS[0])]["kD_face"] for t in sound)) if sound else None,
        "max_abs_fem_level_change": float(max(abs(t["fem_level_change"]) for t in sound
                                              if t["fem_level_change"] is not None)) if sound else None,
        "max_abs_fem_8conn_vs_4conn": float(max(abs(t["fem_8conn_vs_4conn"]) for t in sound
                                                 if t.get("fem_8conn_vs_4conn") is not None)) if sound else None,
        # every sub-window with a face-connected path must come out sound: no silent drops
        "PASS_i": bool(sound and len(sound) == len([t for t in out["tiles"] if t.get("status") in
                                                    ("ok", "fv_failed", "fem_failed")])),
        "PASS_ii": bool(devs and np.median(devs) < 0.10 and devs[-1] < 0.25),
    }
    p = ROOT / f"results/b7_subwindow_pair_{CASE}_tile{TILE}.json"
    json.dump(out, open(p, "w"), indent=1)
    print("\ncriteria:", json.dumps(out["criteria"], indent=1))
    print("->", p)


if __name__ == "__main__":
    main()
