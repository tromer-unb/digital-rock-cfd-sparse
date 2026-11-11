#!/usr/bin/env python3
"""B6 - independent-code benchmark on a real digital-rock window (finding F-008, docs/EXTERNAL-BENCHMARK-PLAN.md step 1):
pressure-driven Stokes flow through one segmented rock window solved with OpenFOAM v1912 (finite volume, simpleFoam/SIMPLEC)
on a PIXEL mesh - one hexahedral cell per pore pixel, written directly as a polyMesh - instead of scikit-fem P2/P1 on a
Gmsh contour mesh. Independent of the authors' Gmsh meshing AND of the 0.25 px contour tolerance, so the comparison tests
the discretised geometry as well as the solver.

Protocol P1 is reproduced exactly (fix round, finding F-001, user decision 2026-10-21):
  mask -> clean_mask -> crop_between_fixed_planes(band) -> select_full_percolating_components   (authors' functions, read-only)
  -> pad_reservoirs(b = 0.10 * L_s) on both flow-direction sides (fully open pore rows/columns)
  pressure inlet / traction-free outlet, no-slip on pore walls AND on the lateral window walls (decision 7)
  k_D = mu * Q * L_s / (W * dP),  Q = integrated boundary flux,  W = 4096 px (full transverse width),
  dP = difference of the FLUX-WEIGHTED plane pressures taken 1 px inside each reservoir ("face" definition, not "outer").
Pixel geometry follows the authors' convention (pipeline extract_contours): pore pixel (row, col) is the unit square
centred on x = col, y = (nrows - 1) - row, so the cell spans [col-0.5, col+0.5] x [y-0.5, y+0.5].

Pass criteria (declared before looking at the result):
  (i)   mass conservation: |Q_in - Q_out| / Q_in < 1e-6, and simpleFoam reports "SIMPLE solution converged";
  (ii)  k_D agrees with the P1 P2/P1 reference of the same case and window within the documented mesh/window spread
        (EXTERNAL-BENCHMARK-PLAN.md step 3 acceptance). The P2/P1 level-1 -> level-5 drift is the mesh spread
        (R3-H -0.77%, R2-V -0.25%); the staircase pixel boundary is a DIFFERENT discretisation of the same geometry,
        so a deviation of a few per cent is expected and is reported, not thresholded. The hard claim under test is the
        ORDER OF MAGNITUDE and the sign of the correction relative to the manuscript: the manuscript reports
        k = 270.3 px^2 for R2-V and 110.9 px^2 for R3-H against the P1 values 60.90 and 18.28 px^2, i.e. factors of
        4.4 and 6.1. Criterion (ii) is met if |k_OF / k_P2P1 - 1| < 0.25, which still separates P1 from the manuscript
        by more than a factor of 3.
  (iii) reported, no threshold: the same comparison for the "outer" pressure definition (imposed inlet/outlet values),
        which must come out LOWER than the face definition, as in the fix round (R3-V -10.7%).
Self-test (B6-S, `--selftest`, run and passed before the rock windows): a fully open channel of N pixels across,
with the same mesh writer, the same pressure BCs and the same probe planes, must reproduce the exact cell-centred
finite-volume result k = (N^2 + 2)/12 px^2 (the analytic h^2/12 plus the known O(1/N^2) wall discretisation).
Measured for N = 10, 20, 40, 80: relative deviation 2.000e-2, 5.000e-3, 1.250e-3, 3.125e-4 = exactly 2/N^2,
error ratio 4.00 per refinement, mass imbalance <= 1.2e-9. So the pixel mesh, the boundary conditions, the flux
integration and the flux-weighted probe-plane pressures carry no error of their own; the residual deviation on a
rock window is the staircase resolution of the throats, which is the quantity criterion (ii) is about.
Note on the normalisation: the two probe planes span L_s + 2*EPS_PROBE, while k (like scripts/fr_reservoir) divides
by L_s. The script reports both; at L_s = 3440 px the difference is 2/3442 = 0.058%.

Requires: the OpenFOAM v1912 environment on PATH (source /usr/share/openfoam/etc/bashrc on this host).
OUTPUT: results/b6_openfoam_window_<case>.json ; case under $B6_WORK.
"""
import json, os, re, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT.parents[1]                                    # .../QC-20261019-v1
sys.path.insert(0, str(QC / "scripts"))

DZ = 1.0                                               # one pixel thick, 2D case (empty front/back)
NU, RHO = 1.0, 1.0                                     # mu = rho*nu = 1, so k comes out in px^2 directly
DP = 1e-4                                              # kinematic pressure drop; keeps Re ~ 1e-3 (creeping flow)
EPS_PROBE = 1.0                                        # distance of the pressure planes inside the reservoirs, in px
END_TIME = int(os.environ.get("B6_END_TIME", "200000"))
RES_P = float(os.environ.get("B6_RES_P", "1e-8"))
RES_U = float(os.environ.get("B6_RES_U", "1e-9"))
COARSEN = int(os.environ.get("B6_COARSEN", "1"))       # block-OR pixel coarsening; 1 = native resolution
RELAX_U = float(os.environ.get("B6_RELAX_U", "0.9"))   # SIMPLEC tolerates 1.0; only changes the path, not the answer

HEADER = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"


# ----------------------------------------------------------------------------------------------- mask (protocol P1)
def build_mask(case, band, pad, level_for_cfg=5, tol=0.25):
    """Authors' pipeline up to the percolating selection, then the P1 reservoirs. Returns the padded boolean mask."""
    import fr_reservoir as fr
    m = fr.load_pipeline()
    patch, direction = fr.CASES[case]
    cfg = fr.build_cfg(m, level_for_cfg, tol, band, Path(os.environ.get("B6_WORK", "/tmp")) / "_cfg")
    sel, info = fr.prepare(m, cfg, patch, direction)
    padded, Ls, b = fr.pad_reservoirs(sel, direction, pad)
    W = sel.shape[0] if direction == "lr" else sel.shape[1]
    return padded, direction, int(Ls), int(b), int(W), int(info["band_px"])


# ----------------------------------------------------------------------------------------------- polyMesh writer
def write_poly_mesh(mask, direction, d):
    """One hex cell per pore pixel, written straight to constant/polyMesh in ascii.

    Cells are numbered row-major over the pore pixels. Internal faces are emitted in upper-triangular order
    (sorted by owner, then neighbour) as OpenFOAM requires; boundary faces follow, grouped per patch.
    Node (ix, iy, k) sits at x = ix - 0.5, y = iy - 0.5, z = k*DZ, so pore pixel (r, c) occupies
    ix in {c, c+1}, iy in {nrows-1-r, nrows-r}: exactly the pixel square of the authors' convention.
    """
    nrows, ncols = mask.shape
    cid = np.full((nrows, ncols), -1, np.int64)
    cid[mask] = np.arange(int(mask.sum()), dtype=np.int64)       # row-major ids
    n_cells = int(mask.sum())

    # --- nodes actually used (corners of pore pixels)
    need = np.zeros((nrows + 1, ncols + 1), bool)
    for dr in (0, 1):
        for dc in (0, 1):
            need[dr:dr + nrows, dc:dc + ncols] |= mask
    nid = np.full((nrows + 1, ncols + 1), -1, np.int64)
    nid[need] = np.arange(int(need.sum()), dtype=np.int64)
    n_pts2 = int(need.sum())
    rr, cc = np.nonzero(need)
    # grid row index r of the node array maps to iy = nrows - r  (y decreases with r)
    xs = cc.astype(np.float64) - 0.5
    ys = (nrows - rr).astype(np.float64) - 0.5
    n_points = 2 * n_pts2                                        # k = 0 then k = 1

    def P(r_, c_, k_):
        return nid[r_, c_] + k_ * n_pts2

    # --- internal faces
    right = mask[:, :-1] & mask[:, 1:]                           # (r,c)-(r,c+1)
    down = mask[:-1, :] & mask[1:, :]                            # (r,c)-(r+1,c)
    r_r, r_c = np.nonzero(right)
    d_r, d_c = np.nonzero(down)
    own = np.concatenate([cid[r_r, r_c], cid[d_r, d_c]])
    nei = np.concatenate([cid[r_r, r_c + 1], cid[d_r + 1, d_c]])
    kind = np.concatenate([np.zeros(r_r.size, np.int8), np.ones(d_r.size, np.int8)])
    fr_, fc_ = np.concatenate([r_r, d_r]), np.concatenate([r_c, d_c])
    order = np.lexsort((nei, own))                               # upper-triangular order
    own, nei, kind, fr_, fc_ = own[order], nei[order], kind[order], fr_[order], fc_[order]
    n_internal = own.size

    # --- boundary faces: (mask shifted outside) per side, plus the two empty z sides
    def side(kind_name):
        if kind_name == "xm":                                    # face at x = c-0.5, outward -x
            sel = mask & ~np.pad(mask, ((0, 0), (1, 0)))[:, :-1]
        elif kind_name == "xp":                                  # face at x = c+0.5, outward +x
            sel = mask & ~np.pad(mask, ((0, 0), (0, 1)))[:, 1:]
        elif kind_name == "yp":                                  # upper y face (row r-1 missing), outward +y
            sel = mask & ~np.pad(mask, ((1, 0), (0, 0)))[:-1, :]
        else:                                                    # "ym": lower y face (row r+1 missing), outward -y
            sel = mask & ~np.pad(mask, ((0, 1), (0, 0)))[1:, :]
        return np.nonzero(sel)

    def classify(kind_name, r_, c_):
        """Split a side's faces into the named patches: inlet / outlet / walls."""
        if direction == "lr":
            if kind_name == "xm":
                return np.where(c_ == 0, "inlet", "walls")
            if kind_name == "xp":
                return np.where(c_ == ncols - 1, "outlet", "walls")
            return np.full(r_.size, "walls")                     # lateral y sides and all internal pore walls
        else:                                                    # "tb": inlet at the top (row 0 = largest y)
            if kind_name == "yp":
                return np.where(r_ == 0, "inlet", "walls")
            if kind_name == "ym":
                return np.where(r_ == nrows - 1, "outlet", "walls")
            return np.full(r_.size, "walls")

    bnd = {"inlet": [], "outlet": [], "walls": []}
    for kind_name in ("xm", "xp", "yp", "ym"):
        r_, c_ = side(kind_name)
        if r_.size == 0:
            continue
        tags = classify(kind_name, r_, c_)
        for name in bnd:
            m_ = tags == name
            if m_.any():
                bnd[name].append((kind_name, r_[m_], c_[m_]))
    r_all, c_all = np.nonzero(mask)
    bnd["frontAndBack"] = [("zm", r_all, c_all), ("zp", r_all, c_all)]

    # --- face vertex order; each list gives the normal stated in the comment (right-hand rule)
    def verts(kind_name, r_, c_):
        iy0, iy1 = nrows - 1 - r_, nrows - r_                     # node rows: nid row index = nrows - iy
        g0, g1 = nrows - iy1, nrows - iy0                         # = r_, r_+1
        c0, c1 = c_, c_ + 1
        if kind_name == "int_x" or kind_name == "xp":             # +x
            return [P(g1, c1, 0), P(g0, c1, 0), P(g0, c1, 1), P(g1, c1, 1)]
        if kind_name == "xm":                                     # -x
            return [P(g1, c0, 0), P(g1, c0, 1), P(g0, c0, 1), P(g0, c0, 0)]
        if kind_name == "int_y" or kind_name == "ym":              # -y (lower face)
            return [P(g1, c0, 0), P(g1, c1, 0), P(g1, c1, 1), P(g1, c0, 1)]
        if kind_name == "yp":                                     # +y (upper face)
            return [P(g0, c0, 0), P(g0, c0, 1), P(g0, c1, 1), P(g0, c1, 0)]
        if kind_name == "zm":                                     # -z
            return [P(g1, c0, 0), P(g0, c0, 0), P(g0, c1, 0), P(g1, c1, 0)]
        return [P(g1, c0, 1), P(g1, c1, 1), P(g0, c1, 1), P(g0, c0, 1)]   # +z

    pm = d / "constant/polyMesh"
    pm.mkdir(parents=True, exist_ok=True)

    def write_list(name, cls, rows_iter, n, pre=""):
        with open(pm / name, "w") as fh:
            fh.write(HEADER.format(cls=cls, obj=name) + pre + f"{n}\n(\n")
            for chunk in rows_iter:
                fh.write(chunk)
            fh.write(")\n")

    # points
    def pts_chunks():
        for k_ in (0, 1):
            z = k_ * DZ
            for i in range(0, xs.size, 200000):
                sl = slice(i, i + 200000)
                yield "".join(f"({x:.1f} {y:.1f} {z:.1f})\n" for x, y in zip(xs[sl], ys[sl]))
    write_list("points", "vectorField", pts_chunks(), n_points)

    # faces: internal (upper-triangular) then each patch
    def face_rows(kind_name, r_, c_):
        """Vertex index array (n, 4) for a block of faces of one kind."""
        vv = verts(kind_name, r_, c_)
        return np.stack([np.broadcast_to(np.asarray(v), r_.shape) for v in vv], axis=1)

    def fmt(arr):
        return "".join("4(%d %d %d %d)\n" % tuple(row) for row in arr.tolist())

    def face_chunks():
        for i in range(0, n_internal, 200000):
            sl = slice(i, i + 200000)
            r_, c_, k_ = fr_[sl], fc_[sl], kind[sl]
            rows = np.empty((r_.size, 4), np.int64)
            for kk, nm in ((0, "int_x"), (1, "int_y")):
                m_ = k_ == kk
                if m_.any():
                    rows[m_] = face_rows(nm, r_[m_], c_[m_])
            yield fmt(rows)
        for name in PATCH_ORDER:
            for kind_name, r_, c_ in bnd[name]:
                for i in range(0, r_.size, 200000):
                    sl = slice(i, i + 200000)
                    yield fmt(face_rows(kind_name, r_[sl], c_[sl]))
    n_bnd = {name: sum(r_.size for _, r_, _ in v) for name, v in bnd.items()}
    n_faces = n_internal + sum(n_bnd.values())
    write_list("faces", "faceList", face_chunks(), n_faces)

    # owner / neighbour
    def own_chunks():
        for i in range(0, n_internal, 500000):
            yield "\n".join(map(str, own[i:i + 500000].tolist())) + "\n"
        for name in PATCH_ORDER:
            for kind_name, r_, c_ in bnd[name]:
                ids = cid[r_, c_]
                for i in range(0, ids.size, 500000):
                    yield "\n".join(map(str, ids[i:i + 500000].tolist())) + "\n"
    write_list("owner", "labelList", own_chunks(), n_faces,
               pre=f"// nCells: {n_cells} nFaces: {n_faces} nInternalFaces: {n_internal}\n")

    def nei_chunks():
        for i in range(0, n_internal, 500000):
            yield "\n".join(map(str, nei[i:i + 500000].tolist())) + "\n"
    write_list("neighbour", "labelList", nei_chunks(), n_internal)

    # boundary
    start = n_internal
    body = f"{len(PATCH_ORDER)}\n(\n"
    for name in PATCH_ORDER:
        t = "empty" if name == "frontAndBack" else ("wall" if name == "walls" else "patch")
        body += f"  {name} {{ type {t}; nFaces {n_bnd[name]}; startFace {start}; }}\n"
        start += n_bnd[name]
    body += ")\n"
    (pm / "boundary").write_text(HEADER.format(cls="polyBoundaryMesh", obj="boundary") + body)
    return dict(n_cells=n_cells, n_points=n_points, n_faces=n_faces, n_internal=n_internal,
                n_bnd={k: int(v) for k, v in n_bnd.items()},
                own=own, nei=nei, kind=kind, fr=fr_, fc=fc_)


PATCH_ORDER = ("inlet", "outlet", "walls", "frontAndBack")


# ----------------------------------------------------------------------------------------------- case files
def write_case(d, n_cells):
    files = {
        "constant/transportProperties": ("dictionary", f"transportModel Newtonian;\nnu {NU};\n"),
        "constant/turbulenceProperties": ("dictionary", "simulationType laminar;\n"),
        "system/controlDict": ("dictionary",
            f"application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime {END_TIME}; deltaT 1;\n"
            f"writeControl timeStep; writeInterval {END_TIME}; writeFormat ascii; writePrecision 14; runTimeModifiable false;\n"
            # NOTE: a surfaceFieldValue functionObject was tried here to log the boundary flux every iteration, but
            # this OpenFOAM build raises "FOAM FATAL IO ERROR: error in IOstream sha1" when it writes its properties
            # under -parallel. The cheaper equivalent used instead is to compare two runs with different
            # residualControl targets, which answers the same question (can we stop earlier?) directly in k.
            ),
        "system/fvSchemes": ("dictionary",
            "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
            "divSchemes { default none; div(phi,U) bounded Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }\n"
            "laplacianSchemes { default Gauss linear corrected; }\ninterpolationSchemes { default linear; }\n"
            "snGradSchemes { default corrected; }\n"),
        "system/fvSolution": ("dictionary",
            "solvers { p { solver GAMG; smoother GaussSeidel; tolerance 1e-11; relTol 0.001;\n"
            "    nCellsInCoarsestLevel 50; mergeLevels 1; agglomerator faceAreaPair; cacheAgglomeration true; }\n"
            "U { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0.001; } }\n"
            "SIMPLE { consistent yes; nNonOrthogonalCorrectors 0;\n"
            f"  residualControl {{ p {RES_P}; U {RES_U}; }} }}\nrelaxationFactors {{ equations {{ U {RELAX_U}; }} fields {{ p 1; }} }}\n"),
        # pressure-driven: fixed kinematic pressure at the inlet, 0 at the outlet, zero-gradient velocity on both
        # (OpenFOAM's traction-free equivalent); no-slip on pore walls and lateral window walls.
        "0/U": ("volVectorField", "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField {\n"
                "  inlet { type zeroGradient; }\n  outlet { type zeroGradient; }\n"
                "  walls { type noSlip; }\n  frontAndBack { type empty; } }\n"),
        "0/p": ("volScalarField", "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n"
                f"  inlet {{ type fixedValue; value uniform {DP}; }}\n  outlet {{ type fixedValue; value uniform 0; }}\n"
                "  walls { type zeroGradient; }\n  frontAndBack { type empty; } }\n"),
    }
    for rel, (cls, body) in files.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(HEADER.format(cls=cls, obj=p.name) + body)


def write_decompose(d, nproc, shape):
    """Geometric (`simple`) decomposition: this OpenFOAM build ships only the dummy scotch/metis stubs.
    The split is chosen as the factorisation of nproc that comes closest to square blocks, so that the
    communication surface stays small; load balance is set by the pore distribution and is reported by decomposePar."""
    nrows, ncols = shape
    best = min(((a, nproc // a) for a in range(1, nproc + 1) if nproc % a == 0),
               key=lambda ab: abs((ncols / ab[0]) / (nrows / ab[1]) - 1.0))
    nx, ny = best
    (d / "system/decomposeParDict").write_text(
        HEADER.format(cls="dictionary", obj="decomposeParDict")
        + f"numberOfSubdomains {nproc};\nmethod simple;\nsimpleCoeffs {{ n ({nx} {ny} 1); delta 0.001; }}\n")
    return nx, ny


def run(cmd, d, log, env=None):
    with open(d / log, "w") as fh:
        r = subprocess.run(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT, env=env)
    if r.returncode:
        tail = "".join(open(d / log).readlines()[-25:])
        raise RuntimeError(f"{cmd} failed in {d} (rc={r.returncode}); tail of {log}:\n{tail}")


# ----------------------------------------------------------------------------------------------- post-processing
def read_scalar_field(path, n_cells):
    """internalField of an ascii volScalarField (uniform or nonuniform)."""
    txt = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    if m:
        vals = np.fromstring(txt[m.end():txt.index(")", m.end())], sep=" ")
        assert vals.size == int(m.group(1)) == n_cells, (vals.size, m.group(1), n_cells)
        return vals
    m = re.search(r"internalField\s+uniform\s+([-\deE.+]+)", txt)
    return np.full(n_cells, float(m.group(1)))


def read_surface_field(path, n_internal):
    """internalField of an ascii surfaceScalarField (phi)."""
    txt = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    vals = np.fromstring(txt[m.end():txt.index(")", m.end())], sep=" ")
    assert vals.size == int(m.group(1)) == n_internal, (vals.size, m.group(1), n_internal)
    return vals


def final_residuals(log):
    """Initial residuals of the last SIMPLE iteration, per field. simpleFoam's own "solution converged" message is
    not usable as the sole criterion: when a velocity component is identically zero (a straight channel), its
    normalised residual stays at O(0.1) for ever and residualControl is never met even though the solution is exact."""
    blocks = re.split(r"^Time = \d+$", log, flags=re.M)
    last = blocks[-1] if blocks[-1].strip() else blocks[-2]
    out = {}
    for m in re.finditer(r"Solving for (\w+), Initial residual = ([\d.eE+-]+)", last):
        out[m.group(1)] = float(m.group(2))
    return out


def patch_sum(path, patch):
    txt = path.read_text()
    blk = txt[re.search(rf"\b{patch}\s*\{{", txt).end():]
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", blk)
    if not m:
        m2 = re.search(r"value\s+uniform\s+([-\deE.+]+)", blk)
        return 0.0 if m2 is None else float(m2.group(1))
    vals = np.fromstring(blk[m.end():blk.index(")", m.end())], sep=" ")
    return float(vals.sum())


def plane_pressures(t, mask, direction, b, Ls, n_internal, geom):
    """Flux-weighted and area-weighted mean pressure on the two probe planes, 1 px inside each reservoir.

    The probe planes coincide exactly with internal cell-face planes, and the reservoirs are fully open pore,
    so each plane is a complete row/column of faces: the sum runs over every internal face on that plane.
    Face pressure is the arithmetic mean of the two adjacent cell values (the mesh is uniform, so this is the
    second-order face interpolation OpenFOAM itself uses).
    """
    p = read_scalar_field(t / "p", geom["n_cells"])
    phi = read_surface_field(t / "phi", n_internal)
    own, nei, kind, fr_, fc_ = geom["own"], geom["nei"], geom["kind"], geom["fr"], geom["fc"]
    nrows, ncols = mask.shape
    pf = 0.5 * (p[own] + p[nei])
    out = {}
    if direction == "lr":
        # sample columns are b .. b+Ls-1; interfaces at x = b-0.5 and x = b+Ls-0.5
        x_in, x_out = (b - 0.5) - EPS_PROBE, (b + Ls - 0.5) + EPS_PROBE
        # internal +x face of column c sits at x = c+0.5
        sel_in = (kind == 0) & (fc_ == int(round(x_in - 0.5)))
        sel_out = (kind == 0) & (fc_ == int(round(x_out - 0.5)))
    else:
        # inlet at row 0 (largest y); sample rows are b .. b+Ls-1. The internal -y face indexed fr_ = r lies
        # between rows r and r+1, at y = nrows - r - 1.5. Setting that equal to the probe planes of
        # scripts/fr_reservoir.interface_pressures (y = nrows - b + 0.5 and y = nrows - b - Ls - 1.5) gives:
        r_in, r_out = b - 2, b + Ls
        sel_in = (kind == 1) & (fr_ == r_in)
        sel_out = (kind == 1) & (fr_ == r_out)
    for name, sel in (("in", sel_in), ("out", sel_out)):
        q = np.abs(phi[sel])
        # the probe plane lies inside a fully open reservoir, so it must be a complete row/column of faces
        assert int(sel.sum()) == (nrows if direction == "lr" else ncols), (name, int(sel.sum()), nrows, ncols)
        out[name] = dict(n_faces=int(sel.sum()), q=float(q.sum()),
                         p_flux=float(np.sum(pf[sel] * q) / np.sum(q)), p_area=float(np.mean(pf[sel])))
    return out


# ----------------------------------------------------------------------------------------------- driver
def main():
    case = os.environ.get("B6_CASE", "R2-V")
    band = float(os.environ.get("B6_BAND", "0.08"))
    pad = float(os.environ.get("B6_PAD", "0.10"))
    nproc = int(os.environ.get("B6_NPROC", "12"))
    work = Path(os.environ["B6_WORK"])
    d = work / (f"{case}_band{band:g}_pad{pad:g}" + ("" if COARSEN == 1 else f"_c{COARSEN}"))
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    mask, direction, Ls, b, W, band_px = build_mask(case, band, pad)
    if COARSEN > 1:
        # scaling probe only: block-OR by COARSEN px. Lengths are rescaled so k stays in native px^2,
        # and connectivity is preserved (a block is pore if ANY of its pixels is pore).
        k_ = COARSEN
        nr, nc = (mask.shape[0] // k_) * k_, (mask.shape[1] // k_) * k_
        mask = mask[:nr, :nc].reshape(nr // k_, k_, nc // k_, k_).any(axis=(1, 3))
        Ls, b, W = Ls // k_, b // k_, W // k_
        print(f"[coarsen x{k_}] mask {mask.shape} Ls={Ls} b={b} W={W} pore={int(mask.sum()):_}", flush=True)
    t_mask = time.time() - t0
    print(f"[{case}] mask {mask.shape} dir={direction} Ls={Ls} b={b} W={W} pore_px={int(mask.sum()):_}", flush=True)

    t1 = time.time()
    geom = write_poly_mesh(mask, direction, d)
    write_case(d, geom["n_cells"])
    t_mesh = time.time() - t1
    print(f"[{case}] mesh {geom['n_cells']:_} cells, {geom['n_faces']:_} faces "
          f"({geom['n_internal']:_} internal) in {t_mesh:.0f}s; patches {geom['n_bnd']}", flush=True)

    run(["checkMesh"], d, "log.checkMesh")
    cm = (d / "log.checkMesh").read_text()
    print(f"[{case}] checkMesh: {'OK' if 'Mesh OK' in cm else 'FAILED'}", flush=True)

    t2 = time.time()
    if nproc > 1:
        nx, ny = write_decompose(d, nproc, mask.shape)
        print(f"[{case}] decomposing into {nx} x {ny} = {nproc} subdomains", flush=True)
        run(["decomposePar"], d, "log.decomposePar")
        run(["mpirun", "-np", str(nproc), "simpleFoam", "-parallel"], d, "log.simpleFoam")
        run(["reconstructPar", "-latestTime"], d, "log.reconstructPar")
    else:
        run(["simpleFoam"], d, "log.simpleFoam")
    t_solve = time.time() - t2
    log = (d / "log.simpleFoam").read_text()
    res = final_residuals(log)
    u_main = "Ux" if direction == "lr" else "Uy"
    res_ok = res.get("p", 1.0) < 1e-8 and res.get(u_main, 1.0) < 1e-8

    t = max((p for p in d.iterdir() if re.fullmatch(r"\d+", p.name) and p.name != "0"), key=lambda p: int(p.name))
    q_in, q_out = abs(patch_sum(t / "phi", "inlet")), abs(patch_sum(t / "phi", "outlet"))
    Q = 0.5 * (q_in + q_out) / DZ                                 # flux per unit depth
    pl = plane_pressures(t, mask, direction, b, Ls, geom["n_internal"], geom)
    mu = RHO * NU
    rec = dict(case=case, solver="OpenFOAM v1912 simpleFoam (SIMPLEC), pixel hex polyMesh, pressure-driven",
               protocol="P1", coarsen=COARSEN, relax_U=RELAX_U, band=band, band_px=band_px, pad_frac=pad, pad_px=b, Ls=Ls, W=W, nu=NU, rho=RHO, dP_imposed=DP,
               n_cells=geom["n_cells"], n_faces=geom["n_faces"], n_internal=geom["n_internal"], n_bnd=geom["n_bnd"],
               mesh_ok="Mesh OK" in cm, iterations=int(t.name), hit_end_time=int(t.name) >= END_TIME,
               converged=bool(res_ok), foam_converged_flag="SIMPLE solution converged" in log,
               final_residuals=res, residual_field_checked=u_main,
               Q=Q, q_in=q_in / DZ, q_out=q_out / DZ, mass_imbalance=abs(q_in - q_out) / q_in,
               planes=pl, t_mask_s=t_mask, t_mesh_s=t_mesh, t_solve_s=t_solve, nproc=nproc, case_dir=str(d))
    dP_face = abs(pl["in"]["p_flux"] - pl["out"]["p_flux"])
    rec["dP_face"] = dP_face
    rec["kD_face"] = mu * Q * Ls / (W * dP_face * RHO)
    # the two planes actually span Ls + 2*EPS_PROBE; reported for completeness (0.058% at Ls = 3440)
    rec["kD_face_Lprobe"] = mu * Q * (Ls + 2 * EPS_PROBE) / (W * dP_face * RHO)
    rec["dP_area"] = abs(pl["in"]["p_area"] - pl["out"]["p_area"])
    rec["kD_area"] = mu * Q * Ls / (W * rec["dP_area"] * RHO)
    rec["dP_outer"] = DP
    rec["kD_outer"] = mu * Q * Ls / (W * DP * RHO)
    rec["Re"] = (Q / W) * Ls / NU
    if COARSEN > 1:
        # k has dimensions of area: in coarse-block units it must be multiplied by COARSEN^2 to read as px^2
        for key in ("kD_face", "kD_face_Lprobe", "kD_area", "kD_outer"):
            rec[key + "_px2"] = rec[key] * COARSEN ** 2
    
    ref = json.load(open(QC / "fix-round/fixround_summary.json"))["convergence"][case]
    rec["k_p2p1_ref"] = ref["k_final_iface"]
    rec["rel_diff_vs_p2p1"] = rec["kD_face"] * COARSEN ** 2 / ref["k_final_iface"] - 1
    rec["PASS_i"] = bool(res_ok and rec["mass_imbalance"] < 1e-6 and rec["mesh_ok"])
    rec["PASS_ii"] = bool(abs(rec["rel_diff_vs_p2p1"]) < 0.25)
    tagc = f"_coarsen{COARSEN}" if COARSEN > 1 else ""
    out = ROOT / f"results/b6_openfoam_window_{case}{tagc}.json"
    json.dump(rec, open(out, "w"), indent=1)
    print(f"[{case}] k_face={rec['kD_face'] * COARSEN ** 2:.4f} px^2  vs P2/P1 {ref['k_final_iface']:.4f}  "
          f"{100 * rec['rel_diff_vs_p2p1']:+.2f}%  | k_outer={rec['kD_outer'] * COARSEN ** 2:.4f}  Re={rec['Re']:.2e}  "
          f"imbalance={rec['mass_imbalance']:.2e}  iters={rec['iterations']}  PASS_i={rec['PASS_i']} PASS_ii={rec['PASS_ii']}")
    print("->", out)


def selftest(root):
    """B6-S: open channel of N px, pixel mesh -> exact cell-centred FV result k = (N^2 + 2)/12 px^2."""
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    Ls, b, ok = 60, 6, True
    print(f"{'N':>4} {'k(Ls+2)':>13} {'(N^2+2)/12':>13} {'rel':>12} {'2/N^2':>10} {'imbalance':>11}")
    for N in (10, 20, 40, 80):
        mask = np.ones((N, Ls + 2 * b), bool)
        d = root / f"N{N}"
        geom = write_poly_mesh(mask, "lr", d)
        write_case(d, geom["n_cells"])
        run(["checkMesh"], d, "log.checkMesh")
        assert "Mesh OK" in (d / "log.checkMesh").read_text(), f"checkMesh failed for N={N}"
        run(["simpleFoam"], d, "log.simpleFoam")
        t = max((q for q in d.iterdir() if re.fullmatch(r"\d+", q.name) and q.name != "0"), key=lambda q: int(q.name))
        qi, qo = abs(patch_sum(t / "phi", "inlet")), abs(patch_sum(t / "phi", "outlet"))
        Q = 0.5 * (qi + qo) / DZ
        pl = plane_pressures(t, mask, "lr", b, Ls, geom["n_internal"], geom)
        dP = abs(pl["in"]["p_flux"] - pl["out"]["p_flux"])
        k = Q * (Ls + 2 * EPS_PROBE) / (N * dP)
        exact = (N ** 2 + 2) / 12
        rel = k / exact - 1
        ok &= abs(rel) < 1e-6 and abs(qi - qo) / qi < 1e-6
        print(f"{N:>4} {k:>13.6f} {exact:>13.6f} {rel:>+12.2e} {2 / N ** 2:>10.3e} {abs(qi - qo) / qi:>11.1e}")
    print("B6-S", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest(os.environ.get("B6_WORK", "/tmp/b6_selftest") + "/selftest"))
    main()
