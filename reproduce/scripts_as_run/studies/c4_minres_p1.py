#!/usr/bin/env python3
"""C4 - the screening solver (Mini + preconditioned MINRES) under protocol P1 (npj revision, 2026-11-03).

The Mini/MINRES path of the pipeline was never run under P1 (the fix round solved Mini with sparse LU). This script
solves the P1 problem (pressure-driven, reservoirs 10%, tolerance 0.25 px, band 8%) on the fix-round meshes of
FR-E1_convergence_tol0.25 with exactly the pipeline's preconditioner (diag of the velocity stiffness, diag of the pressure
mass matrix / mu), its tolerance (rtol = 1e-7) and its per-level iteration limits (10000, 10000, 10000, 15000, 20000),
and repeats the solve with sparse LU on the same matrix. Each (case, level, method) runs in its own process so that
the peak resident set size belongs to that solve alone.

SciPy's minres (1.18.1) stops when ||r|| / (||A|| ||x||) <= rtol in the preconditioned norm (istop = 1) or when
||A r|| / (||A|| ||r||) <= rtol (istop = 2, a least-squares exit), and returns info = 0 for every istop except the
iteration limit. The stop reason is parsed from its verbose output and stored.

Declared criterion (before running): |k_D,face(MINRES) / k_D,face(LU) - 1| < 1e-3 for every case and level.
Also reported: iterations, stop reason, true relative residual ||b - K x|| / ||b||, wall time, peak RSS of both methods.
OUTPUT: results/c4_minres_p1.json
"""
import contextlib, io, json, os, re, resource, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT.parents[1]
FR = QC / "fix-round" / "FR-E1_convergence_tol0.25"
MAXITER = {1: 10000, 2: 10000, 3: 10000, 4: 15000, 5: 20000}
DIRECTION = {"R3H": "lr", "R2V": "tb", "R3V": "tb"}


def worker(case, level, method):
    sys.path.insert(0, str(ROOT / "src"))
    from scipy.sparse.linalg import LinearOperator, minres, splu
    from skfem import (Basis, ElementTriMini, ElementTriP1, ElementVector, FacetBasis, MeshTri, asm, bmat, LinearForm)
    from skfem.helpers import dot
    from skfem.models.general import divergence
    from skfem.models.poisson import vector_laplace, mass
    from drc_audit import stokes
    tag = f"{case}_L{level}_pad0.1_tol0.25_band0.08"
    ref = json.load(open(FR / f"{tag}_mini.json"))
    direction = DIRECTION[case]
    t0 = time.time()
    mesh = MeshTri.load(str(FR / "mesh" / tag / "mesh.msh"))
    ev = ElementVector(ElementTriMini()); bu = Basis(mesh, ev, intorder=4); bp = Basis(mesh, ElementTriP1(), intorder=4)
    A = asm(vector_laplace, bu).tocsr(); B = (-asm(divergence, bu, bp)).tocsr(); Mp = asm(mass, bp).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr"); nu = bu.N
    fixed = np.unique(bu.get_dofs("WALLS").all())
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    fin = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)

    @LinearForm
    def trac(v, w):
        return -1.0 * dot(w.n, v)
    f = np.zeros(K.shape[0]); f[:nu] = trac.assemble(fin)
    rhs = f[free]; Kf = K[free][:, free].tocsr()
    t_asm = time.time() - t0
    out = dict(case=case, level=level, method=method, n_free=int(len(free)), nnz=int(Kf.nnz),
               csr_MB=float((Kf.data.nbytes + Kf.indices.nbytes + Kf.indptr.nbytes) / 1024 ** 2), t_assemble_s=t_asm)
    t1 = time.time()
    if method == "lu":
        xs = splu(Kf.tocsc()).solve(rhs)
    else:
        d = np.concatenate([np.maximum(np.abs(A.diagonal()), 1e-14), np.maximum(np.abs(Mp.diagonal()), 1e-14)])[free]
        P = LinearOperator(Kf.shape, matvec=lambda v: v / d, dtype=float)
        it = [0]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            xs, info = minres(Kf, rhs, M=P, rtol=1e-7, maxiter=MAXITER[level], show=True, check=False,
                              callback=lambda x: it.__setitem__(0, it[0] + 1))
        log = buf.getvalue()
        stop = [l for l in log.splitlines() if l.strip().startswith("minres:") and ("solution" in l.lower() or "iteration" in l.lower()
                                                                                 or "condition" in l.lower() or "eps" in l.lower())]
        out.update(iterations=it[0], info=int(info), stop_message=stop[-1].strip() if stop else log.splitlines()[-1].strip())
    out["t_solve_s"] = time.time() - t1
    out["peak_rss_MB_after_solve"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0   # before post-processing
    out["rel_residual_true"] = float(np.linalg.norm(Kf @ xs - rhs) / np.linalg.norm(rhs))
    x = np.zeros(K.shape[0]); x[free] = xs
    fields = dict(mesh=mesh, u=x[:nu], p=x[nu:], bu=bu)
    face = stokes.interface_pressure_drop(fields, direction, ref["pad_px"], ref["Ls"], ref["W"])
    fb = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)
    from skfem import Functional

    @Functional
    def flux(w):
        return dot(w["u"], w.n)
    Q = abs(float(flux.assemble(fb, u=fb.interpolate(x[:nu]))))
    out.update(Q=Q, kD_face=Q * ref["Ls"] / (ref["W"] * face["flux"]), kD_face_LU_reference=ref["kD_iface_flux"],
               peak_rss_MB_incl_postprocessing=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0)
    out["rel_k_vs_LU_reference"] = out["kD_face"] / ref["kD_iface_flux"] - 1
    print(json.dumps(out), flush=True)


def driver():
    from concurrent.futures import ThreadPoolExecutor
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4", MKL_NUM_THREADS="4")
    rows = []

    def run_case(case):
        res = []
        for level in (1, 2, 3, 4, 5):
            for method in ("lu", "minres"):
                p = subprocess.run(["nice", "-n", "10", sys.executable, __file__, "worker", case, str(level), method],
                                   capture_output=True, text=True, env=env)
                line = [l for l in p.stdout.splitlines() if l.startswith("{")]
                r = json.loads(line[-1]) if line else dict(case=case, level=level, method=method, error=p.stderr[-2000:])
                res.append(r); print(json.dumps({k: r.get(k) for k in ("case", "level", "method", "iterations", "t_solve_s",
                                                                      "peak_rss_MB_after_solve", "rel_residual_true", "rel_k_vs_LU_reference", "error")}), flush=True)
        return res
    with ThreadPoolExecutor(3) as ex:
        for res in ex.map(run_case, ("R3H", "R2V", "R3V")):
            rows += res
    mr = [r for r in rows if r.get("method") == "minres" and "error" not in r]
    worst = max(abs(r["rel_k_vs_LU_reference"]) for r in mr) if mr else None
    crit = dict(declared="|kD_face(MINRES)/kD_face(LU) - 1| < 1e-3 for every case and level", n_minres=len(mr), worst=worst,
                PASS=bool(mr and len(mr) == 15 and worst < 1e-3))
    json.dump(dict(rows=rows, criteria=crit), open(ROOT / "results" / "c4_minres_p1.json", "w"), indent=1)
    print(json.dumps(crit, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    else:
        driver()
