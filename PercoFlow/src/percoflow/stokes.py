"""Steady incompressible Stokes flow with stable mixed finite elements (scikit-fem).

Strong form (Omega = fluid domain, Gamma_w walls, Gamma_in / Gamma_out inlet and outlet curves):

    -mu Laplace(u) + grad(p) = 0,   div(u) = 0                in Omega
    u = 0                                                     on Gamma_w
    mu du/dn - p n = -p_in n   (pressure-driven inlet)         on Gamma_in
    mu du/dn - p n = 0         (traction-free outlet)          on Gamma_out

The open-boundary conditions are the natural ("do-nothing") conditions of the Laplacian form of the viscous term; for the
uniform parallel flow of an open reservoir they coincide with an imposed pressure. No pressure value is pinned: the outlet
condition fixes the pressure level.

Weak form: find (u, p) with u = 0 on Gamma_w such that for all (v, q)

    mu (grad u, grad v) - (p, div v) = -p_in <v.n, 1>_Gamma_in,      (q, div u) = 0,

which gives the symmetric indefinite saddle-point system [[A, B^T], [B, 0]] [u; p] = [f; 0].

Elements:
    "p2"   Taylor-Hood P2/P1 (reference). Contains the Poiseuille parabola exactly, so a straight channel is reproduced
           to round-off on any mesh (benchmarks/poiseuille.py checks this).
    "mini" P1 + cubic bubble / P1 (screening). Fewer unknowns; converges from below.
Both satisfy the inf-sup (LBB) condition with a mesh-independent constant.

Linear solvers:
    "lu"      sparse LU of the reduced system (wall unknowns eliminated). Algebraic error at round-off.
    "minres"  MINRES with the block-diagonal preconditioner diag(diag(A), diag(M_p)/mu) (diagonal version of the
              Silvester-Wathen preconditioner). MINRES stops on a backward error of the PRECONDITIONED system, which can
              be far below the true residual. The true relative residual ||K x - b|| / ||b|| is therefore always computed
              and a solution above ``accept_residual`` raises ``SolverNotConverged``.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, minres, splu
from skfem import (Basis, ElementTriMini, ElementTriP1, ElementTriP2, ElementVector, FacetBasis, Functional, LinearForm,
                   MeshTri, asm, bmat)
from skfem.helpers import dot
from skfem.models.general import divergence
from skfem.models.poisson import mass, vector_laplace

from .config import SolverOptions


class SolverNotConverged(RuntimeError):
    """Raised when the true residual of an iterative solution is above the acceptance threshold."""


@dataclass
class StokesSolution:
    """Solution of one Stokes problem: mesh, finite-element bases, degrees of freedom and diagnostics."""
    mesh: MeshTri
    basis_u: Basis
    basis_p: Basis
    u: np.ndarray            # velocity degrees of freedom (vector element)
    p: np.ndarray            # pressure degrees of freedom (P1, one per vertex)
    info: dict               # sizes, residuals, timings, boundary integrals

    def vertex_velocity(self):
        """(u_x, u_y) at the mesh vertices."""
        nd = self.basis_u.nodal_dofs
        return self.u[nd[0]], self.u[nd[1]]

    def vertex_pressure(self):
        """Pressure at the mesh vertices (the first N_vertices P1 degrees of freedom)."""
        return self.p[: self.mesh.p.shape[1]]


def _flow_axis(direction: str):
    """(component index, sign) of the main flow direction."""
    return (0, 1.0) if direction == "lr" else (1, -1.0)


def solve(mesh: MeshTri, direction: str = "lr", opts: SolverOptions | None = None) -> StokesSolution:
    """Assemble and solve the Stokes problem on ``mesh`` (physical groups INLET, OUTLET, WALLS)."""
    opts = opts or SolverOptions()
    ev = ElementVector(ElementTriMini() if opts.element == "mini" else ElementTriP2())
    bu = Basis(mesh, ev, intorder=4)
    bp = Basis(mesh, ElementTriP1(), intorder=4)
    n_total = bu.N + bp.N
    if n_total > opts.max_dofs:
        raise MemoryError(f"{n_total} unknowns exceed solver.max_dofs = {opts.max_dofs}; coarsen the mesh or raise the limit")
    t0 = time.time()
    # saddle-point matrix K = [[A, B^T], [B, 0]]: A viscous (vector Laplacian), B = -divergence
    A = (opts.mu * asm(vector_laplace, bu)).tocsr()
    B = (-asm(divergence, bu, bp)).tocsr()
    K = bmat([[A, B.T], [B, None]], format="csr")
    n_u = bu.N
    f_in = FacetBasis(mesh, ev, facets=mesh.boundaries["INLET"], intorder=4)

    x = np.zeros(K.shape[0])
    d_walls = bu.get_dofs("WALLS").all()
    if opts.inlet == "velocity":
        comp, sign = _flow_axis(direction)

        def u_in(xx):
            v = np.zeros((2,) + xx.shape[1:])
            v[comp] = sign * opts.u0
            return v
        d_in = bu.get_dofs("INLET").all()
        x[d_in] = f_in.project(u_in)[d_in]
        x[d_walls] = 0.0
        fixed = np.unique(np.concatenate([d_in, d_walls]))
        rhs_full = np.zeros(K.shape[0])
    else:
        fixed = np.unique(d_walls)

        @LinearForm
        def traction(v, w):
            return -opts.p_in * dot(w.n, v)
        rhs_full = np.zeros(K.shape[0])
        rhs_full[:n_u] = traction.assemble(f_in)
    # eliminate the prescribed unknowns by index slicing (keeps the reduced matrix symmetric)
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    rhs = rhs_full[free] - K[free][:, fixed] @ x[fixed]
    Kf = K[free][:, free].tocsc()
    t_asm = time.time() - t0

    t0 = time.time()
    info = {"element": opts.element, "method": opts.method, "inlet": opts.inlet}
    if opts.method == "lu":
        xs = splu(Kf).solve(rhs)
    else:
        Mp = asm(mass, bp)
        diag = np.concatenate([A.diagonal(), Mp.diagonal() / opts.mu])[free]
        P = LinearOperator(Kf.shape, matvec=lambda r: r / diag, dtype=float)
        iters = {"n": 0}

        def count(_):
            iters["n"] += 1
        try:
            xs, flag = minres(Kf, rhs, M=P, rtol=opts.minres_tol, maxiter=opts.minres_maxiter, callback=count)
        except TypeError:   # SciPy < 1.12 uses `tol`
            xs, flag = minres(Kf, rhs, M=P, tol=opts.minres_tol, maxiter=opts.minres_maxiter, callback=count)
        info.update(minres_flag=int(flag), minres_iterations=iters["n"])
    t_solve = time.time() - t0
    # TRUE relative residual of the reduced system, computed for every solver and stored with the result
    residual = float(np.linalg.norm(Kf @ xs - rhs) / np.linalg.norm(rhs))
    info.update(n_triangles=int(mesh.t.shape[1]), n_dof=int(K.shape[0]), n_free=int(len(free)), residual=residual,
                t_assembly_s=t_asm, t_solve_s=t_solve, nnz_reduced=int(Kf.nnz))
    if opts.method == "minres" and residual > opts.accept_residual:
        raise SolverNotConverged(f"MINRES true relative residual {residual:.2e} > solver.accept_residual = {opts.accept_residual:.1e}")
    x[free] = xs
    sol = StokesSolution(mesh, bu, bp, x[:n_u], x[n_u:], info)
    info.update(boundary_integrals(sol))
    info["dissipation"] = float(opts.mu * sol.u @ (asm(vector_laplace, bu) @ sol.u))
    if opts.inlet == "pressure":
        # discrete energy identity: viscous dissipation = p_in * Q (follows from the weak form with v = u, q = p)
        info["energy_identity_relerr"] = float(abs(info["dissipation"] - opts.p_in * info["Q"]) / max(abs(opts.p_in * info["Q"]), 1e-300))
    return sol


def boundary_integrals(sol: StokesSolution) -> dict:
    """Flux through INLET and OUTLET, mass imbalance, and three boundary means of the pressure (length, flux, nodal)."""
    mesh, bu, bp, u, p = sol.mesh, sol.basis_u, sol.basis_p, sol.u, sol.p
    out = {}
    for name in ("INLET", "OUTLET"):
        fb = FacetBasis(mesh, bu.elem, facets=mesh.boundaries[name], intorder=4)
        fbp = FacetBasis(mesh, ElementTriP1(), facets=mesh.boundaries[name], intorder=4)

        @Functional
        def one(w):
            return 1.0 + 0.0 * w.x[0]

        @Functional
        def flux(w):
            return dot(w["u"], w.n)

        @Functional
        def pint(w):
            return w["p"]

        @Functional
        def pflux(w):
            return w["p"] * dot(w["u"], w.n)
        uu, pp = fb.interpolate(u), fbp.interpolate(p)
        L = one.assemble(fb)
        Qn = flux.assemble(fb, u=uu)
        out[name] = {"length": float(L), "flux": float(Qn), "p_mean_length": float(pint.assemble(fbp, p=pp) / L),
                     "p_mean_flux": float(pflux.assemble(fb, u=uu, p=pp) / Qn) if Qn != 0 else float("nan"),
                     "p_mean_nodal": float(np.mean(p[bp.get_dofs(name).all()]))}
    Q = abs(out["INLET"]["flux"])
    return {"Q": Q, "mass_imbalance": float((out["INLET"]["flux"] + out["OUTLET"]["flux"]) / Q) if Q else float("nan"),
            "inlet": out["INLET"], "outlet": out["OUTLET"]}
