"""simpleFoam case for a pressure-driven 2D Stokes problem, and post-processing of its result.

Physics: laminar, steady, incompressible; kinematic pressure p (OpenFOAM convention) fixed to p_in on "inlet" and 0 on
"outlet", zero-gradient velocity there, no-slip on walls, "frontAndBack" empty. A large kinematic viscosity keeps the
Reynolds number far below one, so the convective term is negligible and the solution is the Stokes solution:
k = nu * U * L / dp with kinematic quantities is the same permeability as mu * U * L / dP with dynamic ones.

Numerics: SIMPLEC (consistent yes), linear Gauss schemes, GAMG for p, residual targets 1e-9 (p) and 1e-10 (U).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

HEADER = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"


def write_case(case_dir, nu: float, p_in: float = 1.0, max_iterations: int = 200_000, nprocs: int = 1):
    """Write controlDict, fvSchemes, fvSolution, transport/turbulence settings and the 0/U, 0/p fields of the case."""
    files = {
        "constant/transportProperties": ("dictionary", f"transportModel Newtonian;\nnu {nu:.10g};\n"),
        "constant/turbulenceProperties": ("dictionary", "simulationType laminar;\n"),
        "constant/momentumTransport": ("dictionary", "simulationType laminar;\n"),
        "system/controlDict": ("dictionary",
                               f"application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime {max_iterations};\n"
                               f"deltaT 1; writeControl timeStep; writeInterval {max_iterations}; writeFormat ascii; writePrecision 14;\n"
                               "runTimeModifiable false;\n"),
        "system/fvSchemes": ("dictionary",
                             "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
                             "divSchemes { default none; div(phi,U) bounded Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }\n"
                             "laplacianSchemes { default Gauss linear corrected; }\ninterpolationSchemes { default linear; }\n"
                             "snGradSchemes { default corrected; }\n"),
        "system/fvSolution": ("dictionary",
                              "solvers { p { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0.01; }\n"
                              "  U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0.01; } }\n"
                              "SIMPLE { consistent yes; nNonOrthogonalCorrectors 1;\n"
                              "  residualControl { p 1e-9; U 1e-10; } }\n"
                              "relaxationFactors { equations { U 0.9; } fields { p 1; } }\n"),
        "system/decomposeParDict": ("dictionary", f"numberOfSubdomains {max(1, nprocs)}; method scotch;\n"),
        "0/U": ("volVectorField", "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField {\n"
                "  inlet { type zeroGradient; }\n  outlet { type zeroGradient; }\n  walls { type noSlip; }\n"
                "  frontAndBack { type empty; } }\n"),
        "0/p": ("volScalarField", "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n"
                f"  inlet {{ type fixedValue; value uniform {p_in:.10g}; }}\n  outlet {{ type fixedValue; value uniform 0; }}\n"
                "  walls { type zeroGradient; }\n  frontAndBack { type empty; } }\n"),
    }
    for rel, (cls, body) in files.items():
        p = Path(case_dir) / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(HEADER.format(cls=cls, obj=p.name) + body)


def run(of, case_dir, nprocs: int = 1, check_mesh: bool = True) -> dict:
    """checkMesh (must report "Mesh OK"), then simpleFoam (serial or MPI). Returns run diagnostics."""
    case_dir = Path(case_dir)
    out = {}
    if check_mesh:
        of.command(["checkMesh"], case_dir, "log.checkMesh")
        cm = (case_dir / "log.checkMesh").read_text()
        out["mesh_ok"] = "Mesh OK" in cm
        m = re.search(r"non-orthogonality Max:\s+([\d.eE+-]+)", cm)
        out["max_non_orthogonality"] = float(m.group(1)) if m else None
        if not out["mesh_ok"]:
            raise RuntimeError(f"checkMesh failed for {case_dir}, see log.checkMesh")
    if nprocs > 1:
        of.command(["decomposePar", "-force"], case_dir, "log.decomposePar")
        of.command(["simpleFoam"], case_dir, "log.simpleFoam", nprocs=nprocs)
        of.command(["reconstructPar", "-latestTime"], case_dir, "log.reconstructPar")
    else:
        of.command(["simpleFoam"], case_dir, "log.simpleFoam")
    log = (case_dir / "log.simpleFoam").read_text()
    out.update(residual_history(log))
    out["iterations"] = int(latest_time(case_dir).name)
    return out


def residual_history(log: str, p_target: float = 1e-9) -> dict:
    """Convergence assessment. simpleFoam's residualControl stops the run when all initial residuals fall below their
    targets. A velocity component that is identically zero in the exact solution (e.g. u_y in a straight channel) keeps a
    normalized residual of order one made of round-off noise, so the run then ends at endTime although it has converged.
    The run is therefore accepted as converged when simpleFoam reports convergence OR the final pressure initial residual
    is below ``p_target``; both residuals are recorded."""
    p = [float(x) for x in re.findall(r"Solving for p, Initial residual = ([\d.eE+-]+)", log)]
    ux = [float(x) for x in re.findall(r"Solving for Ux, Initial residual = ([\d.eE+-]+)", log)]
    uy = [float(x) for x in re.findall(r"Solving for Uy, Initial residual = ([\d.eE+-]+)", log)]
    reported = "SIMPLE solution converged" in log
    p_last = p[-1] if p else float("nan")
    return {"converged": bool(reported or (p and p_last < p_target)), "converged_by": "simpleFoam residualControl" if reported
            else (f"pressure initial residual < {p_target:g}" if p and p_last < p_target else "not converged"),
            "p_initial_residual_last": p_last, "Ux_initial_residual_last": ux[-1] if ux else None,
            "Uy_initial_residual_last": uy[-1] if uy else None}


def latest_time(case_dir) -> Path:
    """The latest result time directory of a case (the converged solution of a steady run)."""
    dirs = [p for p in Path(case_dir).iterdir() if p.is_dir() and re.fullmatch(r"\d+(\.\d+)?", p.name) and float(p.name) > 0]
    if not dirs:
        raise RuntimeError(f"no result time directory in {case_dir}")
    return max(dirs, key=lambda p: float(p.name))


def _list_after(text: str, start: int, dtype=float):
    """Parse the OpenFOAM list "N ( v1 v2 ... )" that follows position ``start`` in ``text``."""
    m = re.compile(r"(\d+)\s*\(").search(text, start)
    n = int(m.group(1))
    body = text[m.end(): text.index(")", m.end())]
    vals = np.array(body.split(), dtype=float)
    if len(vals) != n:
        raise ValueError("unexpected OpenFOAM list length")
    return vals


def read_internal_scalar(path) -> np.ndarray:
    """Values of the internal field of an ASCII scalar field file (one per cell, or per internal face for phi)."""
    text = Path(path).read_text()
    i = text.index("internalField")
    if re.match(r"internalField\s+uniform", text[i:]):
        raise ValueError(f"{path}: uniform internal field (solver did not run?)")
    return _list_after(text, i)


def read_patch_values(path, patch: str) -> np.ndarray:
    """Face values of one boundary patch of an ASCII scalar field file."""
    text = Path(path).read_text()
    blk = re.search(rf"\b{re.escape(patch)}\s*\{{", text[text.index("boundaryField"):])
    start = text.index("boundaryField") + blk.end()
    m = re.compile(r"value\s+(nonuniform\s+List<scalar>|uniform)\s*").search(text, start)
    if m.group(1) == "uniform":
        return np.array([float(text[m.end():].split(";")[0])])
    return _list_after(text, m.end())


def line_pressures(case_dir, internal: dict, axis: int, coords: dict, dz: float = 1.0, tol: float = 1e-6) -> dict:
    """Flux-weighted and plain mean pressures on lines of internal faces normal to ``axis`` at the given coordinates.
    Face pressure = mean of the two adjacent cell pressures (exact for the uniform pixel grid)."""
    t = latest_time(case_dir)
    p = read_internal_scalar(t / "p")
    phi = read_internal_scalar(t / "phi")
    out = {}
    for key, c in coords.items():
        sel = (np.abs(internal["centre"][:, axis] - c) < tol) & (np.abs(internal["normal"][:, axis]) > 0)
        if not sel.any():
            raise ValueError(f"no internal faces on the probe line {key} at {c}")
        pf = 0.5 * (p[internal["owner"][sel]] + p[internal["neighbour"][sel]])
        sign = np.sign(internal["normal"][sel, axis])
        fl = phi[sel] * sign
        out[key] = {"p_flux": float(np.sum(pf * fl) / np.sum(fl)), "p_mean": float(np.mean(pf)), "Q": float(np.sum(fl) / dz),
                    "n_faces": int(sel.sum())}
    return out
