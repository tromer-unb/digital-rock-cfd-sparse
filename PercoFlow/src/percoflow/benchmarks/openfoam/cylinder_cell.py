"""Periodic cylinder array solved with OpenFOAM (finite volumes) and compared with the finite-element value of this package.

Case: unit cell, cylinder of solid fraction c, structured O-grid of 4 transfinite hexahedral blocks
between the cylinder and the cell edges (Gmsh), extruded by one layer (empty front and back), translational cyclic patches
built by createPatch (gmshToFoam alone leaves the cyclic faces unordered), body force through fvOptions
(vectorSemiImplicitSource), simpleFoam with SIMPLEC, nu = 100 and f = 1 so that Re ~ 1e-6.
k / L^2 = nu <u> / f, with <u> from the flux through a cyclic patch.

Pass criteria (declared before running):
  (i)  FV mesh convergence: k changes by < 1 % between the two finest FV meshes;
  (ii) agreement with the finite-element (Taylor-Hood) value within 2 % at c = 0.1 and 0.2 on the finest FV mesh.
The full suite (``percoflow benchmark --full``) uses four FV meshes (1/h = 25 ... 200) and c = 0.05, 0.1, 0.2, 0.4.
"""
from __future__ import annotations

import re
from pathlib import Path

import gmsh
import numpy as np

from ...config import OpenFOAMOptions
from .. import cylinder_array
from .environment import find_openfoam

NU, F, DZ = 100.0, 1.0, 0.01
SEP = {"left": "(1 0 0)", "right": "(-1 0 0)", "bottom": "(0 1 0)", "top": "(0 -1 0)"}
PAIR = {"left": "right", "right": "left", "bottom": "top", "top": "bottom"}
HEADER = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
FILES = {
    "constant/transportProperties": ("dictionary", f"transportModel Newtonian;\nnu {NU};\n"),
    "constant/turbulenceProperties": ("dictionary", "simulationType laminar;\n"),
    "constant/momentumTransport": ("dictionary", "simulationType laminar;\n"),
    "system/controlDict": ("dictionary", "application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime 200000;\n"
                           "deltaT 1; writeControl timeStep; writeInterval 200000; writeFormat ascii; writePrecision 14; runTimeModifiable false;\n"),
    "system/fvSchemes": ("dictionary", "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
                         "divSchemes { default none; div(phi,U) bounded Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }\n"
                         "laplacianSchemes { default Gauss linear corrected; }\ninterpolationSchemes { default linear; }\n"
                         "snGradSchemes { default corrected; }\n"),
    "system/fvSolution": ("dictionary", "solvers { p { solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0.01; }\n"
                          "U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0.01; } }\n"
                          "SIMPLE { consistent yes; nNonOrthogonalCorrectors 1; pRefCell 0; pRefValue 0;\n"
                          "  residualControl { p 1e-10; U 1e-11; } }\nrelaxationFactors { equations { U 0.9; } fields { p 1; } }\n"),
    "system/fvOptions": ("dictionary", "bodyForce { type vectorSemiImplicitSource; active true; selectionMode all; volumeMode specific;\n"
                         f"  injectionRateSuSp {{ U (({F} 0 0) 0); }} }}\n"),
    "0/U": ("volVectorField", "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField {\n"
            "  cylinder { type noSlip; }\n  \".*_c\" { type cyclic; }\n  frontAndBack { type empty; } }\n"),
    "0/p": ("volScalarField", "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n"
            "  cylinder { type zeroGradient; }\n  \".*_c\" { type cyclic; }\n  frontAndBack { type empty; } }\n"),
}


def ogrid_mesh(c: float, h: float, path) -> tuple:
    """Structured O-grid, msh 2.2 for gmshToFoam. Returns (cells per edge, cells along the block diagonals)."""
    R = float(np.sqrt(c / np.pi))
    nt = max(4, int(round(1 / h)))
    nr = max(4, int(round((np.sqrt(0.5) - R) / h)))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        g = gmsh.model.geo
        gmsh.model.add("cell")
        C = [g.addPoint(x, y, 0) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))]
        cc = g.addPoint(0.5, 0.5, 0)
        P = [g.addPoint(0.5 + R * np.cos(a), 0.5 + R * np.sin(a), 0) for a in np.radians((225, 315, 45, 135))]
        S = [g.addLine(C[i], C[(i + 1) % 4]) for i in range(4)]
        A = [g.addCircleArc(P[i], cc, P[(i + 1) % 4]) for i in range(4)]
        D = [g.addLine(P[i], C[i]) for i in range(4)]
        for t in S + A:
            g.mesh.setTransfiniteCurve(t, nt + 1)
        for t in D:
            g.mesh.setTransfiniteCurve(t, nr + 1)
        surfs = []
        for i in range(4):
            s = g.addPlaneSurface([g.addCurveLoop([A[i], D[(i + 1) % 4], -S[i], -D[i]])])
            g.mesh.setTransfiniteSurface(s)
            g.mesh.setRecombine(2, s)
            surfs.append(s)
        g.extrude([(2, s) for s in surfs], 0, 0, DZ, numElements=[1], recombine=True)
        g.synchronize()
        vols = [t for _, t in gmsh.model.getEntities(3)]
        groups = {k: [] for k in ("left", "right", "bottom", "top", "cylinder", "frontAndBack")}
        for _, t in gmsh.model.getBoundary([(3, v) for v in vols], combined=True, oriented=False):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, t)
            if z1 - z0 < 1e-9:
                groups["frontAndBack"].append(t)
            elif x1 < 1e-9:
                groups["left"].append(t)
            elif x0 > 1 - 1e-9:
                groups["right"].append(t)
            elif y1 < 1e-9:
                groups["bottom"].append(t)
            elif y0 > 1 - 1e-9:
                groups["top"].append(t)
            else:
                groups["cylinder"].append(t)
        for name, tags in groups.items():
            gmsh.model.addPhysicalGroup(2, tags, name=name)
        gmsh.model.addPhysicalGroup(3, vols, name="fluid")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return nt, nr


def _patch_sum(path, patch):
    """Sum of the face values of a patch in an ASCII field file (used for the flux through a cyclic patch)."""
    txt = Path(path).read_text()
    blk = txt[re.search(rf"\b{patch}\s*\{{", txt).end():]
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", blk)
    vals = np.array(blk[m.end():].split(")")[0].split(), float)
    return float(vals.sum())


def solve_fv(of, c: float, h: float, workdir) -> dict:
    """One finite-volume solution of the cell: O-grid mesh, gmshToFoam, cyclic patches, checkMesh, simpleFoam."""
    d = Path(workdir) / f"of_c{c}_h{h}"
    d.mkdir(parents=True, exist_ok=True)
    for rel, (cls, body) in FILES.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(HEADER.format(cls=cls, obj=p.name) + body)
    nt, nr = ogrid_mesh(c, h, d / "cell.msh")
    of.command(["gmshToFoam", "cell.msh"], d, "log.gmshToFoam")
    b = d / "constant/polyMesh/boundary"
    txt = b.read_text()
    txt = re.sub(r"(cylinder\s*\{[^}]*?type\s+)\w+", r"\1wall", txt)
    txt = re.sub(r"(frontAndBack\s*\{[^}]*?type\s+)\w+", r"\1empty", txt)
    b.write_text(txt)
    body = "pointSync false;\npatches (\n" + "".join(
        f"  {{ name {a}_c; patchInfo {{ type cyclic; neighbourPatch {PAIR[a]}_c; transform translational; separationVector {SEP[a]}; "
        f"matchTolerance 1e-6; }} constructFrom patches; patches ({a}); }}\n" for a in PAIR) + ");\n"
    (d / "system/createPatchDict").write_text(HEADER.format(cls="dictionary", obj="createPatchDict") + body)
    of.command(["createPatch", "-overwrite"], d, "log.createPatch")
    of.command(["checkMesh"], d, "log.checkMesh")
    cm = (d / "log.checkMesh").read_text()
    of.command(["simpleFoam"], d, "log.simpleFoam")
    t = max((p for p in d.iterdir() if re.fullmatch(r"\d+", p.name) and p.name != "0"), key=lambda p: int(p.name))
    ql, qr = _patch_sum(t / "phi", "left_c") / DZ, _patch_sum(t / "phi", "right_c") / DZ
    u = 0.5 * (abs(ql) + abs(qr))
    return {"c": c, "h": h, "cells_per_edge": nt, "cells_radial": nr, "mesh_ok": "Mesh OK" in cm,
            "converged": "SIMPLE solution converged" in (d / "log.simpleFoam").read_text(), "iterations": int(t.name),
            "k": NU * u / F, "Re": u / NU, "flux_mismatch": abs(abs(ql) - abs(qr)) / u}


def run(workdir, cs=(0.1, 0.2), hs=(0.04, 0.02, 0.01), h_fe: float = 0.01, opts: OpenFOAMOptions | None = None) -> dict:
    """B4. Default sizes keep the run short (minutes); ``--full`` uses hs = (0.04, 0.02, 0.01, 0.005) and
    cs = (0.05, 0.1, 0.2, 0.4)."""
    opts = opts or OpenFOAMOptions()
    of = find_openfoam(opts.env_script)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in cs:
        fv = [solve_fv(of, c, h, workdir) for h in hs]
        cylinder_array.cell_mesh(c, h_fe, workdir / f"fe_c{c}.msh")
        fe = cylinder_array.solve_cell(workdir / f"fe_c{c}.msh")
        k_fv = fv[-1]["k"]
        rows.append({"c": c, "fv_meshes": fv, "k_fv_finest": k_fv, "conv_last_step": k_fv / fv[-2]["k"] - 1.0,
                     "k_fe": fe["k"], "fe_mesh_h": h_fe, "rel_diff_fv_vs_fe": k_fv / fe["k"] - 1.0,
                     "k_hasimoto": float(cylinder_array.hasimoto(c))})
    crit_i = max(abs(r["conv_last_step"]) for r in rows)
    crit_ii = max(abs(r["rel_diff_fv_vs_fe"]) for r in rows if r["c"] in (0.1, 0.2)) if any(r["c"] in (0.1, 0.2) for r in rows) else None
    ok_runs = all(m["mesh_ok"] and m["converged"] for r in rows for m in r["fv_meshes"])
    return {"benchmark": "openfoam_cylinder_cell", "openfoam": of.describe(), "nu": NU, "f": F, "rows": rows,
            "criteria": {"i_fv_convergence_lt_1pct": crit_i, "ii_fv_vs_fe_lt_2pct": crit_ii, "all_mesh_ok_and_converged": ok_runs},
            "pass": bool(ok_runs and crit_i < 0.01 and (crit_ii is None or crit_ii < 0.02))}
