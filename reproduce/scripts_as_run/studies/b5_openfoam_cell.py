#!/usr/bin/env python3
"""B5 - independent-code benchmark (finding F-008, docs/EXTERNAL-BENCHMARK-PLAN.md step 2): periodic square array of cylinders
solved with OpenFOAM v1912 (finite volume, SIMPLEC, simpleFoam) instead of scikit-fem P2/P1, same geometry and same definition of k
(k/L^2 = nu <u> / f, <u> superficial velocity, f body force per unit mass, fully periodic cell, no-slip cylinder).
Creeping flow: nu = 100, f = 1 -> Re = <u> L / nu ~ 4e-6.
Pass criteria (declared before looking at the result):
  (i)  OpenFOAM mesh convergence: k changes < 1% between the two finest meshes, for every c;
  (ii) agreement with the P2/P1 reference (results/b2_periodic_cell.json, finest mesh) within 2% for c = 0.1 and c = 0.2
       (acceptance of the plan), using the finest OpenFOAM mesh;
  (iii) reported, no threshold: the same comparison for c = 0.05 and 0.4 and against the Hasimoto series.
Requires: source ~/opt/openfoam-sysroot/env_openfoam.sh (simpleFoam, gmshToFoam, createPatch, checkMesh on PATH).
OUTPUT: results/b5_openfoam_cell.json ; cases under $B5_WORK (default: a temporary directory).
"""
import json, os, re, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import gmsh

ROOT = Path(__file__).resolve().parents[1]

NU, F = 100.0, 1.0
CS = (0.05, 0.1, 0.2, 0.4)
HS = (0.04, 0.02, 0.01, 0.005)
DZ = 0.01
SEP = {"left": "(1 0 0)", "right": "(-1 0 0)", "bottom": "(0 1 0)", "top": "(0 -1 0)"}
NEIGHBOUR = {"left": "right", "right": "left", "bottom": "top", "top": "bottom"}

HEADER = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
FILES = {
    "constant/transportProperties": ("dictionary", f"transportModel Newtonian;\nnu {NU};\n"),
    "constant/turbulenceProperties": ("dictionary", "simulationType laminar;\n"),
    "system/controlDict": ("dictionary", "application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime 200000; deltaT 1;\n"
                           "writeControl timeStep; writeInterval 200000; writeFormat ascii; writePrecision 14; runTimeModifiable false;\n"),
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


def make_mesh(c, h, path):
    """Structured O-grid of the unit cell (4 transfinite hex blocks between the cylinder and the cell edges), extruded one layer,
    msh 2.2 for gmshToFoam. Opposite cell edges carry identical uniform node distributions, so the cyclic faces match."""
    R = float(np.sqrt(c / np.pi))
    nt = max(4, int(round(1 / h)))                                   # cells per cell edge (and per quarter arc)
    nr = max(4, int(round((np.sqrt(0.5) - R) / h)))                  # cells along the block diagonals
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    try:
        g = gmsh.model.geo; gmsh.model.add("cell")
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
            g.mesh.setTransfiniteSurface(s); g.mesh.setRecombine(2, s); surfs.append(s)
        g.extrude([(2, s) for s in surfs], 0, 0, DZ, numElements=[1], recombine=True)
        g.synchronize()
        vols = [t for _, t in gmsh.model.getEntities(3)]
        groups = {k: [] for k in ("left", "right", "bottom", "top", "cylinder", "frontAndBack")}
        for _, t in gmsh.model.getBoundary([(3, v) for v in vols], combined=True, oriented=False):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, t)
            if z1 - z0 < 1e-9: groups["frontAndBack"].append(t)
            elif x1 < 1e-9: groups["left"].append(t)
            elif x0 > 1 - 1e-9: groups["right"].append(t)
            elif y1 < 1e-9: groups["bottom"].append(t)
            elif y0 > 1 - 1e-9: groups["top"].append(t)
            else: groups["cylinder"].append(t)
        assert [len(v) for v in groups.values()] == [1, 1, 1, 1, 4, 8], groups
        for name, tags in groups.items():
            gmsh.model.addPhysicalGroup(2, tags, name=name)
        gmsh.model.addPhysicalGroup(3, vols, name="fluid")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(3); gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return nt, nr


def write_case(d):
    for rel, (cls, body) in FILES.items():
        p = d / rel; p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(HEADER.format(cls=cls, obj=p.name) + body)


def make_cyclic(d):
    """Set wall/empty types, then let createPatch build (and face-order) the translational cyclic pairs left_c/right_c, bottom_c/top_c."""
    b = d / "constant/polyMesh/boundary"; txt = b.read_text()
    txt = re.sub(r"(cylinder\s*\{[^}]*?type\s+)\w+", r"\1wall", txt)
    txt = re.sub(r"(frontAndBack\s*\{[^}]*?type\s+)\w+", r"\1empty", txt)
    b.write_text(txt)
    body = "pointSync false;\npatches (\n" + "".join(
        f"  {{ name {a}_c; patchInfo {{ type cyclic; neighbourPatch {NEIGHBOUR[a]}_c; transform translational; separationVector {SEP[a]}; "
        f"matchTolerance 1e-6; }} constructFrom patches; patches ({a}); }}\n" for a in NEIGHBOUR) + ");\n"
    (d / "system/createPatchDict").write_text(HEADER.format(cls="dictionary", obj="createPatchDict") + body)
    run(["createPatch", "-overwrite"], d, "log.createPatch")


def patch_sum(path, patch):
    """Sum of a surfaceScalarField (ascii) over one boundary patch."""
    txt = path.read_text()
    blk = txt[re.search(rf"\b{patch}\s*\{{", txt).end():]
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", blk)
    vals = np.array(blk[m.end():].split(")")[0].split(), float)
    assert len(vals) == int(m.group(1))
    return float(vals.sum())


def run(cmd, d, log):
    with open(d / log, "w") as fh:
        r = subprocess.run(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"{cmd} failed in {d}, see {log}")


def solve(args):
    c, h, work = args
    d = Path(work) / f"c{c}_h{h}"; d.mkdir(parents=True, exist_ok=True)
    write_case(d); nt, nr = make_mesh(c, h, d / "cell.msh")
    run(["gmshToFoam", "cell.msh"], d, "log.gmshToFoam"); make_cyclic(d)
    run(["checkMesh"], d, "log.checkMesh")
    cm = (d / "log.checkMesh").read_text()
    run(["simpleFoam"], d, "log.simpleFoam")
    t = max((p for p in d.iterdir() if re.fullmatch(r"\d+", p.name) and p.name != "0"), key=lambda p: int(p.name))
    # periodic incompressible flow: the flux per unit depth through any x = const line is the same, so <u> = Q'/L with L = 1
    q_left, q_right = patch_sum(t / "phi", "left_c") / DZ, patch_sum(t / "phi", "right_c") / DZ
    u_sup = 0.5 * (abs(q_left) + abs(q_right))
    log = (d / "log.simpleFoam").read_text()
    return dict(c=c, h=h, nt=nt, nr=nr, n_cells=int(re.search(r"cells:\s+(\d+)", cm).group(1)),
                max_non_orth=float(re.search(r"non-orthogonality Max:\s+([\d.eE+-]+)", cm).group(1)),
                mesh_ok="Mesh OK" in cm, iterations=int(t.name), converged="SIMPLE solution converged" in log,
                flux_mismatch=abs(abs(q_left) - abs(q_right)) / u_sup, k=NU * u_sup / F, Re=u_sup / NU, case=str(d))


if __name__ == "__main__":
    work = Path(os.environ.get("B5_WORK") or tempfile.mkdtemp(prefix="b5_"))
    nproc = int(os.environ.get("B5_NPROC", "8"))
    with ProcessPoolExecutor(nproc) as ex:
        runs = list(ex.map(solve, [(c, h, str(work)) for c in CS for h in HS]))
    b2 = {x["c"]: x for x in json.load(open(ROOT / "results/b2_periodic_cell.json"))["cases"]}
    out = {"solver": "OpenFOAM v1912 simpleFoam (SIMPLEC), structured O-grid hex mesh (gmsh transfinite), 1 layer, cyclic x/y",
           "nu": NU, "f": F, "work_dir": str(work), "cases": [], "criteria": {}}
    for c in CS:
        rr = sorted((r for r in runs if r["c"] == c), key=lambda r: -r["h"])
        kof = rr[-1]["k"]
        out["cases"].append(dict(c=c, mesh=rr, k_of=kof, conv_last_step=kof / rr[-2]["k"] - 1, k_p2p1=b2[c]["k_num"],
                                 rel_diff_vs_p2p1=kof / b2[c]["k_num"] - 1, k_hasimoto=b2[c]["k_hasimoto"],
                                 rel_diff_vs_hasimoto=kof / b2[c]["k_hasimoto"] - 1))
    cs = out["cases"]
    out["criteria"]["conv_max_last_step"] = max(abs(x["conv_last_step"]) for x in cs)
    out["criteria"]["max_rel_diff_vs_p2p1(c=0.1,0.2)"] = max(abs(x["rel_diff_vs_p2p1"]) for x in cs if x["c"] in (0.1, 0.2))
    out["criteria"]["all_converged"] = all(r["converged"] for r in runs)
    out["criteria"]["all_mesh_ok"] = all(r["mesh_ok"] for r in runs)
    out["criteria"]["PASS_i"] = out["criteria"]["conv_max_last_step"] < 0.01
    out["criteria"]["PASS_ii"] = out["criteria"]["max_rel_diff_vs_p2p1(c=0.1,0.2)"] < 0.02
    json.dump(out, open(ROOT / "results/b5_openfoam_cell.json", "w"), indent=1)
    for x in cs:
        print(x["c"], [round(m["k"], 6) for m in x["mesh"]], "vs P2/P1", round(x["k_p2p1"], 6),
              f"{100 * x['rel_diff_vs_p2p1']:+.2f}%", f"conv {100 * x['conv_last_step']:+.3f}%")
    print(out["criteria"])
