"""Command line interface: ``percoflow <command> ...`` (``percoflow <command> --help`` for the options).

    percoflow image  ROCK.png  -o runs/rock_lr  [--direction tb] [--level 3] [--config run.yaml]
    percoflow mesh   part.msh  -o runs/part     [--direction lr] [--length L --width W]
    percoflow benchmark -o runs/benchmarks      [--full]
    percoflow openfoam-check                     (find OpenFOAM and print its version)
    percoflow verify runs/rock_lr                (recompute the output hashes of a run)
    percoflow levels                             (print the mesh-level presets)

The OpenFOAM cross-check runs by default; ``--no-openfoam`` disables it and the manifest records that it was disabled.
"""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .config import MESH_LEVELS, RunConfig


def _config(args) -> RunConfig:
    """Configuration file (if any) first, then the command-line options on top of it."""
    cfg = RunConfig.load(args.config) if getattr(args, "config", None) else RunConfig()
    if getattr(args, "direction", None):
        cfg.domain.direction = args.direction
    if getattr(args, "level", None):
        cfg.mesh.level = args.level
    if getattr(args, "element", None):
        cfg.solver.element = args.element
    if getattr(args, "solver", None):
        cfg.solver.method = args.solver
    if getattr(args, "tolerance", None) is not None:
        cfg.mesh.contour_tolerance = args.tolerance
    if getattr(args, "band", None) is not None:
        cfg.domain.band_fraction = args.band
    if getattr(args, "no_openfoam", False):
        cfg.openfoam.enabled = False
    if getattr(args, "openfoam_env", None):
        cfg.openfoam.env_script = args.openfoam_env
    if getattr(args, "openfoam_nprocs", None):
        cfg.openfoam.nprocs = args.openfoam_nprocs
    cfg.validate()
    return cfg


def _common(p):
    """Options shared by the image and mesh commands."""
    p.add_argument("-o", "--out", required=True, help="output directory (created)")
    p.add_argument("--config", help="YAML or JSON configuration (see examples/config_default.yaml)")
    p.add_argument("--direction", choices=["lr", "tb"], help="flow left->right (lr) or top->bottom (tb)")
    p.add_argument("--element", choices=["p2", "mini"], help="p2 = Taylor-Hood (reference), mini = screening")
    p.add_argument("--solver", choices=["lu", "minres"], help="sparse LU (default) or MINRES")
    p.add_argument("--no-openfoam", action="store_true", help="disable the OpenFOAM cross-check (recorded in the manifest)")
    p.add_argument("--openfoam-env", help="script that sets up OpenFOAM (etc/bashrc)")
    p.add_argument("--openfoam-nprocs", type=int, help="MPI ranks for simpleFoam")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="percoflow", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"percoflow {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("image", help="run the workflow on a segmented image")
    p.add_argument("image")
    _common(p)
    p.add_argument("--level", type=int, choices=sorted(MESH_LEVELS), help="mesh level preset (1 coarse ... 5 fine)")
    p.add_argument("--tolerance", type=float, help="contour tolerance in pixels (default 0.25)")
    p.add_argument("--band", type=float, help="boundary band fraction (default 0.08)")
    p = sub.add_parser("mesh", help="run the solver on a Gmsh mesh with groups INLET, OUTLET, WALLS")
    p.add_argument("mesh")
    _common(p)
    p.add_argument("--length", type=float, help="sample length L for k (default: bounding box along the flow)")
    p.add_argument("--width", type=float, help="transverse width W for k (default: bounding box across the flow)")
    p = sub.add_parser("benchmark", help="verification and validation suite (OpenFOAM included)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--full", action="store_true", help="finer meshes and more cases (slower)")
    p.add_argument("--no-openfoam", action="store_true", help="skip the OpenFOAM benchmarks (report marked INCOMPLETE)")
    p.add_argument("--openfoam-env")
    sub.add_parser("openfoam-check", help="locate OpenFOAM and print its version").add_argument("--openfoam-env")
    sub.add_parser("verify", help="recompute the output hashes of a run directory").add_argument("run_dir")
    sub.add_parser("levels", help="print the mesh-level presets")
    args = ap.parse_args(argv)

    if args.cmd == "image":
        from .pipeline import run_image
        r = run_image(args.image, args.out, _config(args))
        _report(r)
    elif args.cmd == "mesh":
        from .pipeline import run_mesh
        r = run_mesh(args.mesh, args.out, _config(args), args.length, args.width)
        _report(r)
    elif args.cmd == "benchmark":
        from .benchmarks import run_all
        from .config import OpenFOAMOptions
        r = run_all(args.out, skip_openfoam=args.no_openfoam, of_opts=OpenFOAMOptions(env_script=args.openfoam_env), quick=not args.full)
        print(json.dumps({k: v.get("pass") for k, v in r["results"].items()}, indent=1))
        print("STATUS:", r["status"])
        return 0 if r["status"] == "PASS" else 1
    elif args.cmd == "openfoam-check":
        from .benchmarks.openfoam import find_openfoam
        of = find_openfoam(args.openfoam_env)
        print(json.dumps(of.describe(), indent=1))
    elif args.cmd == "verify":
        from .provenance import verify_manifest
        r = verify_manifest(args.run_dir)
        print(json.dumps(r, indent=1))
        return 0 if r["ok"] else 1
    elif args.cmd == "levels":
        for k, v in MESH_LEVELS.items():
            print(f"level {k}: h_min={v[0]} h_max={v[1]} d_min={v[2]} d_max={v[3]} px")
    return 0


def _report(r: dict):
    """Short console summary of a run (the full result is in result.json)."""
    if r.get("status") != "ok":
        print("status:", r.get("status"))
        return
    k = r["permeability"]
    print(f"k_D = {k['k_D']:.6g}   (Q = {k['Q']:.6g}, dP = {k['dP']:.6g}, L = {k['L_s']}, W = {k['W']})")
    of = r.get("openfoam_crosscheck", {})
    if of.get("status") == "completed":
        if "summary" in of:
            print(f"OpenFOAM cross-check: {of['scope']}, median |FV/FE - 1| = {of['summary']['median_abs_rel_diff']:.3%}")
        else:
            print(f"OpenFOAM cross-check: FV/FE - 1 = {of['rel_diff_fv_vs_fe']:+.3%}")
    else:
        print("OpenFOAM cross-check:", of.get("status"))


if __name__ == "__main__":
    sys.exit(main())
