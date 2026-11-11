"""Triangle meshes: generated from a pore mask with Gmsh, or read from a user file.

Mesh contract (for both routes): a 2D triangle mesh with the physical groups

    FLUID   (surfaces)  the fluid domain
    INLET   (curves)    where the pressure p_in (or the velocity u0) is imposed
    OUTLET  (curves)    traction-free outlet
    WALLS   (curves)    no-slip boundaries (pore walls and lateral edges)

Size field of generated meshes (docs/THEORY.md, Section 4): with d the distance to the nearest wall,
r = clip((d - d_min)/(d_max - d_min), 0, 1) and s = 1/(1 + exp(-(12 r - 6))), the target size is
h = h_min (1 - s) + h_max s. Small elements sit in throats and next to walls, large ones in pore bodies.
"""
from __future__ import annotations

from pathlib import Path

import gmsh
import numpy as np
from skfem import MeshTri

from .config import MeshOptions
from .geometry import contour_hierarchy, extract_contours

REQUIRED_GROUPS = ("INLET", "OUTLET", "WALLS")


def _loop(points, point_cache, endpoints, h):
    """Add one closed polygon to the Gmsh model as a curve loop of straight lines. Points shared by two polygons are
    created once (``point_cache``); the end points of every line are kept in ``endpoints`` for the boundary labelling."""
    pts = points[:-1]
    tags = []
    for x, y in pts:
        key = (int(round(x * 1e6)), int(round(y * 1e6)))
        if key not in point_cache:
            point_cache[key] = gmsh.model.geo.addPoint(float(x), float(y), 0.0, h)
        tags.append(point_cache[key])
    curves = []
    for k, a in enumerate(tags):
        b = tags[(k + 1) % len(tags)]
        if a == b:
            continue
        c = gmsh.model.geo.addLine(a, b)
        curves.append(c)
        endpoints[c] = (pts[k].copy(), pts[(k + 1) % len(pts)].copy())
    if len(curves) < 3:
        raise RuntimeError("degenerate contour after simplification")
    return gmsh.model.geo.addCurveLoop(curves)


def _classify(endpoints, width, height, direction, tol):
    """Label each line INLET (on the inlet face), OUTLET (on the outlet face) or WALLS (everything else)."""
    groups = {"INLET": [], "OUTLET": [], "WALLS": []}
    xl, xr, yb, yt = -0.5, width - 0.5, -0.5, height - 0.5
    for tag, (p0, p1) in endpoints.items():
        if direction == "lr":
            inlet = abs(p0[0] - xl) <= tol and abs(p1[0] - xl) <= tol
            outlet = abs(p0[0] - xr) <= tol and abs(p1[0] - xr) <= tol
        else:
            inlet = abs(p0[1] - yt) <= tol and abs(p1[1] - yt) <= tol
            outlet = abs(p0[1] - yb) <= tol and abs(p1[1] - yb) <= tol
        groups["INLET" if inlet else "OUTLET" if outlet else "WALLS"].append(tag)
    return groups


def mesh_from_mask(mask: np.ndarray, direction: str, outdir, opts: MeshOptions | None = None, verbose: bool = False) -> dict:
    """Mesh the pore space of ``mask`` (True = fluid). The first/last column (lr) or row (tb) must be the inlet/outlet.

    Writes ``outdir/mesh.msh`` (Gmsh 4.1) and ``outdir/geometry.geo_unrolled`` and returns a summary dict.
    """
    opts = opts or MeshOptions()
    h_min, h_max, d_min, d_max = opts.sizes()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    contours = extract_contours(mask, opts.contour_tolerance, opts.min_contour_points, opts.min_contour_area)
    if not contours:
        raise RuntimeError("no valid contour: the mask has no pore region large enough to mesh")
    depth, children = contour_hierarchy(contours)
    height, width = mask.shape
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    try:
        gmsh.model.add("pore_domain")
        cache, endpoints, loops = {}, {}, {}
        for i, c in enumerate(contours):
            loops[i] = _loop(c, cache, endpoints, h_max)
        surfaces = []
        # even depth = outer boundary of a fluid region; its children of the next depth are holes (solid inclusions)
        for i in range(len(contours)):
            if depth[i] % 2:
                continue                                    # odd depth: a solid inclusion, used as a hole
            holes = [loops[k] for k in children[i] if depth[k] == depth[i] + 1]
            surfaces.append(gmsh.model.geo.addPlaneSurface([loops[i], *holes]))
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, surfaces, name="FLUID")
        groups = _classify(endpoints, width, height, direction, opts.face_tolerance)
        for name, tags in groups.items():
            if tags:
                gmsh.model.addPhysicalGroup(1, tags, name=name)
        if not groups["INLET"] or not groups["OUTLET"]:
            raise RuntimeError("no INLET or OUTLET curve found: does the domain touch both planes?")
        # background size field: h(d) from the distance d to the nearest wall, sigmoid between d_min and d_max
        f = gmsh.model.mesh.field
        dist = f.add("Distance")
        f.setNumbers(dist, "CurvesList", sorted(endpoints))
        f.setNumber(dist, "Sampling", 100)
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        for key, val in (("SizeMin", h_min), ("SizeMax", h_max), ("DistMin", d_min), ("DistMax", d_max), ("Sigmoid", 1)):
            f.setNumber(thr, key, val)
        f.setAsBackgroundMesh(thr)
        for key, val in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                         ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", opts.algorithm),
                         ("Mesh.Smoothing", opts.smoothing), ("Mesh.MshFileVersion", 4.1)):
            gmsh.option.setNumber(key, val)
        gmsh.model.mesh.generate(2)
        mesh_path = outdir / "mesh.msh"
        gmsh.write(str(mesh_path))
        gmsh.write(str(outdir / "geometry.geo_unrolled"))
        n_nodes = len(gmsh.model.mesh.getNodes()[0])
    finally:
        gmsh.finalize()
    return {"mesh_path": str(mesh_path), "contours": len(contours), "surfaces": len(surfaces), "nodes": int(n_nodes),
            "inlet_curves": len(groups["INLET"]), "outlet_curves": len(groups["OUTLET"]), "wall_curves": len(groups["WALLS"]),
            "sizes": {"h_min": h_min, "h_max": h_max, "d_min": d_min, "d_max": d_max}}


def load_mesh(path) -> MeshTri:
    """Read a triangle mesh and check the mesh contract (physical groups INLET, OUTLET, WALLS)."""
    mesh = MeshTri.load(str(path))
    missing = [g for g in REQUIRED_GROUPS if g not in (mesh.boundaries or {})]
    if missing:
        raise ValueError(f"{path}: missing physical group(s) {missing}. Name the boundary curves INLET, OUTLET and WALLS "
                         "(see docs/INPUTS_OUTPUTS.md); the fluid surface may be named FLUID.")
    return mesh


def mesh_summary(mesh: MeshTri) -> dict:
    """Vertices, triangles, fluid area, bounding box and number of boundary facets per group."""
    a = mesh.p[:, mesh.t]
    area = 0.5 * np.abs((a[0, 1] - a[0, 0]) * (a[1, 2] - a[1, 0]) - (a[1, 1] - a[1, 0]) * (a[0, 2] - a[0, 0]))
    return {"n_vertices": int(mesh.p.shape[1]), "n_triangles": int(mesh.t.shape[1]), "area": float(area.sum()),
            "bbox": [float(mesh.p[0].min()), float(mesh.p[0].max()), float(mesh.p[1].min()), float(mesh.p[1].max())],
            "boundary_facets": {k: int(len(v)) for k, v in (mesh.boundaries or {}).items()}}
