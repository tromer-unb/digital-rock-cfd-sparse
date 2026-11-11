"""Write an OpenFOAM polyMesh for a 2D domain extruded by one cell layer (thickness dz).

Input: 2D vertices (N, 2), cells as an (n_cells, k) array of vertex indices (k = 3 triangles or k = 4 quads), and a function
that names the patch of every boundary edge. Output: constant/polyMesh/{points,faces,owner,neighbour,boundary}.

OpenFOAM conventions implemented here:
* every face lists its points so that the right-hand normal points out of the owner cell (into the neighbour);
* internal faces come first, in upper-triangular order (sorted by owner, then neighbour), owner < neighbour;
* boundary faces follow, grouped by patch; the two layer faces of each cell form the "frontAndBack" patch of type empty,
  which makes the problem two-dimensional.
For a counter-clockwise polygon the outward normal of edge a -> b is to its right, so the quad (a0, b0, b1, a1) built from the
owner's own edge orientation already points outward; the bottom layer face is the reversed polygon (normal -z), the top
layer face the polygon itself (+z).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HEADER = ("FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {cls};\n    location \"constant/polyMesh\";\n"
          "    object {obj};\n}}\n")


def _ccw(points2d: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Reorder each cell counter-clockwise (positive signed area)."""
    xy = points2d[cells]                                    # (n, k, 2)
    x, y = xy[..., 0], xy[..., 1]
    area2 = np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
    out = cells.copy()
    out[area2 < 0] = out[area2 < 0, ::-1]
    return out


def write_polymesh(case_dir, points2d, cells, edge_patch, patch_types, dz: float = 1.0) -> dict:
    """Write the polyMesh. ``edge_patch(midpoints (m,2)) -> array of patch names`` labels boundary edges;
    ``patch_types`` maps patch name -> OpenFOAM type ("patch", "wall"). Returns counts per patch."""
    points2d = np.asarray(points2d, float)
    cells = _ccw(points2d, np.asarray(cells, np.int64))
    n_cells, k = cells.shape
    n_pts = len(points2d)
    # all oriented edges (a -> b) of every cell, in the cell's CCW order
    a = cells.reshape(-1)
    b = np.roll(cells, -1, axis=1).reshape(-1)
    owner_cell = np.repeat(np.arange(n_cells), k)
    key = np.minimum(a, b) * n_pts + np.maximum(a, b)
    order = np.argsort(key, kind="stable")
    key_s = key[order]
    first = np.r_[True, key_s[1:] != key_s[:-1]]
    group = np.cumsum(first) - 1
    counts = np.bincount(group)
    if counts.max() > 2:
        raise ValueError("non-manifold mesh: an edge is shared by more than two cells")
    starts = np.flatnonzero(first)
    # internal edges: two oriented copies; the owner is the smaller cell index and keeps its own orientation
    internal_groups = np.flatnonzero(counts == 2)
    e1, e2 = order[starts[internal_groups]], order[starts[internal_groups] + 1]
    c1, c2 = owner_cell[e1], owner_cell[e2]
    own_edge = np.where(c1 < c2, e1, e2)
    own = np.minimum(c1, c2)
    nei = np.maximum(c1, c2)
    srt = np.lexsort((nei, own))
    own, nei, own_edge = own[srt], nei[srt], own_edge[srt]
    # boundary edges
    bnd_edge = order[starts[counts == 1]]
    mid = 0.5 * (points2d[a[bnd_edge]] + points2d[b[bnd_edge]])
    names = np.asarray(edge_patch(mid))
    # points: layer 0 then layer 1
    pts3 = np.vstack([np.column_stack([points2d, np.zeros(n_pts)]), np.column_stack([points2d, np.full(n_pts, dz)])])

    faces, owners = [], []

    def quad(e):
        return np.column_stack([a[e], b[e], b[e] + n_pts, a[e] + n_pts])
    faces.append(quad(own_edge)); owners.append(own)
    patches = []
    start = len(own)
    for name in sorted(set(names.tolist())):
        sel = bnd_edge[names == name]
        faces.append(quad(sel)); owners.append(owner_cell[sel])
        patches.append((name, patch_types.get(name, "wall"), len(sel), start))
        start += len(sel)
    # frontAndBack: bottom (reversed polygon, normal -z) and top (polygon, normal +z)
    bottom = cells[:, ::-1]
    top = cells + n_pts
    fb = np.vstack([bottom, top]) if k == 4 else np.vstack([bottom, top])
    patches.append(("frontAndBack", "empty", 2 * n_cells, start))
    owners.append(np.r_[np.arange(n_cells), np.arange(n_cells)])

    d = Path(case_dir) / "constant" / "polyMesh"
    d.mkdir(parents=True, exist_ok=True)
    n_faces = sum(len(f) for f in faces) + len(fb)
    with open(d / "points", "w") as fh:
        fh.write(HEADER.format(cls="vectorField", obj="points") + f"\n{len(pts3)}\n(\n")
        np.savetxt(fh, pts3, fmt="(%.10g %.10g %.10g)")
        fh.write(")\n")
    with open(d / "faces", "w") as fh:
        fh.write(HEADER.format(cls="faceList", obj="faces") + f"\n{n_faces}\n(\n")
        for f in faces:
            np.savetxt(fh, f, fmt="4(%d %d %d %d)")
        np.savetxt(fh, fb, fmt=f"{k}(" + " ".join(["%d"] * k) + ")")
        fh.write(")\n")
    note = f"nPoints:{len(pts3)}  nCells:{n_cells}  nFaces:{n_faces}  nInternalFaces:{len(own)}"
    with open(d / "owner", "w") as fh:
        fh.write(HEADER.format(cls="labelList", obj="owner").replace("    object", f"    note \"{note}\";\n    object")
                 + f"\n{n_faces}\n(\n")
        np.savetxt(fh, np.concatenate(owners), fmt="%d")
        fh.write(")\n")
    with open(d / "neighbour", "w") as fh:
        fh.write(HEADER.format(cls="labelList", obj="neighbour").replace("    object", f"    note \"{note}\";\n    object")
                 + f"\n{len(own)}\n(\n")
        np.savetxt(fh, nei, fmt="%d")
        fh.write(")\n")
    with open(d / "boundary", "w") as fh:
        fh.write(HEADER.format(cls="polyBoundaryMesh", obj="boundary") + f"\n{len(patches)}\n(\n")
        for name, typ, n, s in patches:
            extra = "    inGroups List<word> 1(wall);\n" if typ == "wall" else ""
            fh.write(f"    {name}\n    {{\n        type {typ};\n{extra}        nFaces {n};\n        startFace {s};\n    }}\n")
        fh.write(")\n")
    # geometry of the internal faces, used later to find probe lines
    fc = 0.5 * (points2d[a[own_edge]] + points2d[b[own_edge]])
    t = points2d[b[own_edge]] - points2d[a[own_edge]]
    normal = np.column_stack([t[:, 1], -t[:, 0]])
    return {"n_cells": int(n_cells), "n_internal_faces": int(len(own)), "patches": {p[0]: p[2] for p in patches},
            "_internal": {"owner": own, "neighbour": nei, "centre": fc, "normal": normal}}


def pixel_cells(mask: np.ndarray):
    """Vertices and quads of the pore pixels of ``mask`` (one cell per True pixel), in the package's coordinates
    (pixel (i, j) centred at (j, N_y - 1 - i)). Returns (points2d, cells, pixel_index_of_cell)."""
    ny, nx = mask.shape
    ii, jj = np.nonzero(mask)
    node = lambda r, c: r * (nx + 1) + c                      # noqa: E731 - grid node index (row edge r, column edge c)
    cells = np.column_stack([node(ii + 1, jj), node(ii + 1, jj + 1), node(ii, jj + 1), node(ii, jj)])
    used, inv = np.unique(cells, return_inverse=True)
    r, c = used // (nx + 1), used % (nx + 1)
    pts = np.column_stack([c - 0.5, (ny - 0.5) - r]).astype(float)
    return pts, inv.reshape(cells.shape), np.column_stack([ii, jj])
