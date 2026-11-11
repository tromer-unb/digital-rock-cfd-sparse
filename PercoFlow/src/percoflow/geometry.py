"""Reservoirs and polygonal geometry of the pore domain.

Steps (docs/THEORY.md, Sections 3 and 4):

1. ``pad_reservoirs``: append fully open strips of width b = round(f * L_s) before the inlet and after the outlet. The
   pressure is imposed on their outer faces and measured one pixel inside them (``permeability.py``), so the inlet
   singularities of a velocity condition never touch the sample.
2. ``extract_contours``: marching squares at level 1/2 on the zero-padded mask, pores fully connected (consistent with
   8-neighbour percolation). Straight pixel walls are reproduced exactly and pixel corners are cut diagonally.
3. Douglas-Peucker simplification with tolerance eps: no wall moves by more than eps. A throat of width w has a cubic
   conductance, so its relative error is up to ~6 eps / w; keep eps well below one pixel (0.25 px is the converged default).
4. ``contour_hierarchy``: nested contours are organized by an even-odd depth rule (depth 0 = outer boundary of a fluid
   region, depth 1 = a solid inclusion inside it, depth 2 = fluid inside that inclusion, ...).
"""
from __future__ import annotations

import numpy as np
from matplotlib.path import Path as MplPath
from skimage import measure


def pad_reservoirs(selected: np.ndarray, direction: str, fraction: float):
    """Return (padded_mask, L_s, b): L_s = sample length along the flow, b = reservoir width in pixels."""
    Ls = selected.shape[1] if direction == "lr" else selected.shape[0]
    b = int(round(fraction * Ls))
    if b == 0:
        return selected.copy(), Ls, 0
    if direction == "lr":
        pad = np.ones((selected.shape[0], b), dtype=bool)
        return np.concatenate([pad, selected, pad], axis=1), Ls, b
    pad = np.ones((b, selected.shape[1]), dtype=bool)
    return np.concatenate([pad, selected, pad], axis=0), Ls, b


def signed_area(points: np.ndarray) -> float:
    """Shoelace formula for a closed polygon (first point repeated); positive for counter-clockwise order."""
    x, y = points[:, 0], points[:, 1]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def _closed(points: np.ndarray) -> np.ndarray:
    """Append the first point at the end if the polygon is not closed yet."""
    if len(points) >= 2 and not np.allclose(points[0], points[-1]):
        points = np.vstack([points, points[0]])
    return points


def _dedup(points: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    """Drop consecutive duplicate points (closer than ``tol``) and close the polygon again."""
    if len(points) <= 2:
        return points
    keep = [0]
    for k in range(1, len(points)):
        if np.linalg.norm(points[k] - points[keep[-1]]) > tol:
            keep.append(k)
    return _closed(points[keep])


def extract_contours(mask: np.ndarray, tolerance: float, min_points: int = 8, min_area: float = 20.0):
    """Closed polygons (x, y in pixel-centre coordinates, first point repeated at the end) of the pore/solid interface.

    ``min_points`` counts the repeated closing point (8 points = 7 distinct vertices); ``min_area`` is in px^2. Contours
    below either limit are dropped; they are tiny solid specks or pore fragments that cannot be meshed sensibly.
    """
    height = mask.shape[0]
    # pad with solid so that every contour is closed, also where the pore space touches the image edge
    padded = np.pad(mask.astype(np.float32), 1, constant_values=0.0)
    out = []
    for rc in measure.find_contours(padded, level=0.5, fully_connected="high", positive_orientation="low"):
        row, col = rc[:, 0] - 1.0, rc[:, 1] - 1.0                  # undo the padding
        pts = np.column_stack([col, (height - 1.0) - row])         # (row, col) -> (x, y), y pointing up
        if tolerance > 0:
            pts = measure.approximate_polygon(pts, tolerance=tolerance)   # Douglas-Peucker
        pts = _closed(_dedup(pts))
        if len(pts) < min_points or abs(signed_area(pts)) < min_area:
            continue
        out.append(pts.astype(np.float64))
    return out


def _interior_point(points: np.ndarray) -> np.ndarray:
    """A point strictly inside a closed polygon (used to test nesting)."""
    path = MplPath(points)
    core = points[:-1] if np.allclose(points[0], points[-1]) else points
    centroid = core.mean(axis=0)
    if path.contains_point(centroid):
        return centroid
    for v in core[:: max(1, len(core) // 25)]:
        cand = 0.75 * v + 0.25 * centroid
        if path.contains_point(cand):
            return cand
    p0, p1 = points[0], points[1]
    t = p1 - p0
    n = np.array([-t[1], t[0]]) / max(np.linalg.norm(t), 1e-300)
    for eps in (0.1, 0.5, 1.0, 2.0):
        for s in (-1.0, 1.0):
            cand = 0.5 * (p0 + p1) + s * eps * n
            if path.contains_point(cand):
                return cand
    return points[0]


def contour_hierarchy(contours):
    """Return (depth, children): depth[i] = number of contours that enclose contour i; children[i] = contours whose
    smallest enclosing contour is i."""
    n = len(contours)
    paths = [MplPath(c) for c in contours]
    areas = [abs(signed_area(c)) for c in contours]
    probes = [_interior_point(c) for c in contours]
    containers = [[j for j in range(n) if j != i and areas[j] > areas[i] and paths[j].contains_point(probes[i])] for i in range(n)]
    depth = [len(c) for c in containers]
    children = [[] for _ in range(n)]
    for i, cont in enumerate(containers):
        if cont:
            children[min(cont, key=lambda j: areas[j])].append(i)
    return depth, children
