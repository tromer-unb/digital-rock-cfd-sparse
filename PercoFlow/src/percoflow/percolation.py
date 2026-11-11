"""Percolation analysis: which pore clusters can carry flow between the inlet and outlet planes.

Why this step is not optional (docs/THEORY.md, Section 2):

* An isolated cluster (touching neither plane) admits "constant pressure, zero velocity" as a non-trivial solution of the
  homogeneous Stokes problem, so every isolated cluster adds one vector to the kernel of the saddle-point matrix and the
  discrete system is singular.
* A dead end (touching the inlet but not the outlet) receives a prescribed inflow under a velocity inlet with nowhere to
  discharge it, so the velocity-driven problem has no solution.

Keeping only clusters that join both planes therefore makes the problem well posed and removes unknowns that cannot carry
flow. Connectivity is 8-neighbour, matching the contour extraction (marching squares with fully connected pores).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

_EIGHT = np.ones((3, 3), dtype=np.uint8)


def _planes(labels: np.ndarray, direction: str):
    """Labels on the inlet plane and on the outlet plane (first/last column for "lr", first/last row for "tb")."""
    if direction == "lr":
        return labels[:, 0], labels[:, -1]
    return labels[0, :], labels[-1, :]


def analyse(window: np.ndarray, direction: str, keep: str = "all"):
    """Label the pore clusters of ``window`` and keep the percolating ones.

    Returns (selected_mask, report). The report counts percolating, dead-end and isolated clusters, gives their areas,
    and the connected pore fraction chi_c = (area kept) / (pore area of the window).
    """
    labels, count = ndi.label(window, structure=_EIGHT)                   # 0 = solid, 1..count = pore clusters
    inlet, outlet = _planes(labels, direction)
    at_in, at_out = set(np.unique(inlet)) - {0}, set(np.unique(outlet)) - {0}
    areas = np.bincount(labels.ravel(), minlength=count + 1)              # pixels per cluster
    perc = sorted(at_in & at_out)                                         # touch both planes
    dead = sorted((at_in ^ at_out))                                       # touch exactly one plane
    isolated = sorted(set(range(1, count + 1)) - at_in - at_out)          # touch neither
    if keep == "largest" and perc:
        perc = [max(perc, key=lambda k: areas[k])]
    selected = np.isin(labels, perc) if perc else np.zeros_like(window, dtype=bool)
    pore_area = int(window.sum())
    report = {
        "clusters_total": int(count),
        "clusters_percolating": len(perc),
        "clusters_dead_end": len(dead),
        "clusters_isolated": len(isolated),
        "area_percolating_px": int(selected.sum()),
        "area_dead_end_px": int(areas[dead].sum()) if dead else 0,
        "area_isolated_px": int(areas[isolated].sum()) if isolated else 0,
        "pore_area_window_px": pore_area,
        "connected_fraction": float(selected.sum() / pore_area) if pore_area else 0.0,
        "percolates": bool(perc),
    }
    return selected, report


def percolates(mask: np.ndarray, direction: str) -> bool:
    """True if at least one 8-connected pore cluster joins the inlet and the outlet planes."""
    labels, _ = ndi.label(mask, structure=_EIGHT)
    a, b = _planes(labels, direction)
    return bool((set(np.unique(a)) - {0}) & (set(np.unique(b)) - {0}))


def critical_throat_width(selected: np.ndarray, direction: str, tol: float = 0.01) -> float:
    """Critical throat width w_c = 2 r_c, with r_c the largest wall distance r for which the pixels farther than r from
    the solid still percolate. Equivalently, twice the largest (over all paths joining the planes) of the smallest wall
    distance along the path: the widest constriction that every path must cross (Katz-Thompson picture).

    The lateral window edges count as solid; the inlet and outlet planes are open (padded with fluid before the distance
    transform so that they do not create artificial walls). Found by bisection to ``tol`` pixels.
    """
    if not percolates(selected, direction):
        return 0.0
    pad = ((0, 0), (1, 1)) if direction == "lr" else ((1, 1), (0, 0))
    open_ends = np.pad(selected, pad, constant_values=True)
    lateral = ((1, 1), (0, 0)) if direction == "lr" else ((0, 0), (1, 1))
    closed = np.pad(open_ends, lateral, constant_values=False)
    dist = ndi.distance_transform_edt(closed)
    core = (slice(1, -1), slice(1, -1))
    dist = dist[core]
    lo, hi = 0.0, float(dist[selected].max())
    # bisection on r: "pixels farther than r from the walls still percolate" is true for r < r_c and false above
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if percolates(selected & (dist > mid), direction):
            lo = mid
        else:
            hi = mid
    return 2.0 * lo
