"""Descriptors and sanity indicators for a 2D digital-rock window (pores True)."""
from __future__ import annotations
import numpy as np
from skimage import measure

UM2_PER_DARCY = 0.986923          # 1 darcy = 0.986923 um^2


def specific_surface(mask):
    """Fluid-solid interface length per unit window area [1/px]. The window edges are not counted (edge-padded before contouring)."""
    pad = np.pad(mask.astype(float), 1, mode="edge")
    L = 0.0
    for c in measure.find_contours(pad, 0.5, fully_connected="high"):
        d = np.diff(c, axis=0); L += float(np.hypot(d[:, 0], d[:, 1]).sum())
    return L / mask.size


def kozeny_carman_K(k_px2, mask):
    """2D Kozeny-Carman number K = phi^3 / (k S^2). Straight parallel channels give exactly K = 3 (k = phi h^2/12 with h = 2 phi / S).
    Tortuous or throat-limited media give K > 3. A value clearly below 3 is a warning sign (heuristic indicator, not a theorem)."""
    phi = float(mask.mean()); S = specific_surface(mask)
    return phi ** 3 / (k_px2 * S ** 2), phi, S


def to_si(k_px2, pixel_size_m):
    """Convert k [px^2] to m^2, um^2/darcy and millidarcy. pixel_size_m MUST be a measured value for real samples."""
    k_m2 = k_px2 * pixel_size_m ** 2
    k_um2 = k_m2 * 1e12
    return dict(k_m2=k_m2, k_darcy=k_um2 / UM2_PER_DARCY, k_mD=1e3 * k_um2 / UM2_PER_DARCY)
