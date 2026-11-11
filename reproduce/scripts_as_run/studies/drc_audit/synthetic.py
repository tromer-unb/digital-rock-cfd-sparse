"""Synthetic digital-rock masks (pores are True). Everything here is SYNTHETIC and must be labelled as such in any output.

cylinder_array_mask: pixelised square array of circular solid cylinders (benchmark with a literature reference).
boolean_disks_mask:  Boolean model of overlapping solid grains of log-normal radius (proof-of-concept rock-like windows with known porosity
                     and known pixel size; the physical scale is an ASSUMPTION of the user of the function).
"""
from __future__ import annotations
import numpy as np
from PIL import Image


def cylinder_array_mask(a=160, c=0.1, M=2, N=4):
    """Pores True. Cell side a [px], solid fraction c (nominal), M cells across the width, N cells along the flow (x). Disk centres at cell centres,
    so both lateral edges are symmetry planes of the array. Returns (mask[H=M*a, L=N*a], R)."""
    R = a * np.sqrt(c / np.pi)
    yy, xx = np.mgrid[0:M * a, 0:N * a] + 0.5
    cx = (np.floor(xx / a) + 0.5) * a; cy = (np.floor(yy / a) + 0.5) * a
    solid = (xx - cx) ** 2 + (yy - cy) ** 2 <= R ** 2
    return ~solid, R


def boolean_disks_mask(shape=(1024, 1024), porosity=0.3, r_mean=14.0, r_sigma=0.35, seed=0, tol=2e-3, max_iter=60):
    """Overlapping solid disks (log-normal radius, median r_mean px) placed as a Poisson process; the intensity is tuned (bisection on a
    scalar multiplier) until the pore fraction is within `tol` of `porosity`. Returns (mask, achieved_porosity)."""
    rng = np.random.default_rng(seed)
    H, W = shape
    area_mean = np.pi * r_mean ** 2 * np.exp(2 * r_sigma ** 2)

    def draw(n):
        cx = rng.uniform(0, W, n); cy = rng.uniform(0, H, n); r = r_mean * np.exp(rng.normal(0, r_sigma, n))
        return cx, cy, r

    def render(n, state):
        solid = np.zeros(shape, bool)
        cx, cy, r = state
        for i in range(n):
            x0, x1 = int(max(cx[i] - r[i] - 1, 0)), int(min(cx[i] + r[i] + 2, W)); y0, y1 = int(max(cy[i] - r[i] - 1, 0)), int(min(cy[i] + r[i] + 2, H))
            yy, xx = np.mgrid[y0:y1, x0:x1] + 0.5
            solid[y0:y1, x0:x1] |= (xx - cx[i]) ** 2 + (yy - cy[i]) ** 2 <= r[i] ** 2
        return solid
    n_max = int(4 * (-np.log(max(porosity, 0.02))) * H * W / area_mean) + 10      # Poisson: phi = exp(-lambda A)
    state = draw(n_max)
    lo, hi = 0, n_max
    best = None
    for _ in range(max_iter):
        mid = (lo + hi) // 2
        phi = 1.0 - render(mid, state).mean()
        best = (mid, phi)
        if abs(phi - porosity) < tol or hi - lo <= 1:
            break
        if phi > porosity: lo = mid
        else: hi = mid
    solid = render(best[0], state)
    return ~solid, float(1.0 - solid.mean())


def save_green_png(mask, path):
    """Author's mask convention: pores pure green (0,255,0), solids black, RGBA."""
    rgba = np.zeros(mask.shape + (4,), np.uint8); rgba[..., 3] = 255; rgba[mask, 1] = 255
    Image.fromarray(rgba, "RGBA").save(path)
