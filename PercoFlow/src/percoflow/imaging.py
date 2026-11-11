"""From an image file to a boolean pore mask, and from the mask to the analysed window.

Coordinate convention used everywhere in the package: the pixel in row i and column j of an array with N_y rows has its
centre at (x, y) = (j, N_y - 1 - i). Pixel edges therefore lie on half-integers, and the outer faces of an
N_y x N_x window are x = -0.5, x = N_x - 0.5, y = -0.5 and y = N_y - 0.5.

Flow directions:
    "lr"  inlet = first column (left),  outlet = last column (right),  flow along +x
    "tb"  inlet = first row (top),       outlet = last row (bottom),    flow along -y
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image
from skimage.morphology import remove_small_objects

from .config import DomainOptions, ImageOptions


def load_mask(path, opts: ImageOptions | None = None) -> np.ndarray:
    """Read an image (PNG, TIFF, JPEG, BMP) or a NumPy array (.npy, boolean or 0/1) and return the pore mask.

    The rule is chosen by ``opts.pore_rule`` (see ImageOptions). For the "green" rule an image with no green pixel at
    all falls back to the mid-range threshold, which is reported by ``describe_mask``; check it when using new data.
    """
    opts = opts or ImageOptions()
    path = Path(path)
    if path.suffix.lower() == ".npy":
        mask = np.load(path).astype(bool)
    else:
        image = Image.open(path)
        if opts.pore_rule == "green" and image.mode not in ("1", "L", "I", "F", "I;16"):     # colour image
            rgb = np.asarray(image.convert("RGB"), dtype=np.int16)
            r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
            mask = (g >= opts.green_min) & (g >= r + opts.green_delta) & (g >= b + opts.green_delta)
            if not mask.any():
                mask = _threshold(np.asarray(image.convert("L"), dtype=np.float64), None)
        elif opts.pore_rule == "nonzero":
            mask = np.asarray(image.convert("L")) > 0
        else:
            mask = _threshold(np.asarray(image.convert("L"), dtype=np.float64), opts.threshold)
    if opts.invert:
        mask = ~mask
    return np.ascontiguousarray(mask, dtype=bool)


def _threshold(gray: np.ndarray, value) -> np.ndarray:
    """Pore = grey value above ``value`` (None: the mid-range of the image)."""
    t = 0.5 * (gray.min() + gray.max()) if value is None else float(value)
    return gray > t


def clean_mask(mask: np.ndarray, opts: ImageOptions | None = None) -> np.ndarray:
    """Remove pore objects smaller than ``min_object_size`` pixels (8-connectivity). No closing or hole filling is
    applied, so the pore geometry is otherwise preserved."""
    opts = opts or ImageOptions()
    if opts.min_object_size > 1:
        try:   # scikit-image >= 0.26 removes objects of size <= max_size
            return remove_small_objects(mask.astype(bool), max_size=opts.min_object_size - 1, connectivity=2)
        except TypeError:   # older scikit-image: objects of size < min_size
            return remove_small_objects(mask.astype(bool), min_size=opts.min_object_size, connectivity=2)
    return mask.astype(bool).copy()


def band_width(mask: np.ndarray, opts: DomainOptions) -> int:
    """Width in pixels of the band removed at each flow-direction side."""
    side = mask.shape[1] if opts.direction == "lr" else mask.shape[0]
    band = opts.band_px if opts.band_px is not None else int(round(side * opts.band_fraction))
    if band <= 0:
        return 0
    return max(1, min(int(band), side // 2 - 1))


def crop_window(mask: np.ndarray, opts: DomainOptions):
    """Remove the inlet- and outlet-side bands. The new edges are the inlet and outlet planes; the bands avoid artificial
    connections along the image border and make the planes identical for every sample. Returns (window, info)."""
    band = band_width(mask, opts)
    if opts.direction == "lr":
        window = mask[:, band: mask.shape[1] - band]
    else:
        window = mask[band: mask.shape[0] - band, :]
    return np.ascontiguousarray(window), {"band_px": int(band), "window_shape": list(window.shape)}


def describe_mask(mask: np.ndarray) -> dict:
    """Shape, number of pore pixels and image porosity (pore fraction of the whole image)."""
    return {"shape": list(mask.shape), "pore_pixels": int(mask.sum()), "porosity": float(mask.mean())}
