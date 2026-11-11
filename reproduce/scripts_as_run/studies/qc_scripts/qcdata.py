"""
Data access for the QC figures. Reads ONLY the author's raw outputs (global_summary.csv, solution.vtk, full_percolating_domain.npy,
structures/*.png) and the QC evidence JSON files. No result is typed by hand.
"""
from __future__ import annotations
import csv, json
from functools import lru_cache
from pathlib import Path
import numpy as np

QC = Path(__file__).resolve().parents[1]
PROJ = QC.parents[1] / "convergencia"
EVID = QC / "evidence"
LEVELS = ["convergence_mesh_01_coarse", "convergence_mesh_02_intermediate", "convergence_mesh_03_fine",
          "convergence_mesh_04_very_fine", "convergence_mesh_05_extra_fine"]
CASES = {"R3-H": dict(direction="lr", patch="patch_y7600_x53200", tag="R3H", sample="R3"),
         "R2-V": dict(direction="tb", patch="patch_y7600_x19000", tag="R2V", sample="R2"),
         "R3-V": dict(direction="tb", patch="patch_y7600_x53200", tag="R3V", sample="R3")}
SAMPLES = {"R1": "patch_y3800_x3800", "R2": "patch_y7600_x19000", "R3": "patch_y7600_x53200"}
U0 = 1.0e-4


@lru_cache(None)
def summary(case: str, level: int) -> dict:
    c = CASES[case]
    for r in csv.DictReader(open(PROJ / LEVELS[level] / "global_summary.csv")):
        if r["filename"].startswith(c["patch"]) and r["direction"] == c["direction"]:
            return {k: (float(v) if _isnum(v) else v) for k, v in r.items()}
    raise KeyError((case, level))


def _isnum(s: str) -> bool:
    try:
        float(s); return True
    except ValueError:
        return False


def series(case: str, key: str) -> np.ndarray:
    return np.array([summary(case, i)[key] for i in range(5)], float)


@lru_cache(None)
def field(case: str, level: int = 4) -> dict:
    """Vertex fields of the extra-fine solution (author's VTK)."""
    import meshio
    c = CASES[case]
    m = meshio.read(PROJ / LEVELS[level] / c["direction"] / (c["patch"] + "_c0_mask") / "solution.vtk")
    s = summary(case, level)
    P = m.points[:, :2]; T = m.cells[0].data
    p = m.point_data["pressure"]; sp = m.point_data["speed"]
    a = 0.5 * np.abs((P[T[:, 1], 0] - P[T[:, 0], 0]) * (P[T[:, 2], 1] - P[T[:, 0], 1])
                     - (P[T[:, 1], 1] - P[T[:, 0], 1]) * (P[T[:, 2], 0] - P[T[:, 0], 0]))
    return dict(P=P, T=T, p=p, speed=sp, vel=m.point_data["velocity"][:, :2], area=a,
                pstar=(p - s["pressure_outlet_mean"]) / s["pressure_drop"], ustar=sp / U0, summary=s)


@lru_cache(None)
def structure(sample: str) -> np.ndarray:
    from PIL import Image
    im = np.asarray(Image.open(PROJ / "github/structures" / (SAMPLES[sample] + "_c0_mask.png")).convert("RGB"))
    return (im[..., 1] > 128) & (im[..., 0] < 100)      # pores are pure green in the author's masks


@lru_cache(None)
def domain(sample: str, direction: str) -> np.ndarray | None:
    p = PROJ / LEVELS[4] / direction / (SAMPLES[sample] + "_c0_mask") / "full_percolating_domain.npy"
    return np.load(p) if p.exists() else None


def block_mode(cls: np.ndarray, b: int, n: int) -> np.ndarray:
    """Down-sample an integer class image by b with majority vote (display only)."""
    h, w = cls.shape
    c = cls[: h // b * b, : w // b * b].reshape(h // b, b, w // b, b).transpose(0, 2, 1, 3).reshape(h // b, w // b, b * b)
    counts = np.stack([(c == k).sum(-1) for k in range(n)], -1)
    return counts.argmax(-1)


def load_json(name: str):
    return json.load(open(EVID / name))
