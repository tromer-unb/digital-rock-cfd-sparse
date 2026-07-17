#!/usr/bin/env python3
"""
Generate lightweight CFD visualizations from the output produced by
``digital_rock_cfd_sparse.py``.

Expected directory layout
-------------------------
The CFD solver and this visualization script may be stored in the same folder:

    digital_rock_cfd_sparse.py
    visualize_cfd_outputs.py
    structures/
    results_complete_lr/
        lr/
            sample_name/
                solution.vtk
                summary.csv

The script automatically searches the current directory for folders whose names
start with ``results_complete``. Every ``solution.vtk`` found below those folders
is processed. A specific result folder can also be supplied with
``--results-root``.

Generated images
----------------
For every CFD case, the script creates:

    normalized_speed.png
    pressure.png
    velocity_vectors.png
    streamlines.png
    flow_zones.png
    overview.png

Memory and file-size control
----------------------------
The original CFD mesh is not plotted directly when it is too large. Instead,
the script applies topology-preserving vertex clustering and then limits the
number of triangles. Streamlines use a regular grid with a configurable maximum
resolution, while velocity vectors use spatially distributed arrow sampling.

Dependencies
------------
    pip install numpy matplotlib meshio

Examples
--------
Automatic discovery in the current directory:

    python3 visualize_cfd_outputs.py

Read one result folder explicitly:

    python3 visualize_cfd_outputs.py \
        --results-root results_complete_lr

Use stronger reduction for very large meshes:

    python3 visualize_cfd_outputs.py \
        --max-points 30000 \
        --max-triangles 60000 \
        --grid-resolution 180 \
        --max-arrows 900
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Patch
import meshio
import numpy as np


@dataclass(frozen=True)
class CasePaths:
    """Paths and identifiers for one CFD result case."""

    results_root: Path
    direction: str
    sample_name: str
    sample_dir: Path
    vtk_path: Path
    summary_path: Path


@dataclass
class CaseData:
    """Reduced mesh and fields used for visualization."""

    case: CasePaths
    points: np.ndarray
    triangles: np.ndarray
    velocity: np.ndarray
    speed: np.ndarray
    normalized_velocity: np.ndarray
    normalized_speed: np.ndarray
    pressure: np.ndarray | None
    active_threshold: float
    dead_threshold: float
    original_point_count: int
    original_triangle_count: int
    reduced_point_count: int
    reduced_triangle_count: int
    inlet_velocity: float


IMAGE_NAMES = (
    "normalized_speed",
    "pressure",
    "velocity_vectors",
    "streamlines",
    "flow_zones",
    "overview",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read digital-rock CFD output folders and generate lightweight "
            "velocity, pressure, streamline, and flow-zone images."
        )
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("."),
        help=(
            "Directory containing results_complete* folders. "
            "The default is the current directory."
        ),
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        action="append",
        default=None,
        help=(
            "Result folder to process. This option may be repeated. "
            "When omitted, results_complete* folders are discovered automatically."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("cfd_visualizations"),
        help="Directory where all generated images are written.",
    )
    parser.add_argument(
        "--direction",
        choices=("auto", "lr", "tb", "both"),
        default="auto",
        help="Optional flow-direction filter.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=None,
        help=(
            "Optional sample-name filter. This option may be repeated and "
            "matches complete sample directory names."
        ),
    )
    parser.add_argument(
        "--format",
        choices=("png", "svg", "pdf"),
        default="png",
        help="Output image format.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="Image resolution.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=60000,
        help=(
            "Maximum number of mesh vertices retained for plotting. "
            "Vertex clustering is used when the original mesh is larger."
        ),
    )
    parser.add_argument(
        "--max-triangles",
        type=int,
        default=120000,
        help="Maximum number of triangles retained for plotting.",
    )
    parser.add_argument(
        "--grid-resolution",
        type=int,
        default=220,
        help=(
            "Number of regular-grid points along the largest dimension for "
            "streamline interpolation."
        ),
    )
    parser.add_argument(
        "--seed-density",
        type=int,
        default=24,
        help="Base number of streamline seeds placed near the inlet.",
    )
    parser.add_argument(
        "--streamline-density",
        type=float,
        default=1.05,
        help="Density passed to matplotlib.streamplot.",
    )
    parser.add_argument(
        "--max-arrows",
        type=int,
        default=1600,
        help="Maximum number of velocity arrows per vector figure.",
    )
    parser.add_argument(
        "--speed-percentile",
        type=float,
        default=99.5,
        help=(
            "Percentile used as the upper color limit for normalized speed. "
            "This avoids a few extreme values dominating the color scale."
        ),
    )
    parser.add_argument(
        "--no-mesh-lines",
        action="store_true",
        help="Do not draw the reduced triangular mesh over scalar fields.",
    )

    args = parser.parse_args()

    if args.dpi < 50:
        parser.error("--dpi must be at least 50.")
    if args.max_points < 100:
        parser.error("--max-points must be at least 100.")
    if args.max_triangles < 100:
        parser.error("--max-triangles must be at least 100.")
    if args.grid_resolution < 30:
        parser.error("--grid-resolution must be at least 30.")
    if args.seed_density < 1:
        parser.error("--seed-density must be positive.")
    if args.streamline_density <= 0:
        parser.error("--streamline-density must be positive.")
    if args.max_arrows < 1:
        parser.error("--max-arrows must be positive.")
    if not (50.0 <= args.speed_percentile <= 100.0):
        parser.error("--speed-percentile must be between 50 and 100.")

    return args


def resolve_path(path: Path, base_dir: Path) -> Path:
    """Resolve a user path relative to the selected base directory."""
    if path.is_absolute():
        return path.resolve()
    return (base_dir / path).resolve()


def discover_results_roots(
    base_dir: Path,
    requested_roots: list[Path] | None,
) -> list[Path]:
    """Return explicit or automatically discovered CFD result roots."""
    base_dir = base_dir.resolve()

    if requested_roots:
        roots = [resolve_path(path, base_dir) for path in requested_roots]
    else:
        roots = sorted(
            path.resolve()
            for path in base_dir.glob("results_complete*")
            if path.is_dir()
        )

    valid_roots: list[Path] = []

    for root in roots:
        if not root.is_dir():
            print(f"Warning: result folder not found: {root}", file=sys.stderr)
            continue

        if not any(root.rglob("solution.vtk")):
            print(
                f"Warning: no solution.vtk files were found below {root}",
                file=sys.stderr,
            )
            continue

        valid_roots.append(root)

    return valid_roots


def read_summary_row(summary_path: Path) -> dict[str, str]:
    """Read the first row of summary.csv, returning an empty mapping if absent."""
    if not summary_path.is_file():
        return {}

    with summary_path.open("r", newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        return next(reader, {})


def get_summary_float(
    summary: dict[str, str],
    key: str,
    default: float,
) -> float:
    """Convert one summary value to float with a safe fallback."""
    raw_value = summary.get(key)

    if raw_value is None or raw_value == "":
        return default

    try:
        value = float(raw_value)
    except (TypeError, ValueError):
        return default

    return value if np.isfinite(value) else default


def infer_direction(sample_dir: Path, results_root: Path) -> str:
    """Infer lr or tb from the path, falling back to summary.csv."""
    try:
        relative_parts = sample_dir.relative_to(results_root).parts
    except ValueError:
        relative_parts = sample_dir.parts

    for part in reversed(relative_parts):
        if part in {"lr", "tb"}:
            return part

    summary = read_summary_row(sample_dir / "summary.csv")
    direction = summary.get("direction", "unknown").strip().lower()

    return direction if direction in {"lr", "tb"} else "unknown"


def discover_cases(
    results_roots: Iterable[Path],
    direction_filter: str,
    sample_filter: list[str] | None,
) -> list[CasePaths]:
    """Discover every CFD sample containing a solution.vtk file."""
    allowed_directions: set[str] | None

    if direction_filter == "lr":
        allowed_directions = {"lr"}
    elif direction_filter == "tb":
        allowed_directions = {"tb"}
    elif direction_filter == "both":
        allowed_directions = {"lr", "tb"}
    else:
        allowed_directions = None

    selected_samples = set(sample_filter or [])
    cases: list[CasePaths] = []

    for results_root in results_roots:
        for vtk_path in sorted(results_root.rglob("solution.vtk")):
            sample_dir = vtk_path.parent
            sample_name = sample_dir.name
            direction = infer_direction(sample_dir, results_root)

            if allowed_directions is not None and direction not in allowed_directions:
                continue
            if selected_samples and sample_name not in selected_samples:
                continue

            cases.append(
                CasePaths(
                    results_root=results_root,
                    direction=direction,
                    sample_name=sample_name,
                    sample_dir=sample_dir,
                    vtk_path=vtk_path,
                    summary_path=sample_dir / "summary.csv",
                )
            )

    return sorted(
        cases,
        key=lambda case: (
            case.results_root.name.lower(),
            case.direction,
            case.sample_name.lower(),
        ),
    )


def get_triangles(mesh: meshio.Mesh) -> np.ndarray:
    """Collect linear triangular cells from a meshio mesh."""
    blocks = [
        block.data
        for block in mesh.cells
        if block.type.startswith("triangle")
    ]

    if not blocks:
        raise RuntimeError("No triangular elements were found in the VTK file.")

    triangles = np.vstack(blocks)

    if triangles.shape[1] > 3:
        triangles = triangles[:, :3]

    return np.asarray(triangles, dtype=np.int64)


def read_velocity(mesh: meshio.Mesh) -> np.ndarray:
    """Read the two in-plane velocity components."""
    if "velocity" not in mesh.point_data:
        raise KeyError("The VTK file does not contain point field 'velocity'.")

    values = np.asarray(mesh.point_data["velocity"], dtype=float)

    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(f"Unexpected velocity shape: {values.shape}")

    return values[:, :2]


def read_optional_scalar(mesh: meshio.Mesh, name: str) -> np.ndarray | None:
    """Read an optional scalar point field."""
    if name not in mesh.point_data:
        return None

    values = np.asarray(mesh.point_data[name], dtype=float)
    values = np.squeeze(values).reshape(-1)
    return values


def aggregate_field(
    inverse: np.ndarray,
    counts: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    """Average one scalar or vector field inside clustering bins."""
    values = np.asarray(values, dtype=float)

    if values.ndim == 1:
        sums = np.bincount(inverse, weights=values, minlength=len(counts))
        return sums / counts

    columns = []

    for component in range(values.shape[1]):
        sums = np.bincount(
            inverse,
            weights=values[:, component],
            minlength=len(counts),
        )
        columns.append(sums / counts)

    return np.column_stack(columns)


def compact_mesh(
    points: np.ndarray,
    triangles: np.ndarray,
    fields: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Remove vertices that are no longer referenced by any triangle."""
    used = np.unique(triangles.reshape(-1))
    mapping = np.full(points.shape[0], -1, dtype=np.int64)
    mapping[used] = np.arange(len(used), dtype=np.int64)

    compact_triangles = mapping[triangles]
    compact_points = points[used]
    compact_fields = {
        name: np.asarray(values)[used]
        for name, values in fields.items()
    }

    return compact_points, compact_triangles, compact_fields


def cluster_mesh_once(
    points: np.ndarray,
    triangles: np.ndarray,
    fields: dict[str, np.ndarray],
    target_points: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Cluster nearby vertices once without deleting isolated triangles."""
    x = points[:, 0]
    y = points[:, 1]
    width = max(float(np.ptp(x)), np.finfo(float).eps)
    height = max(float(np.ptp(y)), np.finfo(float).eps)
    aspect = width / height

    bins_x = max(2, int(round(math.sqrt(target_points * aspect))))
    bins_y = max(2, int(math.floor(target_points / bins_x)))

    x_scaled = (x - x.min()) / width
    y_scaled = (y - y.min()) / height

    ix = np.clip((x_scaled * bins_x).astype(int), 0, bins_x - 1)
    iy = np.clip((y_scaled * bins_y).astype(int), 0, bins_y - 1)
    cluster_key = iy * bins_x + ix

    _, inverse = np.unique(cluster_key, return_inverse=True)
    counts = np.bincount(inverse).astype(float)

    clustered_x = np.bincount(inverse, weights=x) / counts
    clustered_y = np.bincount(inverse, weights=y) / counts
    clustered_points = np.column_stack([clustered_x, clustered_y])

    clustered_fields = {
        name: aggregate_field(inverse, counts, values)
        for name, values in fields.items()
    }

    clustered_triangles = inverse[triangles]
    nondegenerate = (
        (clustered_triangles[:, 0] != clustered_triangles[:, 1])
        & (clustered_triangles[:, 1] != clustered_triangles[:, 2])
        & (clustered_triangles[:, 0] != clustered_triangles[:, 2])
    )
    clustered_triangles = clustered_triangles[nondegenerate]

    triangle_keys = np.sort(clustered_triangles, axis=1)
    _, unique_indices = np.unique(
        triangle_keys,
        axis=0,
        return_index=True,
    )
    clustered_triangles = clustered_triangles[np.sort(unique_indices)]

    return compact_mesh(
        clustered_points,
        clustered_triangles,
        clustered_fields,
    )


def coarsen_mesh(
    points: np.ndarray,
    triangles: np.ndarray,
    fields: dict[str, np.ndarray],
    max_points: int,
    max_triangles: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """
    Reduce the mesh using topology-preserving vertex clustering.

    The method never removes isolated triangles from an otherwise valid mesh,
    because doing so would create white gaps in the rendered porous domain.
    Instead, clustering is repeated with progressively larger spatial bins until
    both visualization limits are satisfied or no further safe reduction is
    possible.
    """
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int64)
    reduced_fields = {
        name: np.asarray(values)
        for name, values in fields.items()
    }

    minimum_target = 40
    previous_size = (points.shape[0], triangles.shape[0])

    for _ in range(12):
        point_count = points.shape[0]
        triangle_count = triangles.shape[0]

        if point_count <= max_points and triangle_count <= max_triangles:
            break

        target_from_points = min(max_points, point_count - 1)

        if triangle_count > max_triangles:
            ratio = max_triangles / max(triangle_count, 1)
            target_from_triangles = int(point_count * ratio * 0.72)
        else:
            target_from_triangles = target_from_points

        target_points = max(
            minimum_target,
            min(target_from_points, target_from_triangles),
        )

        if target_points >= point_count:
            target_points = max(minimum_target, int(point_count * 0.75))

        points, triangles, reduced_fields = cluster_mesh_once(
            points,
            triangles,
            reduced_fields,
            target_points,
        )

        current_size = (points.shape[0], triangles.shape[0])

        if current_size == previous_size:
            if points.shape[0] <= minimum_target:
                break

            points, triangles, reduced_fields = cluster_mesh_once(
                points,
                triangles,
                reduced_fields,
                max(minimum_target, int(points.shape[0] * 0.60)),
            )
            current_size = (points.shape[0], triangles.shape[0])

            if current_size == previous_size:
                break

        previous_size = current_size

    points, triangles, reduced_fields = compact_mesh(
        points,
        triangles,
        reduced_fields,
    )

    return points, triangles, reduced_fields


def load_case_data(
    case: CasePaths,
    max_points: int,
    max_triangles: int,
) -> CaseData:
    """Load and reduce one CFD result case."""
    mesh = meshio.read(case.vtk_path)
    points = np.asarray(mesh.points[:, :2], dtype=float)
    triangles = get_triangles(mesh)
    velocity = read_velocity(mesh)
    pressure = read_optional_scalar(mesh, "pressure")
    speed_from_file = read_optional_scalar(mesh, "speed")

    if velocity.shape[0] != points.shape[0]:
        raise ValueError(
            "Velocity length does not match the number of mesh points."
        )

    speed = (
        speed_from_file
        if speed_from_file is not None
        else np.linalg.norm(velocity, axis=1)
    )

    if speed.shape[0] != points.shape[0]:
        raise ValueError("Speed length does not match the number of mesh points.")
    if pressure is not None and pressure.shape[0] != points.shape[0]:
        raise ValueError(
            "Pressure length does not match the number of mesh points."
        )

    summary = read_summary_row(case.summary_path)
    default_normalizer = max(
        float(np.nanpercentile(np.abs(speed), 95.0)),
        np.finfo(float).eps,
    )
    inlet_velocity = abs(
        get_summary_float(summary, "inlet_velocity", default_normalizer)
    )

    if inlet_velocity <= np.finfo(float).eps:
        inlet_velocity = default_normalizer

    speed_max = max(float(np.nanmax(speed)), np.finfo(float).eps)
    active_threshold = get_summary_float(
        summary,
        "active_threshold",
        0.01 * speed_max,
    )
    dead_threshold = get_summary_float(
        summary,
        "dead_threshold",
        0.001 * speed_max,
    )

    fields: dict[str, np.ndarray] = {
        "velocity": velocity,
        "speed": speed,
    }

    if pressure is not None:
        fields["pressure"] = pressure

    original_point_count = points.shape[0]
    original_triangle_count = triangles.shape[0]

    points, triangles, fields = coarsen_mesh(
        points=points,
        triangles=triangles,
        fields=fields,
        max_points=max_points,
        max_triangles=max_triangles,
    )

    velocity = fields["velocity"]
    speed = fields["speed"]
    pressure = fields.get("pressure")

    return CaseData(
        case=case,
        points=points,
        triangles=triangles,
        velocity=velocity,
        speed=speed,
        normalized_velocity=velocity / inlet_velocity,
        normalized_speed=speed / inlet_velocity,
        pressure=pressure,
        active_threshold=active_threshold,
        dead_threshold=dead_threshold,
        original_point_count=original_point_count,
        original_triangle_count=original_triangle_count,
        reduced_point_count=points.shape[0],
        reduced_triangle_count=triangles.shape[0],
        inlet_velocity=inlet_velocity,
    )


def make_triangulation(data: CaseData) -> mtri.Triangulation:
    """Build a Matplotlib triangulation for the reduced mesh."""
    return mtri.Triangulation(
        data.points[:, 0],
        data.points[:, 1],
        data.triangles,
    )


def set_clean_axis(
    ax: plt.Axes,
    data: CaseData,
    title: str,
) -> None:
    """Apply common plot formatting."""
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=13, fontweight="bold", pad=7)

    direction_label = {
        "lr": "left to right",
        "tb": "top to bottom",
    }.get(data.case.direction, "unknown direction")

    ax.text(
        0.012,
        0.018,
        direction_label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.78,
        },
        zorder=30,
    )


def draw_mesh_lines(
    ax: plt.Axes,
    triangulation: mtri.Triangulation,
    enabled: bool,
) -> None:
    """Draw faint reduced-mesh lines when requested."""
    if enabled:
        ax.triplot(
            triangulation,
            linewidth=0.08,
            alpha=0.13,
            color="black",
            zorder=3,
        )


def normalized_speed_limit(data: CaseData, percentile: float) -> float:
    """Return a robust upper color limit for normalized speed."""
    finite = data.normalized_speed[np.isfinite(data.normalized_speed)]

    if finite.size == 0:
        return 1.0

    return max(float(np.percentile(finite, percentile)), 1.0)


def pressure_limits(pressure: np.ndarray) -> tuple[float, float]:
    """Return robust pressure color limits."""
    finite = pressure[np.isfinite(pressure)]

    if finite.size == 0:
        return 0.0, 1.0

    vmin, vmax = np.percentile(finite, [1.0, 99.0])

    if np.isclose(vmin, vmax):
        delta = max(abs(float(vmin)) * 0.05, 1.0)
        return float(vmin - delta), float(vmax + delta)

    return float(vmin), float(vmax)


def choose_spatial_arrow_indices(
    points: np.ndarray,
    speed: np.ndarray,
    max_arrows: int,
) -> np.ndarray:
    """Select spatially distributed vector locations, favoring larger speeds."""
    count = points.shape[0]

    if count <= max_arrows:
        return np.arange(count, dtype=np.int64)

    x = points[:, 0]
    y = points[:, 1]
    width = max(float(np.ptp(x)), np.finfo(float).eps)
    height = max(float(np.ptp(y)), np.finfo(float).eps)
    aspect = width / height

    bins_x = max(2, int(round(math.sqrt(max_arrows * aspect))))
    bins_y = max(2, int(math.floor(max_arrows / bins_x)))

    ix = np.clip(
        (((x - x.min()) / width) * bins_x).astype(int),
        0,
        bins_x - 1,
    )
    iy = np.clip(
        (((y - y.min()) / height) * bins_y).astype(int),
        0,
        bins_y - 1,
    )
    keys = iy * bins_x + ix

    order = np.lexsort((-speed, keys))
    ordered_keys = keys[order]
    first_in_bin = np.r_[True, ordered_keys[1:] != ordered_keys[:-1]]
    selected = order[first_in_bin]

    if selected.size > max_arrows:
        selected = selected[:max_arrows]

    return np.sort(selected)


def make_regular_grid(
    data: CaseData,
    grid_resolution: int,
) -> dict[str, np.ndarray | np.ma.MaskedArray]:
    """Interpolate reduced velocity onto a regular grid for streamlines."""
    x = data.points[:, 0]
    y = data.points[:, 1]
    width = max(float(np.ptp(x)), np.finfo(float).eps)
    height = max(float(np.ptp(y)), np.finfo(float).eps)

    if width >= height:
        nx = grid_resolution
        ny = max(35, int(round(grid_resolution * height / width)))
    else:
        ny = grid_resolution
        nx = max(35, int(round(grid_resolution * width / height)))

    x_grid = np.linspace(float(x.min()), float(x.max()), nx)
    y_grid = np.linspace(float(y.min()), float(y.max()), ny)
    xx, yy = np.meshgrid(x_grid, y_grid)

    triangulation = make_triangulation(data)
    interp_u = mtri.LinearTriInterpolator(
        triangulation,
        data.normalized_velocity[:, 0],
    )
    interp_v = mtri.LinearTriInterpolator(
        triangulation,
        data.normalized_velocity[:, 1],
    )

    u_grid = np.ma.masked_invalid(interp_u(xx, yy))
    v_grid = np.ma.masked_invalid(interp_v(xx, yy))
    combined_mask = np.ma.getmaskarray(u_grid) | np.ma.getmaskarray(v_grid)
    u_grid = np.ma.array(u_grid, mask=combined_mask)
    v_grid = np.ma.array(v_grid, mask=combined_mask)
    speed_grid = np.ma.sqrt(u_grid**2 + v_grid**2)

    return {
        "x_grid": x_grid,
        "y_grid": y_grid,
        "xx": xx,
        "yy": yy,
        "u_grid": u_grid,
        "v_grid": v_grid,
        "speed_grid": speed_grid,
    }


def generate_streamline_seeds(
    data: CaseData,
    grid: dict[str, np.ndarray | np.ma.MaskedArray],
    seed_density: int,
) -> np.ndarray | None:
    """Generate valid streamline seeds near the CFD inlet."""
    x_grid = np.asarray(grid["x_grid"])
    y_grid = np.asarray(grid["y_grid"])
    u_grid = np.ma.asarray(grid["u_grid"])
    v_grid = np.ma.asarray(grid["v_grid"])
    seeds: list[list[float]] = []

    if data.case.direction == "lr":
        candidate_columns = range(min(4, len(x_grid)))
        y_candidates = np.linspace(y_grid.min(), y_grid.max(), seed_density)

        for y_value in y_candidates:
            y_index = int(np.argmin(np.abs(y_grid - y_value)))

            for x_index in candidate_columns:
                if np.ma.is_masked(u_grid[y_index, x_index]):
                    continue
                if np.ma.is_masked(v_grid[y_index, x_index]):
                    continue

                magnitude = float(
                    np.hypot(
                        u_grid[y_index, x_index],
                        v_grid[y_index, x_index],
                    )
                )

                if magnitude > 1e-9:
                    seeds.append([x_grid[x_index], y_grid[y_index]])
                    break

    elif data.case.direction == "tb":
        candidate_rows = range(len(y_grid) - 1, max(-1, len(y_grid) - 5), -1)
        x_candidates = np.linspace(x_grid.min(), x_grid.max(), seed_density)

        for x_value in x_candidates:
            x_index = int(np.argmin(np.abs(x_grid - x_value)))

            for y_index in candidate_rows:
                if np.ma.is_masked(u_grid[y_index, x_index]):
                    continue
                if np.ma.is_masked(v_grid[y_index, x_index]):
                    continue

                magnitude = float(
                    np.hypot(
                        u_grid[y_index, x_index],
                        v_grid[y_index, x_index],
                    )
                )

                if magnitude > 1e-9:
                    seeds.append([x_grid[x_index], y_grid[y_index]])
                    break

    if not seeds:
        return None

    return np.asarray(seeds, dtype=float)


def draw_streamlines(
    ax: plt.Axes,
    data: CaseData,
    grid_resolution: int,
    seed_density: int,
    streamline_density: float,
    speed_vmax: float,
) -> None:
    """Draw an interpolated speed background and lightweight streamlines."""
    grid = make_regular_grid(data, grid_resolution)
    x_grid = np.asarray(grid["x_grid"])
    y_grid = np.asarray(grid["y_grid"])
    speed_grid = np.ma.asarray(grid["speed_grid"])
    u_grid = np.ma.asarray(grid["u_grid"])
    v_grid = np.ma.asarray(grid["v_grid"])

    ax.pcolormesh(
        x_grid,
        y_grid,
        speed_grid,
        shading="auto",
        vmin=0.0,
        vmax=speed_vmax,
        rasterized=True,
        zorder=1,
    )

    seeds = generate_streamline_seeds(data, grid, seed_density)
    streamplot_kwargs = {
        "density": streamline_density,
        "color": "white",
        "linewidth": 0.85,
        "arrowsize": 0.85,
        "maxlength": 7.0,
        "integration_direction": "forward",
        "zorder": 8,
    }

    if seeds is not None:
        streamplot_kwargs["start_points"] = seeds

    try:
        ax.streamplot(
            x_grid,
            y_grid,
            u_grid,
            v_grid,
            **streamplot_kwargs,
        )
    except ValueError as exc:
        ax.text(
            0.5,
            0.5,
            f"Streamlines unavailable\n{exc}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            zorder=20,
        )


def plot_speed(
    ax: plt.Axes,
    data: CaseData,
    triangulation: mtri.Triangulation,
    speed_vmax: float,
    show_mesh_lines: bool,
):
    artist = ax.tripcolor(
        triangulation,
        data.normalized_speed,
        shading="gouraud",
        vmin=0.0,
        vmax=speed_vmax,
        rasterized=True,
        zorder=1,
    )
    draw_mesh_lines(ax, triangulation, show_mesh_lines)
    set_clean_axis(ax, data, r"Normalized speed, $|\mathbf{u}|/U_0$")
    return artist


def plot_pressure(
    ax: plt.Axes,
    data: CaseData,
    triangulation: mtri.Triangulation,
    show_mesh_lines: bool,
):
    if data.pressure is None:
        ax.text(
            0.5,
            0.5,
            "Pressure field not found",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
        )
        set_clean_axis(ax, data, "Pressure")
        return None

    vmin, vmax = pressure_limits(data.pressure)
    artist = ax.tripcolor(
        triangulation,
        data.pressure,
        shading="gouraud",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
        zorder=1,
    )
    draw_mesh_lines(ax, triangulation, show_mesh_lines)
    set_clean_axis(ax, data, "Pressure")
    return artist


def plot_vectors(
    ax: plt.Axes,
    data: CaseData,
    triangulation: mtri.Triangulation,
    speed_vmax: float,
    max_arrows: int,
    show_mesh_lines: bool,
):
    background = ax.tripcolor(
        triangulation,
        data.normalized_speed,
        shading="gouraud",
        vmin=0.0,
        vmax=speed_vmax,
        rasterized=True,
        alpha=0.80,
        zorder=1,
    )
    indices = choose_spatial_arrow_indices(
        data.points,
        data.normalized_speed,
        max_arrows,
    )

    ax.quiver(
        data.points[indices, 0],
        data.points[indices, 1],
        data.normalized_velocity[indices, 0],
        data.normalized_velocity[indices, 1],
        color="black",
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0018,
        headwidth=3.2,
        headlength=4.2,
        zorder=10,
    )
    draw_mesh_lines(ax, triangulation, show_mesh_lines)
    set_clean_axis(ax, data, "Velocity vectors")
    return background


def plot_zones(
    ax: plt.Axes,
    data: CaseData,
    triangulation: mtri.Triangulation,
    show_mesh_lines: bool,
):
    triangle_speed = np.mean(data.speed[data.triangles], axis=1)
    zones = np.zeros(triangle_speed.shape[0], dtype=float)
    zones[triangle_speed <= data.dead_threshold] = -1.0
    zones[triangle_speed >= data.active_threshold] = 1.0

    artist = ax.tripcolor(
        triangulation,
        facecolors=zones,
        shading="flat",
        cmap="coolwarm",
        vmin=-1.0,
        vmax=1.0,
        rasterized=True,
        zorder=1,
    )
    draw_mesh_lines(ax, triangulation, show_mesh_lines)
    set_clean_axis(ax, data, "Active, intermediate, and dead flow zones")

    ax.legend(
        handles=[
            Patch(facecolor=plt.cm.coolwarm(1.0), label="Active"),
            Patch(facecolor=plt.cm.coolwarm(0.5), label="Intermediate"),
            Patch(facecolor=plt.cm.coolwarm(0.0), label="Dead"),
        ],
        loc="upper right",
        framealpha=0.85,
        fontsize=8,
    )
    return artist


def add_case_heading(fig: plt.Figure, data: CaseData) -> None:
    """Add sample and reduction information above a figure."""
    direction = data.case.direction.upper()
    fig.suptitle(
        f"{data.case.sample_name} | direction={direction}",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )


def save_single_plot(
    output_path: Path,
    data: CaseData,
    plot_callback,
    colorbar_label: str | None,
    dpi: int,
) -> None:
    """Save one single-panel figure."""
    fig, ax = plt.subplots(figsize=(10.5, 8.5), constrained_layout=True)
    artist = plot_callback(ax)
    add_case_heading(fig, data)

    if artist is not None and colorbar_label:
        colorbar = fig.colorbar(artist, ax=ax, fraction=0.045, pad=0.025)
        colorbar.set_label(colorbar_label, fontsize=11, fontweight="bold")

    fig.savefig(
        output_path,
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)


def create_case_images(
    data: CaseData,
    output_dir: Path,
    image_format: str,
    dpi: int,
    grid_resolution: int,
    seed_density: int,
    streamline_density: float,
    max_arrows: int,
    speed_percentile: float,
    show_mesh_lines: bool,
) -> list[Path]:
    """Create all requested visualizations for one CFD case."""
    output_dir.mkdir(parents=True, exist_ok=True)
    triangulation = make_triangulation(data)
    speed_vmax = normalized_speed_limit(data, speed_percentile)
    outputs: list[Path] = []

    speed_path = output_dir / f"normalized_speed.{image_format}"
    save_single_plot(
        speed_path,
        data,
        lambda ax: plot_speed(
            ax,
            data,
            triangulation,
            speed_vmax,
            show_mesh_lines,
        ),
        r"$|\mathbf{u}|/U_0$",
        dpi,
    )
    outputs.append(speed_path)

    pressure_path = output_dir / f"pressure.{image_format}"
    save_single_plot(
        pressure_path,
        data,
        lambda ax: plot_pressure(
            ax,
            data,
            triangulation,
            show_mesh_lines,
        ),
        "Pressure",
        dpi,
    )
    outputs.append(pressure_path)

    vectors_path = output_dir / f"velocity_vectors.{image_format}"
    save_single_plot(
        vectors_path,
        data,
        lambda ax: plot_vectors(
            ax,
            data,
            triangulation,
            speed_vmax,
            max_arrows,
            show_mesh_lines,
        ),
        r"$|\mathbf{u}|/U_0$",
        dpi,
    )
    outputs.append(vectors_path)

    streamlines_path = output_dir / f"streamlines.{image_format}"
    fig, ax = plt.subplots(figsize=(10.5, 8.5), constrained_layout=True)
    draw_streamlines(
        ax,
        data,
        grid_resolution,
        seed_density,
        streamline_density,
        speed_vmax,
    )
    draw_mesh_lines(ax, triangulation, show_mesh_lines)
    set_clean_axis(ax, data, "Velocity field and streamlines")
    add_case_heading(fig, data)
    fig.savefig(
        streamlines_path,
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)
    outputs.append(streamlines_path)

    zones_path = output_dir / f"flow_zones.{image_format}"
    save_single_plot(
        zones_path,
        data,
        lambda ax: plot_zones(
            ax,
            data,
            triangulation,
            show_mesh_lines,
        ),
        None,
        dpi,
    )
    outputs.append(zones_path)

    overview_path = output_dir / f"overview.{image_format}"
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15.0, 12.0),
        constrained_layout=True,
    )

    speed_artist = plot_speed(
        axes[0, 0],
        data,
        triangulation,
        speed_vmax,
        show_mesh_lines,
    )
    pressure_artist = plot_pressure(
        axes[0, 1],
        data,
        triangulation,
        show_mesh_lines,
    )
    plot_vectors(
        axes[1, 0],
        data,
        triangulation,
        speed_vmax,
        max_arrows,
        show_mesh_lines,
    )
    draw_streamlines(
        axes[1, 1],
        data,
        grid_resolution,
        seed_density,
        streamline_density,
        speed_vmax,
    )
    draw_mesh_lines(axes[1, 1], triangulation, show_mesh_lines)
    set_clean_axis(axes[1, 1], data, "Velocity field and streamlines")

    speed_colorbar = fig.colorbar(
        speed_artist,
        ax=[axes[0, 0], axes[1, 0], axes[1, 1]],
        fraction=0.025,
        pad=0.015,
    )
    speed_colorbar.set_label(r"$|\mathbf{u}|/U_0$", fontweight="bold")

    if pressure_artist is not None:
        pressure_colorbar = fig.colorbar(
            pressure_artist,
            ax=axes[0, 1],
            fraction=0.045,
            pad=0.02,
        )
        pressure_colorbar.set_label("Pressure", fontweight="bold")

    add_case_heading(fig, data)
    fig.savefig(
        overview_path,
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)
    outputs.append(overview_path)

    return outputs


def output_directory_for_case(
    main_output_dir: Path,
    case: CasePaths,
) -> Path:
    """Mirror the result-root, direction, and sample hierarchy."""
    direction = case.direction if case.direction != "unknown" else "unknown_direction"
    return main_output_dir / case.results_root.name / direction / case.sample_name


def write_manifest(rows: list[dict[str, object]], output_path: Path) -> None:
    """Write a CSV summary of all processed visualization cases."""
    if not rows:
        return

    fieldnames: list[str] = []

    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def configure_matplotlib() -> None:
    """Apply compact, publication-friendly defaults."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        }
    )


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    output_dir = resolve_path(args.output_dir, base_dir)

    roots = discover_results_roots(base_dir, args.results_root)

    if not roots:
        print(
            "No valid results_complete* folders were found. "
            "Use --results-root to provide one explicitly.",
            file=sys.stderr,
        )
        return 1

    cases = discover_cases(
        roots,
        direction_filter=args.direction,
        sample_filter=args.sample,
    )

    if not cases:
        print("No CFD cases matched the selected filters.", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    print(f"Result roots: {len(roots)}")
    print(f"CFD cases found: {len(cases)}")
    print(f"Visualization output: {output_dir}")
    print(
        "Reduction limits: "
        f"points={args.max_points:,}, "
        f"triangles={args.max_triangles:,}, "
        f"grid={args.grid_resolution}, "
        f"arrows={args.max_arrows:,}"
    )

    manifest_rows: list[dict[str, object]] = []
    success_count = 0

    for index, case in enumerate(cases, start=1):
        print(
            f"\n[{index}/{len(cases)}] "
            f"{case.results_root.name}/{case.direction}/{case.sample_name}"
        )

        case_output_dir = output_directory_for_case(output_dir, case)

        try:
            data = load_case_data(
                case,
                max_points=args.max_points,
                max_triangles=args.max_triangles,
            )

            print(
                "  Mesh reduction: "
                f"{data.original_point_count:,} -> {data.reduced_point_count:,} points, "
                f"{data.original_triangle_count:,} -> "
                f"{data.reduced_triangle_count:,} triangles"
            )

            outputs = create_case_images(
                data=data,
                output_dir=case_output_dir,
                image_format=args.format,
                dpi=args.dpi,
                grid_resolution=args.grid_resolution,
                seed_density=args.seed_density,
                streamline_density=args.streamline_density,
                max_arrows=args.max_arrows,
                speed_percentile=args.speed_percentile,
                show_mesh_lines=not args.no_mesh_lines,
            )

            success_count += 1
            print(f"  Images written: {len(outputs)}")
            print(f"  Output folder: {case_output_dir}")

            manifest_rows.append(
                {
                    "results_root": str(case.results_root),
                    "direction": case.direction,
                    "sample": case.sample_name,
                    "status": "success",
                    "vtk_path": str(case.vtk_path),
                    "output_dir": str(case_output_dir),
                    "original_points": data.original_point_count,
                    "reduced_points": data.reduced_point_count,
                    "original_triangles": data.original_triangle_count,
                    "reduced_triangles": data.reduced_triangle_count,
                    "image_count": len(outputs),
                }
            )

        except Exception as exc:
            print(f"  Failed: {exc}", file=sys.stderr)
            manifest_rows.append(
                {
                    "results_root": str(case.results_root),
                    "direction": case.direction,
                    "sample": case.sample_name,
                    "status": "failed",
                    "vtk_path": str(case.vtk_path),
                    "output_dir": str(case_output_dir),
                    "error": str(exc),
                }
            )

    manifest_path = output_dir / "visualization_manifest.csv"
    write_manifest(manifest_rows, manifest_path)

    print("\nVisualization pipeline finished.")
    print(f"Successful cases: {success_count}/{len(cases)}")
    print(f"Manifest: {manifest_path}")

    return 0 if success_count > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
