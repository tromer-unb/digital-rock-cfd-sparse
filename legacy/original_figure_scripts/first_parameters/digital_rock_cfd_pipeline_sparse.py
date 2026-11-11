#!/usr/bin/env python3
"""
digital_rock_cfd_pipeline_sparse.py

Pipeline completo, independente e de memória controlada para CFD 2D em rochas digitais:

    máscaras em structures/
        -> leitura e limpeza
        -> recorte entre planos fixos
        -> seleção do domínio poroso percolante completo
        -> extração dos contornos
        -> geração de malha adaptativa no Gmsh
        -> solução das equações de Stokes com scikit-fem
        -> cálculo de vazão, queda de pressão e permeabilidade
        -> mapas de pressão, velocidade e zonas ativas/mortas
        -> exportação VTK e CSV

Não depende dos scripts Step 01, Step 02, Step 03, Step 04 ou Step 05.

Dependências:
    pip install numpy scipy matplotlib pillow scikit-image gmsh meshio "scikit-fem[all]"

Exemplo horizontal:
    python3 digital_rock_cfd_pipeline.py \
        --structures ../../../poros/poros/ \
        --output results_complete_lr \
        --direction lr \
        --boundary-band-fraction 0.08 \
        --closing-radius 0 \
        --h-min 4 \
        --h-max 30 \
        --distance-min 4 \
        --distance-max 60 \
        --viscosity 1.0 \
        --density 1.0 \
        --inlet-velocity 0.001

Exemplo vertical:
    python3 digital_rock_cfd_pipeline.py \
        --structures ../../../poros/poros/ \
        --output results_complete_tb \
        --direction tb \
        --boundary-band-fraction 0.08 \
        --closing-radius 0 \
        --h-min 4 \
        --h-max 30 \
        --distance-min 4 \
        --distance-max 60 \
        --viscosity 1.0 \
        --density 1.0 \
        --inlet-velocity 0.001

Para testar as duas direções:
    --direction both

Observações:
    1. True/branco representa fluido; False/preto representa matriz.
    2. O fechamento morfológico altera a geometria. O padrão é raio zero.
    3. Com pixel_size=1, a permeabilidade é reportada em pixel².
    4. A condição de entrada é velocidade uniforme prescrita.
    5. O modelo é Stokes estacionário, adequado para Reynolds baixo.
    6. O sistema é resolvido por MINRES esparso, sem fatoração LU densa.
    7. A interpolação da solução é feita em lotes para limitar a memória.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import gmsh
import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.path import Path as MplPath
from scipy import ndimage as ndi
from scipy.sparse.linalg import LinearOperator, minres
from PIL import Image
from skimage import measure
from skimage.morphology import closing, disk, remove_small_objects

from skfem import (
    Basis,
    ElementTriMini,
    ElementTriP1,
    ElementVector,
    FacetBasis,
    MeshTri,
    asm,
    bmat,
)
from skfem.models.general import divergence
from skfem.models.poisson import mass, vector_laplace


IMAGE_EXTENSIONS = {
    ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"
}


@dataclass
class Config:
    structures_dir: Path
    output_dir: Path

    directions: tuple[str, ...]

    green_min: int = 80
    green_delta: int = 30

    min_object_size: int = 20
    fill_holes_smaller_than: int = 0

    boundary_band_px: int = 0
    boundary_band_fraction: float = 0.08

    closing_radius: int = 0
    closing_iterations: int = 1

    keep_percolating: str = "all"
    min_percolating_area_px: int = 0

    pixel_size: float = 1.0
    unit_name: str = "pixel"
    thickness: float = 1.0

    h_min: float = 4.0
    h_max: float = 30.0
    distance_min: float = 4.0
    distance_max: float = 60.0

    contour_tolerance_px: float = 1.5
    min_contour_points: int = 8
    min_contour_area_px2: float = 20.0

    mesh_algorithm: int = 6
    smoothing_steps: int = 5
    msh_version: float = 4.1

    viscosity: float = 1.0
    density: float = 1.0
    inlet_velocity: float = 0.001

    active_threshold_relative: float = 0.01
    dead_threshold_relative: float = 0.001

    solver_tolerance: float = 1e-8
    solver_max_iterations: int = 5000
    interpolation_batch_size: int = 2048
    max_total_dofs: int = 2_000_000

    dpi: int = 220
    max_preview_triangles: int = 250000


def parse_args() -> Config:
    parser = argparse.ArgumentParser(
        description=(
            "Pipeline completo: máscara -> malha Gmsh -> Stokes -> "
            "permeabilidade e diagnóstico."
        )
    )

    parser.add_argument(
        "--structures",
        type=Path,
        default=Path("structures"),
        help="Pasta com as máscaras.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results_complete"),
        help="Pasta principal de saída.",
    )
    parser.add_argument(
        "--direction",
        choices=["lr", "tb", "both"],
        default="lr",
        help="Direção do escoamento.",
    )

    parser.add_argument(
        "--boundary-band-px",
        type=int,
        default=0,
        help=(
            "Largura das faixas removidas em pixels. "
            "Se zero, usa --boundary-band-fraction."
        ),
    )
    parser.add_argument(
        "--boundary-band-fraction",
        type=float,
        default=0.08,
        help="Fração removida em cada extremidade.",
    )

    parser.add_argument(
        "--closing-radius",
        type=int,
        default=0,
        help="Raio do fechamento morfológico; zero desativa.",
    )
    parser.add_argument(
        "--closing-iterations",
        type=int,
        default=1,
        help="Número de aplicações do fechamento.",
    )

    parser.add_argument(
        "--keep-percolating",
        choices=["all", "largest"],
        default="all",
        help="Mantém todos os componentes percolantes ou apenas o maior.",
    )
    parser.add_argument(
        "--min-percolating-area",
        type=int,
        default=0,
        help="Área mínima de cada componente percolante em pixels.",
    )

    parser.add_argument(
        "--min-object-size",
        type=int,
        default=20,
        help="Remove componentes menores que este valor.",
    )
    parser.add_argument(
        "--fill-holes-smaller-than",
        type=int,
        default=0,
        help="Preenche cavidades menores que este valor.",
    )

    parser.add_argument(
        "--pixel-size",
        type=float,
        default=1.0,
        help="Comprimento físico correspondente a um pixel.",
    )
    parser.add_argument(
        "--unit-name",
        type=str,
        default="pixel",
        help="Nome da unidade espacial.",
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=1.0,
        help="Espessura fora do plano.",
    )

    parser.add_argument("--h-min", type=float, default=4.0)
    parser.add_argument("--h-max", type=float, default=30.0)
    parser.add_argument("--distance-min", type=float, default=4.0)
    parser.add_argument("--distance-max", type=float, default=60.0)

    parser.add_argument(
        "--contour-tolerance",
        type=float,
        default=1.5,
        help="Tolerância de simplificação em pixels.",
    )
    parser.add_argument(
        "--min-contour-area",
        type=float,
        default=20.0,
        help="Área mínima dos anéis em pixels².",
    )
    parser.add_argument(
        "--mesh-algorithm",
        type=int,
        default=6,
        help="Algoritmo 2D do Gmsh.",
    )
    parser.add_argument(
        "--smoothing",
        type=int,
        default=5,
        help="Passos de suavização da malha.",
    )

    parser.add_argument(
        "--viscosity",
        type=float,
        default=1.0,
        help="Viscosidade dinâmica.",
    )
    parser.add_argument(
        "--density",
        type=float,
        default=1.0,
        help="Densidade.",
    )
    parser.add_argument(
        "--inlet-velocity",
        type=float,
        default=0.001,
        help="Velocidade uniforme prescrita na entrada.",
    )

    parser.add_argument(
        "--active-threshold-relative",
        type=float,
        default=0.01,
        help="Limiar relativo para zona ativa.",
    )
    parser.add_argument(
        "--dead-threshold-relative",
        type=float,
        default=0.001,
        help="Limiar relativo para zona morta.",
    )
    parser.add_argument(
        "--solver-tolerance",
        type=float,
        default=1e-8,
        help="Tolerância relativa do MINRES esparso.",
    )
    parser.add_argument(
        "--solver-max-iterations",
        type=int,
        default=5000,
        help="Número máximo de iterações do MINRES.",
    )
    parser.add_argument(
        "--interpolation-batch-size",
        type=int,
        default=2048,
        help=(
            "Número de vértices interpolados por lote. "
            "Valores menores reduzem o pico de memória."
        ),
    )
    parser.add_argument(
        "--max-total-dofs",
        type=int,
        default=2000000,
        help=(
            "Limite de segurança para o total de graus de liberdade. "
            "Use 0 para desativar."
        ),
    )
    parser.add_argument("--dpi", type=int, default=220)

    args = parser.parse_args()

    if args.direction == "both":
        directions = ("lr", "tb")
    else:
        directions = (args.direction,)

    if args.boundary_band_px < 0:
        parser.error("--boundary-band-px não pode ser negativo.")
    if not (0.0 < args.boundary_band_fraction < 0.5):
        parser.error("--boundary-band-fraction deve estar entre 0 e 0.5.")
    if args.closing_radius < 0:
        parser.error("--closing-radius não pode ser negativo.")
    if args.closing_iterations < 1:
        parser.error("--closing-iterations deve ser pelo menos 1.")
    if args.min_object_size < 0:
        parser.error("--min-object-size não pode ser negativo.")
    if args.fill_holes_smaller_than < 0:
        parser.error("--fill-holes-smaller-than não pode ser negativo.")
    if args.min_percolating_area < 0:
        parser.error("--min-percolating-area não pode ser negativo.")
    if args.pixel_size <= 0:
        parser.error("--pixel-size deve ser positivo.")
    if args.thickness <= 0:
        parser.error("--thickness deve ser positiva.")
    if args.h_min <= 0 or args.h_max <= args.h_min:
        parser.error("É necessário 0 < h-min < h-max.")
    if args.distance_min < 0 or args.distance_max <= args.distance_min:
        parser.error("É necessário 0 <= distance-min < distance-max.")
    if args.viscosity <= 0:
        parser.error("--viscosity deve ser positiva.")
    if args.density <= 0:
        parser.error("--density deve ser positiva.")
    if args.inlet_velocity == 0:
        parser.error("--inlet-velocity não pode ser zero.")
    if not (0 < args.active_threshold_relative < 1):
        parser.error("--active-threshold-relative deve estar entre 0 e 1.")
    if not (0 <= args.dead_threshold_relative < 1):
        parser.error("--dead-threshold-relative deve estar entre 0 e 1.")
    if args.solver_tolerance <= 0:
        parser.error("--solver-tolerance deve ser positiva.")
    if args.solver_max_iterations < 1:
        parser.error("--solver-max-iterations deve ser positivo.")
    if args.interpolation_batch_size < 1:
        parser.error("--interpolation-batch-size deve ser positivo.")
    if args.max_total_dofs < 0:
        parser.error("--max-total-dofs não pode ser negativo.")

    return Config(
        structures_dir=args.structures,
        output_dir=args.output,
        directions=directions,
        min_object_size=args.min_object_size,
        fill_holes_smaller_than=args.fill_holes_smaller_than,
        boundary_band_px=args.boundary_band_px,
        boundary_band_fraction=args.boundary_band_fraction,
        closing_radius=args.closing_radius,
        closing_iterations=args.closing_iterations,
        keep_percolating=args.keep_percolating,
        min_percolating_area_px=args.min_percolating_area,
        pixel_size=args.pixel_size,
        unit_name=args.unit_name,
        thickness=args.thickness,
        h_min=args.h_min,
        h_max=args.h_max,
        distance_min=args.distance_min,
        distance_max=args.distance_max,
        contour_tolerance_px=args.contour_tolerance,
        min_contour_area_px2=args.min_contour_area,
        mesh_algorithm=args.mesh_algorithm,
        smoothing_steps=args.smoothing,
        viscosity=args.viscosity,
        density=args.density,
        inlet_velocity=args.inlet_velocity,
        active_threshold_relative=args.active_threshold_relative,
        dead_threshold_relative=args.dead_threshold_relative,
        solver_tolerance=args.solver_tolerance,
        solver_max_iterations=args.solver_max_iterations,
        interpolation_batch_size=args.interpolation_batch_size,
        max_total_dofs=args.max_total_dofs,
        dpi=args.dpi,
    )


# ---------------------------------------------------------------------------
# LEITURA E PREPARAÇÃO DA MÁSCARA
# ---------------------------------------------------------------------------

def find_images(folder: Path) -> list[Path]:
    if not folder.exists():
        raise FileNotFoundError(f"A pasta '{folder}' não existe.")

    return sorted(
        [
            path
            for path in folder.iterdir()
            if path.is_file() and path.suffix.lower() in IMAGE_EXTENSIONS
        ],
        key=lambda path: path.name.lower(),
    )


def load_pore_mask(path: Path, cfg: Config) -> np.ndarray:
    image = Image.open(path)

    if image.mode in {"1", "L", "I", "F"}:
        gray = np.asarray(image.convert("L"), dtype=np.float32)
        threshold = 0.5 * (float(gray.min()) + float(gray.max()))
        return (gray > threshold).astype(bool)

    rgb = np.asarray(image.convert("RGB"), dtype=np.int16)

    red = rgb[..., 0]
    green = rgb[..., 1]
    blue = rgb[..., 2]

    pore = (
        (green >= cfg.green_min)
        & (green >= red + cfg.green_delta)
        & (green >= blue + cfg.green_delta)
    )

    if not np.any(pore):
        gray = np.asarray(image.convert("L"), dtype=np.float32)
        threshold = 0.5 * (float(gray.min()) + float(gray.max()))
        pore = gray > threshold

    return pore.astype(bool)


def remove_small_holes_binary(
    mask: np.ndarray,
    max_hole_size: int,
) -> np.ndarray:
    if max_hole_size <= 0:
        return mask.copy()

    labels, count = ndi.label(
        ~mask,
        structure=np.ones((3, 3), dtype=np.uint8),
    )

    border_labels = set(
        np.unique(
            np.concatenate(
                [
                    labels[0, :],
                    labels[-1, :],
                    labels[:, 0],
                    labels[:, -1],
                ]
            )
        )
    )

    counts = np.bincount(labels.ravel())
    result = mask.copy()

    for label_id in range(1, count + 1):
        if label_id in border_labels:
            continue

        if counts[label_id] <= max_hole_size:
            result[labels == label_id] = True

    return result


def clean_mask(mask: np.ndarray, cfg: Config) -> np.ndarray:
    result = mask.astype(bool)

    if cfg.min_object_size > 1:
        try:
            result = remove_small_objects(
                result,
                max_size=cfg.min_object_size - 1,
                connectivity=2,
            )
        except TypeError:
            result = remove_small_objects(
                result,
                min_size=cfg.min_object_size,
                connectivity=2,
            )

    if cfg.fill_holes_smaller_than > 0:
        result = remove_small_holes_binary(
            result,
            cfg.fill_holes_smaller_than,
        )

    return result


def apply_optional_closing(
    mask: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    if cfg.closing_radius <= 0:
        return mask.copy()

    result = mask.copy()
    footprint = disk(cfg.closing_radius)

    for _ in range(cfg.closing_iterations):
        result = closing(result, footprint=footprint)

    return np.asarray(result, dtype=bool)


def get_band_size(
    mask: np.ndarray,
    direction: str,
    cfg: Config,
) -> int:
    dimension = mask.shape[1] if direction == "lr" else mask.shape[0]

    if cfg.boundary_band_px > 0:
        band = cfg.boundary_band_px
    else:
        band = int(round(dimension * cfg.boundary_band_fraction))

    return max(1, min(band, dimension // 2 - 1))


def crop_between_fixed_planes(
    mask: np.ndarray,
    direction: str,
    cfg: Config,
) -> tuple[np.ndarray, dict[str, int]]:
    band = get_band_size(mask, direction, cfg)

    if direction == "lr":
        start = band
        end = mask.shape[1] - band - 1

        if end <= start:
            raise RuntimeError("Faixas horizontais grandes demais.")

        cropped = mask[:, start : end + 1]

        info = {
            "band_px": band,
            "x_offset_px": start,
            "y_offset_px": 0,
            "plane_start_original_px": start,
            "plane_end_original_px": end,
        }

    else:
        start = band
        end = mask.shape[0] - band - 1

        if end <= start:
            raise RuntimeError("Faixas verticais grandes demais.")

        cropped = mask[start : end + 1, :]

        info = {
            "band_px": band,
            "x_offset_px": 0,
            "y_offset_px": start,
            "plane_start_original_px": start,
            "plane_end_original_px": end,
        }

    return cropped, info


def select_full_percolating_components(
    cropped_mask: np.ndarray,
    direction: str,
    cfg: Config,
) -> tuple[np.ndarray, dict[str, object]]:
    labels, count = ndi.label(
        cropped_mask,
        structure=np.ones((3, 3), dtype=np.uint8),
    )

    if count == 0:
        return np.zeros_like(cropped_mask), {
            "candidate_count": 0,
            "selected_count": 0,
        }

    counts = np.bincount(labels.ravel())

    if direction == "lr":
        inlet_labels = set(np.unique(labels[:, 0])) - {0}
        outlet_labels = set(np.unique(labels[:, -1])) - {0}
    else:
        inlet_labels = set(np.unique(labels[0, :])) - {0}
        outlet_labels = set(np.unique(labels[-1, :])) - {0}

    candidates = sorted(inlet_labels.intersection(outlet_labels))

    if cfg.min_percolating_area_px > 0:
        candidates = [
            label_id
            for label_id in candidates
            if counts[label_id] >= cfg.min_percolating_area_px
        ]

    if not candidates:
        return np.zeros_like(cropped_mask), {
            "candidate_count": 0,
            "selected_count": 0,
        }

    if cfg.keep_percolating == "largest":
        selected_labels = [
            max(candidates, key=lambda label_id: counts[label_id])
        ]
    else:
        selected_labels = candidates

    selected = np.isin(labels, selected_labels)

    return selected, {
        "candidate_count": len(candidates),
        "selected_count": len(selected_labels),
        "selected_areas_px": [
            int(counts[label_id])
            for label_id in selected_labels
        ],
    }


# ---------------------------------------------------------------------------
# CONTORNOS E GEOMETRIA GMSH
# ---------------------------------------------------------------------------

def polygon_signed_area(points: np.ndarray) -> float:
    x = points[:, 0]
    y = points[:, 1]

    return 0.5 * float(
        np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    )


def ensure_closed(points: np.ndarray) -> np.ndarray:
    if points.shape[0] >= 2 and not np.allclose(points[0], points[-1]):
        points = np.vstack([points, points[0]])

    return points


def remove_near_duplicate_points(
    points: np.ndarray,
    tolerance: float = 1e-8,
) -> np.ndarray:
    if points.shape[0] <= 2:
        return points

    keep = [0]

    for index in range(1, points.shape[0]):
        if np.linalg.norm(points[index] - points[keep[-1]]) > tolerance:
            keep.append(index)

    return ensure_closed(points[keep])


def extract_contours(
    mask: np.ndarray,
    cfg: Config,
) -> list[np.ndarray]:
    height, _ = mask.shape

    padded = np.pad(
        mask.astype(np.uint8),
        pad_width=1,
        mode="constant",
        constant_values=0,
    )

    raw_contours = measure.find_contours(
        padded.astype(np.float32),
        level=0.5,
        fully_connected="high",
        positive_orientation="low",
    )

    contours: list[np.ndarray] = []

    for contour_rc in raw_contours:
        row = contour_rc[:, 0] - 1.0
        col = contour_rc[:, 1] - 1.0

        x = col
        y = (height - 1.0) - row

        points = np.column_stack([x, y])

        if cfg.contour_tolerance_px > 0:
            points = measure.approximate_polygon(
                points,
                tolerance=cfg.contour_tolerance_px,
            )

        points = remove_near_duplicate_points(points)
        points = ensure_closed(points)

        if points.shape[0] < cfg.min_contour_points:
            continue

        if abs(polygon_signed_area(points)) < cfg.min_contour_area_px2:
            continue

        contours.append(points.astype(np.float64))

    return contours


def find_interior_test_point(points: np.ndarray) -> np.ndarray:
    path = MplPath(points)
    core = points[:-1] if np.allclose(points[0], points[-1]) else points
    centroid = np.mean(core, axis=0)

    if path.contains_point(centroid):
        return centroid

    step = max(1, len(core) // 25)

    for vertex in core[::step]:
        candidate = 0.75 * vertex + 0.25 * centroid

        if path.contains_point(candidate):
            return candidate

    p0 = points[0]
    p1 = points[1]
    midpoint = 0.5 * (p0 + p1)
    tangent = p1 - p0
    norm = np.linalg.norm(tangent)

    if norm > 0:
        normal = np.array([-tangent[1], tangent[0]]) / norm

        for epsilon in [0.1, 0.5, 1.0, 2.0]:
            for sign in [-1.0, 1.0]:
                candidate = midpoint + sign * epsilon * normal

                if path.contains_point(candidate):
                    return candidate

    return points[0]


def classify_contour_hierarchy(
    contours: list[np.ndarray],
) -> tuple[list[int], list[int | None], list[list[int]]]:
    count = len(contours)
    paths = [MplPath(contour) for contour in contours]
    areas = [abs(polygon_signed_area(contour)) for contour in contours]
    test_points = [
        find_interior_test_point(contour)
        for contour in contours
    ]

    containers: list[list[int]] = [[] for _ in range(count)]

    for i in range(count):
        for j in range(count):
            if i == j:
                continue

            if areas[j] <= areas[i]:
                continue

            if paths[j].contains_point(test_points[i]):
                containers[i].append(j)

    depth = [len(values) for values in containers]
    parent: list[int | None] = [None] * count

    for i in range(count):
        if containers[i]:
            parent[i] = min(
                containers[i],
                key=lambda j: areas[j],
            )

    children: list[list[int]] = [[] for _ in range(count)]

    for index, parent_id in enumerate(parent):
        if parent_id is not None:
            children[parent_id].append(index)

    return depth, parent, children


def add_contour_to_gmsh(
    points_px: np.ndarray,
    cfg: Config,
    point_cache: dict[tuple[int, int], int],
    curve_endpoint_map: dict[int, tuple[np.ndarray, np.ndarray]],
) -> int:
    points = points_px[:-1]
    point_tags: list[int] = []
    quantization = 1e6

    for x_px, y_px in points:
        x = float(x_px * cfg.pixel_size)
        y = float(y_px * cfg.pixel_size)

        key = (
            int(round(x * quantization)),
            int(round(y * quantization)),
        )

        if key not in point_cache:
            point_cache[key] = gmsh.model.geo.addPoint(
                x,
                y,
                0.0,
                cfg.h_max,
            )

        point_tags.append(point_cache[key])

    if len(point_tags) < 3:
        raise RuntimeError("Contorno degenerado.")

    curve_tags: list[int] = []

    for index, start_tag in enumerate(point_tags):
        end_tag = point_tags[(index + 1) % len(point_tags)]

        if start_tag == end_tag:
            continue

        curve_tag = gmsh.model.geo.addLine(start_tag, end_tag)
        curve_tags.append(curve_tag)

        curve_endpoint_map[curve_tag] = (
            points[index].copy(),
            points[(index + 1) % len(points)].copy(),
        )

    if len(curve_tags) < 3:
        raise RuntimeError("Contorno degenerado após simplificação.")

    return gmsh.model.geo.addCurveLoop(curve_tags)


def classify_boundary_curves(
    curve_endpoint_map: dict[int, tuple[np.ndarray, np.ndarray]],
    width: int,
    height: int,
    direction: str,
    tolerance_px: float = 0.75,
) -> dict[str, list[int]]:
    groups = {
        "INLET": [],
        "OUTLET": [],
        "WALLS": [],
    }

    x_left = -0.5
    x_right = width - 0.5
    y_bottom = -0.5
    y_top = height - 0.5

    for curve_tag, (p0, p1) in curve_endpoint_map.items():
        x0, y0 = p0
        x1, y1 = p1

        if direction == "lr":
            is_inlet = (
                abs(x0 - x_left) <= tolerance_px
                and abs(x1 - x_left) <= tolerance_px
            )
            is_outlet = (
                abs(x0 - x_right) <= tolerance_px
                and abs(x1 - x_right) <= tolerance_px
            )
        else:
            is_inlet = (
                abs(y0 - y_top) <= tolerance_px
                and abs(y1 - y_top) <= tolerance_px
            )
            is_outlet = (
                abs(y0 - y_bottom) <= tolerance_px
                and abs(y1 - y_bottom) <= tolerance_px
            )

        if is_inlet:
            groups["INLET"].append(curve_tag)
        elif is_outlet:
            groups["OUTLET"].append(curve_tag)
        else:
            groups["WALLS"].append(curve_tag)

    return groups


def create_gmsh_mesh(
    mask: np.ndarray,
    direction: str,
    sample_dir: Path,
    cfg: Config,
) -> dict[str, object]:
    contours = extract_contours(mask, cfg)

    if not contours:
        raise RuntimeError("Nenhum contorno válido foi extraído.")

    depth, _, children = classify_contour_hierarchy(contours)

    height, width = mask.shape

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    try:
        gmsh.model.add("digital_rock_fluid")

        point_cache: dict[tuple[int, int], int] = {}
        curve_endpoint_map: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        loop_tags: dict[int, int] = {}

        for index, contour in enumerate(contours):
            loop_tags[index] = add_contour_to_gmsh(
                contour,
                cfg,
                point_cache,
                curve_endpoint_map,
            )

        surface_tags: list[int] = []

        for index in range(len(contours)):
            if depth[index] % 2 != 0:
                continue

            hole_loops = [
                loop_tags[child]
                for child in children[index]
                if (
                    depth[child] == depth[index] + 1
                    and depth[child] % 2 == 1
                )
            ]

            surface_tags.append(
                gmsh.model.geo.addPlaneSurface(
                    [loop_tags[index], *hole_loops]
                )
            )

        if not surface_tags:
            raise RuntimeError("Nenhuma superfície válida foi criada.")

        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(
            2,
            surface_tags,
            name="FLUID",
        )

        boundary_groups = classify_boundary_curves(
            curve_endpoint_map,
            width,
            height,
            direction,
        )

        for name, tags in boundary_groups.items():
            if tags:
                gmsh.model.addPhysicalGroup(1, tags, name=name)

        if not boundary_groups["INLET"]:
            raise RuntimeError("Nenhuma curva INLET foi detectada.")
        if not boundary_groups["OUTLET"]:
            raise RuntimeError("Nenhuma curva OUTLET foi detectada.")

        all_curves = sorted(curve_endpoint_map)

        distance_field = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(
            distance_field,
            "CurvesList",
            all_curves,
        )
        gmsh.model.mesh.field.setNumber(
            distance_field,
            "Sampling",
            100,
        )

        threshold_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "InField",
            distance_field,
        )
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "SizeMin",
            cfg.h_min,
        )
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "SizeMax",
            cfg.h_max,
        )
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "DistMin",
            cfg.distance_min,
        )
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "DistMax",
            cfg.distance_max,
        )
        gmsh.model.mesh.field.setNumber(
            threshold_field,
            "Sigmoid",
            1,
        )
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", cfg.mesh_algorithm)
        gmsh.option.setNumber("Mesh.Smoothing", cfg.smoothing_steps)
        gmsh.option.setNumber("Mesh.MshFileVersion", cfg.msh_version)

        gmsh.model.mesh.generate(2)

        mesh_path = sample_dir / "mesh.msh"
        geo_path = sample_dir / "geometry.geo_unrolled"

        gmsh.write(str(mesh_path))
        gmsh.write(str(geo_path))

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        element_types, element_tags, _ = gmsh.model.mesh.getElements(dim=2)

        triangle_count = 0

        for element_type, tags in zip(element_types, element_tags):
            name, dim, _, _, _, _ = (
                gmsh.model.mesh.getElementProperties(element_type)
            )

            if dim == 2 and "Triangle" in name:
                triangle_count += len(tags)

        return {
            "mesh_path": mesh_path,
            "geometry_path": geo_path,
            "contour_count": len(contours),
            "surface_count": len(surface_tags),
            "node_count": len(node_tags),
            "triangle_count": triangle_count,
            "inlet_curve_count": len(boundary_groups["INLET"]),
            "outlet_curve_count": len(boundary_groups["OUTLET"]),
            "wall_curve_count": len(boundary_groups["WALLS"]),
        }

    finally:
        gmsh.finalize()


# ---------------------------------------------------------------------------
# SOLUÇÃO DE STOKES
# ---------------------------------------------------------------------------

def inlet_velocity_function(direction: str, magnitude: float):
    def velocity(x: np.ndarray) -> np.ndarray:
        values = np.zeros(
            (2, x.shape[1], x.shape[2]),
            dtype=float,
        )

        if direction == "lr":
            values[0, :, :] = magnitude
        else:
            values[1, :, :] = -magnitude

        return values

    return velocity


def interpolate_basis_in_batches(
    basis: Basis,
    coefficients: np.ndarray,
    points: np.ndarray,
    batch_size: int,
) -> np.ndarray:
    """
    Avalia uma função de elementos finitos em pequenos lotes.

    A chamada única basis.interpolator(coefficients)(points) pode construir
    arrays temporários enormes em malhas grandes. O processamento por lotes
    limita o pico de memória.
    """
    evaluator = basis.interpolator(coefficients)
    batches = []

    for start in range(0, points.shape[1], batch_size):
        end = min(start + batch_size, points.shape[1])
        values = np.asarray(evaluator(points[:, start:end]))
        batches.append(values)

    if not batches:
        return np.empty((0,), dtype=float)

    first = batches[0]

    if first.ndim == 1:
        return np.concatenate(batches, axis=0)

    return np.concatenate(batches, axis=-1)


def run_minres(
    matrix,
    rhs: np.ndarray,
    preconditioner: LinearOperator,
    tolerance: float,
    max_iterations: int,
) -> tuple[np.ndarray, int]:
    """
    Compatibilidade com versões do SciPy que usam rtol ou tol.
    """
    try:
        return minres(
            matrix,
            rhs,
            M=preconditioner,
            rtol=tolerance,
            maxiter=max_iterations,
            show=False,
            check=False,
        )
    except TypeError:
        return minres(
            matrix,
            rhs,
            M=preconditioner,
            tol=tolerance,
            maxiter=max_iterations,
            show=False,
            check=False,
        )


def solve_stokes(
    mesh_path: Path,
    direction: str,
    cfg: Config,
) -> dict[str, object]:
    """
    Solução de Stokes totalmente esparsa.

    Diferenças em relação à versão anterior:
        - não usa fatoração LU direta;
        - não cria matrizes densas;
        - elimina Dirichlet por indexação esparsa;
        - usa MINRES com precondicionador diagonal em blocos;
        - interpola os resultados nos vértices em lotes.
    """
    mesh = MeshTri.load(str(mesh_path))

    if mesh.boundaries is None:
        raise RuntimeError("A malha foi carregada sem grupos físicos.")

    missing = {"INLET", "OUTLET", "WALLS"} - set(mesh.boundaries)

    if missing:
        raise RuntimeError(
            f"Grupos físicos ausentes: {sorted(missing)}."
        )

    velocity_element = ElementVector(ElementTriMini())
    pressure_element = ElementTriP1()

    basis_u = Basis(mesh, velocity_element, intorder=4)
    basis_p = Basis(mesh, pressure_element, intorder=4)

    total_dofs = basis_u.N + basis_p.N

    if cfg.max_total_dofs > 0 and total_dofs > cfg.max_total_dofs:
        raise RuntimeError(
            f"A malha produz {total_dofs:,} graus de liberdade, acima do "
            f"limite de segurança {cfg.max_total_dofs:,}. "
            "Gere uma malha mais grossa, por exemplo aumentando "
            "--h-min e --h-max, ou ajuste --max-total-dofs."
        )

    print(
        f"  DOFs velocidade: {basis_u.N:,} | "
        f"pressão: {basis_p.N:,} | total: {total_dofs:,}"
    )

    # Todas as matrizes abaixo são esparsas.
    A = (cfg.viscosity * asm(vector_laplace, basis_u)).tocsr()
    B = (-asm(divergence, basis_u, basis_p)).tocsr()
    pressure_mass = asm(mass, basis_p).tocsr()

    K = bmat(
        [
            [A, B.T],
            [B, None],
        ],
        format="csr",
    )

    rhs = np.zeros(total_dofs, dtype=float)

    x_u = np.zeros(basis_u.N, dtype=float)
    x_p = np.zeros(basis_p.N, dtype=float)

    inlet_fbasis = FacetBasis(
        mesh,
        velocity_element,
        facets=mesh.boundaries["INLET"],
        intorder=4,
    )

    projected_inlet = inlet_fbasis.project(
        inlet_velocity_function(direction, cfg.inlet_velocity)
    )

    inlet_dofs = basis_u.get_dofs("INLET").all()
    wall_dofs = basis_u.get_dofs("WALLS").all()

    # Nos cantos INLET-WALLS, a condição de parede prevalece.
    x_u[inlet_dofs] = projected_inlet[inlet_dofs]
    x_u[wall_dofs] = 0.0

    outlet_pressure_dofs = basis_p.get_dofs("OUTLET").all()

    if outlet_pressure_dofs.size == 0:
        raise RuntimeError("OUTLET não possui DOFs de pressão.")

    pressure_reference = int(outlet_pressure_dofs[0])
    x_p[pressure_reference] = 0.0

    prescribed_values = np.concatenate([x_u, x_p])

    prescribed_dofs = np.unique(
        np.concatenate(
            [
                inlet_dofs,
                wall_dofs,
                np.array(
                    [basis_u.N + pressure_reference],
                    dtype=np.int64,
                ),
            ]
        )
    )

    all_dofs = np.arange(total_dofs, dtype=np.int64)
    free_mask = np.ones(total_dofs, dtype=bool)
    free_mask[prescribed_dofs] = False
    free_dofs = all_dofs[free_mask]

    # Eliminação de Dirichlet mantendo o sistema esparso e simétrico.
    rhs_free = (
        rhs[free_dofs]
        - K[free_dofs][:, prescribed_dofs]
        @ prescribed_values[prescribed_dofs]
    )

    K_free = K[free_dofs][:, free_dofs].tocsr()

    # Precondicionador SPD em blocos:
    #   velocidade ~ diagonal da rigidez;
    #   pressão   ~ diagonal da matriz de massa / viscosidade.
    velocity_diag = np.maximum(np.abs(A.diagonal()), 1e-14)
    pressure_diag = np.maximum(
        np.abs(pressure_mass.diagonal()) / cfg.viscosity,
        1e-14,
    )

    preconditioner_diag_full = np.concatenate(
        [velocity_diag, pressure_diag]
    )

    preconditioner_diag = preconditioner_diag_full[free_dofs]

    preconditioner = LinearOperator(
        shape=K_free.shape,
        matvec=lambda vector: vector / preconditioner_diag,
        dtype=float,
    )

    estimated_sparse_mb = (
        K_free.data.nbytes
        + K_free.indices.nbytes
        + K_free.indptr.nbytes
    ) / 1024**2

    print(
        f"  Sistema livre: {K_free.shape[0]:,} incógnitas | "
        f"nnz={K_free.nnz:,} | armazenamento CSR≈{estimated_sparse_mb:.1f} MB"
    )
    print(
        f"  Resolvendo por MINRES: tol={cfg.solver_tolerance:g}, "
        f"maxiter={cfg.solver_max_iterations}"
    )

    solution_free, info = run_minres(
        K_free,
        rhs_free,
        preconditioner,
        cfg.solver_tolerance,
        cfg.solver_max_iterations,
    )

    if info < 0:
        raise RuntimeError(
            f"MINRES terminou com erro interno, código {info}."
        )

    if info > 0:
        raise RuntimeError(
            f"MINRES não convergiu em {info} iterações. "
            "Tente uma malha mais grossa, aumente "
            "--solver-max-iterations ou relaxe --solver-tolerance "
            "para 1e-7."
        )

    solution = prescribed_values.copy()
    solution[free_dofs] = solution_free

    residual = K_free @ solution_free - rhs_free
    relative_residual = float(
        np.linalg.norm(residual)
        / max(np.linalg.norm(rhs_free), 1e-30)
    )

    velocity = solution[: basis_u.N]
    pressure = solution[basis_u.N :]

    print(
        "  Interpolando solução nos vértices em lotes de "
        f"{cfg.interpolation_batch_size:,}..."
    )

    velocity_vertices = interpolate_basis_in_batches(
        basis_u,
        velocity,
        mesh.p,
        cfg.interpolation_batch_size,
    )

    pressure_vertices = interpolate_basis_in_batches(
        basis_p,
        pressure,
        mesh.p,
        cfg.interpolation_batch_size,
    ).reshape(-1)

    if velocity_vertices.shape != (2, mesh.p.shape[1]):
        raise RuntimeError(
            "Formato inesperado da velocidade interpolada: "
            f"{velocity_vertices.shape}."
        )

    speed_vertices = np.linalg.norm(
        velocity_vertices,
        axis=0,
    )

    return {
        "mesh": mesh,
        "basis_u": basis_u,
        "basis_p": basis_p,
        "velocity": velocity,
        "pressure": pressure,
        "velocity_vertices": velocity_vertices,
        "pressure_vertices": pressure_vertices,
        "speed_vertices": speed_vertices,
        "solver_name": "MINRES sparse",
        "solver_info": int(info),
        "solver_relative_residual": relative_residual,
        "velocity_dof_count": int(basis_u.N),
        "pressure_dof_count": int(basis_p.N),
        "total_dof_count": int(total_dofs),
        "free_dof_count": int(len(free_dofs)),
        "system_nnz": int(K_free.nnz),
        "system_csr_megabytes": float(estimated_sparse_mb),
    }


# ---------------------------------------------------------------------------
# MÉTRICAS E EXPORTAÇÃO
# ---------------------------------------------------------------------------

def boundary_length(mesh: MeshTri, name: str) -> float:
    facets = np.asarray(
        mesh.boundaries[name],
        dtype=np.int64,
    )

    edge_nodes = mesh.facets[:, facets]
    p0 = mesh.p[:, edge_nodes[0]]
    p1 = mesh.p[:, edge_nodes[1]]

    return float(
        np.linalg.norm(p1 - p0, axis=0).sum()
    )


def boundary_pressure_mean(
    basis_p: Basis,
    pressure: np.ndarray,
    name: str,
) -> float:
    dofs = basis_p.get_dofs(name).all()

    if dofs.size == 0:
        raise RuntimeError(f"Fronteira {name} sem DOFs de pressão.")

    return float(np.mean(pressure[dofs]))


def triangle_areas(
    points: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    p0 = points[:, triangles[0]].T
    p1 = points[:, triangles[1]].T
    p2 = points[:, triangles[2]].T

    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0])
        * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1])
        * (p2[:, 0] - p0[:, 0])
    )


def calculate_metrics(
    result: dict[str, object],
    direction: str,
    cfg: Config,
) -> dict[str, object]:
    mesh: MeshTri = result["mesh"]
    basis_p: Basis = result["basis_p"]
    pressure = result["pressure"]
    speed = result["speed_vertices"]

    inlet_length = boundary_length(mesh, "INLET")
    outlet_length = boundary_length(mesh, "OUTLET")
    wall_length = boundary_length(mesh, "WALLS")

    axis = 0 if direction == "lr" else 1
    transverse_axis = 1 - axis

    flow_length = float(
        mesh.p[axis].max() - mesh.p[axis].min()
    )

    transverse_length = float(
        mesh.p[transverse_axis].max()
        - mesh.p[transverse_axis].min()
    )

    p_in = boundary_pressure_mean(
        basis_p,
        pressure,
        "INLET",
    )

    p_out = boundary_pressure_mean(
        basis_p,
        pressure,
        "OUTLET",
    )

    delta_p = abs(p_in - p_out)

    flow_rate = (
        abs(cfg.inlet_velocity)
        * inlet_length
        * cfg.thickness
    )

    cross_section = inlet_length * cfg.thickness
    darcy_velocity = flow_rate / cross_section

    permeability = (
        cfg.viscosity
        * darcy_velocity
        * flow_length
        / delta_p
        if delta_p > 1e-30
        else float("nan")
    )

    areas = triangle_areas(mesh.p, mesh.t)
    fluid_area = float(areas.sum())

    bounding_area = flow_length * transverse_length

    porosity_domain = (
        fluid_area / bounding_area
        if bounding_area > 0
        else float("nan")
    )

    triangle_speed = np.mean(speed[mesh.t], axis=0)

    speed_max = float(np.max(speed))
    active_threshold = (
        cfg.active_threshold_relative * speed_max
    )
    dead_threshold = (
        cfg.dead_threshold_relative * speed_max
    )

    active_area = float(
        areas[triangle_speed >= active_threshold].sum()
    )

    dead_area = float(
        areas[triangle_speed <= dead_threshold].sum()
    )

    active_fraction = active_area / fluid_area
    dead_fraction = dead_area / fluid_area

    hydraulic_diameter = (
        4.0 * fluid_area / wall_length
        if wall_length > 0
        else float("nan")
    )

    reynolds = (
        cfg.density
        * abs(darcy_velocity)
        * hydraulic_diameter
        / cfg.viscosity
    )

    return {
        "direction": direction,
        "length_unit": cfg.unit_name,
        "viscosity": cfg.viscosity,
        "density": cfg.density,
        "inlet_velocity": cfg.inlet_velocity,
        "thickness": cfg.thickness,
        "vertex_count": int(mesh.p.shape[1]),
        "triangle_count": int(mesh.t.shape[1]),
        "fluid_area": fluid_area,
        "bounding_area": bounding_area,
        "domain_porosity": porosity_domain,
        "flow_length": flow_length,
        "transverse_length": transverse_length,
        "inlet_length": inlet_length,
        "outlet_length": outlet_length,
        "wall_length": wall_length,
        "pressure_inlet_mean": p_in,
        "pressure_outlet_mean": p_out,
        "pressure_drop": delta_p,
        "flow_rate_prescribed": flow_rate,
        "darcy_velocity": darcy_velocity,
        "apparent_permeability": permeability,
        "permeability_unit": f"{cfg.unit_name}^2",
        "speed_max": speed_max,
        "active_threshold": active_threshold,
        "dead_threshold": dead_threshold,
        "active_area": active_area,
        "active_area_fraction": active_fraction,
        "dead_area": dead_area,
        "dead_area_fraction": dead_fraction,
        "hydraulic_diameter_2d": hydraulic_diameter,
        "reynolds_number": reynolds,
    }


def save_vtk(
    result: dict[str, object],
    metrics: dict[str, object],
    output_path: Path,
) -> None:
    mesh: MeshTri = result["mesh"]
    velocity_vertices = result["velocity_vertices"]
    pressure_vertices = result["pressure_vertices"]
    speed_vertices = result["speed_vertices"]

    points_3d = np.column_stack(
        [
            mesh.p[0],
            mesh.p[1],
            np.zeros(mesh.p.shape[1]),
        ]
    )

    velocity_3d = np.column_stack(
        [
            velocity_vertices[0],
            velocity_vertices[1],
            np.zeros(mesh.p.shape[1]),
        ]
    )

    active = (
        speed_vertices >= metrics["active_threshold"]
    ).astype(np.int8)

    dead = (
        speed_vertices <= metrics["dead_threshold"]
    ).astype(np.int8)

    vtk_mesh = meshio.Mesh(
        points=points_3d,
        cells=[("triangle", mesh.t.T)],
        point_data={
            "velocity": velocity_3d,
            "speed": speed_vertices,
            "pressure": pressure_vertices,
            "active_zone": active,
            "dead_zone": dead,
        },
    )

    vtk_mesh.write(output_path)


def save_solution_figures(
    result: dict[str, object],
    metrics: dict[str, object],
    output_dir: Path,
    cfg: Config,
) -> None:
    mesh: MeshTri = result["mesh"]
    pressure = result["pressure_vertices"]
    speed = result["speed_vertices"]
    velocity = result["velocity_vertices"]

    x = mesh.p[0]
    y = mesh.p[1]
    triangles = mesh.t.T

    fig, ax = plt.subplots(figsize=(12, 10))
    image = ax.tripcolor(
        x,
        y,
        triangles,
        pressure,
        shading="gouraud",
    )
    ax.set_aspect("equal")
    ax.set_title("Pressão")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="p")
    fig.savefig(
        output_dir / "pressure.png",
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 10))
    image = ax.tripcolor(
        x,
        y,
        triangles,
        speed,
        shading="gouraud",
    )
    ax.set_aspect("equal")
    ax.set_title("Módulo da velocidade")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="|u|")
    fig.savefig(
        output_dir / "velocity_magnitude.png",
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    active = (
        speed >= metrics["active_threshold"]
    ).astype(float)

    fig, ax = plt.subplots(figsize=(12, 10))
    image = ax.tripcolor(
        x,
        y,
        triangles,
        active,
        shading="flat",
        vmin=0,
        vmax=1,
    )
    ax.set_aspect("equal")
    ax.set_title("Zona ativa de escoamento")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="ativa")
    fig.savefig(
        output_dir / "active_zone.png",
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    dead = (
        speed <= metrics["dead_threshold"]
    ).astype(float)

    fig, ax = plt.subplots(figsize=(12, 10))
    image = ax.tripcolor(
        x,
        y,
        triangles,
        dead,
        shading="flat",
        vmin=0,
        vmax=1,
    )
    ax.set_aspect("equal")
    ax.set_title("Zona morta ou quase estagnada")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="morta")
    fig.savefig(
        output_dir / "dead_zone.png",
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    count = mesh.p.shape[1]
    max_arrows = 2500
    stride = max(1, int(math.ceil(count / max_arrows)))
    indices = np.arange(0, count, stride)

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.triplot(x, y, triangles, linewidth=0.05)
    ax.quiver(
        x[indices],
        y[indices],
        velocity[0, indices],
        velocity[1, indices],
        speed[indices],
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0015,
    )
    ax.set_aspect("equal")
    ax.set_title("Campo de velocidade")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.savefig(
        output_dir / "velocity_vectors.png",
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_mask_diagnostic(
    original: np.ndarray,
    clean: np.ndarray,
    processed: np.ndarray,
    cropped: np.ndarray,
    selected: np.ndarray,
    added: np.ndarray,
    output_path: Path,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(18, 12),
        constrained_layout=True,
    )

    images = [
        (original, "Máscara original"),
        (clean, "Máscara limpa"),
        (processed, "Após fechamento opcional"),
        (added, "Pixels adicionados"),
        (cropped, "Todos os poros entre os planos"),
        (selected, "Domínio percolante completo"),
    ]

    for ax, (image, title) in zip(axes.ravel(), images):
        ax.imshow(image, cmap="gray", origin="upper")
        ax.set_title(title)
        ax.axis("off")

    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_mesh_preview(
    mesh_path: Path,
    output_path: Path,
    cfg: Config,
) -> None:
    mesh = meshio.read(mesh_path)
    points = mesh.points[:, :2]

    triangles = []

    for block in mesh.cells:
        if block.type == "triangle":
            triangles.append(block.data)

    if not triangles:
        return

    triangles = np.vstack(triangles)

    if len(triangles) > cfg.max_preview_triangles:
        step = int(
            math.ceil(
                len(triangles)
                / cfg.max_preview_triangles
            )
        )
        triangles = triangles[::step]

    fig, ax = plt.subplots(figsize=(11, 11))
    ax.triplot(
        points[:, 0],
        points[:, 1],
        triangles,
        linewidth=0.2,
    )
    ax.set_aspect("equal")
    ax.set_title("Malha triangular adaptativa")
    ax.set_xlabel(f"x [{cfg.unit_name}]")
    ax.set_ylabel(f"y [{cfg.unit_name}]")
    fig.savefig(
        output_path,
        dpi=cfg.dpi,
        bbox_inches="tight",
    )
    plt.close(fig)


def write_rows(
    rows: Iterable[dict[str, object]],
    output_path: Path,
) -> None:
    rows = list(rows)

    if not rows:
        return

    fieldnames: list[str] = []

    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def write_single_summary(
    summary: dict[str, object],
    output_dir: Path,
) -> None:
    write_rows(
        [summary],
        output_dir / "summary.csv",
    )

    with (
        output_dir / "summary.json"
    ).open("w", encoding="utf-8") as file:
        json.dump(
            summary,
            file,
            indent=2,
            ensure_ascii=False,
            default=str,
        )


# ---------------------------------------------------------------------------
# PROCESSAMENTO DE UMA AMOSTRA E DIREÇÃO
# ---------------------------------------------------------------------------

def process_sample_direction(
    image_path: Path,
    direction: str,
    cfg: Config,
) -> dict[str, object]:
    sample_dir = (
        cfg.output_dir
        / direction
        / image_path.stem
    )
    sample_dir.mkdir(parents=True, exist_ok=True)

    pore_original = load_pore_mask(image_path, cfg)
    pore_clean = clean_mask(pore_original, cfg)
    pore_processed = apply_optional_closing(pore_clean, cfg)

    pixels_added = pore_processed & ~pore_clean

    cropped_all, crop_info = crop_between_fixed_planes(
        pore_processed,
        direction,
        cfg,
    )

    selected, selection_info = select_full_percolating_components(
        cropped_all,
        direction,
        cfg,
    )

    if not np.any(selected):
        raise RuntimeError(
            "Nenhum componente conecta os dois planos fixos."
        )

    np.save(
        sample_dir / "full_percolating_domain.npy",
        selected,
    )

    Image.fromarray(
        selected.astype(np.uint8) * 255,
        mode="L",
    ).save(
        sample_dir / "full_percolating_domain.png"
    )

    save_mask_diagnostic(
        original=pore_original,
        clean=pore_clean,
        processed=pore_processed,
        cropped=cropped_all,
        selected=selected,
        added=pixels_added,
        output_path=sample_dir / "domain_diagnostic.png",
        dpi=cfg.dpi,
    )

    mesh_info = create_gmsh_mesh(
        selected,
        direction,
        sample_dir,
        cfg,
    )

    save_mesh_preview(
        mesh_info["mesh_path"],
        sample_dir / "mesh_preview.png",
        cfg,
    )

    stokes_result = solve_stokes(
        mesh_info["mesh_path"],
        direction,
        cfg,
    )

    metrics = calculate_metrics(
        stokes_result,
        direction,
        cfg,
    )

    save_vtk(
        stokes_result,
        metrics,
        sample_dir / "solution.vtk",
    )

    save_solution_figures(
        stokes_result,
        metrics,
        sample_dir,
        cfg,
    )

    clean_pixels = int(pore_clean.sum())
    selected_pixels = int(selected.sum())
    cropped_pixels = int(cropped_all.sum())
    added_pixels = int(pixels_added.sum())

    summary = {
        "filename": image_path.name,
        "direction": direction,
        "status": "success",
        "band_px": crop_info["band_px"],
        "x_offset_px": crop_info["x_offset_px"],
        "y_offset_px": crop_info["y_offset_px"],
        "plane_start_original_px": (
            crop_info["plane_start_original_px"]
        ),
        "plane_end_original_px": (
            crop_info["plane_end_original_px"]
        ),
        "candidate_percolating_count": (
            selection_info["candidate_count"]
        ),
        "selected_percolating_count": (
            selection_info["selected_count"]
        ),
        "clean_pore_pixels": clean_pixels,
        "cropped_pore_pixels": cropped_pixels,
        "selected_pore_pixels": selected_pixels,
        "selected_fraction_of_cropped_pores": (
            selected_pixels / cropped_pixels
            if cropped_pixels
            else 0.0
        ),
        "closing_radius_px": cfg.closing_radius,
        "pixels_added_by_closing": added_pixels,
        "pixels_added_fraction": (
            added_pixels / clean_pixels
            if clean_pixels
            else 0.0
        ),
        "contour_count": mesh_info["contour_count"],
        "surface_count": mesh_info["surface_count"],
        "gmsh_node_count": mesh_info["node_count"],
        "gmsh_triangle_count": mesh_info["triangle_count"],
        "inlet_curve_count": mesh_info["inlet_curve_count"],
        "outlet_curve_count": mesh_info["outlet_curve_count"],
        "wall_curve_count": mesh_info["wall_curve_count"],
        "mesh_path": str(mesh_info["mesh_path"]),
        "solution_vtk": str(sample_dir / "solution.vtk"),
        "solver_name": stokes_result["solver_name"],
        "solver_info": stokes_result["solver_info"],
        "solver_relative_residual": (
            stokes_result["solver_relative_residual"]
        ),
        "velocity_dof_count": stokes_result["velocity_dof_count"],
        "pressure_dof_count": stokes_result["pressure_dof_count"],
        "total_dof_count": stokes_result["total_dof_count"],
        "free_dof_count": stokes_result["free_dof_count"],
        "system_nnz": stokes_result["system_nnz"],
        "system_csr_megabytes": (
            stokes_result["system_csr_megabytes"]
        ),
        **metrics,
    }

    write_single_summary(summary, sample_dir)

    return summary


def main() -> int:
    cfg = parse_args()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        images = find_images(cfg.structures_dir)
    except Exception as exc:
        print(f"ERRO: {exc}", file=sys.stderr)
        return 1

    if not images:
        print("ERRO: nenhuma imagem encontrada.", file=sys.stderr)
        return 1

    print(f"Imagens encontradas: {len(images)}")
    print(f"Direções: {', '.join(cfg.directions)}")
    print(f"Entrada: {cfg.structures_dir.resolve()}")
    print(f"Saída: {cfg.output_dir.resolve()}")
    print(
        "Pipeline: máscara -> domínio percolante -> malha -> Stokes -> análise."
    )

    rows: list[dict[str, object]] = []

    total_jobs = len(images) * len(cfg.directions)
    job_index = 0

    for direction in cfg.directions:
        for image_path in images:
            job_index += 1

            print(
                f"\n[{job_index}/{total_jobs}] "
                f"{image_path.name} | direção={direction}"
            )

            try:
                summary = process_sample_direction(
                    image_path,
                    direction,
                    cfg,
                )

                rows.append(summary)

                print(
                    f"  Nós: {summary['gmsh_node_count']}"
                )
                print(
                    f"  Triângulos: {summary['gmsh_triangle_count']}"
                )
                print(
                    "  ΔP: "
                    f"{summary['pressure_drop']:.8e}"
                )
                print(
                    "  k: "
                    f"{summary['apparent_permeability']:.8e} "
                    f"{summary['permeability_unit']}"
                )
                print(
                    "  Zona ativa: "
                    f"{100.0 * summary['active_area_fraction']:.3f}%"
                )
                print(
                    "  Zona morta: "
                    f"{100.0 * summary['dead_area_fraction']:.3f}%"
                )
                print(
                    "  Reynolds: "
                    f"{summary['reynolds_number']:.8e}"
                )

            except Exception as exc:
                error_message = str(exc)

                if "Nenhum componente conecta os dois planos fixos" in error_message:
                    print(
                        "  NÃO PERCOLANTE: nenhum caminho contínuo entre "
                        "entrada e saída."
                    )
                    rows.append(
                        {
                            "filename": image_path.name,
                            "direction": direction,
                            "status": "not_percolating",
                            "apparent_permeability": 0.0,
                            "permeability_unit": f"{cfg.unit_name}^2",
                            "error": error_message,
                        }
                    )
                else:
                    print(f"  FALHA: {error_message}", file=sys.stderr)
                    rows.append(
                        {
                            "filename": image_path.name,
                            "direction": direction,
                            "status": "failed",
                            "error": error_message,
                        }
                    )

    write_rows(
        rows,
        cfg.output_dir / "global_summary.csv",
    )

    with (
        cfg.output_dir / "run_configuration.json"
    ).open("w", encoding="utf-8") as file:
        json.dump(
            cfg.__dict__,
            file,
            indent=2,
            ensure_ascii=False,
            default=str,
        )

    success_count = sum(
        row.get("status") == "success"
        for row in rows
    )
    nonpercolating_count = sum(
        row.get("status") == "not_percolating"
        for row in rows
    )

    print("\nPipeline concluído.")
    print(f"Casos CFD concluídos: {success_count}/{len(rows)}")
    print(f"Casos não percolantes: {nonpercolating_count}")
    print(
        f"Resumo global: "
        f"{cfg.output_dir / 'global_summary.csv'}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
