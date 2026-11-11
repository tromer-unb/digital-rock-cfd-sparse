#!/usr/bin/env python3
"""
Figura 2 — versão leve para SVG híbrido.

Os campos de pressão/velocidade e as máscaras são incorporados como imagens
raster dentro do SVG. Letras, títulos, setas, molduras e barras de cor
permanecem vetoriais e editáveis.

Isso evita SVGs com centenas de milhares de polígonos triangulares, que podem
consumir vários GB de memória ao serem abertos no Inkscape.

Dependências:
    pip install numpy pandas matplotlib meshio pillow

Uso:
    python3 figure2_light.py \
        --results-root ../convergence_mesh_05_extra_fine \
        --output Figure_2_flow_fields.svg
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import meshio
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PIL import Image


@dataclass(frozen=True)
class FlowCase:
    row_name: str
    filename_stem: str
    direction: str


CASES = (
    FlowCase(
        row_name="R3 horizontal",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="lr",
    ),
    FlowCase(
        row_name="R2 vertical",
        filename_stem="patch_y7600_x19000_c0_mask",
        direction="tb",
    ),
    FlowCase(
        row_name="R3 vertical",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="tb",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera uma Figura 2 leve em SVG híbrido."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Pasta da malha escolhida, contendo os diretórios lr/ e tb/.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_2_flow_fields.svg"),
        help="Arquivo SVG de saída.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=170,
        help="Resolução dos mapas raster incorporados no SVG.",
    )
    return parser.parse_args()


def get_triangles(mesh: meshio.Mesh) -> np.ndarray:
    blocks = [
        block.data
        for block in mesh.cells
        if block.type.startswith("triangle")
    ]
    if not blocks:
        raise RuntimeError("Nenhum triângulo encontrado no VTK.")

    triangles = np.vstack(blocks)
    return np.asarray(triangles[:, :3], dtype=np.int64)


def read_scalar(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.point_data:
        raise KeyError(f"Campo ausente no VTK: {name}")

    return np.asarray(mesh.point_data[name], dtype=float).squeeze().reshape(-1)


def read_speed(mesh: meshio.Mesh) -> np.ndarray:
    if "speed" in mesh.point_data:
        return read_scalar(mesh, "speed")

    if "velocity" not in mesh.point_data:
        raise KeyError("O VTK não contém 'speed' nem 'velocity'.")

    velocity = np.asarray(mesh.point_data["velocity"], dtype=float)
    return np.linalg.norm(velocity[:, :2], axis=1)


def read_mask(sample_dir: Path) -> np.ndarray:
    npy_path = sample_dir / "full_percolating_domain.npy"
    png_path = sample_dir / "full_percolating_domain.png"

    if npy_path.is_file():
        return np.asarray(np.load(npy_path), dtype=bool)

    if png_path.is_file():
        return np.asarray(Image.open(png_path).convert("L")) > 0

    raise FileNotFoundError(
        "Máscara percolante não encontrada em "
        f"{npy_path} ou {png_path}."
    )


def read_summary(sample_dir: Path) -> dict[str, float]:
    path = sample_dir / "summary.csv"

    if not path.is_file():
        raise FileNotFoundError(f"Resumo não encontrado: {path}")

    frame = pd.read_csv(path)
    if frame.empty:
        raise RuntimeError(f"Resumo vazio: {path}")

    row = frame.iloc[0]

    return {
        "pressure_outlet_mean": float(row["pressure_outlet_mean"]),
        "pressure_drop": float(row["pressure_drop"]),
        "inlet_velocity": abs(float(row["inlet_velocity"])),
    }


def load_case(results_root: Path, case: FlowCase) -> dict[str, object]:
    sample_dir = results_root / case.direction / case.filename_stem
    vtk_path = sample_dir / "solution.vtk"

    if not vtk_path.is_file():
        raise FileNotFoundError(f"VTK não encontrado: {vtk_path}")

    mesh = meshio.read(vtk_path)
    points = np.asarray(mesh.points, dtype=float)
    triangles = get_triangles(mesh)

    pressure = read_scalar(mesh, "pressure")
    speed = read_speed(mesh)
    summary = read_summary(sample_dir)
    mask = read_mask(sample_dir)

    delta_p = summary["pressure_drop"]
    p_out = summary["pressure_outlet_mean"]
    inlet_velocity = summary["inlet_velocity"]

    if delta_p <= np.finfo(float).eps:
        raise RuntimeError(f"Queda de pressão inválida em {case.row_name}.")
    if inlet_velocity <= np.finfo(float).eps:
        raise RuntimeError(f"Velocidade de entrada inválida em {case.row_name}.")

    p_star = (pressure - p_out) / delta_p
    p_star = np.clip(p_star, 0.0, 1.0)

    speed_star = speed / inlet_velocity

    triangulation = mtri.Triangulation(
        points[:, 0],
        points[:, 1],
        triangles,
    )

    return {
        "case": case,
        "sample_dir": sample_dir,
        "mask": mask,
        "triangulation": triangulation,
        "pressure_star": p_star,
        "speed_star": speed_star,
    }


def clean_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.018,
        0.982,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        bbox={
            "boxstyle": "square,pad=0.10",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.88,
        },
        zorder=20,
    )


def add_flow_arrow(ax: plt.Axes, direction: str) -> None:
    if direction == "lr":
        start = (0.14, 0.10)
        end = (0.35, 0.10)
    else:
        start = (0.10, 0.84)
        end = (0.10, 0.63)

    ax.annotate(
        "",
        xy=end,
        xytext=start,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "linewidth": 2.1,
            "mutation_scale": 16,
        },
        zorder=20,
    )


def plot_domain(
    ax: plt.Axes,
    data: dict[str, object],
    panel_label: str,
) -> None:
    # imshow cria somente uma imagem raster, em vez de centenas de milhares
    # de polígonos vetoriais.
    ax.imshow(
        data["mask"],
        origin="upper",
        cmap="gray",
        interpolation="nearest",
        vmin=0,
        vmax=1,
        rasterized=True,
    )

    add_flow_arrow(ax, data["case"].direction)
    add_panel_label(ax, panel_label)
    clean_axis(ax)


def plot_field(
    ax: plt.Axes,
    triangulation: mtri.Triangulation,
    values: np.ndarray,
    norm: Normalize,
    panel_label: str,
):
    # rasterized=True é essencial. Sem isso, cada triângulo vira um caminho
    # vetorial dentro do SVG.
    artist = ax.tripcolor(
        triangulation,
        values,
        shading="gouraud",
        norm=norm,
        edgecolors="none",
        antialiased=False,
        rasterized=True,
        zorder=0,
    )

    add_panel_label(ax, panel_label)
    clean_axis(ax)
    return artist


def main() -> int:
    args = parse_args()

    if not args.results_root.is_dir():
        raise NotADirectoryError(
            f"Pasta de resultados inexistente: {args.results_root}"
        )

    datasets = [load_case(args.results_root, case) for case in CASES]

    speed_max = max(
        float(np.nanmax(data["speed_star"]))
        for data in datasets
    )

    pressure_norm = Normalize(vmin=0.0, vmax=1.0)
    speed_norm = Normalize(vmin=0.0, vmax=speed_max)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "svg.fonttype": "none",
            "svg.image_inline": True,
            "image.composite_image": True,
            "savefig.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(
        3,
        3,
        figsize=(12.4, 12.8),
        dpi=args.dpi,
    )

    fig.subplots_adjust(
        left=0.105,
        right=0.965,
        bottom=0.105,
        top=0.935,
        wspace=0.055,
        hspace=0.075,
    )

    column_titles = (
        "Percolating domain",
        r"Normalized pressure, $p^{*}$",
        r"Normalized speed, $|\mathbf{u}|/U_0$",
    )

    for column, title in enumerate(column_titles):
        axes[0, column].set_title(
            title,
            fontsize=16,
            fontweight="bold",
            pad=11,
        )

    labels = iter(
        ("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)")
    )

    pressure_artist = None
    speed_artist = None

    for row, data in enumerate(datasets):
        plot_domain(axes[row, 0], data, next(labels))

        pressure_artist = plot_field(
            axes[row, 1],
            data["triangulation"],
            data["pressure_star"],
            pressure_norm,
            next(labels),
        )

        speed_artist = plot_field(
            axes[row, 2],
            data["triangulation"],
            data["speed_star"],
            speed_norm,
            next(labels),
        )

        axes[row, 0].text(
            -0.105,
            0.50,
            data["case"].row_name,
            transform=axes[row, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=15,
            fontweight="bold",
        )

    pressure_map = ScalarMappable(
        norm=pressure_norm,
        cmap=pressure_artist.cmap,
    )
    pressure_map.set_array([])

    speed_map = ScalarMappable(
        norm=speed_norm,
        cmap=speed_artist.cmap,
    )
    speed_map.set_array([])

    cbar_p = fig.colorbar(
        pressure_map,
        ax=axes[:, 1],
        orientation="horizontal",
        fraction=0.032,
        pad=0.035,
        aspect=34,
    )
    cbar_p.set_label(
        r"$p^{*}=(p-\overline{p}_{\mathrm{out}})/\Delta p$",
        fontsize=13,
        fontweight="bold",
    )

    cbar_u = fig.colorbar(
        speed_map,
        ax=axes[:, 2],
        orientation="horizontal",
        fraction=0.032,
        pad=0.035,
        aspect=34,
    )
    cbar_u.set_label(
        r"$|\mathbf{u}|/U_0$",
        fontsize=13,
        fontweight="bold",
    )

    for colorbar in (cbar_p, cbar_u):
        colorbar.ax.tick_params(
            labelsize=11,
            width=1.2,
            length=4,
        )
        for tick in colorbar.ax.get_xticklabels():
            tick.set_fontweight("bold")
        colorbar.outline.set_linewidth(1.1)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Evitar bbox_inches="tight": em alguns casos ele produz dimensões
    # excessivas devido aos rótulos posicionados fora dos eixos.
    fig.savefig(
        args.output,
        format="svg",
        dpi=args.dpi,
        facecolor="white",
        metadata={
            "Title": "Pore-scale pressure and velocity fields",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    size_mb = args.output.stat().st_size / 1024**2
    print(f"Figura salva em: {args.output.resolve()}")
    print(f"Tamanho do SVG: {size_mb:.2f} MB")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
