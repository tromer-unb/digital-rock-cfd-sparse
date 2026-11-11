#!/usr/bin/env python3
"""
Gera a Figura 3 do artigo:

(a) mapas das zonas hidraulicamente ativas;
(b) mapas das zonas mortas ou quase estagnadas;
(c) partição da área porosa em ativa, intermediária e morta;
(d) permeabilidade aparente e anisotropia direcional.

Os dados são lidos diretamente da pasta de resultados da malha selecionada.

Estrutura esperada:
    RESULTS_ROOT/
        lr/
            patch_y7600_x53200_c0_mask/
                solution.vtk
                summary.csv
        tb/
            patch_y7600_x19000_c0_mask/
                solution.vtk
                summary.csv
            patch_y7600_x53200_c0_mask/
                solution.vtk
                summary.csv

Dependências:
    pip install numpy pandas matplotlib meshio

Exemplo:
    python3 figure_3_flow_partition.py \
        --results-root convergence_mesh_05_extra_fine \
        --output Figure_3_flow_partition.svg

Por padrão, somente os mapas com centenas de milhares de triângulos são
rasterizados dentro do SVG. Textos, eixos, barras e anotações permanecem
vetoriais. Use --fully-vectorial apenas se aceitar um SVG potencialmente
muito grande.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import meshio
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


@dataclass(frozen=True)
class FlowCase:
    short_name: str
    display_name: str
    filename_stem: str
    direction: str


CASES = (
    FlowCase(
        short_name="R3-H",
        display_name="R3 horizontal",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="lr",
    ),
    FlowCase(
        short_name="R2-V",
        display_name="R2 vertical",
        filename_stem="patch_y7600_x19000_c0_mask",
        direction="tb",
    ),
    FlowCase(
        short_name="R3-V",
        display_name="R3 vertical",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="tb",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera a Figura 3 sobre partição hidráulica e anisotropia."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Pasta de resultados da malha selecionada.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_3_flow_partition.svg"),
        help="Arquivo SVG de saída.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=350,
        help="Resolução dos mapas rasterizados no interior do SVG.",
    )
    parser.add_argument(
        "--fully-vectorial",
        action="store_true",
        help="Mantém também os mapas de triângulos em formato vetorial.",
    )
    return parser.parse_args()


def get_triangles(mesh: meshio.Mesh) -> np.ndarray:
    blocks = [
        block.data
        for block in mesh.cells
        if block.type.startswith("triangle")
    ]

    if not blocks:
        raise RuntimeError("Nenhum elemento triangular foi encontrado.")

    triangles = np.vstack(blocks)

    if triangles.shape[1] > 3:
        triangles = triangles[:, :3]

    return np.asarray(triangles, dtype=np.int64)


def read_speed(mesh: meshio.Mesh) -> np.ndarray:
    if "speed" in mesh.point_data:
        values = np.asarray(mesh.point_data["speed"], dtype=float)
        return np.squeeze(values).reshape(-1)

    if "velocity" not in mesh.point_data:
        raise KeyError("O VTK não possui os campos 'speed' ou 'velocity'.")

    velocity = np.asarray(mesh.point_data["velocity"], dtype=float)

    if velocity.ndim != 2 or velocity.shape[1] < 2:
        raise ValueError(
            f"Formato inesperado do campo velocity: {velocity.shape}"
        )

    return np.linalg.norm(velocity[:, :2], axis=1)


def read_summary(sample_dir: Path) -> pd.Series:
    summary_path = sample_dir / "summary.csv"

    if not summary_path.is_file():
        raise FileNotFoundError(f"Resumo não encontrado: {summary_path}")

    frame = pd.read_csv(summary_path)

    if frame.empty:
        raise RuntimeError(f"Resumo vazio: {summary_path}")

    required = {
        "active_threshold",
        "dead_threshold",
        "active_area_fraction",
        "dead_area_fraction",
        "apparent_permeability",
    }

    missing = required.difference(frame.columns)

    if missing:
        raise ValueError(
            f"{summary_path} não possui as colunas: {sorted(missing)}"
        )

    return frame.iloc[0]


def load_case(results_root: Path, case: FlowCase) -> dict[str, object]:
    sample_dir = (
        results_root
        / case.direction
        / case.filename_stem
    )

    vtk_path = sample_dir / "solution.vtk"

    if not vtk_path.is_file():
        raise FileNotFoundError(
            f"VTK não encontrado para {case.display_name}: {vtk_path}"
        )

    mesh = meshio.read(vtk_path)
    triangles = get_triangles(mesh)

    points = np.asarray(mesh.points, dtype=float)
    x = points[:, 0]
    y = points[:, 1]

    speed = read_speed(mesh)
    summary = read_summary(sample_dir)

    triangle_speed = np.mean(speed[triangles], axis=1)

    active_threshold = float(summary["active_threshold"])
    dead_threshold = float(summary["dead_threshold"])

    active_cells = (triangle_speed >= active_threshold).astype(float)
    dead_cells = (triangle_speed <= dead_threshold).astype(float)

    active_fraction = float(summary["active_area_fraction"])
    dead_fraction = float(summary["dead_area_fraction"])
    intermediate_fraction = max(
        0.0,
        1.0 - active_fraction - dead_fraction,
    )

    triangulation = mtri.Triangulation(x, y, triangles)

    return {
        "case": case,
        "triangulation": triangulation,
        "active_cells": active_cells,
        "dead_cells": dead_cells,
        "active_fraction": active_fraction,
        "intermediate_fraction": intermediate_fraction,
        "dead_fraction": dead_fraction,
        "permeability": float(summary["apparent_permeability"]),
    }


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        clip_on=False,
    )


def clean_map_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(1.3)


def add_flow_arrow(ax: plt.Axes, direction: str) -> None:
    if direction == "lr":
        xytext = (0.14, 0.09)
        xy = (0.36, 0.09)
    else:
        xytext = (0.10, 0.86)
        xy = (0.10, 0.64)

    ax.annotate(
        "",
        xy=xy,
        xytext=xytext,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "linewidth": 2.1,
            "mutation_scale": 15,
        },
        zorder=20,
    )


def plot_map_group(
    fig: plt.Figure,
    parent_spec,
    datasets: list[dict[str, object]],
    field_name: str,
    title: str,
    panel_label: str,
    rasterized: bool,
) -> None:
    group = parent_spec.subgridspec(
        2,
        3,
        height_ratios=[0.16, 0.84],
        hspace=0.02,
        wspace=0.055,
    )

    title_ax = fig.add_subplot(group[0, :])
    title_ax.axis("off")
    title_ax.text(
        0.5,
        0.66,
        title,
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
    )
    title_ax.text(
        -0.03,
        0.95,
        panel_label,
        transform=title_ax.transAxes,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        clip_on=False,
    )

    for column, data in enumerate(datasets):
        ax = fig.add_subplot(group[1, column])

        artist = ax.tripcolor(
            data["triangulation"],
            data[field_name],
            shading="flat",
            vmin=0.0,
            vmax=1.0,
            edgecolors="none",
        )
        artist.set_rasterized(rasterized)

        ax.set_title(
            data["case"].short_name,
            fontsize=15,
            fontweight="bold",
            pad=5,
        )

        add_flow_arrow(ax, data["case"].direction)
        clean_map_axis(ax)


def style_bar_axis(ax: plt.Axes) -> None:
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=13,
        width=1.5,
        length=6,
        direction="out",
    )

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

    for spine in ax.spines.values():
        spine.set_linewidth(1.4)

    ax.grid(
        axis="y",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
    )


def plot_area_partition(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
) -> None:
    names = [data["case"].short_name for data in datasets]
    active = np.array(
        [data["active_fraction"] for data in datasets],
        dtype=float,
    )
    intermediate = np.array(
        [data["intermediate_fraction"] for data in datasets],
        dtype=float,
    )
    dead = np.array(
        [data["dead_fraction"] for data in datasets],
        dtype=float,
    )

    positions = np.arange(len(names))

    active_bars = ax.bar(
        positions,
        100.0 * active,
        label="Active",
    )
    intermediate_bars = ax.bar(
        positions,
        100.0 * intermediate,
        bottom=100.0 * active,
        label="Intermediate",
    )
    dead_bars = ax.bar(
        positions,
        100.0 * dead,
        bottom=100.0 * (active + intermediate),
        label="Dead",
    )

    ax.set_xticks(positions, names)
    ax.set_ylim(0.0, 100.0)
    ax.set_ylabel(
        "Pore-domain area (%)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Hydraulic partition",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    # Rótulos apenas para frações relevantes, evitando excesso de texto.
    for bars, fractions in (
        (active_bars, active),
        (dead_bars, dead),
    ):
        for bar, fraction in zip(bars, fractions):
            if fraction < 0.035:
                continue

            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_y() + bar.get_height() / 2.0,
                f"{100.0 * fraction:.1f}",
                ha="center",
                va="center",
                fontsize=11.5,
                fontweight="bold",
            )

    legend = ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=12,
    )
    for text in legend.get_texts():
        text.set_fontweight("bold")

    style_bar_axis(ax)
    add_panel_label(ax, "c)")


def format_permeability(value: float) -> str:
    if abs(value) >= 10_000:
        return f"{value / 1000:.1f}k"
    if abs(value) >= 1000:
        return f"{value / 1000:.2f}k"
    if abs(value) >= 100:
        return f"{value:.1f}"
    return f"{value:.2f}"


def plot_permeability(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
) -> None:
    names = [data["case"].short_name for data in datasets]
    values = np.array(
        [data["permeability"] for data in datasets],
        dtype=float,
    )

    positions = np.arange(len(names))
    bars = []

    # Chamadas separadas fazem o Matplotlib usar automaticamente
    # cores sucessivas do ciclo padrão, sem fixar uma paleta.
    for position, value in zip(positions, values):
        bar = ax.bar(position, value)
        bars.append(bar[0])

    ax.set_yscale("log")
    ax.set_xticks(positions, names)
    ax.set_ylabel(
        r"$k_{\mathrm{app}}$ (pixel$^2$)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Directional permeability",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    for bar, value in zip(bars, values):
        ax.annotate(
            format_permeability(value),
            xy=(
                bar.get_x() + bar.get_width() / 2.0,
                value,
            ),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11.5,
            fontweight="bold",
        )

    # R3-H é o primeiro caso e R3-V é o terceiro.
    r3_horizontal = values[0]
    r3_vertical = values[2]
    anisotropy = (
        max(r3_horizontal, r3_vertical)
        / min(r3_horizontal, r3_vertical)
    )

    ax.text(
        0.96,
        0.08,
        rf"$A_k={anisotropy:.1f}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=14,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.92,
        },
    )

    style_bar_axis(ax)
    add_panel_label(ax, "d)")


def main() -> int:
    args = parse_args()

    if not args.results_root.is_dir():
        raise NotADirectoryError(
            f"Pasta de resultados não encontrada: {args.results_root}"
        )

    datasets = [
        load_case(args.results_root, case)
        for case in CASES
    ]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.linewidth": 1.4,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.07,
        }
    )

    fig = plt.figure(figsize=(15.0, 10.8))

    grid = GridSpec(
        nrows=2,
        ncols=2,
        figure=fig,
        left=0.065,
        right=0.975,
        bottom=0.105,
        top=0.965,
        wspace=0.25,
        hspace=0.31,
    )

    rasterized = not args.fully_vectorial

    plot_map_group(
        fig=fig,
        parent_spec=grid[0, 0],
        datasets=datasets,
        field_name="active_cells",
        title="Active-flow regions",
        panel_label="a)",
        rasterized=rasterized,
    )

    plot_map_group(
        fig=fig,
        parent_spec=grid[0, 1],
        datasets=datasets,
        field_name="dead_cells",
        title="Dead-flow regions",
        panel_label="b)",
        rasterized=rasterized,
    )

    area_ax = fig.add_subplot(grid[1, 0])
    plot_area_partition(area_ax, datasets)

    permeability_ax = fig.add_subplot(grid[1, 1])
    plot_permeability(permeability_ax, datasets)

    fig.canvas.draw()

    # Reforça o peso tipográfico após a criação dos ticks logarítmicos.
    for ax in (area_ax, permeability_ax):
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        args.output,
        format="svg",
        dpi=args.dpi,
        facecolor="white",
        metadata={
            "Title": "Hydraulic partition and permeability anisotropy",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    print(f"Figura salva em: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
