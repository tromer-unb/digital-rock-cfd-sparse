#!/usr/bin/env python3
"""
Gera a Figura 1 do artigo:

(a) imagens das três rochas digitais;
(b) convergência de malha de R3 na direção horizontal;
(c) convergência de malha de R2 na direção vertical;
(d) convergência de malha de R3 na direção vertical.

A saída é salva em SVG, mantendo textos e linhas editáveis como vetores.
As imagens das rochas permanecem rasterizadas dentro do arquivo SVG.

Dependências:
    pip install pandas numpy matplotlib pillow

Exemplo:
    python3 figure_1_convergence.py \
        --structures ../../../poros/poros \
        --summaries \
            convergence_mesh_01_coarse/global_summary.csv \
            convergence_mesh_02_intermediate/global_summary.csv \
            convergence_mesh_03_fine/global_summary.csv \
            convergence_mesh_04_very_fine/global_summary.csv \
            convergence_mesh_05_extra_fine/global_summary.csv \
        --output Figure_1_convergence.svg
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter
from PIL import Image, ImageOps


# ---------------------------------------------------------------------------
# IDENTIFICAÇÃO DAS AMOSTRAS
# ---------------------------------------------------------------------------

ROCKS = [
    {
        "label": "R1",
        "filename": "patch_y3800_x3800_c0_mask.png",
        "subtitle": "non-percolating",
    },
    {
        "label": "R2",
        "filename": "patch_y7600_x19000_c0_mask.png",
        "subtitle": "vertical",
    },
    {
        "label": "R3",
        "filename": "patch_y7600_x53200_c0_mask.png",
        "subtitle": "horizontal / vertical",
    },
]

CONVERGENCE_CASES = [
    {
        "panel": "b)",
        "title": "R3 — horizontal",
        "filename": "patch_y7600_x53200_c0_mask.png",
        "direction": "lr",
        "marker": "o",
    },
    {
        "panel": "c)",
        "title": "R2 — vertical",
        "filename": "patch_y7600_x19000_c0_mask.png",
        "direction": "tb",
        "marker": "s",
    },
    {
        "panel": "d)",
        "title": "R3 — vertical",
        "filename": "patch_y7600_x53200_c0_mask.png",
        "direction": "tb",
        "marker": "^",
    },
]


# ---------------------------------------------------------------------------
# LEITURA DOS DADOS
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera a Figura 1 com amostras e convergência de malha."
    )
    parser.add_argument(
        "--structures",
        type=Path,
        required=True,
        help="Pasta contendo as três imagens das rochas.",
    )
    parser.add_argument(
        "--summaries",
        type=Path,
        nargs="+",
        required=True,
        help="Arquivos global_summary.csv dos níveis de refinamento.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_1_convergence.svg"),
        help="Arquivo SVG de saída.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolução usada para os elementos rasterizados no SVG.",
    )
    return parser.parse_args()


def load_summary_files(summary_paths: list[Path]) -> pd.DataFrame:
    """Lê e combina os resultados dos diferentes níveis de malha."""
    frames: list[pd.DataFrame] = []

    required_columns = {
        "filename",
        "direction",
        "status",
        "gmsh_triangle_count",
        "apparent_permeability",
    }

    for level, path in enumerate(summary_paths, start=1):
        if not path.is_file():
            raise FileNotFoundError(f"Resumo não encontrado: {path}")

        frame = pd.read_csv(path)

        missing = required_columns.difference(frame.columns)
        if missing:
            raise ValueError(
                f"O arquivo {path} não possui as colunas: {sorted(missing)}"
            )

        frame = frame.copy()
        frame["mesh_level"] = level
        frame["summary_path"] = str(path)
        frames.append(frame)

    data = pd.concat(frames, ignore_index=True)

    data["gmsh_triangle_count"] = pd.to_numeric(
        data["gmsh_triangle_count"], errors="coerce"
    )
    data["apparent_permeability"] = pd.to_numeric(
        data["apparent_permeability"], errors="coerce"
    )

    return data


def extract_curve(
    data: pd.DataFrame,
    filename: str,
    direction: str,
) -> pd.DataFrame:
    """Seleciona e ordena uma curva de convergência."""
    curve = data.loc[
        (data["filename"] == filename)
        & (data["direction"] == direction)
        & (data["status"] == "success"),
        [
            "mesh_level",
            "gmsh_triangle_count",
            "apparent_permeability",
        ],
    ].dropna()

    curve = curve.sort_values("gmsh_triangle_count").reset_index(drop=True)

    if len(curve) < 2:
        raise RuntimeError(
            f"Dados insuficientes para {filename}, direção={direction}."
        )

    return curve


def load_display_image(path: Path, max_side: int = 1600) -> np.ndarray:
    """
    Lê a imagem e limita seu tamanho para evitar um SVG excessivamente grande.
    O conteúdo visual e a proporção são preservados.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Imagem não encontrada: {path}")

    image = Image.open(path).convert("L")
    image = ImageOps.autocontrast(image)

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    image.thumbnail((max_side, max_side), resampling)

    return np.asarray(image)


# ---------------------------------------------------------------------------
# FORMATAÇÃO GRÁFICA
# ---------------------------------------------------------------------------

def format_permeability(value: float) -> str:
    """Formatação compacta para os valores mostrados sobre os pontos."""
    magnitude = abs(value)

    if magnitude >= 10_000:
        return f"{value / 1000:.1f}k"
    if magnitude >= 1000:
        return f"{value / 1000:.2f}k"
    if magnitude >= 100:
        return f"{value:.1f}"
    return f"{value:.2f}"


def set_bold_tick_labels(ax: plt.Axes, fontsize: int = 13) -> None:
    """Aplica fonte grande e em negrito aos números dos eixos."""
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=fontsize,
        width=1.6,
        length=6,
        direction="out",
    )
    ax.tick_params(
        axis="both",
        which="minor",
        width=1.2,
        length=3,
        direction="out",
    )

    # Necessário após o canvas ser desenhado para garantir os rótulos atuais.
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.10,
        label,
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def style_convergence_axis(ax: plt.Axes) -> None:
    ax.set_xscale("log")

    ax.set_xlabel(
        "Number of triangles",
        fontsize=15,
        fontweight="bold",
        labelpad=8,
    )
    ax.set_ylabel(
        r"Apparent permeability, $k_{\mathrm{app}}$ (pixel$^2$)",
        fontsize=15,
        fontweight="bold",
        labelpad=8,
    )

    ax.grid(
        True,
        which="major",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
    )

    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=5))
    ax.xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=20)
    )
    ax.xaxis.set_minor_formatter(NullFormatter())

    y_formatter = ScalarFormatter(useMathText=True)
    y_formatter.set_scientific(False)
    y_formatter.set_useOffset(False)
    ax.yaxis.set_major_formatter(y_formatter)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    set_bold_tick_labels(ax)


def plot_convergence(
    ax: plt.Axes,
    curve: pd.DataFrame,
    panel: str,
    title: str,
    marker: str,
) -> None:
    triangles = curve["gmsh_triangle_count"].to_numpy(dtype=float)
    permeability = curve["apparent_permeability"].to_numpy(dtype=float)

    # A cor é obtida automaticamente do ciclo padrão do Matplotlib.
    line = ax.plot(
        triangles,
        permeability,
        marker=marker,
        markersize=8.5,
        markeredgewidth=1.3,
        linewidth=2.4,
        zorder=3,
    )[0]

    line_color = line.get_color()

    # Margens para impedir que rótulos e marcadores encostem nas bordas.
    x_min = triangles.min() / 1.35
    x_max = triangles.max() * 1.35
    ax.set_xlim(x_min, x_max)

    y_range = np.ptp(permeability)
    if y_range <= 0:
        y_range = max(abs(permeability.mean()) * 0.08, 1.0)

    lower = permeability.min() - 0.30 * y_range
    upper = permeability.max() + 0.42 * y_range
    ax.set_ylim(lower, upper)

    for x_value, y_value in zip(triangles, permeability):
        ax.annotate(
            format_permeability(y_value),
            xy=(x_value, y_value),
            xytext=(0, 11),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11.5,
            fontweight="bold",
            color=line_color,
            clip_on=False,
        )

    relative_change = (
        abs(permeability[-1] - permeability[-2])
        / max(abs(permeability[-1]), np.finfo(float).eps)
        * 100.0
    )

    ax.text(
        0.97,
        0.08,
        rf"$\Delta k_{{\mathrm{{last}}}}={relative_change:.2f}\%$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=13,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.92,
        },
    )

    ax.set_title(
        title,
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    style_convergence_axis(ax)
    add_panel_label(ax, panel)


def plot_rock_panel(
    fig: plt.Figure,
    parent_spec,
    structures_dir: Path,
) -> None:
    """
    Cria o painel (a) com as três rochas em uma única região da figura.
    """
    panel_grid = parent_spec.subgridspec(
        2,
        3,
        height_ratios=[0.16, 0.84],
        hspace=0.03,
        wspace=0.06,
    )

    title_ax = fig.add_subplot(panel_grid[0, :])
    title_ax.axis("off")

    title_ax.text(
        0.5,
        0.72,
        "Digital rock samples",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
    )
    title_ax.text(
        -0.04,
        0.95,
        "a)",
        transform=title_ax.transAxes,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        clip_on=False,
    )

    for column, rock in enumerate(ROCKS):
        ax = fig.add_subplot(panel_grid[1, column])

        image_path = structures_dir / rock["filename"]
        image = load_display_image(image_path)

        ax.imshow(
            image,
            cmap="gray",
            interpolation="nearest",
            vmin=0,
            vmax=255,
        )

        ax.set_title(
            rock["label"],
            fontsize=17,
            fontweight="bold",
            pad=6,
        )

        ax.text(
            0.5,
            -0.045,
            rock["subtitle"],
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11.5,
            fontweight="bold",
        )

        ax.set_xticks([])
        ax.set_yticks([])

        for spine in ax.spines.values():
            spine.set_linewidth(1.4)


# ---------------------------------------------------------------------------
# FIGURA PRINCIPAL
# ---------------------------------------------------------------------------

def create_figure(
    structures_dir: Path,
    summary_paths: list[Path],
    output_path: Path,
    dpi: int,
) -> None:
    if not structures_dir.is_dir():
        raise NotADirectoryError(
            f"Pasta de estruturas não encontrada: {structures_dir}"
        )

    data = load_summary_files(summary_paths)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.linewidth": 1.5,
            "svg.fonttype": "none",  # texto editável no Illustrator/Inkscape
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )

    fig = plt.figure(figsize=(15.5, 10.8))

    grid = GridSpec(
        nrows=2,
        ncols=2,
        figure=fig,
        width_ratios=[1.0, 1.14],
        height_ratios=[1.0, 1.0],
        left=0.055,
        right=0.985,
        bottom=0.075,
        top=0.965,
        wspace=0.27,
        hspace=0.34,
    )

    plot_rock_panel(
        fig=fig,
        parent_spec=grid[0, 0],
        structures_dir=structures_dir,
    )

    axes = [
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]

    for ax, case in zip(axes, CONVERGENCE_CASES):
        curve = extract_curve(
            data,
            filename=case["filename"],
            direction=case["direction"],
        )
        plot_convergence(
            ax=ax,
            curve=curve,
            panel=case["panel"],
            title=case["title"],
            marker=case["marker"],
        )

    # Atualiza os tick labels antes de forçar o peso em negrito.
    fig.canvas.draw()
    for ax in axes:
        set_bold_tick_labels(ax)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_path,
        format="svg",
        dpi=dpi,
        facecolor="white",
        metadata={
            "Title": "Digital rock samples and mesh-convergence analysis",
            "Creator": "Python/Matplotlib",
        },
    )
    plt.close(fig)

    print(f"Figura salva em: {output_path.resolve()}")


def main() -> int:
    args = parse_args()

    create_figure(
        structures_dir=args.structures,
        summary_paths=args.summaries,
        output_path=args.output,
        dpi=args.dpi,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
