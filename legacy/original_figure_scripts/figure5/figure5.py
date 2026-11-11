#!/usr/bin/env python3
"""
Gera a Figura 5 do artigo: custo computacional, escalabilidade e compromisso
entre precisão e memória para os três casos percolantes.

Painéis:
    (a) graus de liberdade versus número de triângulos;
    (b) memória CSR do sistema versus graus de liberdade;
    (c) erro relativo da permeabilidade em relação à malha extrafina;
    (d) custo incremental da malha extrafina versus alteração em permeabilidade.

A figura foi concebida para sustentar a contribuição de baixo custo
computacional do pipeline esparso.

Dependências:
    pip install numpy pandas matplotlib

Exemplo:
    python3 figure_5_computational_efficiency.py \
        --summaries \
            convergence_mesh_01_coarse/global_summary.csv \
            convergence_mesh_02_intermediate/global_summary.csv \
            convergence_mesh_03_fine/global_summary.csv \
            convergence_mesh_04_very_fine/global_summary.csv \
            convergence_mesh_05_extra_fine/global_summary.csv \
        --output Figure_5_computational_efficiency.svg
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator, NullFormatter


@dataclass(frozen=True)
class FlowCase:
    short_name: str
    filename: str
    direction: str
    marker: str


CASES = (
    FlowCase(
        short_name="R3-H",
        filename="patch_y7600_x53200_c0_mask.png",
        direction="lr",
        marker="o",
    ),
    FlowCase(
        short_name="R2-V",
        filename="patch_y7600_x19000_c0_mask.png",
        direction="tb",
        marker="s",
    ),
    FlowCase(
        short_name="R3-V",
        filename="patch_y7600_x53200_c0_mask.png",
        direction="tb",
        marker="^",
    ),
)

MESH_NAMES = (
    "Coarse",
    "Intermediate",
    "Fine",
    "Very fine",
    "Extra fine",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera a Figura 5 sobre eficiência computacional."
    )
    parser.add_argument(
        "--summaries",
        type=Path,
        nargs="+",
        required=True,
        help=(
            "Arquivos global_summary.csv em ordem crescente de refinamento."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_5_computational_efficiency.svg"),
        help="Arquivo SVG de saída.",
    )
    return parser.parse_args()


def load_summaries(paths: list[Path]) -> pd.DataFrame:
    if len(paths) < 3:
        raise ValueError(
            "Forneça pelo menos três níveis sucessivos de refinamento."
        )

    frames: list[pd.DataFrame] = []

    required = {
        "filename",
        "direction",
        "status",
        "gmsh_triangle_count",
        "total_dof_count",
        "system_csr_megabytes",
        "apparent_permeability",
    }

    for level, path in enumerate(paths, start=1):
        if not path.is_file():
            raise FileNotFoundError(f"Resumo não encontrado: {path}")

        frame = pd.read_csv(path)
        missing = required.difference(frame.columns)

        if missing:
            raise ValueError(
                f"{path} não contém as colunas: {sorted(missing)}"
            )

        frame = frame.copy()
        frame["mesh_level"] = level
        frame["mesh_name"] = (
            MESH_NAMES[level - 1]
            if level <= len(MESH_NAMES)
            else f"Level {level}"
        )
        frames.append(frame)

    data = pd.concat(frames, ignore_index=True)

    numeric_columns = (
        "gmsh_triangle_count",
        "total_dof_count",
        "system_csr_megabytes",
        "apparent_permeability",
    )

    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    return data


def extract_case(data: pd.DataFrame, case: FlowCase) -> pd.DataFrame:
    curve = data.loc[
        (data["filename"] == case.filename)
        & (data["direction"] == case.direction)
        & (data["status"] == "success"),
        [
            "mesh_level",
            "mesh_name",
            "gmsh_triangle_count",
            "total_dof_count",
            "system_csr_megabytes",
            "apparent_permeability",
        ],
    ].dropna()

    curve = curve.sort_values("mesh_level").reset_index(drop=True)

    if len(curve) < 3:
        raise RuntimeError(
            f"Dados insuficientes para {case.short_name}."
        )

    return curve


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.10,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        clip_on=False,
    )


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=13,
        width=1.5,
        length=6,
        direction="out",
    )
    ax.tick_params(
        axis="both",
        which="minor",
        width=1.1,
        length=3,
        direction="out",
    )

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

    for spine in ax.spines.values():
        spine.set_linewidth(1.4)

    ax.grid(
        True,
        which="major",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
    )


def configure_log_axis(ax: plt.Axes, axis: str) -> None:
    locator_major = LogLocator(base=10.0)
    locator_minor = LogLocator(
        base=10.0,
        subs=np.arange(2, 10) * 0.1,
    )

    if axis == "x":
        ax.xaxis.set_major_locator(locator_major)
        ax.xaxis.set_minor_locator(locator_minor)
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_major_locator(locator_major)
        ax.yaxis.set_minor_locator(locator_minor)
        ax.yaxis.set_minor_formatter(NullFormatter())


def add_legend(ax: plt.Axes, location: str = "best") -> None:
    legend = ax.legend(
        loc=location,
        frameon=True,
        fontsize=11.5,
        handlelength=2.3,
        borderpad=0.45,
    )

    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(0.92)

    for text in legend.get_texts():
        text.set_fontweight("bold")


def plot_dofs_vs_triangles(
    ax: plt.Axes,
    datasets: list[tuple[FlowCase, pd.DataFrame]],
) -> None:
    exponents = []

    for case, curve in datasets:
        triangles = curve["gmsh_triangle_count"].to_numpy(dtype=float)
        dofs = curve["total_dof_count"].to_numpy(dtype=float)

        ax.plot(
            triangles,
            dofs,
            linewidth=2.4,
            marker=case.marker,
            markersize=7.2,
            markeredgewidth=1.2,
            label=case.short_name,
        )

        exponent = np.polyfit(
            np.log10(triangles),
            np.log10(dofs),
            deg=1,
        )[0]
        exponents.append(exponent)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(
        "Number of triangles",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        "Total degrees of freedom",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Discretization scaling",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    ax.text(
        0.96,
        0.08,
        rf"median slope $={np.median(exponents):.2f}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12.5,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.92,
        },
    )

    configure_log_axis(ax, "x")
    configure_log_axis(ax, "y")
    style_axis(ax)
    add_legend(ax)
    add_panel_label(ax, "a)")


def plot_memory_vs_dofs(
    ax: plt.Axes,
    datasets: list[tuple[FlowCase, pd.DataFrame]],
) -> None:
    memory_rates = []

    for case, curve in datasets:
        dofs = curve["total_dof_count"].to_numpy(dtype=float)
        memory = curve["system_csr_megabytes"].to_numpy(dtype=float)

        ax.plot(
            dofs,
            memory,
            linewidth=2.4,
            marker=case.marker,
            markersize=7.2,
            markeredgewidth=1.2,
            label=case.short_name,
        )

        memory_rates.extend(memory / dofs * 1.0e6)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(
        "Total degrees of freedom",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        "Sparse-system storage (MB)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Sparse-memory scaling",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    ax.text(
        0.96,
        0.08,
        rf"median $\approx {np.median(memory_rates):.0f}$ MB/MDOF",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12.5,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.92,
        },
    )

    configure_log_axis(ax, "x")
    configure_log_axis(ax, "y")
    style_axis(ax)
    add_legend(ax)
    add_panel_label(ax, "b)")


def plot_error_vs_memory(
    ax: plt.Axes,
    datasets: list[tuple[FlowCase, pd.DataFrame]],
) -> None:
    for case, curve in datasets:
        memory = curve["system_csr_megabytes"].to_numpy(dtype=float)
        permeability = curve["apparent_permeability"].to_numpy(dtype=float)

        reference = permeability[-1]
        error = (
            np.abs(permeability - reference)
            / max(abs(reference), np.finfo(float).eps)
            * 100.0
        )

        ax.plot(
            memory,
            error,
            linewidth=2.4,
            marker=case.marker,
            markersize=7.2,
            markeredgewidth=1.2,
            label=case.short_name,
        )

    ax.set_xscale("log")
    ax.set_ylim(bottom=0.0)

    ax.set_xlabel(
        "Sparse-system storage (MB)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        r"Relative error in $k_{\mathrm{app}}$ (\%)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Accuracy–memory trade-off",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    ax.axhline(
        1.0,
        linestyle=":",
        linewidth=1.7,
        label="1% criterion",
    )
    ax.axhline(
        2.0,
        linestyle="--",
        linewidth=1.5,
        label="2% criterion",
    )

    configure_log_axis(ax, "x")
    style_axis(ax)
    add_legend(ax)
    add_panel_label(ax, "c)")


def plot_incremental_cost(
    ax: plt.Axes,
    datasets: list[tuple[FlowCase, pd.DataFrame]],
) -> None:
    x_values = []
    y_values = []

    for case, curve in datasets:
        if len(curve) < 2:
            continue

        memory = curve["system_csr_megabytes"].to_numpy(dtype=float)
        permeability = curve["apparent_permeability"].to_numpy(dtype=float)

        memory_increase = (
            (memory[-1] - memory[-2])
            / memory[-2]
            * 100.0
        )
        permeability_change = (
            abs(permeability[-1] - permeability[-2])
            / max(abs(permeability[-1]), np.finfo(float).eps)
            * 100.0
        )

        x_values.append(memory_increase)
        y_values.append(permeability_change)

        ax.plot(
            memory_increase,
            permeability_change,
            linestyle="none",
            marker=case.marker,
            markersize=11,
            markeredgewidth=1.5,
            label=case.short_name,
        )

        ax.annotate(
            case.short_name,
            xy=(memory_increase, permeability_change),
            xytext=(7, 7),
            textcoords="offset points",
            fontsize=12,
            fontweight="bold",
        )

    ax.axhline(
        1.0,
        linestyle=":",
        linewidth=1.7,
        label="1% criterion",
    )
    ax.axhline(
        2.0,
        linestyle="--",
        linewidth=1.5,
        label="2% criterion",
    )

    if x_values:
        x_margin = max(np.ptp(x_values), 1.0) * 0.45
        ax.set_xlim(
            min(x_values) - x_margin,
            max(x_values) + x_margin,
        )

    ax.set_ylim(
        0.0,
        max(2.25, max(y_values, default=0.0) * 1.35),
    )

    ax.set_xlabel(
        "Additional CSR memory from\nvery fine to extra fine (%)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        r"Change in $k_{\mathrm{app}}$ (\%)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Diminishing return of refinement",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    style_axis(ax)
    add_legend(ax, "upper right")
    add_panel_label(ax, "d)")


def main() -> int:
    args = parse_args()

    data = load_summaries(args.summaries)

    datasets = [
        (case, extract_case(data, case))
        for case in CASES
    ]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.linewidth": 1.4,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14.5, 10.8),
    )

    fig.subplots_adjust(
        left=0.09,
        right=0.975,
        bottom=0.09,
        top=0.96,
        wspace=0.29,
        hspace=0.31,
    )

    plot_dofs_vs_triangles(axes[0, 0], datasets)
    plot_memory_vs_dofs(axes[0, 1], datasets)
    plot_error_vs_memory(axes[1, 0], datasets)
    plot_incremental_cost(axes[1, 1], datasets)

    fig.canvas.draw()

    for ax in axes.ravel():
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        args.output,
        format="svg",
        facecolor="white",
        metadata={
            "Title": "Computational efficiency and accuracy of the sparse pipeline",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    print(f"Figura salva em: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
