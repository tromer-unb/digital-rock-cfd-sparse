#!/usr/bin/env python3
"""
Gera a Figura 6 do artigo: síntese das propriedades hidráulicas obtidas
com a malha final convergida.

Painéis:
    (a) permeabilidade aparente;
    (b) gradiente médio de pressão;
    (c) porosidade do domínio e fração porosa hidraulicamente conectada;
    (d) número de Reynolds e amplificação máxima da velocidade.

A figura utiliza apenas o global_summary.csv da malha escolhida, normalmente
a malha extrafina.

Dependências:
    pip install numpy pandas matplotlib

Exemplo:
    python3 figure_6_hydraulic_properties.py \
        --summary convergence_mesh_05_extra_fine/global_summary.csv \
        --output Figure_6_hydraulic_properties.svg
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
    display_name: str
    filename: str
    direction: str
    marker: str


CASES = (
    FlowCase(
        short_name="R3-H",
        display_name="R3 horizontal",
        filename="patch_y7600_x53200_c0_mask.png",
        direction="lr",
        marker="o",
    ),
    FlowCase(
        short_name="R2-V",
        display_name="R2 vertical",
        filename="patch_y7600_x19000_c0_mask.png",
        direction="tb",
        marker="s",
    ),
    FlowCase(
        short_name="R3-V",
        display_name="R3 vertical",
        filename="patch_y7600_x53200_c0_mask.png",
        direction="tb",
        marker="^",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera a Figura 6 com a síntese das propriedades hidráulicas."
    )
    parser.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="global_summary.csv correspondente à malha final.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_6_hydraulic_properties.svg"),
        help="Arquivo SVG de saída.",
    )
    return parser.parse_args()


def load_final_results(summary_path: Path) -> pd.DataFrame:
    if not summary_path.is_file():
        raise FileNotFoundError(f"Resumo não encontrado: {summary_path}")

    data = pd.read_csv(summary_path)

    required = {
        "filename",
        "direction",
        "status",
        "apparent_permeability",
        "pressure_drop",
        "flow_length",
        "domain_porosity",
        "selected_fraction_of_cropped_pores",
        "hydraulic_diameter_2d",
        "reynolds_number",
        "speed_max",
        "inlet_velocity",
        "active_area_fraction",
        "dead_area_fraction",
    }

    missing = required.difference(data.columns)

    if missing:
        raise ValueError(
            f"{summary_path} não possui as colunas: {sorted(missing)}"
        )

    rows = []

    for case in CASES:
        selected = data.loc[
            (data["filename"] == case.filename)
            & (data["direction"] == case.direction)
            & (data["status"] == "success")
        ]

        if len(selected) != 1:
            raise RuntimeError(
                f"Esperava uma linha para {case.display_name}, "
                f"mas encontrei {len(selected)}."
            )

        row = selected.iloc[0].copy()
        row["case_short_name"] = case.short_name
        row["case_display_name"] = case.display_name
        row["case_marker"] = case.marker
        rows.append(row)

    result = pd.DataFrame(rows).reset_index(drop=True)

    numeric_columns = (
        "apparent_permeability",
        "pressure_drop",
        "flow_length",
        "domain_porosity",
        "selected_fraction_of_cropped_pores",
        "hydraulic_diameter_2d",
        "reynolds_number",
        "speed_max",
        "inlet_velocity",
        "active_area_fraction",
        "dead_area_fraction",
    )

    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="raise")

    result["pressure_gradient"] = (
        result["pressure_drop"] / result["flow_length"]
    )

    result["velocity_amplification"] = (
        result["speed_max"] / result["inlet_velocity"].abs()
    )

    return result


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


def style_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
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
        axis=grid_axis,
        which="major",
        linestyle="--",
        linewidth=0.8,
        alpha=0.35,
    )


def configure_log_y(ax: plt.Axes) -> None:
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    ax.yaxis.set_minor_formatter(NullFormatter())


def format_value(value: float) -> str:
    magnitude = abs(value)

    if magnitude >= 10_000:
        return f"{value / 1000:.1f}k"
    if magnitude >= 1000:
        return f"{value / 1000:.2f}k"
    if magnitude >= 100:
        return f"{value:.1f}"
    if magnitude >= 1:
        return f"{value:.2f}"
    if magnitude >= 0.01:
        return f"{value:.3f}"
    return f"{value:.2e}"


def add_bar_values(
    ax: plt.Axes,
    bars,
    values: np.ndarray,
    logarithmic: bool = False,
    fontsize: float = 11.5,
) -> None:
    for bar, value in zip(bars, values):
        if logarithmic:
            offset = 7
        else:
            offset = 6

        ax.annotate(
            format_value(float(value)),
            xy=(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height(),
            ),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            fontweight="bold",
        )


def plot_permeability(ax: plt.Axes, data: pd.DataFrame) -> None:
    positions = np.arange(len(data))
    values = data["apparent_permeability"].to_numpy(dtype=float)

    bars = []

    # Chamadas separadas usam automaticamente cores consecutivas
    # do ciclo padrão do Matplotlib.
    for position, value in zip(positions, values):
        bars.append(ax.bar(position, value, width=0.65)[0])

    configure_log_y(ax)

    ax.set_xticks(
        positions,
        data["case_short_name"],
    )
    ax.set_ylabel(
        r"Apparent permeability, $k_{\mathrm{app}}$ (pixel$^2$)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Hydraulic conductivity",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    add_bar_values(ax, bars, values, logarithmic=True)

    r3_h = values[0]
    r3_v = values[2]
    anisotropy = max(r3_h, r3_v) / min(r3_h, r3_v)

    ax.text(
        0.96,
        0.08,
        rf"$A_k={anisotropy:.1f}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=13.5,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.93,
        },
    )

    style_axis(ax)
    add_panel_label(ax, "a)")


def plot_pressure_gradient(ax: plt.Axes, data: pd.DataFrame) -> None:
    positions = np.arange(len(data))
    values = data["pressure_gradient"].to_numpy(dtype=float)

    bars = []

    for position, value in zip(positions, values):
        bars.append(ax.bar(position, value, width=0.65)[0])

    configure_log_y(ax)

    ax.set_xticks(
        positions,
        data["case_short_name"],
    )
    ax.set_ylabel(
        r"Mean pressure gradient, $\Delta p/L$",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Hydraulic resistance",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    add_bar_values(
        ax,
        bars,
        values,
        logarithmic=True,
        fontsize=10.8,
    )

    style_axis(ax)
    add_panel_label(ax, "b)")


def plot_connectivity(ax: plt.Axes, data: pd.DataFrame) -> None:
    positions = np.arange(len(data))
    width = 0.34

    porosity = (
        100.0
        * data["domain_porosity"].to_numpy(dtype=float)
    )
    connected = (
        100.0
        * data["selected_fraction_of_cropped_pores"].to_numpy(dtype=float)
    )

    bars_porosity = ax.bar(
        positions - width / 2.0,
        porosity,
        width=width,
        label=r"Domain porosity, $\phi$",
    )
    bars_connected = ax.bar(
        positions + width / 2.0,
        connected,
        width=width,
        label="Connected pore fraction",
    )

    ax.set_xticks(
        positions,
        data["case_short_name"],
    )
    ax.set_ylim(0.0, 110.0)
    ax.set_ylabel(
        "Geometric fraction (%)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Pore connectivity",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    add_bar_values(
        ax,
        bars_porosity,
        porosity,
        fontsize=10.8,
    )
    add_bar_values(
        ax,
        bars_connected,
        connected,
        fontsize=10.8,
    )

    diameters = data["hydraulic_diameter_2d"].to_numpy(dtype=float)

    diameter_text = "   ".join(
        rf"{name}: $D_h={value:.0f}$ px"
        for name, value in zip(
            data["case_short_name"],
            diameters,
        )
    )

    ax.text(
        0.5,
        0.04,
        diameter_text,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10.8,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.93,
        },
    )

    legend = ax.legend(
        loc="upper center",
        ncol=2,
        frameon=True,
        fontsize=11.2,
    )
    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(0.92)

    for text in legend.get_texts():
        text.set_fontweight("bold")

    style_axis(ax)
    add_panel_label(ax, "c)")


def plot_flow_regime(ax: plt.Axes, data: pd.DataFrame) -> None:
    positions = np.arange(len(data))

    reynolds = data["reynolds_number"].to_numpy(dtype=float)
    amplification = data["velocity_amplification"].to_numpy(dtype=float)

    for position, case, re_value, amp_value in zip(
        positions,
        CASES,
        reynolds,
        amplification,
    ):
        ax.plot(
            position,
            re_value,
            linestyle="none",
            marker=case.marker,
            markersize=10.5,
            markeredgewidth=1.4,
            label=(
                r"$Re$"
                if position == 0
                else None
            ),
        )

        ax.plot(
            position,
            amp_value,
            linestyle="none",
            marker=case.marker,
            markersize=10.5,
            markerfacecolor="none",
            markeredgewidth=2.0,
            label=(
                r"$|\mathbf{u}|_{\max}/U_0$"
                if position == 0
                else None
            ),
        )

        ax.annotate(
            f"{re_value:.3f}",
            xy=(position, re_value),
            xytext=(-8, -15),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=10.8,
            fontweight="bold",
        )

        ax.annotate(
            f"{amp_value:.1f}",
            xy=(position, amp_value),
            xytext=(7, 6),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=10.8,
            fontweight="bold",
        )

    configure_log_y(ax)

    ax.axhline(
        1.0,
        linestyle="--",
        linewidth=1.6,
        label=r"$Re=1$",
    )

    ax.set_xticks(
        positions,
        data["case_short_name"],
    )
    ax.set_xlim(-0.55, len(data) - 0.45)
    ax.set_ylabel(
        "Dimensionless magnitude",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Flow regime and localization",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    legend = ax.legend(
        loc="best",
        frameon=True,
        fontsize=11.2,
    )
    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(0.92)

    for text in legend.get_texts():
        text.set_fontweight("bold")

    ax.text(
        0.97,
        0.07,
        r"$Re<1$ for all cases",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12.5,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.27",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.93,
        },
    )

    style_axis(ax)
    add_panel_label(ax, "d)")


def main() -> int:
    args = parse_args()
    data = load_final_results(args.summary)

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
        left=0.095,
        right=0.975,
        bottom=0.09,
        top=0.96,
        wspace=0.30,
        hspace=0.31,
    )

    plot_permeability(axes[0, 0], data)
    plot_pressure_gradient(axes[0, 1], data)
    plot_connectivity(axes[1, 0], data)
    plot_flow_regime(axes[1, 1], data)

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
            "Title": "Hydraulic properties of the converged digital-rock models",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    print(f"Figura salva em: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
