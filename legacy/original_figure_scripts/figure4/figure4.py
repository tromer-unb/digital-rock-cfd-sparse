#!/usr/bin/env python3
"""
Gera a Figura 4 do artigo: heterogeneidade hidrodinâmica e localização
do transporte nos três casos percolantes.

Painéis:
    (a) perfil axial médio da pressão normalizada;
    (b) distribuição ponderada por área da velocidade normalizada;
    (c) fração acumulada da área porosa abaixo de uma dada velocidade;
    (d) localização do transporte ponderado por velocidade.

A figura é construída diretamente a partir de solution.vtk e summary.csv
da malha selecionada.

Dependências:
    pip install numpy pandas matplotlib meshio

Exemplo:
    python3 figure_4_flow_heterogeneity.py \
        --results-root convergence_mesh_05_extra_fine \
        --output Figure_4_flow_heterogeneity.svg
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter


@dataclass(frozen=True)
class FlowCase:
    short_name: str
    display_name: str
    filename_stem: str
    direction: str
    marker: str


CASES = (
    FlowCase(
        short_name="R3-H",
        display_name="R3 horizontal",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="lr",
        marker="o",
    ),
    FlowCase(
        short_name="R2-V",
        display_name="R2 vertical",
        filename_stem="patch_y7600_x19000_c0_mask",
        direction="tb",
        marker="s",
    ),
    FlowCase(
        short_name="R3-V",
        display_name="R3 vertical",
        filename_stem="patch_y7600_x53200_c0_mask",
        direction="tb",
        marker="^",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera a Figura 4 sobre heterogeneidade do escoamento."
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
        default=Path("Figure_4_flow_heterogeneity.svg"),
        help="Arquivo SVG de saída.",
    )
    parser.add_argument(
        "--axial-bins",
        type=int,
        default=60,
        help="Número de faixas para o perfil axial médio.",
    )
    parser.add_argument(
        "--speed-bins",
        type=int,
        default=55,
        help="Número de classes logarítmicas para a distribuição de velocidade.",
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


def triangle_areas(points_xy: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    p0 = points_xy[triangles[:, 0]]
    p1 = points_xy[triangles[:, 1]]
    p2 = points_xy[triangles[:, 2]]

    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def read_scalar(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.point_data:
        raise KeyError(f"Campo ausente no VTK: {name}")

    values = np.asarray(mesh.point_data[name], dtype=float)
    return np.squeeze(values).reshape(-1)


def read_speed(mesh: meshio.Mesh) -> np.ndarray:
    if "speed" in mesh.point_data:
        return read_scalar(mesh, "speed")

    if "velocity" not in mesh.point_data:
        raise KeyError("O VTK não contém 'speed' nem 'velocity'.")

    velocity = np.asarray(mesh.point_data["velocity"], dtype=float)

    if velocity.ndim != 2 or velocity.shape[1] < 2:
        raise ValueError(
            f"Formato inesperado para velocity: {velocity.shape}"
        )

    return np.linalg.norm(velocity[:, :2], axis=1)


def read_summary(sample_dir: Path) -> pd.Series:
    path = sample_dir / "summary.csv"

    if not path.is_file():
        raise FileNotFoundError(f"Resumo não encontrado: {path}")

    frame = pd.read_csv(path)

    if frame.empty:
        raise RuntimeError(f"Resumo vazio: {path}")

    required = {
        "pressure_outlet_mean",
        "pressure_drop",
        "inlet_velocity",
    }

    missing = required.difference(frame.columns)

    if missing:
        raise ValueError(
            f"{path} não possui as colunas: {sorted(missing)}"
        )

    return frame.iloc[0]


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: np.ndarray,
) -> np.ndarray:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]

    cumulative = np.cumsum(sorted_weights)
    cumulative /= cumulative[-1]

    return np.interp(quantiles, cumulative, sorted_values)


def weighted_histogram_density(
    values: np.ndarray,
    weights: np.ndarray,
    bins: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    counts, edges = np.histogram(
        values,
        bins=bins,
        weights=weights,
    )

    total_weight = np.sum(weights)
    widths = np.diff(edges)

    density = counts / np.maximum(
        total_weight * widths,
        np.finfo(float).eps,
    )

    centers = np.sqrt(edges[:-1] * edges[1:])
    return centers, density


def axial_pressure_profile(
    coordinate_star: np.ndarray,
    pressure_star: np.ndarray,
    weights: np.ndarray,
    number_of_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(0.0, 1.0, number_of_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    bin_ids = np.digitize(coordinate_star, edges) - 1
    bin_ids = np.clip(bin_ids, 0, number_of_bins - 1)

    numerator = np.bincount(
        bin_ids,
        weights=pressure_star * weights,
        minlength=number_of_bins,
    )
    denominator = np.bincount(
        bin_ids,
        weights=weights,
        minlength=number_of_bins,
    )

    profile = np.full(number_of_bins, np.nan, dtype=float)
    valid = denominator > 0.0
    profile[valid] = numerator[valid] / denominator[valid]

    return centers[valid], profile[valid]


def load_case(
    results_root: Path,
    case: FlowCase,
    axial_bins: int,
) -> dict[str, object]:
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

    points = np.asarray(mesh.points[:, :2], dtype=float)
    pressure_nodes = read_scalar(mesh, "pressure")
    speed_nodes = read_speed(mesh)
    summary = read_summary(sample_dir)

    areas = triangle_areas(points, triangles)

    valid = np.isfinite(areas) & (areas > 0.0)
    triangles = triangles[valid]
    areas = areas[valid]

    centroids = np.mean(points[triangles], axis=1)
    triangle_pressure = np.mean(pressure_nodes[triangles], axis=1)
    triangle_speed = np.mean(speed_nodes[triangles], axis=1)

    delta_p = float(summary["pressure_drop"])
    p_out = float(summary["pressure_outlet_mean"])
    inlet_velocity = abs(float(summary["inlet_velocity"]))

    if delta_p <= np.finfo(float).eps:
        raise RuntimeError(
            f"Queda de pressão inválida para {case.display_name}."
        )

    if inlet_velocity <= np.finfo(float).eps:
        raise RuntimeError(
            f"Velocidade de entrada inválida para {case.display_name}."
        )

    if case.direction == "lr":
        axial_coordinate = centroids[:, 0]
        coordinate_star = (
            axial_coordinate - np.min(axial_coordinate)
        ) / np.ptp(axial_coordinate)
    else:
        # Entrada no topo e saída na parte inferior:
        # x*=0 na entrada e x*=1 na saída.
        axial_coordinate = centroids[:, 1]
        coordinate_star = (
            np.max(axial_coordinate) - axial_coordinate
        ) / np.ptp(axial_coordinate)

    pressure_star = (triangle_pressure - p_out) / delta_p
    speed_star = triangle_speed / inlet_velocity

    pressure_star = np.clip(pressure_star, -0.05, 1.05)
    speed_star = np.maximum(speed_star, np.finfo(float).tiny)

    axial_x, axial_p = axial_pressure_profile(
        coordinate_star=coordinate_star,
        pressure_star=pressure_star,
        weights=areas,
        number_of_bins=axial_bins,
    )

    order_ascending = np.argsort(speed_star)
    speed_sorted_ascending = speed_star[order_ascending]
    area_sorted_ascending = areas[order_ascending]
    cumulative_area = np.cumsum(area_sorted_ascending)
    cumulative_area /= cumulative_area[-1]

    order_descending = np.argsort(speed_star)[::-1]
    speed_sorted_descending = speed_star[order_descending]
    area_sorted_descending = areas[order_descending]

    cumulative_ranked_area = np.cumsum(area_sorted_descending)
    cumulative_ranked_area /= cumulative_ranked_area[-1]

    transport_proxy = speed_sorted_descending * area_sorted_descending
    cumulative_transport = np.cumsum(transport_proxy)
    cumulative_transport /= cumulative_transport[-1]

    # Área necessária para carregar 50% e 90% do transporte ponderado.
    area_for_50 = float(
        np.interp(
            0.50,
            cumulative_transport,
            cumulative_ranked_area,
        )
    )
    area_for_90 = float(
        np.interp(
            0.90,
            cumulative_transport,
            cumulative_ranked_area,
        )
    )

    return {
        "case": case,
        "areas": areas,
        "speed_star": speed_star,
        "axial_x": axial_x,
        "axial_p": axial_p,
        "speed_sorted": speed_sorted_ascending,
        "cumulative_area": cumulative_area,
        "ranked_area": cumulative_ranked_area,
        "cumulative_transport": cumulative_transport,
        "area_for_50": area_for_50,
        "area_for_90": area_for_90,
    }


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


def add_legend(ax: plt.Axes, location: str = "best") -> None:
    legend = ax.legend(
        loc=location,
        frameon=True,
        fontsize=11.5,
        borderpad=0.45,
        handlelength=2.2,
    )

    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(0.92)

    for text in legend.get_texts():
        text.set_fontweight("bold")


def plot_pressure_profiles(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
) -> None:
    for data in datasets:
        ax.plot(
            data["axial_x"],
            data["axial_p"],
            linewidth=2.4,
            marker=data["case"].marker,
            markersize=5.5,
            markevery=max(1, len(data["axial_x"]) // 8),
            label=data["case"].short_name,
        )

    ax.plot(
        [0.0, 1.0],
        [1.0, 0.0],
        linestyle=":",
        linewidth=1.8,
        label="Linear reference",
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xlabel(
        r"Normalized axial coordinate, $s/L$",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        r"Area-averaged pressure, $\langle p^{*}\rangle$",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Axial pressure profile",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    style_axis(ax)
    add_legend(ax, "best")
    add_panel_label(ax, "a)")


def plot_speed_pdf(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
    speed_bins: int,
) -> None:
    global_min = min(
        float(weighted_quantile(
            data["speed_star"],
            data["areas"],
            np.array([0.001]),
        )[0])
        for data in datasets
    )
    global_max = max(
        float(weighted_quantile(
            data["speed_star"],
            data["areas"],
            np.array([0.999]),
        )[0])
        for data in datasets
    )

    global_min = max(global_min, 1e-5)
    global_max = max(global_max, global_min * 10.0)

    bins = np.logspace(
        np.log10(global_min),
        np.log10(global_max),
        speed_bins + 1,
    )

    for data in datasets:
        centers, density = weighted_histogram_density(
            values=data["speed_star"],
            weights=data["areas"],
            bins=bins,
        )

        positive = density > 0.0

        ax.plot(
            centers[positive],
            density[positive],
            linewidth=2.4,
            marker=data["case"].marker,
            markersize=5.0,
            markevery=max(1, np.count_nonzero(positive) // 8),
            label=data["case"].short_name,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(
        r"Normalized speed, $|\mathbf{u}|/U_0$",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        "Area-weighted probability density",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Velocity distribution",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    ax.xaxis.set_minor_formatter(NullFormatter())

    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    ax.yaxis.set_minor_formatter(NullFormatter())

    style_axis(ax)
    add_legend(ax, "best")
    add_panel_label(ax, "b)")


def plot_area_cdf(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
) -> None:
    for data in datasets:
        ax.plot(
            data["speed_sorted"],
            data["cumulative_area"],
            linewidth=2.4,
            marker=data["case"].marker,
            markersize=5.2,
            markevery=max(1, len(data["speed_sorted"]) // 10),
            label=data["case"].short_name,
        )

    ax.set_xscale("log")
    ax.set_xlim(left=1e-5)
    ax.set_ylim(0.0, 1.0)

    ax.set_xlabel(
        r"Normalized speed, $|\mathbf{u}|/U_0$",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        "Cumulative pore-area fraction",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Low-velocity pore occupancy",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    ax.xaxis.set_minor_formatter(NullFormatter())

    style_axis(ax)
    add_legend(ax, "lower right")
    add_panel_label(ax, "c)")


def plot_transport_localization(
    ax: plt.Axes,
    datasets: list[dict[str, object]],
) -> None:
    for data in datasets:
        ax.plot(
            data["ranked_area"],
            data["cumulative_transport"],
            linewidth=2.5,
            marker=data["case"].marker,
            markersize=5.3,
            markevery=max(1, len(data["ranked_area"]) // 10),
            label=data["case"].short_name,
        )

    ax.plot(
        [0.0, 1.0],
        [0.0, 1.0],
        linestyle=":",
        linewidth=1.8,
        label="Uniform transport",
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)

    ax.set_xlabel(
        "Cumulative pore-area fraction\n(sorted by decreasing speed)",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_ylabel(
        "Cumulative velocity-weighted transport",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_title(
        "Transport localization",
        fontsize=18,
        fontweight="bold",
        pad=10,
    )

    # Apenas dois números por caso: área que concentra 50% e 90%.
    annotation_lines = []
    for data in datasets:
        annotation_lines.append(
            rf"{data['case'].short_name}: "
            rf"$A_{{50}}={100.0 * data['area_for_50']:.1f}\%$, "
            rf"$A_{{90}}={100.0 * data['area_for_90']:.1f}\%$"
        )

    ax.text(
        0.97,
        0.07,
        "\n".join(annotation_lines),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10.8,
        fontweight="bold",
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 1.0,
            "alpha": 0.93,
        },
    )

    style_axis(ax)
    add_legend(ax, "lower right")
    add_panel_label(ax, "d)")


def main() -> int:
    args = parse_args()

    if not args.results_root.is_dir():
        raise NotADirectoryError(
            f"Pasta de resultados não encontrada: {args.results_root}"
        )

    datasets = [
        load_case(
            results_root=args.results_root,
            case=case,
            axial_bins=args.axial_bins,
        )
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
        left=0.085,
        right=0.975,
        bottom=0.09,
        top=0.96,
        wspace=0.28,
        hspace=0.31,
    )

    plot_pressure_profiles(axes[0, 0], datasets)
    plot_speed_pdf(axes[0, 1], datasets, args.speed_bins)
    plot_area_cdf(axes[1, 0], datasets)
    plot_transport_localization(axes[1, 1], datasets)

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
            "Title": "Pore-scale flow heterogeneity and transport localization",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    print(f"Figura salva em: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
