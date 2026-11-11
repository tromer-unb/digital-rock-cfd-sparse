#!/usr/bin/env python3
"""
Gera a Figura 7 do artigo: visualização qualitativa do escoamento em estilo
CFD/OpenFOAM-like, com domínio poroso, campo de velocidade e streamlines.

Layout:
    três painéis horizontais:
        (a) R3 horizontal
        (b) R2 vertical
        (c) R3 vertical

Cada painel mostra:
    - domínio poroso percolante em fundo neutro;
    - mapa do módulo de velocidade normalizado, |u|/U0;
    - linhas de corrente (streamlines) sobrepostas;
    - seta indicando a direção média do escoamento.

Leitura:
    O script usa diretamente:
        solution.vtk
        summary.csv
    da malha selecionada (tipicamente a malha extrafina).

Dependências:
    pip install numpy pandas matplotlib meshio scipy

Exemplo:
    python3 figure_7_streamlines_visual.py \
        --results-root convergence_mesh_05_extra_fine \
        --output Figure_7_streamlines_visual.svg

Opcionalmente:
    --seed-density 26
    --streamline-density 1.15
    --rasterize-field

Observação:
    Por padrão, o mapa colorido é mantido rasterizado dentro do SVG para
    evitar arquivos enormes; linhas, eixos, textos e streamlines permanecem
    vetoriais.
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
from scipy.interpolate import LinearNDInterpolator


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
        description="Gera a Figura 7 com streamlines e campo de velocidade."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Pasta de resultados da malha escolhida.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Figure_7_streamlines_visual.svg"),
        help="Arquivo SVG de saída.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=350,
        help="Resolução do campo rasterizado dentro do SVG.",
    )
    parser.add_argument(
        "--grid-resolution",
        type=int,
        default=260,
        help=(
            "Número de pontos da maior dimensão do grid regular usado "
            "para interpolar o campo de velocidade."
        ),
    )
    parser.add_argument(
        "--seed-density",
        type=int,
        default=24,
        help=(
            "Número-base de sementes para streamlines ao longo da entrada."
        ),
    )
    parser.add_argument(
        "--streamline-density",
        type=float,
        default=1.05,
        help=(
            "Densidade usada pelo matplotlib.streamplot."
        ),
    )
    parser.add_argument(
        "--rasterize-field",
        action="store_true",
        help=(
            "Rasteriza apenas o mapa colorido; streamlines, textos e "
            "contornos permanecem vetoriais."
        ),
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


def read_scalar(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.point_data:
        raise KeyError(f"Campo '{name}' ausente no VTK.")

    values = np.asarray(mesh.point_data[name], dtype=float)
    return np.squeeze(values).reshape(-1)


def read_velocity(mesh: meshio.Mesh) -> np.ndarray:
    if "velocity" not in mesh.point_data:
        raise KeyError("O arquivo VTK não contém o campo 'velocity'.")

    values = np.asarray(mesh.point_data["velocity"], dtype=float)

    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(
            f"Formato inesperado para velocity: {values.shape}"
        )

    return values[:, :2]


def read_summary(sample_dir: Path) -> pd.Series:
    summary_path = sample_dir / "summary.csv"

    if not summary_path.is_file():
        raise FileNotFoundError(f"Resumo não encontrado: {summary_path}")

    frame = pd.read_csv(summary_path)

    if frame.empty:
        raise RuntimeError(f"Resumo vazio: {summary_path}")

    required = {
        "inlet_velocity",
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

    points = np.asarray(mesh.points[:, :2], dtype=float)
    x = points[:, 0]
    y = points[:, 1]

    velocity = read_velocity(mesh)
    speed = np.linalg.norm(velocity, axis=1)

    summary = read_summary(sample_dir)
    inlet_velocity = abs(float(summary["inlet_velocity"]))

    if inlet_velocity <= np.finfo(float).eps:
        raise RuntimeError(
            f"Velocidade de entrada inválida para {case.display_name}."
        )

    speed_star = speed / inlet_velocity
    velocity_star = velocity / inlet_velocity

    triangulation = mtri.Triangulation(x, y, triangles)

    return {
        "case": case,
        "sample_dir": sample_dir,
        "points": points,
        "triangles": triangles,
        "triangulation": triangulation,
        "velocity_star": velocity_star,
        "speed_star": speed_star,
    }


def make_regular_grid(
    points: np.ndarray,
    velocity: np.ndarray,
    grid_resolution: int,
) -> dict[str, np.ndarray]:
    x = points[:, 0]
    y = points[:, 1]

    width = float(np.max(x) - np.min(x))
    height = float(np.max(y) - np.min(y))

    if width >= height:
        nx = grid_resolution
        ny = max(60, int(round(grid_resolution * height / width)))
    else:
        ny = grid_resolution
        nx = max(60, int(round(grid_resolution * width / height)))

    x_grid = np.linspace(np.min(x), np.max(x), nx)
    y_grid = np.linspace(np.min(y), np.max(y), ny)

    xx, yy = np.meshgrid(x_grid, y_grid)

    interp_u = LinearNDInterpolator(points, velocity[:, 0], fill_value=np.nan)
    interp_v = LinearNDInterpolator(points, velocity[:, 1], fill_value=np.nan)

    uu = interp_u(xx, yy)
    vv = interp_v(xx, yy)

    speed = np.sqrt(uu**2 + vv**2)

    return {
        "x_grid": x_grid,
        "y_grid": y_grid,
        "xx": xx,
        "yy": yy,
        "u_grid": uu,
        "v_grid": vv,
        "speed_grid": speed,
    }


def generate_streamline_seeds(
    case: FlowCase,
    grid: dict[str, np.ndarray],
    seed_density: int,
) -> np.ndarray:
    x_grid = grid["x_grid"]
    y_grid = grid["y_grid"]
    u_grid = grid["u_grid"]
    v_grid = grid["v_grid"]

    seeds = []

    if case.direction == "lr":
        x_seed = x_grid[1] if len(x_grid) > 1 else x_grid[0]
        y_candidates = np.linspace(y_grid.min(), y_grid.max(), seed_density)

        x_index = 1 if len(x_grid) > 1 else 0
        for y_value in y_candidates:
            y_index = int(np.argmin(np.abs(y_grid - y_value)))

            u_val = u_grid[y_index, x_index]
            v_val = v_grid[y_index, x_index]

            if np.isfinite(u_val) and np.isfinite(v_val):
                if np.hypot(u_val, v_val) > 1e-6:
                    seeds.append([x_seed, y_grid[y_index]])
    else:
        y_seed = y_grid[-2] if len(y_grid) > 1 else y_grid[-1]
        x_candidates = np.linspace(x_grid.min(), x_grid.max(), seed_density)

        y_index = len(y_grid) - 2 if len(y_grid) > 1 else len(y_grid) - 1
        for x_value in x_candidates:
            x_index = int(np.argmin(np.abs(x_grid - x_value)))

            u_val = u_grid[y_index, x_index]
            v_val = v_grid[y_index, x_index]

            if np.isfinite(u_val) and np.isfinite(v_val):
                if np.hypot(u_val, v_val) > 1e-6:
                    seeds.append([x_grid[x_index], y_seed])

    if not seeds:
        # fallback seguro: uma semente central
        seeds = [[np.mean(x_grid), np.mean(y_grid)]]

    return np.asarray(seeds, dtype=float)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=21,
        fontweight="bold",
        bbox={
            "boxstyle": "square,pad=0.12",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.88,
        },
        zorder=20,
    )


def add_flow_arrow(ax: plt.Axes, direction: str) -> None:
    if direction == "lr":
        start = (0.12, 0.09)
        end = (0.33, 0.09)
    else:
        start = (0.10, 0.88)
        end = (0.10, 0.67)

    ax.annotate(
        "",
        xy=end,
        xytext=start,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "linewidth": 2.3,
            "mutation_scale": 17,
        },
        zorder=25,
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(1.35)


def plot_case_panel(
    ax: plt.Axes,
    data: dict[str, object],
    panel_label: str,
    title: str,
    grid_resolution: int,
    seed_density: int,
    streamline_density: float,
    rasterize_field: bool,
    speed_norm_max: float,
):
    triangulation = data["triangulation"]
    points = data["points"]
    velocity_star = data["velocity_star"]
    speed_star = data["speed_star"]

    field_artist = ax.tripcolor(
        triangulation,
        speed_star,
        shading="gouraud",
        vmin=0.0,
        vmax=speed_norm_max,
        edgecolors="none",
        zorder=1,
    )
    field_artist.set_rasterized(rasterize_field)

    grid = make_regular_grid(
        points=points,
        velocity=velocity_star,
        grid_resolution=grid_resolution,
    )

    seeds = generate_streamline_seeds(
        case=data["case"],
        grid=grid,
        seed_density=seed_density,
    )

    mask = ~np.isfinite(grid["u_grid"]) | ~np.isfinite(grid["v_grid"])
    u_masked = np.ma.array(grid["u_grid"], mask=mask)
    v_masked = np.ma.array(grid["v_grid"], mask=mask)

    ax.streamplot(
        grid["x_grid"],
        grid["y_grid"],
        u_masked,
        v_masked,
        start_points=seeds,
        density=streamline_density,
        color="white",
        linewidth=0.95,
        arrowsize=0.95,
        maxlength=7.0,
        integration_direction="forward",
        zorder=8,
    )

    # Contorno do domínio para melhorar o contraste visual.
    ax.triplot(
        triangulation,
        linewidth=0.14,
        alpha=0.18,
        zorder=3,
    )

    ax.set_title(
        title,
        fontsize=17,
        fontweight="bold",
        pad=8,
    )

    add_panel_label(ax, panel_label)
    add_flow_arrow(ax, data["case"].direction)
    clean_axis(ax)

    return field_artist


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

    global_speed_max = max(
        float(np.nanpercentile(data["speed_star"], 99.6))
        for data in datasets
    )
    global_speed_max = max(global_speed_max, 1.0)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.linewidth": 1.35,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.06,
        }
    )

    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(16.0, 5.7),
        constrained_layout=False,
    )

    fig.subplots_adjust(
        left=0.045,
        right=0.985,
        bottom=0.13,
        top=0.90,
        wspace=0.065,
    )

    panel_labels = ("a)", "b)", "c)")
    field_artist = None

    for ax, data, label in zip(axes, datasets, panel_labels):
        field_artist = plot_case_panel(
            ax=ax,
            data=data,
            panel_label=label,
            title=data["case"].display_name,
            grid_resolution=args.grid_resolution,
            seed_density=args.seed_density,
            streamline_density=args.streamline_density,
            rasterize_field=args.rasterize_field,
            speed_norm_max=global_speed_max,
        )

    colorbar = fig.colorbar(
        field_artist,
        ax=axes,
        orientation="horizontal",
        fraction=0.055,
        pad=0.09,
        aspect=42,
    )
    colorbar.set_label(
        r"Normalized speed, $|\mathbf{u}|/U_0$",
        fontsize=15,
        fontweight="bold",
        labelpad=6,
    )
    colorbar.ax.tick_params(
        labelsize=12,
        width=1.3,
        length=5,
    )
    for label in colorbar.ax.get_xticklabels():
        label.set_fontweight("bold")
    colorbar.outline.set_linewidth(1.2)

    fig.suptitle(
        "Velocity fields and streamlines in the converged pore-scale simulations",
        fontsize=18,
        fontweight="bold",
        y=0.97,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        args.output,
        format="svg",
        dpi=args.dpi,
        facecolor="white",
        metadata={
            "Title": "Velocity fields and streamlines in digital rocks",
            "Creator": "Python/Matplotlib",
        },
    )

    plt.close(fig)

    print(f"Figura salva em: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
