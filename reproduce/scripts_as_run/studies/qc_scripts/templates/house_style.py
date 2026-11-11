"""
TEMPLATE-ID:   FS-house-style-01
FIGURE TYPE:   style module (import from every other template)
DESCRIPTION:   LCCMat house figure style: LaTeX text in Times (mathptmx), 10 pt, 0.5 pt axes and ticks,
               frameless legends, units in square brackets, tab: colours plus the lab's custom colours,
               SVG + PDF + high-dpi PNG output.
INPUT FILES:   none
OUTPUT FILES:  none (provides rcParams, colours, panel_label(), unit(), save())
PARAMETERS:    none
DEPENDENCIES:  matplotlib >= 3.7; LaTeX with psnfss/mathptmx (optional, see USE_TEX fallback)
STYLE:         House Style Contract v1 (derived from plot_psiavg.py, phase_diagram_plot.py,
               plot_alpha-plus-spectrum.py, exciton_lifetimes_2D_systems_v2.py in FIGURES-SCRIPTS/20128410)
PHYSICS NOTES: none. Templates are style references, not physics authorities: audit constants and units
               of every script that is adapted from them.
ORIGIN:        ga-reaxff paper_style.py (project projeto03-26092026), generalized during the quality-check
               round QC-20260930-v1; reviewer model claude-sonnet-5-5, high effort
REVISION LOG:  2026-09-30  first deposit; adds the mathtext fallback for machines without LaTeX
               2026-10-17  QC-20261017-v1 (claude-sonnet-5-5, high): _latex_ok() also requires dvipng; without it
                           usetex raised "No pages of output" on a machine that has latex and mathptmx.sty.
                           Figures made with the fallback are "style-approximate".
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")   # batch/headless: figures are written to files, never shown
import matplotlib.pyplot as plt  # noqa: E402

SINGLE = 3.25   # in, one column (ACS/JCTC). Check the target journal and change if needed.
DOUBLE = 7.0    # in, two columns


def _latex_ok() -> bool:
    """Can matplotlib actually render text through LaTeX here?

    QC-20261017-v1 required dvipng, after usetex raised "No pages of output" on a machine that had latex and
    mathptmx. That requirement is not true of every matplotlib: 3.11 renders usetex PNG and PDF through its own
    dvi reader and needs no dvipng. Rather than guess from the presence of binaries, render one throwaway label
    and see whether it works, which is the only question that matters.
    """
    if not shutil.which("latex"):
        return False
    try:
        out = subprocess.run(["kpsewhich", "mathptmx.sty"], capture_output=True, text=True, timeout=10)
        if not out.stdout.strip():
            return False
    except Exception:
        return False
    try:
        import matplotlib.pyplot as _plt
        with mpl.rc_context({"text.usetex": True, "font.family": "serif", "font.serif": ["Times"],
                             "text.latex.preamble": r"\usepackage{mathptmx}"}):
            _f = _plt.figure(figsize=(1, 1))
            _f.text(0.1, 0.5, r"$3\times10^{4}\,\%$")
            _f.canvas.draw()
            _plt.close(_f)
        return True
    except Exception:
        return False

USE_TEX = _latex_ok()

_base = {
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}
if USE_TEX:
    _base.update({"text.usetex": True, "font.serif": ["Times"],
                  "text.latex.preamble": r"\usepackage{mathptmx}"})
else:   # deviation from the house style: report it ("style-approximate") when using this branch
    _base.update({"text.usetex": False, "mathtext.fontset": "stix",
                  "font.serif": ["Times New Roman", "Times", "STIXGeneral"]})
mpl.rcParams.update(_base)

# Colours seen in the lab scripts
C = {
    "purple": "tab:purple", "orange": "tab:orange", "blue": "tab:blue", "red": "tab:red",
    "green": "tab:green", "chocolate": "chocolate", "steelblue": "steelblue",
    "darkred": "#8b0000", "grey": "#c5d0d5", "deep_red": "#B2182B", "strong_orange": "#EF8A00",
    "sand": "#E6BE6A", "pale_yellow": "#F2E8B6", "light_green": "#A6C96A",
    "mid_green": "#5AA14F", "dark_green": "#0B6B2E", "pale_blue": "#b5d8ef",
}
# Default cycle for categories (keep one colour per entity across all figures of a paper)
CYCLE = [C["darkred"], C["steelblue"], C["sand"], C["grey"], C["mid_green"], C["purple"]]
MARKERS = ["o", "s", "^", "D", "v", "P"]   # second cue for grayscale readers


def unit(quantity: str, u: str) -> str:
    """Axis label in the house notation: quantity followed by the unit in square brackets."""
    return rf"{quantity} [{u}]"


def panel_label(ax, letter: str, x: float = -0.18, y: float = 1.02):
    ax.text(x, y, (r"\textbf{(%s)}" % letter) if USE_TEX else f"({letter})", transform=ax.transAxes,
            fontsize=10, ha="left", va="bottom", fontweight=None if USE_TEX else "bold")


def demo_banner(fig, text: str = "DEMO DATA, NOT RESULTS"):
    fig.text(0.5, 0.5, text, ha="center", va="center", rotation=25, color="0.6", fontsize=14,
             alpha=0.5, zorder=100)


def save(fig, name: str, outdir: str | Path = "figures", dpi_png: int = 500):
    base = Path(outdir) / name            # `name` may contain sub-directories
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_name(base.name + ".svg"), dpi=300)
    fig.savefig(base.with_name(base.name + ".pdf"))
    fig.savefig(base.with_name(base.name + ".png"), dpi=dpi_png)
    plt.close(fig)
    print(f"saved {base}.[svg|pdf|png]  (usetex={USE_TEX})")
