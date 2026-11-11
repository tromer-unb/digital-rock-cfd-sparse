#!/usr/bin/env python3
"""
TEMPLATE-ID: FS-house-style-01 (ADAPT: house_style and qcdata of QC-20261019-v1 used as-is; layout new)
FIGURE TYPE: binary image mosaic (Fig. 1), mesh-convergence curves (Fig. 2)
DESCRIPTION: Splits the former Figure 1 (QC-F1 of make_figures_qc.py) in two, at the author's request (2026-11-04).
  Figure 1  former panel (a) alone, drawing unchanged, panel letter removed, larger type: the three samples R1-R3 with the
            percolating domain of each direction.
  Figure 2  former panels (b)-(d) stacked in one column (3x1), relabelled (a)-(c): k_app against the number of
            triangles for R3-H, R2-V and R3-V; shared x axis, x label on the bottom panel only; larger type.
INPUT FILES: same as QC-F1 (convergencia/convergence_mesh_0*/global_summary.csv, full_percolating_domain.npy,
             structures/*.png), read through qcdata; nothing typed by hand
OUTPUT FILES: <out>/figure1.{svg,pdf,png}, <out>/figure2.{svg,pdf,png}
STYLE: House Style Contract v1 (Times through LaTeX/mathptmx), Okabe-Ito colours as in QC-F1; closed axes boxes
ORIGIN: CFD digital-rock manuscript, npj Computational Materials draft rev3, 2026-11-04
"""
import argparse, sys
from pathlib import Path
import numpy as np

QCS = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1/scripts")
sys.path.insert(0, str(QCS / "templates")); sys.path.insert(0, str(QCS))
import house_style as hs                                        # noqa: E402
from house_style import plt, panel_label, save                  # noqa: E402
import matplotlib as mpl                                        # noqa: E402
from matplotlib.colors import ListedColormap                    # noqa: E402
import make_figures_qc as Q                                     # noqa: E402  (helpers and colours of QC-F1)
import qcdata as D                                              # noqa: E402

COL, MK, MKW, tex, lab = Q.COL, Q.MK, Q.MKW, Q.tex, Q.lab


def figure1(out):
    """Former panel (a) of QC-F1, same drawing code, no panel letter; type enlarged to the Fig. 2 standard (89 mm)."""
    mpl.rcParams.update({"xtick.top": False, "ytick.right": False})
    fig = plt.figure(figsize=(3.5, 2.75))
    sub = fig.add_gridspec(2, 3, left=0.08, right=0.99, bottom=0.08, top=0.92, wspace=0.07, hspace=0.20)
    cmap = ListedColormap([Q.SOLID, "white", Q.SKY, Q.BAND])
    for j, smp in enumerate(["R1", "R2", "R3"]):
        pore = D.structure(smp)
        for i, (d, nm) in enumerate((("lr", "horizontal"), ("tb", "vertical"))):
            ax = fig.add_subplot(sub[i, j])
            cls = np.where(pore, 1, 0).astype(np.uint8)
            dom = D.domain(smp, d)
            band = 328
            if dom is not None:
                if d == "lr": cls[:, band:band + dom.shape[1]][dom] = 2
                else: cls[band:band + dom.shape[0], :][dom] = 2
            sm = D.block_mode(cls, 8, 3)
            nb = band // 8
            if d == "lr": sm[:, :nb] = 3; sm[:, -nb:] = 3
            else: sm[:nb, :] = 3; sm[-nb:, :] = 3
            ax.imshow(sm, cmap=cmap, vmin=0, vmax=3, interpolation="nearest", aspect="equal"); Q.imgax(ax)
            if i == 0: ax.set_title(smp, fontsize=11, pad=3)
            chi = None
            try:
                chi = D.summary({"R2": "R2-V", "R3": "R3-H" if d == "lr" else "R3-V"}[smp], 4)["selected_fraction_of_cropped_pores"] if (smp != "R2" or d == "tb") and smp != "R1" else None
            except KeyError:
                chi = None
            ax.set_xlabel(tex(r"$\chi_c=%.1f\,\%%$" % (100 * chi)) if chi else "non-percolating", fontsize=9.5, labelpad=2)
            if j == 0: ax.set_ylabel(nm, fontsize=10.5, labelpad=3)
    save(fig, "figure1", out)


def figure2(out):
    """Former panels (b)-(d) of QC-F1 in one column; larger type for a single-column (89 mm) figure."""
    FS = dict(label=11, tick=10, title=11, note=10, letter=12)
    mpl.rcParams.update({"xtick.top": False, "ytick.right": False, "axes.labelsize": FS["label"],
                         "xtick.labelsize": FS["tick"], "ytick.labelsize": FS["tick"]})
    fig, axes = plt.subplots(3, 1, figsize=(3.5, 7.4), sharex=True, gridspec_kw=dict(hspace=0.28))
    for ax, case, letter in zip(axes, ("R3-H", "R2-V", "R3-V"), "abc"):
        n = D.series(case, "triangle_count"); k = D.series(case, "apparent_permeability")
        ax.axhspan(k[-1] * 0.99, k[-1] * 1.01, color="0.88", lw=0, zorder=0)
        ax.axhline(k[-1], color="black", lw=0.6, ls="--", zorder=1)
        ax.plot(n, k, "-", color=COL[case], lw=1.2, zorder=2)
        ax.plot(n, k, MK[case], color=COL[case], zorder=3, **dict(MKW, ms=6))
        Q.ntri_axis(ax)
        ax.set_ylabel(lab(r"$k_{\mathrm{app}}$", r"px$^2$"))
        dk = 100 * abs(k[-1] - k[-2]) / k[-1]
        up = case != "R3-H"
        ax.text(0.97, 0.94 if up else 0.06, tex(r"$\Delta k_{\mathrm{last}}=%.2f\,\%%$" % dk), transform=ax.transAxes,
                ha="right", va="top" if up else "bottom", fontsize=FS["note"])
        ax.set_title(case, fontsize=FS["title"], color=COL[case], pad=3)
        ax.text(-0.22, 1.02, (r"\textbf{(%s)}" % letter) if hs.USE_TEX else f"({letter})", transform=ax.transAxes,
                fontsize=FS["letter"], ha="left", va="bottom", fontweight=None if hs.USE_TEX else "bold")
    axes[-1].set_xlabel(lab(r"Number of triangles $N_{\mathrm{tri}}$", "-"))
    fig.align_ylabels(axes)
    save(fig, "figure2", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    figure1(a.out); figure2(a.out)
