#!/usr/bin/env python3
"""
FIGURE TYPE: verification and validation curves
DESCRIPTION: Figure 5 (verification and validation; formerly drawn by make_figures_npj.py fig3) redrawn in the figure
  standard of figstd.py (2026-11-04): one column (3x1; the author compared it with a 1x3 row and kept 3x1), no panel titles (content in the caption), closed boxes, ticks
  inward, no text over data. Same data:
  (a) plane Poiseuille flow under P1 (b9): |k/k_exact - 1| against the mesh level, Taylor-Hood and Mini under the
      pressure-driven definition and Taylor-Hood under the velocity-driven one, four channel widths;
  (b) periodic square array of cylinders: Taylor-Hood, OpenFOAM and the dilute series (b2, b5);
  (c) OpenFOAM mesh convergence relative to Taylor-Hood (b5).
INPUT FILES: RES/b9_poiseuille_p1.json, RES/b2_periodic_cell.json, RES/b5_openfoam_cell.json (nothing typed by hand)
OUTPUT FILES: <out>/figure5.{svg,pdf,png}
ORIGIN: CFD digital-rock manuscript, npj Computational Materials draft rev3, 2026-11-04
"""
import argparse, json
import numpy as np
from pathlib import Path
from matplotlib.lines import Line2D
import figstd as S
from figstd import plt, mpl, GREY, ref, letter

RES = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1/package/drc_audit_pack/results")


def figure5(out):
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(S.WIDTH, 8.15), gridspec_kw=dict(hspace=0.40, height_ratios=[1.4, 1, 1]))

    B9 = json.load(open(RES / "b9_poiseuille_p1.json"))["rows"]
    widths = sorted({r["h"] for r in B9}); shades = dict(zip(widths, ["#000000", "#555555", "#8c8c8c", "#bdbdbd"]))
    floor = 1e-12
    for h in widths:
        rr = sorted((r for r in B9 if r["h"] == h), key=lambda r: r["level"]); lv = [r["level"] for r in rr]
        col = shades[h]
        a.plot(lv, [max(abs(r["P1_p2"]["rel_face"]), floor) for r in rr], "-", color=col, marker="o", mfc=col, mec="black", mew=0.5, ms=5.5)
        a.plot(lv, [abs(r["P1_mini"]["rel_face"]) for r in rr], "-", color=col, marker="s", mfc="white", mec=col, ms=5)
        a.plot(lv, [abs(r["velocity_p2"]["rel_vs_poiseuille"]) for r in rr], "--", color=col, marker="^", mfc="white", mec=col, ms=5)
    hnd = [Line2D([], [], color="k", marker="o", mec="black", mew=0.5, ms=5.5, label="Taylor--Hood, pressure"),
           Line2D([], [], color="k", marker="s", mfc="white", ms=5, label="Mini, pressure"),
           Line2D([], [], color="k", ls="--", marker="^", mfc="white", ms=5, label="Taylor--Hood, velocity")]
    a.set_yscale("log"); a.set_ylim(1e-12, 3); a.set_xticks([1, 2, 3, 4, 5]); a.set_xlim(0.7, 5.3)
    a.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=20, subs=(1.0,)))
    a.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: r"$10^{%d}$" % round(np.log10(v)) if round(np.log10(v)) % 3 == 0 else ""))
    a.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    a.set_xlabel("Mesh level"); a.set_ylabel(r"$|k/k_{\mathrm{exact}}-1|$")
    a.legend(handles=hnd, frameon=False, loc="center right", handlelength=1.8, borderaxespad=0.2, labelspacing=0.25, handletextpad=0.5, bbox_to_anchor=(1.0, 0.585))
    letter(a, "a")

    B2 = json.load(open(RES / "b2_periodic_cell.json"))["cases"]
    B5 = json.load(open(RES / "b5_openfoam_cell.json"))["cases"]
    b.plot([x["c"] for x in B2 if x["c"] <= 0.15], [x["k_hasimoto"] for x in B2 if x["c"] <= 0.15], "-", color=GREY, lw=1.2, label="dilute series")
    b.plot([x["c"] for x in B2 if x["c"] >= 0.15], [x["k_hasimoto"] for x in B2 if x["c"] >= 0.15], "--", color=GREY, lw=1.0)
    b.plot([x["c"] for x in B2], [x["k_num"] for x in B2], "o", color="#0072B2", mfc="white", mew=1.2, ms=6.5, label="FE, Taylor--Hood")
    b.plot([x["c"] for x in B5], [x["k_of"] for x in B5], "+", color="#D55E00", ms=9, mew=1.5, label="FV, OpenFOAM")
    b.set_yscale("log"); b.set_xlabel(r"Solid fraction $c$"); b.set_ylabel(r"$k/L^{2}$")
    b.legend(frameon=False, loc="lower left", handlelength=1.6, borderaxespad=0.3)
    letter(b, "b")

    cols = {0.05: "#000000", 0.1: "#0072B2", 0.2: "#009E73", 0.4: "#CC79A7"}
    for x in B5:
        mm = sorted(x["mesh"], key=lambda m: -m["h"])
        c.plot([1 / m["h"] for m in mm], [100 * (m["k"] / x["k_p2p1"] - 1) for m in mm], "-", marker="o", ms=5.5,
               color=cols[x["c"]], mec="black", mew=0.5, label=r"$c=%g$" % x["c"])
    ref(c.axhline(0, color=GREY, lw=0.6, ls=":"))
    c.set_xscale("log"); c.set_xticks([25, 50, 100, 200]); c.set_xticklabels(["25", "50", "100", "200"])
    c.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    c.set_xlabel(r"Cells per cell side, $1/h$"); c.set_ylabel(r"$k_{\mathrm{FV}}/k_{\mathrm{FE}}-1$ [\%]")
    c.legend(frameon=False, loc="upper right", handlelength=1.6, borderaxespad=0.3)
    letter(c, "c")
    fig.align_ylabels([a, b, c])
    S.save(fig, "figure5", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    figure5(a.out)
