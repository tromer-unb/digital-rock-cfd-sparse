#!/usr/bin/env python3
"""
FIGURE TYPE: solver cost and accuracy curves
DESCRIPTION: Former Figure 6 (accuracy and cost of the reference and screening solvers; formerly drawn by
  make_figures_npj.py fig4) split in two at the author's request (2026-11-04), in the figure standard of figstd.py:
  Figure 6 = former a-c (3x1), Figure 7 = former d-e (2x1, relabelled a-b); no panel titles (content in the
  caption), closed boxes, ticks inward, no text over data. The axis
  quantity N ("Unknowns N") is written as the number of degrees of freedom N_dof, the symbol of Table 1. Same data:
  (a) N_dof against N_tri, Taylor-Hood and Mini; (b) solve time; (c) peak memory; (d) MINRES vs LU and Mini vs
  Taylor-Hood differences in k_D; (e) wall time of contour FE and pixel FV on the 11 matched sub-windows of R2-V.
INPUT FILES: QC/fix-round/FR-E1_convergence_tol0.25/*_{mini,p2}.json, RES/c4_minres_p1.json,
             RES/b7_subwindow_pair_R2-V_tile512.json (nothing typed by hand)
OUTPUT FILES: <out>/figure6.{svg,pdf,png}, <out>/figure7.{svg,pdf,png}
ORIGIN: CFD digital-rock manuscript, npj Computational Materials draft rev3, 2026-11-04
"""
import argparse, json
from pathlib import Path
import numpy as np
from matplotlib.lines import Line2D
import figstd as S
from figstd import plt, mpl, COL, MK, MKW, GREY, CASES, ref, letter

QC = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1")
RES = QC / "package/drc_audit_pack/results"
FRE1 = QC / "fix-round/FR-E1_convergence_tol0.25"
FILE = {"R3-H": "R3H", "R2-V": "R2V", "R3-V": "R3V"}
OPEN = dict(mfc="white", mew=1.0, ms=6)
NDOF = r"Degrees of freedom $N_{\mathrm{dof}}$"


def lab(v):
    e_ = int(np.floor(np.log10(v))); m_ = v / 10 ** e_
    return rf"$10^{{{e_}}}$" if abs(m_ - 1) < 1e-9 else rf"${m_:g}\times10^{{{e_}}}$"


def ticks(ax, axis, vals):
    getattr(ax, f"set_{axis}ticks")(vals); getattr(ax, f"set_{axis}ticklabels")([lab(v) for v in vals])
    getattr(ax, f"{axis}axis").set_minor_formatter(mpl.ticker.NullFormatter())


def figure6(out):
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(S.WIDTH, 7.8), gridspec_kw=dict(hspace=0.42))
    fig7, (d, e) = plt.subplots(2, 1, figsize=(S.WIDTH, 5.4), gridspec_kw=dict(hspace=0.42))
    C4 = json.load(open(RES / "c4_minres_p1.json"))["rows"]
    for cs in CASES:
        mini = [json.load(open(FRE1 / f"{FILE[cs]}_L{l}_pad0.1_tol0.25_band0.08_mini.json")) for l in range(1, 6)]
        p2 = [json.load(open(FRE1 / f"{FILE[cs]}_L{l}_pad0.1_tol0.25_band0.08_p2.json")) for l in range(1, 6)]
        nt = np.array([m["n_tri"] for m in mini]); nm = np.array([m["n_dof"] for m in mini]); n2 = np.array([m["n_dof"] for m in p2])
        a.plot(nt, nm, "-", color=COL[cs], marker=MK[cs], mec=COL[cs], **OPEN)
        a.plot(nt, n2, "-", color=COL[cs], marker=MK[cs], **MKW)
        lu = sorted([r for r in C4 if r["case"] == FILE[cs] and r["method"] == "lu"], key=lambda r: r["level"])
        mr = sorted([r for r in C4 if r["case"] == FILE[cs] and r["method"] == "minres"], key=lambda r: r["level"])
        nf = [r["n_free"] for r in lu]
        b.plot(n2, [p["t_lu_s"] for p in p2], ":", color=COL[cs], marker=MK[cs], mfc=COL[cs], mec="black", mew=0.5, ms=4.5)
        b.plot(nf, [r["t_solve_s"] for r in lu], "-", color=COL[cs], marker=MK[cs], **MKW)
        b.plot([r["n_free"] for r in mr], [r["t_solve_s"] for r in mr], "--", color=COL[cs], marker=MK[cs], mec=COL[cs], **OPEN)
        c.plot(nf, [r["peak_rss_MB_after_solve"] / 1024 for r in lu], "-", color=COL[cs], marker=MK[cs], **MKW)
        c.plot([r["n_free"] for r in mr], [r["peak_rss_MB_after_solve"] / 1024 for r in mr], "--", color=COL[cs], marker=MK[cs], mec=COL[cs], **OPEN)
        d.plot([r["n_free"] for r in mr], [abs(r["rel_k_vs_LU_reference"]) for r in mr], "--", color=COL[cs], marker=MK[cs], mec=COL[cs], **OPEN)
        d.plot(nm, [abs(m["kD_iface_flux"] / p["kD_iface_flux"] - 1) for m, p in zip(mini, p2)], "-", color=COL[cs], marker=MK[cs], **MKW)

    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"Number of triangles $N_{\mathrm{tri}}$"); a.set_ylabel(NDOF)
    a.legend(handles=[Line2D([], [], color="k", marker="o", mec="black", mew=0.5, ms=6, label="Taylor--Hood"),
                      Line2D([], [], color="k", marker="o", **OPEN, label="Mini")] +
             [Line2D([], [], color=COL[cs], marker=MK[cs], ls="none", **MKW, label=cs) for cs in CASES],
             frameon=False, loc="lower right", handlelength=1.4, borderaxespad=0.2, labelspacing=0.2, handletextpad=0.4)
    ticks(a, "x", [6e4, 2e5, 5e5]); ticks(a, "y", [1e4, 1e5, 1e6]); a.set_ylim(8e3, 4e6)
    letter(a, "a", x=-0.26)

    b.set_xscale("log"); b.set_yscale("log"); b.set_xlabel(NDOF); b.set_ylabel("Solve time [s]")
    b.legend(handles=[Line2D([], [], color="k", ls=":", marker="o", ms=4.5, label="Taylor--Hood, LU"),
                      Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=6, label="Mini, LU"),
                      Line2D([], [], color="k", ls="--", marker="o", **OPEN, label="Mini, MINRES")],
             frameon=False, loc="upper left", handlelength=1.8, borderaxespad=0.2, labelspacing=0.2, handletextpad=0.4)
    ticks(b, "x", [2e5, 5e5, 2e6]); ticks(b, "y", [1, 10, 100, 1000]); b.set_ylim(0.5, 2e4)
    letter(b, "b", x=-0.26)

    c.set_xscale("log"); c.set_yscale("log"); c.set_xlabel(NDOF); c.set_ylabel("Peak memory [GB]")
    c.legend(handles=[Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=6, label="Mini, LU"),
                      Line2D([], [], color="k", ls="--", marker="o", **OPEN, label="Mini, MINRES")],
             frameon=False, loc="upper left", handlelength=1.8, borderaxespad=0.2, labelspacing=0.2, handletextpad=0.4)
    ticks(c, "x", [2e5, 5e5, 2e6]); c.set_ylim(0.25, 30)
    c.set_yticks([0.5, 1, 2, 4, 8]); c.set_yticklabels(["0.5", "1", "2", "4", "8"]); c.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    letter(c, "c", x=-0.26)
    fig.align_ylabels([a, b, c])
    S.save(fig, "figure6", out)

    d.set_xscale("log"); d.set_yscale("log"); d.set_xlabel(NDOF + " (Mini)"); d.set_ylabel(r"Relative difference in $k_D$")
    ref(d.axhline(1e-3, color=GREY, lw=0.6, ls=":"))
    d.legend(handles=[Line2D([], [], color="k", ls="--", marker="o", **OPEN, label="MINRES vs LU"),
                      Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=6, label="Mini vs Taylor--Hood")],
             frameon=False, loc="upper right", handlelength=1.8, borderaxespad=0.2, labelspacing=0.2, handletextpad=0.4)
    ticks(d, "x", [2e5, 5e5, 2e6]); d.set_ylim(1e-5, 30); ticks(d, "y", [1e-5, 1e-3, 1e-1, 10])
    letter(d, "a", x=-0.26)

    b7 = json.load(open(RES / "b7_subwindow_pair_R2-V_tile512.json"))
    sound = sorted([t for t in b7["tiles"] if t.get("status") == "ok" and t.get("sound")], key=lambda t: t["fem"]["5"]["kD_face"])
    fe = np.array([t["fem"]["5"]["t_mesh_s"] + t["fem"]["5"]["t_lu_s"] for t in sound])
    fv = np.array([t["fv"]["t_mesh_s"] + t["fv"]["t_solve_s"] for t in sound])
    x = np.arange(1, len(sound) + 1)
    e.plot(x, fe, "s", color=COL["R2-V"], mec="black", mew=0.5, ms=6.5, label="contour FE, 1 process")
    e.plot(x, fv, "D", color="#56B4E9", mec="black", mew=0.5, ms=6, label="pixel FV, 6--10 processes")
    e.set_yscale("log"); e.set_xticks(x); e.set_xlim(0.4, len(sound) + 0.6)
    e.set_xlabel(r"Sub-window, in order of increasing $k_D$"); e.set_ylabel("Wall time [s]")
    e.set_ylim(0.5, 1e4)
    e.legend(frameon=False, loc="center left", handlelength=1.2, borderaxespad=0.2, labelspacing=0.2, bbox_to_anchor=(0.0, 0.36))
    letter(e, "b", x=-0.26)
    print("FE median %.2f s, FV median %.0f s" % (np.median(fe), np.median(fv)))
    fig7.align_ylabels([d, e])
    S.save(fig7, "figure7", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a_ = ap.parse_args()
    figure6(a_.out)
