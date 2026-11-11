#!/usr/bin/env python3
"""
FIGURE TYPE: sensitivity curves (Fig. 2), verification and validation (Fig. 3), solver accuracy and cost (Fig. 4),
             mesh-size field check (Supplementary Fig. 2)
DESCRIPTION: Display items of the npj Computational Materials version that are built from protocol-P1 result files.
  Figure 2  Definitional sensitivities of the permeability (was Figure 7 of the Nature Communications draft).
  Figure 3  Verification and validation: (a) plane Poiseuille flow under P1 (b9), (b) periodic cylinder array,
            Taylor-Hood, OpenFOAM and the dilute series, (c) OpenFOAM mesh convergence relative to Taylor-Hood.
  Figure 4  Screening and reference solvers under P1: (a) unknowns vs triangles, (b) wall time, (c) peak memory,
            (d) accuracy of the screening solver, (e) contour FE vs pixel FV wall time on matched sub-windows.
  Figure S2 Gmsh Threshold size field with Sigmoid = 1 against the sigmoid law and a linear transition (c2).
INPUT FILES (nothing typed by hand):
  QC/evidence/indep_sens/tol{0.5,1.0,1.5}_{R3H,R2V,R3V}_{velocity,pressure}.json, QC/fix-round/fixround_summary.json,
  QC/fix-round/FR-E1_convergence_tol0.25/*_{mini,p2}.json,
  RES/{b2_periodic_cell,b5_openfoam_cell,b7_subwindow_pair_R2-V_tile512,b8_pixel_convergence_R2-V_t1_5,
       b9_poiseuille_p1,c2_gmsh_size_field,c4_minres_p1}.json
OUTPUT FILES: figures/figure{2,3,4}.{pdf,png,svg}, figures/figureS2.{pdf,png,svg}; numbers printed for the text
STYLE: npj figure guidelines: sans-serif lettering (Liberation Sans, metric-compatible with Arial), lower-case bold panel
       letters, Okabe-Ito colours, 180 mm width, no red-green contrast. R3-H blue circle, R2-V orange square,
       R3-V green triangle, as in the other figures.
ORIGIN: CFD digital-rock manuscript, npj Computational Materials revision, 2026-11-03
"""
import json
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

QC = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1")
RES = QC / "package/drc_audit_pack/results"
FRE1 = QC / "fix-round/FR-E1_convergence_tol0.25"
OUT = Path(__file__).resolve().parents[1] / "figures"

mpl.rcParams.update({
    "text.usetex": False, "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular", "axes.unicode_minus": True,
    "figure.dpi": 300, "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "lines.linewidth": 1.0, "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
W = 7.09
CASES = ["R3-H", "R2-V", "R3-V"]
FILE = {"R3-H": "R3H", "R2-V": "R2V", "R3-V": "R3V"}
COL = {"R3-H": "#0072B2", "R2-V": "#E69F00", "R3-V": "#009E73"}
MK = {"R3-H": "o", "R2-V": "s", "R3-V": "^"}
MKW = dict(mec="black", mew=0.5, ms=4.5)
GREY = "0.45"
log = {}


def PL(ax, s, x=-0.2, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom")


def save(fig, name):
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=500 if ext == "png" else 300, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def nolog_minor(ax, axis="y"):
    getattr(ax, f"{axis}axis").set_minor_formatter(mpl.ticker.NullFormatter())


# ======================================================================= Figure 2 (sensitivities)
def fig2():
    FR = json.load(open(QC / "fix-round/fixround_summary.json"))
    fig = plt.figure(figsize=(W, 4.9))
    gs = fig.add_gridspec(2, 6, hspace=0.62, wspace=1.15)
    a = fig.add_subplot(gs[0, 0:2]); b = fig.add_subplot(gs[0, 2:4]); c = fig.add_subplot(gs[0, 4:6])
    d = fig.add_subplot(gs[1, 0:3]); e = fig.add_subplot(gs[1, 3:6])
    tols = [0.5, 1.0, 1.5]
    log["2a"] = {}
    for cs in CASES:
        r = []
        for t in tols:
            v = json.load(open(QC / f"evidence/indep_sens/tol{t}_{FILE[cs]}_velocity.json"))
            p = json.load(open(QC / f"evidence/indep_sens/tol{t}_{FILE[cs]}_pressure.json"))
            r.append(p["k"]["nodal"]["k_author_def"] / v["k"]["nodal"]["k_author_def"])
        log["2a"][cs] = dict(zip(tols, r))
        a.plot(tols, r, "-", color=COL[cs], marker=MK[cs], label=cs, **MKW)
    a.axhline(1, color=GREY, lw=0.6, ls=":")
    a.set_xlabel("Contour tolerance (px)"); a.set_ylabel(r"$\mathit{k}_{\rm pressure}/\mathit{k}_{\rm velocity}$")
    a.set_yscale("log"); a.set_yticks([1, 2, 5]); a.set_yticklabels(["1", "2", "5"]); a.set_ylim(0.8, 7); nolog_minor(a)
    a.set_xticks(tols); a.legend(frameon=False, loc="upper right", handlelength=1.2)
    a.set_title("Inlet condition"); PL(a, "a")
    for cs in CASES:
        dv = FR["reservoir"][cs]["dev_iface_pct"]; bb = sorted(float(x) for x in dv)
        b.plot([100 * x for x in bb], [dv[f"{x:g}"] for x in bb], "-", color=COL[cs], marker=MK[cs], **MKW)
    b.axhline(0, color=GREY, lw=0.6, ls=":"); b.axvline(10, color=GREY, lw=0.6, ls="--")
    b.set_xlabel(r"Reservoir width $\mathit{b}/\mathit{L}_{s}$ (%)"); b.set_ylabel("Deviation (%)")
    b.set_title("Reservoir width"); PL(b, "b")
    for cs in CASES:
        dv = FR["tolerance"][cs]["dev_pct"]; tt = sorted(float(x) for x in dv)
        c.plot(tt, [dv[f"{x:g}" if f"{x:g}" in dv else str(x)] for x in tt], "-", color=COL[cs], marker=MK[cs], **MKW)
    c.axhline(0, color=GREY, lw=0.6, ls=":"); c.axvline(1.5, color=GREY, lw=0.6, ls="--")
    c.set_xlabel("Contour tolerance (px)"); c.set_ylabel("Deviation (%)")
    c.set_title("Geometry simplification"); PL(c, "c")
    for cs in CASES:
        w = FR["window"][cs]; k8 = w["k_iface"][w["band_pct"].index(8.0)]
        d.plot(w["band_pct"], np.array(w["k_iface"]) / k8, "-", color=COL[cs], marker=MK[cs], label=cs, **MKW)
        for bp in w["not_percolating_band_pct"]:
            d.plot([bp], [0.62], marker="x", color=COL[cs], ms=6, mew=1.2, ls="none")
    d.axvline(8, color=GREY, lw=0.6, ls="--"); d.axhline(1, color=GREY, lw=0.6, ls=":")
    d.text(4.6, 0.62, "not percolating", fontsize=6.5, color=COL["R2-V"], ha="left", va="center")
    d.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
    d.set_xlabel("Boundary band (% of the image side)"); d.set_ylabel(r"$\mathit{k}_{D}/\mathit{k}_{D}$(8%)")
    d.set_ylim(0.42, 2.0); d.legend(frameon=False, loc="upper left", ncol=3, handlelength=1.2)
    d.set_title("Cropping window"); PL(d, "d", x=-0.12)
    b7 = json.load(open(RES / "b7_subwindow_pair_R2-V_tile512.json"))
    b8 = json.load(open(RES / "b8_pixel_convergence_R2-V_t1_5.json"))
    sound = [t for t in b7["tiles"] if t.get("status") == "ok" and t.get("sound")]
    kfem = np.array([t["fem"]["5"]["kD_face"] for t in sound]); rd = 100 * np.array([t["rel_diff"] for t in sound])
    e.plot(kfem, rd, "s", color=COL["R2-V"], label="pixel FV vs contour FE, 11 pairs", **MKW)
    t15 = [t for t in sound if t["tag"] == b8["tag"]][0]
    e.plot([t15["fem"]["5"]["kD_face"]], [100 * b8["richardson_vs_fem"]], marker="*", color="black", ms=9, ls="none",
           label=r"same pair, cell size $\rightarrow$ 0")
    e.axhline(0, color=GREY, lw=0.6, ls=":")
    e.set_xscale("log"); e.set_xlabel(r"$\mathit{k}_{D}$ of the sub-window, contour FE (px$^2$)")
    e.set_ylabel(r"$\mathit{k}_{\rm pixel}/\mathit{k}_{\rm contour}-1$ (%)")
    e.legend(frameon=False, loc="lower left", handlelength=1.2)
    e.set_title("Image representation (R2-V, 512 px)"); PL(e, "e", x=-0.14)
    log["2e"] = dict(n=len(sound), median_abs=float(np.median(np.abs(rd))), max_abs=float(np.max(np.abs(rd))))
    save(fig, "figure2")


# ======================================================================= Figure 3 (verification and validation)
def fig3():
    fig = plt.figure(figsize=(W, 2.5))
    gs = fig.add_gridspec(1, 3, wspace=0.55)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[0, 2])
    B9 = json.load(open(RES / "b9_poiseuille_p1.json"))["rows"]
    widths = sorted({r["h"] for r in B9}); shades = dict(zip(widths, ["#000000", "#555555", "#8c8c8c", "#bdbdbd"]))
    floor = 1e-12
    for h in widths:
        rr = sorted((r for r in B9 if r["h"] == h), key=lambda r: r["level"]); lv = [r["level"] for r in rr]
        col = shades[h]
        a.plot(lv, [max(abs(r["P1_p2"]["rel_face"]), floor) for r in rr], "-", color=col, marker="o", mfc=col, mec="black", mew=0.4, ms=4)
        a.plot(lv, [abs(r["P1_mini"]["rel_face"]) for r in rr], "-", color=col, marker="s", mfc="white", mec=col, ms=3.8)
        a.plot(lv, [abs(r["velocity_p2"]["rel_vs_poiseuille"]) for r in rr], "--", color=col, marker="^", mfc="white", mec=col, ms=3.8)
    from matplotlib.lines import Line2D
    hnd = [Line2D([], [], color="k", marker="o", mec="black", mew=0.4, ms=4, label="Taylor-Hood, pressure"),
           Line2D([], [], color="k", marker="s", mfc="white", ms=3.8, label="Mini, pressure"),
           Line2D([], [], color="k", ls="--", marker="^", mfc="white", ms=3.8, label="Taylor-Hood, velocity")]
    a.set_yscale("log"); a.set_ylim(3e-13, 1); a.set_xticks([1, 2, 3, 4, 5])
    a.set_xlabel("Mesh level"); a.set_ylabel(r"$|\mathit{k}/\mathit{k}_{\rm exact}-1|$")
    a.legend(handles=hnd, frameon=False, loc="center right", fontsize=6, handlelength=1.6, bbox_to_anchor=(1.02, 0.665))
    a.text(5.0, 4e-13, r"$\mathit{h}$ = 20, 40, 80, 160 px" + "\n(dark to light)", fontsize=5.8, ha="right", va="bottom", color="0.25")
    a.set_title("Plane Poiseuille flow"); PL(a, "a", x=-0.3)
    log["3a"] = dict(worst_p2=max(abs(r["P1_p2"]["rel_face"]) for r in B9),
                     mini_range=(min(r["P1_mini"]["rel_face"] for r in B9), max(r["P1_mini"]["rel_face"] for r in B9)),
                     vel_range=(min(r["velocity_p2"]["rel_vs_poiseuille"] for r in B9), max(r["velocity_p2"]["rel_vs_poiseuille"] for r in B9)))
    B2 = json.load(open(RES / "b2_periodic_cell.json"))["cases"]
    B5 = json.load(open(RES / "b5_openfoam_cell.json"))["cases"]
    b.plot([x["c"] for x in B2 if x["c"] <= 0.15], [x["k_hasimoto"] for x in B2 if x["c"] <= 0.15], "-", color=GREY, lw=1.0, label="dilute series")
    b.plot([x["c"] for x in B2 if x["c"] >= 0.15], [x["k_hasimoto"] for x in B2 if x["c"] >= 0.15], "--", color=GREY, lw=0.8)
    b.plot([x["c"] for x in B2], [x["k_num"] for x in B2], "o", color="#0072B2", mfc="white", mew=1.0, ms=5, label="FE, Taylor-Hood")
    b.plot([x["c"] for x in B5], [x["k_of"] for x in B5], "+", color="#D55E00", ms=7, mew=1.3, label="FV, OpenFOAM")
    b.set_yscale("log"); b.set_xlabel(r"Solid fraction $\mathit{c}$"); b.set_ylabel(r"$\mathit{k}/\mathit{L}^{2}$")
    b.legend(frameon=False, loc="lower left", fontsize=6.5, handlelength=1.4)
    b.set_title("Periodic cylinder array"); PL(b, "b", x=-0.3)
    cols = {0.05: "#000000", 0.1: "#0072B2", 0.2: "#009E73", 0.4: "#CC79A7"}
    for x in B5:
        mm = sorted(x["mesh"], key=lambda m: -m["h"])
        c.plot([1 / m["h"] for m in mm], [100 * (m["k"] / x["k_p2p1"] - 1) for m in mm], "-", marker="o", ms=4,
               color=cols[x["c"]], mec="black", mew=0.4, label=r"$\mathit{c}$ = " + f"{x['c']:g}")
    c.axhline(0, color=GREY, lw=0.6, ls=":")
    c.set_xscale("log"); c.set_xticks([25, 50, 100, 200]); c.set_xticklabels(["25", "50", "100", "200"]); nolog_minor(c, "x")
    c.set_xlabel(r"Cells per cell side, 1/$\mathit{h}$"); c.set_ylabel(r"$\mathit{k}_{\rm FV}/\mathit{k}_{\rm FE}-1$ (%)")
    c.legend(frameon=False, loc="upper right", fontsize=6.5, handlelength=1.4)
    c.set_title("Finite volume vs finite element"); PL(c, "c", x=-0.3)
    save(fig, "figure3")


# ======================================================================= Figure 4 (solvers under P1)
def fig4():
    C4 = json.load(open(RES / "c4_minres_p1.json"))["rows"]
    fig = plt.figure(figsize=(W, 4.7))
    gs = fig.add_gridspec(2, 6, hspace=0.62, wspace=1.3)
    a = fig.add_subplot(gs[0, 0:2]); b = fig.add_subplot(gs[0, 2:4]); c = fig.add_subplot(gs[0, 4:6])
    d = fig.add_subplot(gs[1, 0:3]); e = fig.add_subplot(gs[1, 3:6])
    fits = {}
    for cs in CASES:
        mini = [json.load(open(FRE1 / f"{FILE[cs]}_L{l}_pad0.1_tol0.25_band0.08_mini.json")) for l in range(1, 6)]
        p2 = [json.load(open(FRE1 / f"{FILE[cs]}_L{l}_pad0.1_tol0.25_band0.08_p2.json")) for l in range(1, 6)]
        nt = np.array([m["n_tri"] for m in mini]); nm = np.array([m["n_dof"] for m in mini]); n2 = np.array([m["n_dof"] for m in p2])
        a.plot(nt, nm, "-", color=COL[cs], marker=MK[cs], mfc="white", mec=COL[cs], ms=4)
        a.plot(nt, n2, "-", color=COL[cs], marker=MK[cs], label=cs, **MKW)
        fits[cs] = dict(mini=float(np.polyfit(np.log(nt), np.log(nm), 1)[0]), p2=float(np.polyfit(np.log(nt), np.log(n2), 1)[0]),
                        p2_over_mini=float(np.median(n2 / nm)))
        lu = sorted([r for r in C4 if r["case"] == FILE[cs] and r["method"] == "lu"], key=lambda r: r["level"])
        mr = sorted([r for r in C4 if r["case"] == FILE[cs] and r["method"] == "minres"], key=lambda r: r["level"])
        nf = [r["n_free"] for r in lu]
        b.plot(n2, [p["t_lu_s"] for p in p2], ":", color=COL[cs], marker=MK[cs], mfc=COL[cs], mec="black", mew=0.4, ms=3.6)
        b.plot(nf, [r["t_solve_s"] for r in lu], "-", color=COL[cs], marker=MK[cs], **MKW)
        b.plot([r["n_free"] for r in mr], [r["t_solve_s"] for r in mr], "--", color=COL[cs], marker=MK[cs], mfc="white", mec=COL[cs], ms=4)
        c.plot(nf, [r["peak_rss_MB_after_solve"] / 1024 for r in lu], "-", color=COL[cs], marker=MK[cs], **MKW)
        c.plot([r["n_free"] for r in mr], [r["peak_rss_MB_after_solve"] / 1024 for r in mr], "--", color=COL[cs], marker=MK[cs], mfc="white", mec=COL[cs], ms=4)
        d.plot([r["n_free"] for r in mr], [abs(r["rel_k_vs_LU_reference"]) for r in mr], "--", color=COL[cs], marker=MK[cs], mfc="white", mec=COL[cs], ms=4)
        d.plot(nm, [abs(m["kD_iface_flux"] / p["kD_iface_flux"] - 1) for m, p in zip(mini, p2)], "-", color=COL[cs], marker=MK[cs], **MKW)
    log["4a"] = fits
    from matplotlib.lines import Line2D
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"Triangles $\mathit{N}_{\rm tri}$"); a.set_ylabel(r"Unknowns $\mathit{N}$")
    a.legend(handles=[Line2D([], [], color="k", marker="o", mec="black", mew=0.5, ms=4, label="Taylor-Hood"),
                      Line2D([], [], color="k", marker="o", mfc="white", ms=4, label="Mini")] +
             [Line2D([], [], color=COL[cs], marker=MK[cs], ls="none", mec="black", mew=0.5, ms=4, label=cs) for cs in CASES],
             frameon=False, loc="upper left", fontsize=6, handlelength=1.4)
    a.set_title("Problem size"); PL(a, "a", x=-0.3)

    b.set_xscale("log"); b.set_yscale("log"); b.set_xlabel(r"Unknowns $\mathit{N}$"); b.set_ylabel("Solve time (s)")
    b.legend(handles=[Line2D([], [], color="k", ls=":", marker="o", ms=3.6, label="Taylor-Hood, LU"),
                      Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=4, label="Mini, LU"),
                      Line2D([], [], color="k", ls="--", marker="o", mfc="white", ms=4, label="Mini, MINRES")],
             frameon=False, loc="lower right", fontsize=6, handlelength=1.8, bbox_to_anchor=(1.02, 0.0))
    b.set_title("Solve time"); PL(b, "b", x=-0.3)
    c.set_xscale("log"); c.set_yscale("log"); c.set_xlabel(r"Unknowns $\mathit{N}$"); c.set_ylabel("Peak memory (GB)")
    c.legend(handles=[Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=4, label="Mini, LU"),
                      Line2D([], [], color="k", ls="--", marker="o", mfc="white", ms=4, label="Mini, MINRES")],
             frameon=False, loc="upper left", fontsize=6, handlelength=1.8)
    c.set_yticks([0.5, 1, 2, 4, 8]); c.set_yticklabels(["0.5", "1", "2", "4", "8"])
    c.set_title("Peak memory"); PL(c, "c", x=-0.3)
    d.set_xscale("log"); d.set_yscale("log"); d.set_xlabel(r"Unknowns $\mathit{N}$ (Mini)"); d.set_ylabel(r"Relative difference in $\mathit{k}_{D}$")
    d.axhline(1e-3, color=GREY, lw=0.6, ls=":")
    d.legend(handles=[Line2D([], [], color="k", ls="--", marker="o", mfc="white", ms=4, label="MINRES vs LU (algebraic)"),
                      Line2D([], [], color="k", ls="-", marker="o", mec="black", mew=0.5, ms=4, label="Mini vs Taylor-Hood (discretization)")],
             frameon=False, loc="upper right", fontsize=6, handlelength=1.8)
    d.set_ylim(1e-5, 3)
    d.set_title("Accuracy of the screening solver"); PL(d, "d", x=-0.13)
    b7 = json.load(open(RES / "b7_subwindow_pair_R2-V_tile512.json"))
    sound = sorted([t for t in b7["tiles"] if t.get("status") == "ok" and t.get("sound")], key=lambda t: t["fem"]["5"]["kD_face"])
    fe = np.array([t["fem"]["5"]["t_mesh_s"] + t["fem"]["5"]["t_lu_s"] for t in sound])
    fv = np.array([t["fv"]["t_mesh_s"] + t["fv"]["t_solve_s"] for t in sound]); npr = np.array([t["fv"]["nproc"] for t in sound])
    x = np.arange(len(sound))
    e.plot(x, fe, "s", color=COL["R2-V"], mec="black", mew=0.5, ms=4.5, label="contour FE, Taylor-Hood + LU, 1 process")
    e.plot(x, fv, "D", color="#56B4E9", mec="black", mew=0.5, ms=4, label="pixel FV, OpenFOAM SIMPLEC, 6-10 processes")
    e.set_yscale("log"); e.set_xticks(x); e.set_xticklabels([str(i + 1) for i in x])
    e.set_xlabel(r"Sub-window, in order of increasing $\mathit{k}_{D}$"); e.set_ylabel("Wall time (s)")
    e.set_ylim(0.5, 1e4); e.legend(frameon=False, loc="center left", fontsize=6, handlelength=1.2, bbox_to_anchor=(0.0, 0.40))
    e.set_title("Matched sub-windows of R2-V (512 px)"); PL(e, "e", x=-0.13)
    log["4e"] = dict(fe_median=float(np.median(fe)), fv_median=float(np.median(fv)), ratio_median=float(np.median(fv / fe)),
                     core_ratio_median=float(np.median(fv * npr / fe)), fe_range=(float(fe.min()), float(fe.max())), fv_range=(float(fv.min()), float(fv.max())))
    mr = [r for r in C4 if r["method"] == "minres"]; lu = [r for r in C4 if r["method"] == "lu"]
    pair = {(r["case"], r["level"]): r for r in lu}
    log["4bcd"] = dict(minres_iter=(min(r["iterations"] for r in mr), max(r["iterations"] for r in mr)),
                       k_err_max=max(abs(r["rel_k_vs_LU_reference"]) for r in mr),
                       true_res=(min(r["rel_residual_true"] for r in mr), max(r["rel_residual_true"] for r in mr)),
                       time_ratio=(min(r["t_solve_s"] / pair[(r["case"], r["level"])]["t_solve_s"] for r in mr),
                                   max(r["t_solve_s"] / pair[(r["case"], r["level"])]["t_solve_s"] for r in mr)),
                       mem_ratio_L45=[(r["case"], r["level"], pair[(r["case"], r["level"])]["peak_rss_MB_after_solve"] / r["peak_rss_MB_after_solve"]) for r in mr if r["level"] >= 4])
    def lab(v):
        e_ = int(np.floor(np.log10(v))); m_ = v / 10 ** e_
        return rf"$10^{{{e_}}}$" if abs(m_ - 1) < 1e-9 else rf"${m_:g}\times10^{{{e_}}}$"
    def ticks(ax_, axis, vals):
        getattr(ax_, f"set_{axis}ticks")(vals); getattr(ax_, f"set_{axis}ticklabels")([lab(v) for v in vals])
        getattr(ax_, f"{axis}axis").set_minor_formatter(mpl.ticker.NullFormatter())
    ticks(a, "x", [6e4, 1e5, 2e5, 5e5]); ticks(a, "y", [2e5, 5e5, 1e6, 2e6])
    for ax_ in (b, c, d):
        ticks(ax_, "x", [2e5, 5e5, 1e6, 2e6])
    ticks(b, "y", [1, 10, 100]); ticks(d, "y", [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    c.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    b.set_ylim(0.25, 600)
    save(fig, "figure4")


# ======================================================================= Supplementary Figure 2 (size field)
def figS2():
    C2 = json.load(open(RES / "c2_gmsh_size_field.json")); P = C2["params"]
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    dd = np.linspace(0, 600, 601)
    r = np.clip((dd - P["d_min"]) / (P["d_max"] - P["d_min"]), 0, 1); s = 1 / (1 + np.exp(-(12 * r - 6)))
    ax.plot(dd, P["h_min"] * (1 - s) + P["h_max"] * s, "-", color="#0072B2", label="sigmoid law")
    ax.plot(dd, P["h_min"] + (P["h_max"] - P["h_min"]) * r, "--", color=GREY, label="linear transition")
    ax.plot([q["d"] for q in C2["rows"]], [q["h_meas"] for q in C2["rows"]], "o", color="#E69F00", mec="black", mew=0.4, ms=3.6, label="measured mean edge")
    for v in (P["d_min"], P["d_max"]):
        ax.axvline(v, color=GREY, lw=0.5, ls=":")
    ax.set_xlabel("Distance to the nearest wall (px)"); ax.set_ylabel("Element size (px)")
    ax.legend(frameon=False, loc="upper left", fontsize=6.5)
    save(fig, "figureS1")


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["2", "3", "4", "S2"]
    for w in which:
        {"2": fig2, "3": fig3, "4": fig4, "S2": figS2}[w]()
    print(json.dumps(log, indent=1, default=float))
