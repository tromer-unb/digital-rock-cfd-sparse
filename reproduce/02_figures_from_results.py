#!/usr/bin/env python3
"""
Step 2 of the reproduction: regenerate the display items of the paper that are drawn from the stored run outputs and
fields, and check the plotted numbers against the archived figure source data.

  file figure2  = Figure 4   mesh convergence    <- data/results/fix_round/FR-E1_convergence_tol0.25/*_{p2,mini}.json
  file figure8  = Figure 8   flow fields         <- data/results/fix_round/FR-H_fields_P1/*_L5_fields.npz + percolating masks
  file figure9  = Figure 9   hydraulic partition <- idem + p1_descriptors.json
  file figure10 = Figure 10  heterogeneity       <- FR-H_fields_P1/*_tri.npz, p1_descriptors.json
  file figure11 = Figure 11  integrated response <- FR-H_fields_P1/p1_descriptors.json (+ connected fraction from the masks)
Each figure writes the numbers it plots to <figure>_source.json. At the end every *_source.json is compared with
data/results/figure_sources/ (acceptance: identical numbers, relative 1e-9), and the result is printed as PASS/FAIL.
Style: Times text via LaTeX, Okabe-Ito colours (R3-H blue circles, R2-V orange squares, R3-V green triangles).

Usage:  python reproduce/02_figures_from_results.py [outdir] [2 8 9 10 11]     (default outdir: runs/02_figures)
Needs: the PercoFlow environment and a LaTeX installation (text.usetex).
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap, Normalize
from matplotlib.tri import Triangulation

# Repository version of scripts_as_run/studies/make_figures_p1.py: identical plotting code; paths point to data/results and the
# percolating masks are rebuilt with PercoFlow (same rules as the audited pipeline; checked: identical pixel counts).
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
E1 = ROOT / "data/results/fix_round/FR-E1_convergence_tol0.25"
FH = ROOT / "data/results/fix_round/FR-H_fields_P1"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "runs/02_figures"
OUT.mkdir(parents=True, exist_ok=True)  # noqa

mpl.rcParams.update({
    "text.usetex": True, "text.latex.preamble": r"\usepackage{times}\usepackage{amsmath}",
    "font.family": "serif", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.2, "lines.markersize": 5,
    "savefig.dpi": 600, "pdf.fonttype": 42})
CASES = ["R3-H", "R2-V", "R3-V"]
COL = {"R3-H": "#0072B2", "R2-V": "#E69F00", "R3-V": "#009E73"}
MK = {"R3-H": "o", "R2-V": "s", "R3-V": "^"}
VERM = "#D55E00"
MM = 1 / 25.4


def plabel(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, r"\textbf{(%s)}" % s, transform=ax.transAxes, fontsize=8, va="bottom", ha="left")


def save(fig, name, src):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", pad_inches=0.02, dpi=200)
    json.dump(src, open(OUT / f"{name}_source.json", "w"), indent=1)
    plt.close(fig)


# ------------------------------------------------------------------ masks (geometry of the percolating domains)
def masks():
    """Cropped pore mask and percolating domain of each configuration, rebuilt with PercoFlow."""
    from percoflow import imaging, percolation
    from percoflow.config import DomainOptions, ImageOptions
    files = {"R3-H": ("patch_y7600_x53200_c0_mask.png", "lr"), "R2-V": ("patch_y7600_x19000_c0_mask.png", "tb"),
             "R3-V": ("patch_y7600_x53200_c0_mask.png", "tb")}
    res = {}
    for c, (img, d) in files.items():
        mask = imaging.clean_mask(imaging.load_mask(ROOT / "data/images" / img, ImageOptions()), ImageOptions())
        crop, _ = imaging.crop_window(mask, DomainOptions(direction=d, band_fraction=0.08))
        sel, _ = percolation.analyse(crop, d)
        res[c] = (crop.astype(bool), sel.astype(bool))
    return res


def fields(c):
    cc = c.replace("-", "")
    z = np.load(FH / f"{cc}_L5_fields.npz"); meta = json.load(open(FH / f"{cc}_L5_fields.json"))
    P, T = z["points"].copy(), z["triangles"]
    b, Ls, d = meta["pad_px"], meta["Ls"], meta["direction"]
    ax = 0 if d == "lr" else 1
    cen = P[:, T].mean(axis=1)
    ins = (cen[ax] > b - 0.5) & (cen[ax] < b + Ls - 0.5)
    P[ax] -= b                                         # cropped-frame coordinates
    sp = np.hypot(z["ux"], z["uy"]) / meta["U_darcy"]
    pst = (z["p"] - meta["p_out_bar"]) / meta["dP_iface_flux"]
    return P, T[:, ins], sp, pst, meta


def bbox(sel):
    ii, jj = np.where(sel); ny = sel.shape[0]
    return jj.min() - 0.5, jj.max() + 0.5, ny - 1 - ii.max() - 0.5, ny - 1 - ii.min() + 0.5


def show_domain(ax, pore, sel):
    rgb = np.zeros(sel.shape + (3,)); rgb[:] = (0.20, 0.20, 0.22)
    rgb[pore] = (1, 1, 1); rgb[sel] = (0.34, 0.71, 0.91)
    ny, nx = sel.shape
    ax.imshow(rgb, extent=(-0.5, nx - 0.5, -0.5, ny - 0.5), origin="upper", interpolation="nearest", rasterized=True)


def frame(ax, bb):
    ax.set_xlim(bb[0], bb[1]); ax.set_ylim(bb[2], bb[3]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)


# ------------------------------------------------------------------ figure 2
def fig2():
    fig, axs = plt.subplots(3, 1, figsize=(85 * MM, 158 * MM))
    src = {}
    for ax, c, lab in zip(axs, CASES, "abc"):
        cc = c.replace("-", "")
        rows = {el: [json.load(open(E1 / f"{cc}_L{L}_pad0.1_tol0.25_band0.08_{el}.json")) for L in range(1, 6)]
                for el in ("p2", "mini")}
        nt = [r["n_tri"] for r in rows["p2"]]
        kth = [r["kD_iface_flux"] for r in rows["p2"]]; kmi = [r["kD_iface_flux"] for r in rows["mini"]]
        kref = kth[-1]; dlast = 100 * (kth[-1] - kth[-2]) / kth[-2]
        ax.axhspan(0.99 * kref, 1.01 * kref, color="0.88", lw=0, zorder=0)
        ax.axhline(kref, color="k", ls=(0, (4, 2)), lw=0.7, zorder=1)
        ax.plot(nt, kth, "-", marker=MK[c], color=COL[c], mec="k", mew=0.5, label=r"Taylor--Hood (reference)", zorder=3)
        ax.plot(nt, kmi, "--", marker=MK[c], color=COL[c], mfc="white", mec=COL[c], mew=0.9, label=r"Mini (screening)", zorder=3)
        ax.set_xscale("log")
        ax.set_ylabel(r"$k_D$ [px$^2$]")
        ax.set_title(c, color=COL[c], fontsize=8, pad=3)
        ax.text(0.97, 0.06, r"$\Delta k_{\mathrm{last}}=%+.3f\,\%%$" % dlast, transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
        if c == "R3-V":
            ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: r"$%d$" % v))
        if c == "R3-H":
            ax.legend(loc="center right", frameon=False, handlelength=2.4)
        ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([6e4, 1e5, 2e5, 4e5]))
        ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter([r"$6\times10^4$", r"$10^5$", r"$2\times10^5$", r"$4\times10^5$"]))
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        plabel(ax, lab)
        src[c] = dict(n_tri=nt, kD_taylor_hood=kth, kD_mini=kmi, delta_last_percent=dlast,
                      mini_vs_TH_finest_percent=100 * (kmi[-1] / kth[-1] - 1))
    axs[-1].set_xlabel(r"Number of triangles $N_{\mathrm{tri}}$ [--]")
    fig.tight_layout(h_pad=0.8)
    save(fig, "figure2", src)


# ------------------------------------------------------------------ figure 8
def fig8(M):
    D = {c: fields(c) for c in CASES}
    vmax = max(float(D[c][2][np.unique(D[c][1])].max()) for c in CASES)
    snorm = LogNorm(vmin=1e-3, vmax=10 ** np.ceil(np.log10(vmax)))
    fig = plt.figure(figsize=(131 * MM, 172 * MM))
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 0.06], hspace=0.28, wspace=0.08)
    src = {"speed_normalisation": "U = Q/W", "speed_colour_range": [1e-3, snorm.vmax]}
    for r, c in enumerate(CASES):
        P, Ts, sp, pst, meta = D[c]
        pore, sel = M[c]
        bb = bbox(sel)
        tri = Triangulation(P[0], P[1], Ts.T)
        a0 = fig.add_subplot(gs[r, 0]); show_domain(a0, pore, sel); frame(a0, bb)
        # inlet / outlet bars
        if meta["direction"] == "lr":
            a0.plot([bb[0], bb[0]], [bb[2], bb[3]], color=VERM, lw=3, clip_on=False, solid_capstyle="butt")
            a0.plot([bb[1], bb[1]], [bb[2], bb[3]], color="k", lw=3, clip_on=False, solid_capstyle="butt")
        else:
            a0.plot([bb[0], bb[1]], [bb[3], bb[3]], color=VERM, lw=3, clip_on=False, solid_capstyle="butt")
            a0.plot([bb[0], bb[1]], [bb[2], bb[2]], color="k", lw=3, clip_on=False, solid_capstyle="butt")
        a0.plot([bb[0], bb[0] + 500], [bb[2] - 0.05 * (bb[3] - bb[2])] * 2, color="k", lw=1.5, clip_on=False)
        a0.text(bb[0] + 250, bb[2] - 0.08 * (bb[3] - bb[2]), r"500 px", ha="center", va="top", fontsize=7)
        a0.set_ylabel(c, fontsize=8, labelpad=8)
        a1 = fig.add_subplot(gs[r, 1])
        pt = pst[Ts].mean(axis=0)
        im1 = a1.tripcolor(tri, facecolors=np.clip(pt, 0, 1), cmap="viridis", norm=Normalize(0, 1), rasterized=True)
        frame(a1, bb)
        a2 = fig.add_subplot(gs[r, 2])
        st = np.clip(sp[Ts].mean(axis=0), 1e-3, None)
        im2 = a2.tripcolor(tri, facecolors=st, cmap="inferno", norm=snorm, rasterized=True)
        frame(a2, bb)
        for ax, lab in zip((a0, a1, a2), "abcdefghi"[3 * r:3 * r + 3]):
            plabel(ax, lab, x=0.0, y=1.02)
        vin = np.unique(Ts)
        src[c] = dict(bbox_cropped_frame=[float(v) for v in bb], speed_max_over_U=float(sp[vin].max()),
                      pstar_vertex_below0=float(np.mean(pst[vin] < 0)), pstar_vertex_above1=float(np.mean(pst[vin] > 1)))
    c1 = fig.add_subplot(gs[3, 1]); c2 = fig.add_subplot(gs[3, 2])
    cb1 = fig.colorbar(im1, cax=c1, orientation="horizontal", extend="both"); cb1.set_label(r"$p^{*}$ (clipped to $[0,1]$)")
    cb2 = fig.colorbar(im2, cax=c2, orientation="horizontal", extend="min"); cb2.set_label(r"$|\mathbf{u}|/U$, $U=Q/W$")
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=7, width=0.5)
    save(fig, "figure8", src)


# ------------------------------------------------------------------ figure 9
def fig9(M, desc):
    cmap = ListedColormap(["#D9D9D9", "#9ECAE1", "#08519C"])
    fig = plt.figure(figsize=(146 * MM, 108 * MM))
    D = {c: fields(c) for c in CASES}
    bbs = {c: bbox(M[c][1]) for c in CASES}
    widths = [(bbs[c][1] - bbs[c][0]) / (bbs[c][3] - bbs[c][2]) for c in CASES]
    gs = fig.add_gridspec(2, 3, width_ratios=widths, height_ratios=[1, 0.62], hspace=0.32, wspace=0.06)
    for i, c in enumerate(CASES):
        P, Ts, sp, pst, meta = D[c]
        ut = sp[Ts].mean(axis=0); umax = sp[np.unique(Ts)].max()
        cls = np.where(ut >= 0.01 * umax, 2, np.where(ut <= 0.001 * umax, 0, 1))
        ax = fig.add_subplot(gs[0, i])
        ax.tripcolor(Triangulation(P[0], P[1], Ts.T), facecolors=cls, cmap=cmap, vmin=-0.5, vmax=2.5, rasterized=True)
        frame(ax, bbs[c]); plabel(ax, "abc"[i], x=0.0, y=1.02)
    ax = fig.add_subplot(gs[1, :])
    y = np.arange(3)[::-1]
    lab = {"active": "#08519C", "intermediate": "#9ECAE1", "dead": "#D9D9D9"}
    left = np.zeros(3)
    for k, col in lab.items():
        v = np.array([100 * desc[c][k] for c in CASES])
        ax.barh(y, v, left=left, color=col, edgecolor="k", lw=0.5, height=0.62, label=k)
        for yy, l, vv in zip(y, left, v):
            ax.text(l + vv / 2, yy, r"%.1f" % vv, ha="center", va="center", fontsize=7, color="white" if k == "active" else "k")
        left += v
    ax.set_yticks(y); ax.set_yticklabels(CASES); ax.set_xlim(0, 100)
    ax.set_xlabel(r"Area share of the percolating domain [\%]")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    ax.tick_params(axis="y", which="both", right=False, left=False)
    plabel(ax, "d", x=-0.06, y=1.02)
    save(fig, "figure9", {c: {k: desc[c][k] for k in ("active", "intermediate", "dead")} for c in CASES})


# ------------------------------------------------------------------ figure 10
def fig10(desc):
    fig, axs = plt.subplots(2, 2, figsize=(146 * MM, 118 * MM))
    (a, b), (c_, d) = axs
    src = {}
    edges = np.logspace(-8, 3, 111)
    for c in CASES:
        z = np.load(FH / f"{c.replace('-', '')}_tri.npz")
        u, A, s, ps = z["ut_over_U"], z["area"], z["s_star"], z["pstar"]
        nb = 100; e = np.linspace(0, 1, nb + 1); idx = np.clip(np.digitize(s, e) - 1, 0, nb - 1)
        pa = np.bincount(idx, ps * A, nb) / np.maximum(np.bincount(idx, A, nb), 1e-300)
        sc = 0.5 * (e[1:] + e[:-1]); ok = np.bincount(idx, A, nb) > 0
        a.plot(sc[ok], pa[ok], color=COL[c], label=c)
        h, _ = np.histogram(u, edges, weights=A)
        dens = h / A.sum() / np.diff(np.log10(edges)); xc = np.sqrt(edges[1:] * edges[:-1]); m = dens > 0
        b.plot(xc[m], dens[m], color=COL[c])
        b.axvline(z["umax_over_U"], color=COL[c], ls=":", lw=0.9)
        o = np.argsort(u); F = np.cumsum(A[o]) / A.sum()
        c_.plot(u[o][::50], F[::50], color=COL[c])
        o = np.argsort(-u); ac = np.cumsum(A[o]) / A.sum(); tw = u[o] * A[o]; tau = np.cumsum(tw) / tw.sum()
        k = max(1, len(ac) // 4000)
        d.plot(np.r_[0, ac[::k], 1], np.r_[0, tau[::k], 1], color=COL[c])
        d.plot([desc[c]["A50"], desc[c]["A90"]], [0.5, 0.9], MK[c], color=COL[c], mec="k", mew=0.5, ms=4.5, ls="none")
        src[c] = dict(umax_over_U=float(z["umax_over_U"]), A50=desc[c]["A50"], A90=desc[c]["A90"],
                      area_below_1em8=float(A[u < 1e-8].sum() / A.sum()), area_zero_speed=float(A[u == 0].sum() / A.sum()),
                      F_A={str(x): float(A[u <= x].sum() / A.sum()) for x in (1e-3, 1e-2, 1e-1, 1, 10)},
                      q05=float(np.sort(u)[np.searchsorted(np.cumsum(A[np.argsort(u)]) / A.sum(), 0.05)]))
    a.plot([0, 1], [1, 0], "k:", lw=0.8, label="uniform")
    a.set_xlim(0, 1); a.set_xlabel(r"Normalized flow coordinate $s/L_s$ [--]"); a.set_ylabel(r"$\langle p^{*}\rangle$ [--]")
    a.legend(frameon=False, loc="upper right")
    b.set_xscale("log"); b.set_yscale("log"); b.set_xlabel(r"$|\mathbf{u}|/U$ [--]"); b.set_ylabel(r"Area-weighted density $f_A$ [--]")
    c_.set_xscale("log"); c_.set_xlabel(r"$|\mathbf{u}|/U$ [--]"); c_.set_ylabel(r"Cumulative area fraction $F_A$ [--]"); c_.set_ylim(0, 1)
    b.set_xlim(1e-8, 1e3); c_.set_xlim(1e-8, 1e3)
    d.plot([0, 1], [0, 1], "k:", lw=0.8); d.set_xlim(0, 1); d.set_ylim(0, 1.02)
    d.set_xlabel(r"Cumulative area fraction (fastest first) [--]"); d.set_ylabel(r"Cumulative transport $T$ [--]")
    tab = r"\begin{tabular}{lrr}& $A_{50}$ [\%] & $A_{90}$ [\%]\\" + "".join(
        r"%s & %.1f & %.1f\\" % (c, 100 * desc[c]["A50"], 100 * desc[c]["A90"]) for c in CASES) + r"\end{tabular}"
    d.text(0.97, 0.05, tab, transform=d.transAxes, ha="right", va="bottom", fontsize=6.5)
    for ax, l in zip((a, b, c_, d), "abcd"):
        plabel(ax, l, x=-0.2)
    fig.tight_layout(w_pad=1.2, h_pad=1.0)
    save(fig, "figure10", src)


# ------------------------------------------------------------------ figure 11
def fig11(M, desc):
    fig, axs = plt.subplots(2, 2, figsize=(146 * MM, 128 * MM))
    (a, b), (c_, d) = axs
    x = np.arange(3)
    chi = {c: float(M[c][1].sum() / M[c][0].sum()) for c in CASES}
    for i, c in enumerate(CASES):
        k = desc[c]["kD"]; U = desc[c]["U"]
        a.plot(i, k, MK[c], color=COL[c], mec="k", mew=0.5, ms=6)
        a.annotate(r"$%s$" % (r"%.2f" % k if k < 1e3 else r"%.0f" % k), (i, k), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7)
        b.plot(i, U, MK[c], color=COL[c], mec="k", mew=0.5, ms=6)
        b.annotate((r"$%.2f$" % U) if U >= 1 else (r"$%.2f\times10^{%d}$" % (U / 10 ** np.floor(np.log10(U)), np.floor(np.log10(U)))), (i, U),
                   xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7)
        d.plot(i, desc[c]["umax_over_U"], MK[c], color=COL[c], mec="k", mew=0.5, ms=6)
        d.annotate((r"$%.0f$" if desc[c]["umax_over_U"] >= 10 else r"$%.1f$") % desc[c]["umax_over_U"], (i, desc[c]["umax_over_U"]), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7)
    a.set_yscale("log"); a.set_ylabel(r"$k_D$ [px$^2$]"); a.set_ylim(3, 3e5)
    a.text(0.03, 0.95, r"$A_k=k_D(\mathrm{R3\mbox{-}V})/k_D(\mathrm{R3\mbox{-}H})=%.0f$" % desc["A_k_R3"], transform=a.transAxes, va="top", fontsize=7)
    b.set_yscale("log"); b.set_ylabel(r"$U=Q/W$ at $p_{\mathrm{in}}=1$ [code units]"); b.set_ylim(1e-3, 1e2)
    w = 0.36
    c_.bar(x - w / 2, [100 * desc[c]["porosity"] for c in CASES], w, color="0.35", edgecolor="k", lw=0.5, label=r"$\phi$")
    c_.bar(x + w / 2, [100 * chi[c] for c in CASES], w, color="0.8", edgecolor="k", lw=0.5, label=r"$\chi_c$")
    c_.set_ylabel(r"Fraction [\%]"); c_.set_ylim(0, 115); c_.legend(ncol=2, frameon=False, loc="upper center")
    c_.set_xticks(x); c_.set_xticklabels([r"%s\\$D_h=%.0f$" % (c, desc[c]["D_h"]) for c in CASES])
    c_.set_xticklabels([r"\shortstack{%s\\$D_h=%.0f$}" % (c, desc[c]["D_h"]) for c in CASES])
    d.set_yscale("log"); d.set_ylabel(r"Velocity amplification $|\mathbf{u}|_{\max}/U$ [--]")
    lo = min(desc[c]["umax_over_U"] for c in CASES); hi = max(desc[c]["umax_over_U"] for c in CASES)
    d.set_ylim(10 ** np.floor(np.log10(lo)) , 10 ** (np.ceil(np.log10(hi)) + 0.3))
    for ax, l in zip((a, b, c_, d), "abcd"):
        ax.set_xticks(x); ax.set_xlim(-0.5, 2.5); plabel(ax, l, x=-0.2)
        if ax is not c_:
            ax.set_xticklabels(CASES)
        ax.tick_params(axis="x", which="both", top=False, bottom=False)
    fig.tight_layout(w_pad=1.2, h_pad=1.2)
    save(fig, "figure11", dict(desc=desc, chi_c=chi))


def _same(a, b, rtol=1e-9):
    """Recursive comparison of two JSON trees; floats within a relative rtol."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k], rtol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y, rtol) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return a == b or abs(a - b) <= rtol * max(abs(a), abs(b))
    return a == b


def compare_sources(which):
    """Compare the freshly written *_source.json files with the archived ones. Returns the exit status."""
    ref_dir = ROOT / "data/results/figure_sources"
    ok = True
    for n in which:
        name = f"figure{n}_source.json"
        same = _same(json.load(open(OUT / name)), json.load(open(ref_dir / name)))
        ok &= same
        print(f"{'IDENTICAL' if same else 'DIFFERS':<10} {name}")
    print("PASS" if ok else "FAIL", f"(figures in {OUT})")
    return 0 if ok else 1


if __name__ == "__main__":
    which = sys.argv[2:] or ["2", "8", "9", "10", "11"]
    if "2" in which:
        fig2()
    if not set(which) - {"2"}:
        sys.exit(compare_sources(which))
    if set(which) - {"2"}:
        M = masks(); desc = json.load(open(FH / "p1_descriptors.json"))
        if "8" in which: fig8(M)
        if "9" in which: fig9(M, desc)
        if "10" in which: fig10(desc)
        if "11" in which: fig11(M, desc)
    sys.exit(compare_sources(which))
