#!/usr/bin/env python3
"""
TEMPLATE-ID: FS-house-style-01 (ADAPT: style module used as-is, plotting logic new; no matching template for FE-field or convergence figures)
FIGURE TYPE: convergence curves, scalar field maps on triangular meshes, area-weighted distributions, bar summaries
DESCRIPTION: Regenerates manuscript Figs 1-7 from the author's RAW outputs (global_summary.csv, solution.vtk,
             full_percolating_domain.npy, structures/*.png) and draws the QC evidence figures E1-E3 from evidence/*.json.
INPUT FILES: PROJETO/convergencia/convergence_mesh_0*/global_summary.csv; .../solution.vtk; .../full_percolating_domain.npy;
             convergencia/github/structures/*.png; QC evidence/*.json (E-figures)
OUTPUT FILES: figures-regenerated/QC-F{1..7}.{svg,pdf,png}, QC-E{1..3}.{svg,pdf,png}
PARAMETERS: --fig {1..7|E1|E2|E3|all} --out OUTDIR
DEPENDENCIES: numpy scipy matplotlib meshio pillow; LaTeX (mathptmx, dvipng) for exact house style
STYLE: House Style Contract v1 (FS-house-style-01). Deviations: Okabe-Ito colours instead of the template cycle (colour-vision safe, FIG-07);
       width 180 mm (npj Computational Materials two-column) = 7.09 in; field maps are rasterised inside the vector files.
PHYSICS NOTES: p* = (p - <p>_out)/dP with dP = |<p>_in - <p>_out| as defined by the author's code (nodal means over boundary DOFs);
               u* = |u|/U0 with U0 = 1e-4; active/dead thresholds 1e-2 and 1e-3 of max |u| (author's definition); every number is read, not typed.
ORIGIN: CFD digital-rock project, QC-20261019-v1, Sonnet 5.5 (declared high effort)
REVISION LOG: 2026-10-20 first version
              2026-10-31 ntri_axis: a literal TAB had replaced the backslash of \times, so every N_tri tick
                         read '3imes10^4' in QC-F1/F5/F6. Fixed.
              2026-10-31 tex() fallback fixed: without dvipng the house style falls back to mathtext, which
                         rejected the literal percent inside $...$ of the chi_c panel labels (ValueError from
                         _mathtext.parse). Install dvipng to reproduce the audited typography exactly.
"""
import argparse, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "templates")); sys.path.insert(0, str(Path(__file__).parent))
import house_style as hs                                        # noqa: E402
from house_style import plt, panel_label, save                  # noqa: E402
import matplotlib as mpl                                        # noqa: E402
from matplotlib.gridspec import GridSpec                        # noqa: E402
from matplotlib.colors import ListedColormap, PowerNorm, Normalize, BoundaryNorm  # noqa: E402
from matplotlib.tri import Triangulation                        # noqa: E402
import qcdata as D                                              # noqa: E402

mpl.rcParams.update({"figure.dpi": 300, "xtick.top": False, "ytick.right": False, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8})
W = 7.09                                                         # in, 180 mm
CASES = list(D.CASES)
COL = {"R3-H": "#0072B2", "R2-V": "#E69F00", "R3-V": "#009E73"}  # Okabe-Ito, one colour per configuration in every figure
MK = {"R3-H": "o", "R2-V": "s", "R3-V": "^"}
ACC = "#D55E00"; SKY = "#56B4E9"; SOLID = "#2b2b2b"; BAND = "#cfcfcf"
MKW = dict(mec="black", mew=0.5, ms=4.5)


def PL(ax, s, x=-0.16, y=1.02): panel_label(ax, s, x=x, y=y)
def tex(s):
    """House style renders text through LaTeX. Without it (no dvipng) matplotlib falls back to mathtext, which
    cannot parse a literal percent sign inside math mode, so close the math first and append the sign as text."""
    if hs.USE_TEX:
        return s
    s = s.replace(r"\,\%$", r"$ \%")          # "$...53.8\,\%$" -> "$...53.8$ \%"
    return s.replace(r"\,", " ").replace(r"\%", "%")
def lab(q, u): return hs.unit(q, u)


def ntri_axis(ax):
    """Log x axis for N_tri with a few readable labels (no overlapping minor labels)."""
    ax.set_xscale("log"); ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([3e4, 1e5, 3e5])); ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: tex(r"$%d\times10^{%d}$" % (round(v / 10 ** int(np.floor(np.log10(v)))), int(np.floor(np.log10(v)))))))
    ax.set_xlim(2.4e4, 4.6e5)


def imgax(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_linewidth(0.5)


# ------------------------------------------------------------------ Figure 1
def fig1(out):
    fig = plt.figure(figsize=(W, 5.3))
    gs = GridSpec(2, 2, figure=fig, left=0.075, right=0.985, bottom=0.095, top=0.975, wspace=0.30, hspace=0.42)
    sub = gs[0, 0].subgridspec(2, 3, wspace=0.06, hspace=0.12)
    cmap = ListedColormap([SOLID, "white", SKY, BAND])
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
            ax.imshow(sm, cmap=cmap, vmin=0, vmax=3, interpolation="nearest", aspect="equal"); imgax(ax)
            if i == 0: ax.set_title(smp, fontsize=9, pad=2)
            chi = None
            try:
                chi = D.summary({"R2": "R2-V", "R3": "R3-H" if d == "lr" else "R3-V"}[smp], 4)["selected_fraction_of_cropped_pores"] if (smp != "R2" or d == "tb") and smp != "R1" else None
            except KeyError:
                chi = None
            ax.set_xlabel(tex(r"$\chi_c=%.1f\,\%%$" % (100 * chi)) if chi else "non-percolating", fontsize=7, labelpad=1.5)
            if j == 0: ax.set_ylabel(nm, fontsize=8, labelpad=2)
    PL(fig.axes[0], "a", x=-0.45, y=1.12)
    for (r, c, case, letter) in ((0, 1, "R3-H", "b"), (1, 0, "R2-V", "c"), (1, 1, "R3-V", "d")):
        ax = fig.add_subplot(gs[r, c])
        n = D.series(case, "triangle_count"); k = D.series(case, "apparent_permeability")
        ax.axhspan(k[-1] * 0.99, k[-1] * 1.01, color="0.88", lw=0, zorder=0)
        ax.axhline(k[-1], color="black", lw=0.6, ls="--", zorder=1)
        ax.plot(n, k, "-", color=COL[case], lw=1.0, zorder=2)
        ax.plot(n, k, MK[case], color=COL[case], zorder=3, **MKW)
        ntri_axis(ax); ax.set_xlabel(lab(r"Number of triangles $N_{\mathrm{tri}}$", "-"))
        ax.set_ylabel(lab(r"$k_{\mathrm{app}}$", r"px$^2$"))
        dk = 100 * abs(k[-1] - k[-2]) / k[-1]
        ax.text(0.97, 0.95 if case != "R3-H" else 0.06, tex(r"$\Delta k_{\mathrm{last}}=%.2f\,\%%$" % dk), transform=ax.transAxes, ha="right", va="top" if case != "R3-H" else "bottom", fontsize=8)
        ax.set_title(case, fontsize=9, color=COL[case], pad=3)
        PL(ax, letter)
    save(fig, "QC-F1", out)


# ------------------------------------------------------------------ Figure 5
def fig5(out):
    fig, ax = plt.subplots(2, 2, figsize=(W, 5.0), gridspec_kw=dict(left=0.09, right=0.985, bottom=0.09, top=0.975, wspace=0.32, hspace=0.38))
    a = ax[0, 0]; ex = []
    for c in CASES:
        n = D.series(c, "triangle_count"); d = D.series(c, "total_dof_count")
        m, b = np.polyfit(np.log(n), np.log(d), 1); ex.append(m)
        a.plot(n, d, MK[c], color=COL[c], label=c, **MKW); a.plot(n, np.exp(b) * n ** m, "--", color=COL[c], lw=0.7)
    ntri_axis(a); a.set_yscale("log"); a.set_xlabel(lab(r"$N_{\mathrm{tri}}$", "-")); a.set_ylabel(lab(r"$N_{\mathrm{dof}}$", "-"))
    a.text(0.97, 0.06, tex(r"median $m=%.2f$" % np.median(ex)), transform=a.transAxes, ha="right", fontsize=8); a.legend(loc="upper left")
    b_ = ax[0, 1]; rate = []
    for c in CASES:
        d = D.series(c, "total_dof_count"); mb = D.series(c, "system_csr_megabytes"); rate += list(mb / d * 1e6)
        b_.plot(d / 1e6, mb, MK[c], color=COL[c], label=c, **MKW)
    b_.set_xlabel(lab(r"$N_{\mathrm{dof}}$", "10$^6$")); b_.set_ylabel(lab(r"$M_{\mathrm{CSR}}$", "MB"))
    b_.text(0.03, 0.93, tex(r"median %.0f MB per $10^6$ DOF" % np.median(rate)), transform=b_.transAxes, fontsize=8, va="top")
    c_ = ax[1, 0]
    for c in CASES:
        k = D.series(c, "apparent_permeability"); mb = D.series(c, "system_csr_megabytes")
        c_.plot(mb[:-1], 100 * np.abs(k[:-1] - k[-1]) / k[-1], "-" + MK[c], color=COL[c], lw=0.8, label=c, **MKW)
    c_.set_yscale("log"); c_.axhline(1, color="black", ls=":", lw=0.7); c_.axhline(2, color="black", ls="--", lw=0.7)
    c_.set_xlabel(lab(r"$M_{\mathrm{CSR}}$", "MB")); c_.set_ylabel(lab(r"$|k-k_5|/k_5$", r"\%"))
    d_ = ax[1, 1]; x = np.arange(3); wd = 0.36
    dM = [100 * (D.series(c, "system_csr_megabytes")[4] / D.series(c, "system_csr_megabytes")[3] - 1) for c in CASES]
    dk = [100 * abs(D.series(c, "apparent_permeability")[4] - D.series(c, "apparent_permeability")[3]) / D.series(c, "apparent_permeability")[4] for c in CASES]
    d_.bar(x - wd / 2, dM, wd, color="0.55", ec="black", lw=0.5, label=r"$\Delta M_{\mathrm{last}}$")
    d_.bar(x + wd / 2, dk, wd, color=[COL[c] for c in CASES], ec="black", lw=0.5, label=r"$\Delta k_{\mathrm{last}}$")
    d_.set_yscale("log"); d_.set_xticks(x); d_.set_xticklabels(CASES); d_.set_ylabel(lab("Change, level 4 to 5", r"\%"))
    for xi, v in zip(x + wd / 2, dk): d_.text(xi, v * 1.15, tex("%.2f" % v), ha="center", fontsize=7)
    for xi, v in zip(x - wd / 2, dM): d_.text(xi, v * 1.15, tex("%.0f" % v), ha="center", fontsize=7)
    from matplotlib.patches import Patch
    d_.set_ylim(5e-3, 4e3)
    d_.legend(handles=[Patch(fc="0.55", ec="black", lw=0.5, label=r"$\Delta M_{\mathrm{last}}$"), Patch(fc="white", ec="black", lw=0.5, hatch="///", label=r"$\Delta k_{\mathrm{last}}$ (colour = configuration)")], loc="upper center", ncol=1, fontsize=7)
    for a_, s in zip(ax.ravel(), "abcd"): PL(a_, s)
    c_.legend(loc="lower left")
    save(fig, "QC-F5", out)


# ------------------------------------------------------------------ Figure 6
def fig6(out):
    fig, ax = plt.subplots(2, 2, figsize=(W, 5.0), gridspec_kw=dict(left=0.09, right=0.985, bottom=0.085, top=0.975, wspace=0.32, hspace=0.38))
    x = np.arange(3)
    S = {c: D.summary(c, 4) for c in CASES}
    k = np.array([S[c]["apparent_permeability"] for c in CASES]); a = ax[0, 0]
    for i, c in enumerate(CASES): a.plot(i, k[i], MK[c], color=COL[c], ms=6, mec="black", mew=0.5)
    a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(CASES); a.set_xlim(-0.5, 2.5); a.set_ylim(40, 2e5)
    a.set_ylabel(lab(r"$k_{\mathrm{app}}$", r"px$^2$"))
    for xi, v in zip(x, k): a.text(xi, v * 1.5, tex(r"$%.1f$" % v if v < 1e3 else r"$%.2f\times10^{%d}$" % (v / 10 ** int(np.log10(v)), int(np.log10(v)))), ha="center", fontsize=7.5)
    a.text(0.5, 0.06, tex(r"$A_k=%.1f$ (R3-V/R3-H)" % (k[2] / k[0])), transform=a.transAxes, ha="center", va="bottom", fontsize=8)
    g = [S[c]["pressure_drop"] / S[c]["flow_length"] for c in CASES]; b = ax[0, 1]
    for i, c in enumerate(CASES): b.plot(i, g[i], MK[c], color=COL[c], ms=6, mec="black", mew=0.5)
    b.set_yscale("log"); b.set_xticks(x); b.set_xticklabels(CASES); b.set_xlim(-0.5, 2.5); b.set_ylim(1e-9, 5e-6)
    for xi, v in zip(x, g): b.text(xi, v * 1.8, tex(r"$%.2f\times10^{%d}$" % (v / 10 ** int(np.floor(np.log10(v))), int(np.floor(np.log10(v))))), ha="center", fontsize=7.5)
    b.set_ylabel(lab(r"$\Delta p/L$", "code units"))
    c_ = ax[1, 0]; w = 0.36
    c_.bar(x - w / 2, [100 * S[c]["domain_porosity"] for c in CASES], w, color="0.35", ec="black", lw=0.5, label=r"$\phi$")
    c_.bar(x + w / 2, [100 * S[c]["selected_fraction_of_cropped_pores"] for c in CASES], w, color="0.8", ec="black", lw=0.5, label=r"$\chi_c$")
    c_.set_xticks(x); c_.set_xticklabels([tex(r"%s" % c + "\n" + r"$D_h=%.0f$" % S[c]["hydraulic_diameter_2d"]) for c in CASES]); c_.set_ylabel(lab("Fraction", r"\%"))
    c_.set_ylim(0, 112); c_.legend(loc="upper right", ncol=2)
    d_ = ax[1, 1]
    for i, c in enumerate(CASES):
        d_.plot(i, S[c]["reynolds_number"], MK[c], color=COL[c], **MKW)
        d_.plot(i, S[c]["speed_max"] / D.U0, MK[c], mfc="white", mec=COL[c], mew=1.0, ms=5)
    d_.set_yscale("log"); d_.set_xticks(x); d_.set_xticklabels(CASES); d_.axhline(1, color="black", ls="--", lw=0.6); d_.set_xlim(-0.5, 2.5)
    d_.set_ylabel(lab(r"$\mathrm{Re}$ (filled), $|\mathbf{u}|_{\max}/U_0$ (open)", "-"))
    for a_, s in zip(ax.ravel(), "abcd"): PL(a_, s)
    save(fig, "QC-F6", out)


# ------------------------------------------------------------------ field figures
def tri_of(c):
    f = D.field(c); return Triangulation(f["P"][:, 0], f["P"][:, 1], f["T"]), f


def scale_bar(ax, length, label, f):
    x0 = f["P"][:, 0].min() + 0.04 * np.ptp(f["P"][:, 0]); y0 = f["P"][:, 1].min() + 0.03 * np.ptp(f["P"][:, 1])
    ax.plot([x0, x0 + length], [y0, y0], color="black", lw=1.4, solid_capstyle="butt"); ax.text(x0 + length / 2, y0 + 0.012 * np.ptp(f["P"][:, 1]), label, ha="center", fontsize=7)


def fig2(out):
    fig = plt.figure(figsize=(W, 8.6))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 0.07], left=0.04, right=0.985, bottom=0.05, top=0.975, wspace=0.04, hspace=0.06)
    umax = max(D.field(c)["ustar"].max() for c in CASES)
    pn = Normalize(0, 1); un = PowerNorm(0.5, vmin=0, vmax=umax)
    meta = {}
    for r, c in enumerate(CASES):
        tri, f = tri_of(c)
        ax = fig.add_subplot(gs[r, 0]); a = ax
        # domain panel: percolating domain in the cropped frame (blue), solid not drawn (white = outside retained domain)
        ax.tripcolor(tri, np.ones(len(f["P"])), cmap=ListedColormap([SKY]), rasterized=True, shading="gouraud"); ax.set_aspect("equal"); imgax(ax)
        ax1 = fig.add_subplot(gs[r, 1]); im1 = ax1.tripcolor(tri, np.clip(f["pstar"], 0, 1), cmap="cividis", norm=pn, rasterized=True, shading="gouraud"); ax1.set_aspect("equal"); imgax(ax1)
        ax2 = fig.add_subplot(gs[r, 2]); im2 = ax2.tripcolor(tri, f["ustar"], cmap="viridis", norm=un, rasterized=True, shading="gouraud"); ax2.set_aspect("equal"); imgax(ax2)
        for axx, s in zip((ax, ax1, ax2), "abcdefghi"[3 * r:3 * r + 3]): axx.text(0.02, 0.985, "(%s)" % s, transform=axx.transAxes, fontsize=9, va="top", ha="left")
        ax.set_ylabel(c, fontsize=9)
        scale_bar(ax, 500, "500 px", f)
        over = float((f["pstar"] > 1).mean() * 100); under = float((f["pstar"] < 0).mean() * 100)
        meta[c] = (over, under)
        if r == 0:
            ax.set_title("Percolating domain", fontsize=9); ax1.set_title(tex(r"$p^*=(p-\langle p\rangle_{\mathrm{out}})/\Delta p$"), fontsize=9); ax2.set_title(tex(r"$|\mathbf{u}|/U_0$"), fontsize=9)
    cb1 = fig.colorbar(im1, cax=fig.add_subplot(gs[3, 1]), orientation="horizontal", extend="both"); cb1.set_label(tex(r"$p^*$ (clipped to [0, 1])"), fontsize=8)
    cb2 = fig.colorbar(im2, cax=fig.add_subplot(gs[3, 2]), orientation="horizontal"); cb2.set_label(tex(r"$|\mathbf{u}|/U_0$ (power-law scale, exponent 0.5)"), fontsize=8)
    print("p* outside [0,1] (% of vertices over/under):", meta)
    save(fig, "QC-F2", out)


def classes(f):
    st = f["speed"][f["T"]].mean(1); um = f["speed"].max()
    return np.where(st >= 0.01 * um, 2, np.where(st <= 0.001 * um, 0, 1)), st


def fig3(out):
    fig = plt.figure(figsize=(W, 4.0))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 0.42], left=0.03, right=0.985, bottom=0.12, top=0.95, wspace=0.05, hspace=0.18)
    cm = ListedColormap(["#d9d9d9", "#9ecae1", "#08519c"]); shares = {}
    for j, c in enumerate(CASES):
        tri, f = tri_of(c); cl, st = classes(f)
        ax = fig.add_subplot(gs[0, j]); ax.tripcolor(tri, facecolors=cl, cmap=cm, vmin=0, vmax=2, rasterized=True); ax.set_aspect("equal"); imgax(ax)
        ax.set_title(c, fontsize=9, color=COL[c]); ax.text(0.02, 0.985, "(%s)" % "abc"[j], transform=ax.transAxes, fontsize=9, va="top")
        a = f["area"]; shares[c] = [a[cl == k].sum() / a.sum() * 100 for k in (2, 1, 0)]
    ax = fig.add_subplot(gs[1, :]); names = ["active", "intermediate", "dead"]; cols = ["#08519c", "#9ecae1", "#d9d9d9"]
    for i, c in enumerate(CASES):
        left = 0
        for v, nme, cc in zip(shares[c], names, cols):
            ax.barh(i, v, left=left, color=cc, ec="black", lw=0.5, height=0.6, label=nme if i == 0 else None)
            ax.text(left + v / 2, i, tex("%.1f" % v), ha="center", va="center", fontsize=7.5, color="white" if cc == "#08519c" else "black"); left += v
    ax.set_yticks(range(3)); ax.set_yticklabels(CASES); ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel(lab("Area share of the percolating domain", r"\%"))
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.28)); ax.text(-0.075, 1.0, "(d)", transform=ax.transAxes, fontsize=9)
    save(fig, "QC-F3", out); print("class shares", shares)


def weighted_pdf(vals, w, edges):
    h, _ = np.histogram(vals, bins=edges, weights=w); return h / (w.sum() * np.diff(edges))


def fig4(out):
    fig, ax = plt.subplots(2, 2, figsize=(W, 5.4), gridspec_kw=dict(left=0.095, right=0.985, bottom=0.085, top=0.975, wspace=0.30, hspace=0.38))
    res = {}
    for c in CASES:
        f = D.field(c); T = f["T"]; a = f["area"]; P = f["P"]; s = D.summary(c, 4)
        ax_i = 0 if s["direction"] == "lr" else 1
        # (a) area-weighted axial pressure profile, flow coordinate normalised, 100 bins
        sc = P[T, ax_i].mean(1); pt = f["pstar"][T].mean(1); lo, hi = P[:, ax_i].min(), P[:, ax_i].max(); frac = (sc - lo) / (hi - lo)
        if s["direction"] == "tb": frac = 1 - frac
        nb = 100; idx = np.clip((frac * nb).astype(int), 0, nb - 1); wsum = np.bincount(idx, a, nb); prof = np.bincount(idx, a * pt, nb) / np.maximum(wsum, 1e-30)
        ax[0, 0].plot((np.arange(nb) + 0.5) / nb, prof, "-", color=COL[c], lw=1.0, label=c)
        # (b) full-range pdf, (c) cdf
        st = f["speed"][T].mean(1) / D.U0; pos = st > 0
        edges = np.logspace(np.log10(max(st[pos].min(), 1e-8)), np.log10(st.max()), 70)
        pdf = weighted_pdf(st[pos], a[pos], edges); xc = np.sqrt(edges[1:] * edges[:-1]); m = pdf > 0
        ax[0, 1].plot(xc[m], pdf[m], "-", color=COL[c], lw=1.0, label=c)
        o = np.argsort(st); cum = np.cumsum(a[o]) / a.sum(); ax[1, 0].plot(st[o][::50], cum[::50], "-", color=COL[c], lw=1.0, label=c)
        # (d) transport localisation
        o = np.argsort(-st); A = np.cumsum(a[o]) / a.sum(); Tn = np.cumsum(st[o] * a[o]); Tn /= Tn[-1]
        ax[1, 1].plot(np.r_[0, A[::50]], np.r_[0, Tn[::50]], "-", color=COL[c], lw=1.0, label=c)
        a50 = A[np.searchsorted(Tn, 0.5)]; a90 = A[np.searchsorted(Tn, 0.9)]
        ax[1, 1].plot([a50, a90], [0.5, 0.9], MK[c], color=COL[c], **MKW); res[c] = dict(A50=a50, A90=a90, umax=st.max())
        ax[0, 1].axvline(st.max(), color=COL[c], lw=0.5, ls=":")
    ax[0, 0].plot([0, 1], [1, 0], ":", color="black", lw=0.7, label="uniform"); ax[0, 0].set_xlim(0, 1)
    ax[0, 0].set_xlabel(lab(r"Normalised flow coordinate $s/L$", "-")); ax[0, 0].set_ylabel(lab(r"$\langle p^*\rangle$", "-")); ax[0, 0].legend(loc="upper right")
    ax[0, 1].set_xscale("log"); ax[0, 1].set_yscale("log"); ax[0, 1].set_xlabel(lab(r"$|\mathbf{u}|/U_0$", "-")); ax[0, 1].set_ylabel(lab(r"Area-weighted density $f_A$", "-"))
    ax[1, 0].set_xscale("log"); ax[1, 0].set_xlim(1e-8, 1e2); ax[0, 1].set_xlim(1e-9, 1e2); ax[1, 0].set_xlabel(lab(r"$|\mathbf{u}|/U_0$", "-")); ax[1, 0].set_ylabel(lab(r"Cumulative area fraction $F_A$", "-")); ax[1, 0].set_ylim(0, 1.02)
    ax[1, 1].plot([0, 1], [0, 1], ":", color="black", lw=0.7); ax[1, 1].set_xlim(0, 1); ax[1, 1].set_ylim(0, 1.02)
    ax[1, 1].set_xlabel(lab("Cumulative area fraction (fastest first)", "-")); ax[1, 1].set_ylabel(lab("Cumulative transport $T$", "-"))
    txt = "\n".join(tex(r"%s: $A_{50}=%.1f\,\%%$, $A_{90}=%.1f\,\%%$" % (c, 100 * res[c]["A50"], 100 * res[c]["A90"])) for c in CASES)
    ax[1, 1].text(0.97, 0.05, txt, transform=ax[1, 1].transAxes, ha="right", fontsize=7.5)
    for a_, s in zip(ax.ravel(), "abcd"): PL(a_, s)
    save(fig, "QC-F4", out); print("fig4", res)


def fig7(out):
    fig = plt.figure(figsize=(W, 3.6))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 0.06], left=0.02, right=0.985, bottom=0.06, top=0.93, wspace=0.05, hspace=0.12)
    umax = max(D.field(c)["ustar"].max() for c in CASES); un = PowerNorm(0.5, vmin=0, vmax=umax)
    from matplotlib.tri import LinearTriInterpolator
    for j, c in enumerate(CASES):
        tri, f = tri_of(c); ax = fig.add_subplot(gs[0, j])
        im = ax.tripcolor(tri, f["ustar"], cmap="viridis", norm=un, rasterized=True, shading="gouraud"); ax.set_aspect("equal"); imgax(ax)
        P = f["P"]; nx = 420; xs = np.linspace(P[:, 0].min(), P[:, 0].max(), nx); ys = np.linspace(P[:, 1].min(), P[:, 1].max(), int(nx * np.ptp(P[:, 1]) / np.ptp(P[:, 0])))
        X, Y = np.meshgrid(xs, ys)
        U = LinearTriInterpolator(tri, f["vel"][:, 0])(X, Y); V = LinearTriInterpolator(tri, f["vel"][:, 1])(X, Y)
        ax.streamplot(X, Y, np.ma.filled(U, 0), np.ma.filled(V, 0), color="white", linewidth=0.45, arrowsize=0.5, density=1.6, broken_streamlines=False)
        ax.set_title(c, fontsize=9, color=COL[c]); ax.text(0.02, 0.985, "(%s)" % "abc"[j], transform=ax.transAxes, fontsize=9, va="top", color="white")
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 1]), orientation="horizontal"); cb.set_label(tex(r"$|\mathbf{u}|/U_0$ (power-law scale, exponent 0.5)"), fontsize=8)
    save(fig, "QC-F7", out)


# ------------------------------------------------------------------ evidence figures (QC-E*)
def E1(out):
    """Mesh convergence of k under the author's definition and under independent definitions."""
    S = D.load_json("indep_summary.json")
    fig, ax = plt.subplots(2, 3, figsize=(W, 5.2), gridspec_kw=dict(left=0.085, right=0.985, bottom=0.09, top=0.94, wspace=0.42, hspace=0.42))
    for j, c in enumerate(CASES):
        tag = D.CASES[c]["tag"]; rows = [r for r in S if r["case"] == tag]
        n = np.array([r["tri"] for r in rows]); ka = np.array([r["k_author"] for r in rows])
        a = ax[0, j]
        a.plot(n, ka / ka[-1], "-o", color="black", label="author (nominal $Q$, Mini)", **{**MKW, "mfc": "black"})
        a.plot(n, [r["kA_mini"] / rows[-1]["kA_mini"] for r in rows], "-s", color="#D55E00", label="actual flux, Mini", **MKW)
        a.plot(n, [r["kA_p2"] / rows[-1]["kA_p2"] for r in rows], "-^", color="#CC79A7", label="actual flux, P2/P1", **MKW)
        pm = []; pp = []
        for r in rows:
            L = r["level"] if "level" in r else None
        for lv in ("01", "02", "03", "04", "05"):
            for el, lst in (("mini", pm), ("p2", pp)):
                p = D.EVID / "indep" / f"pdrive_m{lv}_{tag}_{el}.json"
                import json
                q = json.load(open(p)); lst.append(abs(q["inlet"]["Q"]) * q["flow_len"] / q["k"]["nodal"]["dP"])
        a.plot(n, np.array(pm) / pm[-1], "-v", color="#009E73", label="pressure-driven, Mini", **MKW)
        a.plot(n, np.array(pp) / pp[-1], "-D", color="#0072B2", label="pressure-driven, P2/P1", **MKW)
        ntri_axis(a); a.set_xlabel(lab(r"$N_{\mathrm{tri}}$", "-")); a.set_ylabel(lab(r"$k/k_{\mathrm{level\,5}}$", "-")); a.set_title(c, fontsize=9, color=COL[c])
        a.axhspan(0.99, 1.01, color="0.9", lw=0)
        if j == 0: a.legend(fontsize=6.5, loc="lower right")
        b = ax[1, j]
        W_ = 4096.0
        b.plot(n, [r["k_author"] * r["Lin"] / W_ for r in rows], "-o", color="black", label="author, $Q/W$ with nominal $Q$", **{**MKW, "mfc": "black"})
        b.plot(n, [r["kD_p2"] * (r["W"] / W_) for r in rows], "-^", color="#CC79A7", label="velocity-driven P2/P1", **MKW)
        b.plot(n, np.array(pp) / W_, "-D", color="#0072B2", label="pressure-driven P2/P1", **MKW)
        ntri_axis(b); b.set_yscale("log"); b.set_xlabel(lab(r"$N_{\mathrm{tri}}$", "-")); b.set_ylabel(lab(r"$k$ ($W=4096$ px)", r"px$^2$"))
        b.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: tex(r"$%g$" % v))); b.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        b.set_ylim(*( (10, 40) if j == 0 else (20, 150) if j == 1 else (1e4, 5e4) ))
        if j == 0: b.legend(fontsize=6.5, loc="best")
    for a_, s in zip(ax.ravel(), "abcdef"): PL(a_, s, x=-0.28)
    save(fig, "QC-E1", out)


def E2(out):
    """Sensitivity of k to a purely numerical geometry parameter and to the cropping band."""
    import json
    fig, ax = plt.subplots(1, 3, figsize=(W, 2.9), gridspec_kw=dict(left=0.075, right=0.985, bottom=0.17, top=0.92, wspace=0.40))
    # (a) polygon tolerance
    tols = ["tol1.5", "tol1.0", "tol0.5"]; xv = [1.5, 1.0, 0.5]
    for c in CASES:
        tag = D.CASES[c]["tag"]
        v = [json.load(open(D.EVID / "indep_sens" / f"{t}_{tag}_velocity.json")) for t in tols]; p = [json.load(open(D.EVID / "indep_sens" / f"{t}_{tag}_pressure.json")) for t in tols]
        kv = np.array([abs(q["inlet"]["Q"]) * q["flow_len"] / q["k"]["nodal"]["dP"] for q in v]); kp = np.array([abs(q["inlet"]["Q"]) * q["flow_len"] / q["k"]["nodal"]["dP"] for q in p])
        ax[0].plot(xv, kv / kv[0], "-" + MK[c], color=COL[c], label=c + ", velocity-driven", **MKW)
        ax[0].plot(xv, kp / kp[0], "--" + MK[c], color=COL[c], mfc="white", mec=COL[c], mew=1.0, ms=4.5, lw=0.8, label=c + ", pressure-driven")
    ax[0].set_xlabel(lab("Contour tolerance", "px")); ax[0].set_ylabel(lab(r"$k/k(\mathrm{tol}=1.5)$", "-")); ax[0].invert_xaxis(); ax[0].legend(fontsize=5.8, loc="lower left")
    # (b) band fraction: author's k_app and classification
    import csv
    bands = {0.04: "band0.04", 0.08: None, 0.12: "band0.12", 0.16: "band0.16"}
    base = {(r["filename"][6:18], r["direction"]): r for r in csv.DictReader(open(D.PROJ / D.LEVELS[1] / "global_summary.csv"))}
    for c in CASES:
        pt = D.CASES[c]["patch"][6:]; d = D.CASES[c]["direction"]; xs, ys = [], []
        for bf, tg in bands.items():
            r = base[(pt, d)] if tg is None else next((x for x in csv.DictReader(open(D.QC / "rerun/sens" / f"sens_{tg}" / "global_summary.csv")) if x["filename"][6:18] == pt and x["direction"] == d), None)
            if r and r["status"] == "success": xs.append(100 * bf); ys.append(float(r["apparent_permeability"]))
            else: ax[1].plot(100 * bf, np.nan, "x")
        ax[1].plot(xs, np.array(ys) / ys[xs.index(8.0)], "-" + MK[c], color=COL[c], label=c, **MKW)
    ax[1].set_yscale("log"); ax[1].set_xlabel(lab("Cropping band", r"\% of image side")); ax[1].set_ylabel(lab(r"$k_{\mathrm{app}}/k_{\mathrm{app}}(8\,\%)$", "-")); ax[1].legend(fontsize=7)
    ax[1].text(4, 0.9, "R2-V\nnot percolating", fontsize=6.5, ha="left", color=COL["R2-V"], transform=ax[1].get_xaxis_transform())
    # (c) end loss: pressure profile R3-V etc from evidence
    for c in CASES:
        e = json.load(open(D.EVID / "indep" / f"endeff_m03_{D.CASES[c]['tag']}.json")); ax[2].plot(e["frac"], np.array(e["p"]) / max(e["p"]), "-", color=COL[c], lw=1.0, label=c)
    ax[2].plot([0, 1], [1, 0], ":", color="black", lw=0.7); ax[2].set_xlabel(lab(r"Flow coordinate $s/L$", "-")); ax[2].set_ylabel(lab(r"$\langle p\rangle/p_{\max}$", "-")); ax[2].legend(fontsize=7)
    for a_, s in zip(ax, "abc"): PL(a_, s, x=-0.24)
    save(fig, "QC-E2", out)


def E3(out):
    """Poiseuille verification: pipeline versus independent solvers."""
    import json
    P = json.load(open(D.EVID / "poiseuille" / "poiseuille_indep.json"))
    rows = json.load(open(D.QC.parents[2] / "revised" / "validation" / "validation_rows.json")) if (D.QC.parents[2] / "revised" / "validation" / "validation_rows.json").exists() else None
    import csv, glob
    pipe = {}
    for lv, d in enumerate(["v_01_coarse", "v_02_intermediate", "v_03_fine", "v_04_very_fine", "v_05_extra_fine"], 1):
        for r in csv.DictReader(open(D.QC / "rerun/poiseuille" / d / "global_summary.csv")):
            if r["status"] == "success":
                h = float(r["transverse_length"]); pipe[(lv, int(h))] = 100 * (float(r["apparent_permeability"]) / (h * h / 12) - 1)
    fig, ax = plt.subplots(1, 4, figsize=(W, 2.5), sharey=True, gridspec_kw=dict(left=0.075, right=0.99, bottom=0.2, top=0.88, wspace=0.08))
    for a, h in zip(ax, (20, 40, 80, 160)):
        lv = np.arange(1, 6); sel = [r for r in P if r["h"] == h]
        a.plot(lv, [pipe[(l, h)] for l in lv], "-o", color="black", label="author pipeline", **{**MKW, "mfc": "black"})
        a.plot(lv, [100 * r["k_p2_velocity"] for r in sel], "-^", color="#CC79A7", label="velocity-driven P2/P1", **MKW)
        a.plot(lv, [100 * r["k_p2_pressure"] for r in sel], "-D", color="#0072B2", label="pressure-driven P2/P1", **MKW)
        a.axhline(0, color="black", lw=0.5); a.set_xlabel("Mesh level"); a.set_title(tex(r"$h=%d$ px, $L/h=%d$" % (h, 20 if h <= 40 else 10)), fontsize=8)
        a.set_xticks(lv)
    ax[0].set_ylabel(lab(r"$k_{\mathrm{FE}}/(h^2/12)-1$", r"\%")); ax[0].legend(fontsize=6, loc="upper right")
    for a_, s in zip(ax, "abcd"): a_.text(0.03, 0.04, "(%s)" % s, transform=a_.transAxes, fontsize=9)
    save(fig, "QC-E3", out)


FIGS = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "7": fig7, "E1": E1, "E2": E2, "E3": E3}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--fig", default="all"); ap.add_argument("--out", default=str(D.QC / "figures-regenerated")); a = ap.parse_args()
    for k in (FIGS if a.fig == "all" else a.fig.split(",")):
        print("==", k); FIGS[k](a.out)
