#!/usr/bin/env python3
"""
TEMPLATE-ID: FS-house-style-01 (ADAPT: style module used as-is, plotting logic new)
FIGURE TYPE: sensitivity curves (relative deviation vs control parameter), convergence curves, parameter scan
DESCRIPTION: Fix-round figures for finding F-001 (inlet definition): reservoir-width and contour-tolerance sensitivity,
             mesh convergence of the pressure-driven/reservoir configuration (P2/P1 vs Mini), crop-band (window) scan.
             Also writes fix-round/fixround_summary.json with every number quoted in the fix-round report.
INPUT FILES: fix-round/FR-A_reservoir, FR-B_convergence, FR-C_geometry, FR-E1_convergence_tol0.25, FR-E2_window_scan_tol0.25 (*.json)
OUTPUT FILES: fix-round/figures/FR-F{1,2,3}.{svg,pdf,png}, fix-round/fixround_summary.json
PARAMETERS: none
DEPENDENCIES: numpy matplotlib; LaTeX (mathptmx) for exact house style
STYLE: House Style Contract v1 (FS-house-style-01), Okabe-Ito colours, width 180 mm = 7.09 in. Same colour/marker per configuration as QC-F*.
PHYSICS NOTES: k_D = mu Q L_s / (W dP); L_s = sample length without reservoirs; W = full transverse width (4096 px); mu = 1.
               "iface" = dP between the two sample faces, flux-weighted pressure probed 1 px inside the reservoirs (excludes the access
               resistance of the reservoirs). "total" = dP between the outer reservoir faces.
ORIGIN: CFD digital-rock project, QC-20261019-v1 fix round, Sonnet 5.5 (declared high effort)
REVISION LOG: 2026-10-21 first version
"""
import glob, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "templates"))
import house_style as hs                                        # noqa: E402
from house_style import plt, panel_label, save                  # noqa: E402
import matplotlib as mpl                                        # noqa: E402

mpl.rcParams.update({"figure.dpi": 300, "xtick.top": False, "ytick.right": False, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8})
FR = HERE.parent / "fix-round"
W = 7.09
CASES = ["R3-H", "R2-V", "R3-V"]
COL = {"R3-H": "#0072B2", "R2-V": "#E69F00", "R3-V": "#009E73"}
MK = {"R3-H": "o", "R2-V": "s", "R3-V": "^"}
MKW = dict(mec="black", mew=0.5, ms=4.5)
ACC = "#D55E00"


def tex(s): return s if hs.USE_TEX else s.replace("\\%", "%")
def lab(q, u): return hs.unit(q, u)
def load(d): return [json.load(open(f)) for f in sorted(glob.glob(str(FR / d / "*.json")))]


A = [r for r in load("FR-A_reservoir") if r["status"] == "ok"]
C = [r for r in load("FR-C_geometry") if r["status"] == "ok"]
E1 = load("FR-E1_convergence_tol0.25")
E2 = load("FR-E2_window_scan_tol0.25")
E2_ALL = E2


def pick(rows, **kw):
    out = [r for r in rows if all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v for k, v in kw.items())]
    return out


summary = {}

# ---------------------------------------------------------------- FR-F1: reservoir width and contour tolerance
fig, (a, b) = plt.subplots(1, 2, figsize=(W, 2.9), gridspec_kw=dict(left=0.085, right=0.985, bottom=0.17, top=0.96, wspace=0.28))
for c in CASES:
    rows = sorted(pick(A, case=c), key=lambda r: r["pad_frac"])
    ref = [r for r in rows if abs(r["pad_frac"] - 0.4) < 1e-9][0]["kD_iface_flux"]
    p = [r["pad_frac"] * 100 for r in rows if r["pad_frac"] > 0]
    tot = [(r["kD_flux"] / ref - 1) * 100 for r in rows if r["pad_frac"] > 0]
    ifc = [(r["kD_iface_flux"] / ref - 1) * 100 for r in rows if r["pad_frac"] > 0]
    a.plot(p, ifc, "-", color=COL[c], marker=MK[c], mfc=COL[c], label=c, **MKW)
    a.plot(p, tot, "--", color=COL[c], marker=MK[c], mfc="white", **MKW)
    r0 = [r for r in rows if r["pad_frac"] == 0][0]
    a.plot([0.9], [(r0["kD_flux"] / ref - 1) * 100], marker=MK[c], mfc="white", color=COL[c], ls="none", **MKW)
    summary.setdefault("reservoir", {})[c] = dict(ref_iface_b40=ref, pad0_total=r0["kD_flux"],
        dev_total_pct={f"{r['pad_frac']:g}": (r["kD_flux"] / ref - 1) * 100 for r in rows},
        dev_iface_pct={f"{r['pad_frac']:g}": (r["kD_iface_flux"] / ref - 1) * 100 for r in rows if r["pad_frac"] > 0})
a.axhline(0, color="0.4", lw=0.5)
a.set_xscale("log"); a.set_xticks([0.9, 2.5, 5, 10, 20, 40]); a.set_xticklabels(["0", "2.5", "5", "10", "20", "40"])
a.xaxis.set_minor_locator(mpl.ticker.NullLocator())
a.set_xlabel(lab(r"Reservoir width $b/L_s$", r"\%")); a.set_ylabel(lab(r"$k_D/k_{D,\mathrm{ref}}-1$", r"\%".replace("\\%", "%") if not hs.USE_TEX else r"\%"))
a.set_ylim(-12, 11)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], color="0.2", ls="-", marker="o", mfc="0.2", label=tex("flux-weighted face pressure")),
                  Line2D([], [], color="0.2", ls="--", marker="o", mfc="white", label=tex("outer-face pressure")),
                  *[Line2D([], [], color=COL[c], marker=MK[c], ls="none", mfc=COL[c], mec="black", mew=0.5, ms=4.5, label=c) for c in CASES]],
         frameon=False, loc="upper right", handlelength=2.0, fontsize=7.5)
panel_label(a, "a", x=-0.16)

tols = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0]
for c in CASES:
    ys = []
    for t in tols:
        if abs(t - 1.5) < 1e-9:
            r = pick(A, case=c, pad_frac=0.1)[0]
        else:
            r = pick(C, case=c, tol=t, band=0.08)[0]
        ys.append(r["kD_iface_flux"])
    ref = ys[1]
    b.plot(tols, [(y / ref - 1) * 100 for y in ys], "-", color=COL[c], marker=MK[c], mfc=COL[c], label=c, **MKW)
    summary.setdefault("tolerance", {})[c] = dict(k_iface=dict(zip(map(str, tols), ys)), dev_pct={str(t): (y / ref - 1) * 100 for t, y in zip(tols, ys)})
b.axhline(0, color="0.4", lw=0.5)
b.axvline(1.5, color=ACC, lw=0.8, ls=":")
b.text(1.45, 0.55, tex("manuscript (1.5 px)"), color=ACC, fontsize=7.5, va="center", ha="right", rotation=90, transform=b.get_xaxis_transform())
b.set_xscale("log"); b.set_xticks([0.1, 0.25, 0.5, 1, 2]); b.set_xticklabels(["0.1", "0.25", "0.5", "1", "2"])
b.xaxis.set_minor_locator(mpl.ticker.NullLocator())
b.set_xlabel(lab(r"Contour tolerance", "px")); b.set_ylabel(lab(r"$k_D/k_D(0.25\,\mathrm{px})-1$", r"\%"))
panel_label(b, "b", x=-0.16)
save(fig, "FR-F1", FR / "figures")

# ---------------------------------------------------------------- FR-F2: mesh convergence
fig, axs = plt.subplots(1, 3, figsize=(W, 2.75), gridspec_kw=dict(left=0.075, right=0.99, bottom=0.18, top=0.93, wspace=0.34))
for ax, c, L in zip(axs, CASES, "abc"):
    p2 = sorted(pick(E1, case=c, element="p2"), key=lambda r: r["level"]); mi = sorted(pick(E1, case=c, element="mini"), key=lambda r: r["level"])
    ref = p2[-1]["kD_iface_flux"]
    ax.plot([r["n_tri"] for r in p2], [(r["kD_iface_flux"] / ref - 1) * 100 for r in p2], "-", color=COL[c], marker=MK[c], mfc=COL[c], label="P2/P1", **MKW)
    ax.plot([r["n_tri"] for r in mi], [(r["kD_iface_flux"] / ref - 1) * 100 for r in mi], "--", color=COL[c], marker=MK[c], mfc="white", label="Mini", **MKW)
    ax.axhline(0, color="0.4", lw=0.5)
    ax.set_xscale("log"); ax.set_xticks([6e4, 1.5e5, 4e5]); ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xticklabels([tex(r"$6\times10^{4}$"), tex(r"$1.5\times10^{5}$"), tex(r"$4\times10^{5}$")])
    ax.set_xlabel(lab(r"$N_{\mathrm{tri}}$", "-")); ax.set_title(c, fontsize=9, pad=3)
    if L == "a": ax.set_ylabel(lab(r"$k_D/k_D^{\mathrm{P2,L5}}-1$", r"\%"))
    ax.legend(frameon=False, loc="lower right", fontsize=8, handlelength=2.0)
    panel_label(ax, L, x=-0.2)
    summary.setdefault("convergence", {})[c] = dict(
        p2_iface={str(r["level"]): r["kD_iface_flux"] for r in p2}, mini_iface={str(r["level"]): r["kD_iface_flux"] for r in mi},
        p2_total={str(r["level"]): r["kD_flux"] for r in p2}, n_tri={str(r["level"]): r["n_tri"] for r in p2},
        p2_step_4to5_pct=(p2[4]["kD_iface_flux"] / p2[3]["kD_iface_flux"] - 1) * 100,
        p2_L1_vs_L5_pct=(p2[0]["kD_iface_flux"] / ref - 1) * 100,
        mini_L5_vs_p2_L5_pct=(mi[4]["kD_iface_flux"] / ref - 1) * 100,
        mini_step_4to5_pct=(mi[4]["kD_iface_flux"] / mi[3]["kD_iface_flux"] - 1) * 100,
        k_final_iface=ref, k_final_total=p2[4]["kD_flux"], max_balance=max(abs(r["Q_balance"]) for r in p2 + mi),
        max_residual=max(r["residual_lu"] for r in p2 + mi))
save(fig, "FR-F2", FR / "figures")

# ---------------------------------------------------------------- FR-F3: window (crop band) scan
fig, axs = plt.subplots(1, 3, figsize=(W, 2.75), gridspec_kw=dict(left=0.075, right=0.99, bottom=0.18, top=0.93, wspace=0.34))
for ax, c, L in zip(axs, CASES, "abc"):
    rows = [r for r in E2 if r["case"] == c and r["status"] == "ok"]
    rows += [r for r in pick(E1, case=c, element="p2", level=2)]   # band 0.08 from the convergence series (same settings)
    rows = sorted(rows, key=lambda r: r["band"])
    bands = [r["band"] * 100 for r in rows]; ks = [r["kD_iface_flux"] for r in rows]
    ax.plot(bands, ks, "-", color=COL[c], marker=MK[c], mfc=COL[c], **MKW)
    nop = sorted(r["band"] * 100 for r in E2 if r["case"] == c and r["status"] != "ok")
    if nop:
        ax.plot(nop, [min(ks) * 0.98] * len(nop), "x", color="0.3", ms=5, mew=1.0, ls="none")
        ax.annotate(tex("no\npercolation"), (np.mean(nop), min(ks) * 0.98), xytext=(0, 7), textcoords="offset points", fontsize=7,
                    ha="center", va="bottom", color="0.3", linespacing=0.95)
    ax.axvline(8, color=ACC, lw=0.8, ls=":")
    ax.set_xlabel(lab(r"Crop band", r"\%")); ax.set_title(c, fontsize=9, pad=3)
    if L == "a": ax.set_ylabel(lab(r"$k_D$", r"px$^2$"))
    ax.set_xlim(1, 17)
    panel_label(ax, L, x=-0.2)
    summary.setdefault("window", {})[c] = dict(band_pct=bands, k_iface=ks, not_percolating_band_pct=nop,
                                               min=min(ks), max=max(ks), max_over_min=max(ks) / min(ks))
axs[2].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v/1000:g}k"))
axs[2].set_ylim(27000, 31000)
save(fig, "FR-F3", FR / "figures")

# ---------------------------------------------------------------- headline numbers (before/after)
import qcdata as D   # noqa: E402  (scripts/ is on sys.path via HERE)
head = {}
for c in CASES:
    s = D.summary(c, 4)
    head[c] = dict(author_k_L5=s["apparent_permeability"], fix_k_iface=summary["convergence"][c]["k_final_iface"],
                   fix_k_total=summary["convergence"][c]["k_final_total"])
head["anisotropy"] = dict(author=head["R3-V"]["author_k_L5"] / head["R3-H"]["author_k_L5"],
                          fix_iface=head["R3-V"]["fix_k_iface"] / head["R3-H"]["fix_k_iface"],
                          fix_total=head["R3-V"]["fix_k_total"] / head["R3-H"]["fix_k_total"])
summary["headline"] = head
json.dump(summary, open(FR / "fixround_summary.json", "w"), indent=1)
print(json.dumps(head, indent=1))
