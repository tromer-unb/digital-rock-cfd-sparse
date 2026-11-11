#!/usr/bin/env python3
"""
FIGURE TYPE: sensitivity curves
DESCRIPTION: Splits the former Figure 3 (definitional sensitivities, panels a-e; make_figures_npj.py fig2) in two, at
  the author's request (2026-11-04), in the figure standard of figstd.py.
  Figure 3  former panels (a)-(c) in one column (3x1): inlet condition, reservoir width, contour simplification.
  Figure 4  former panels (d)-(e) in one column (2x1), relabelled (a)-(b): cropping window, image representation.
  Panel titles removed (their content is in the captions); the "not percolating" note of the cropping-window panel no
  longer covers data.
INPUT FILES (same as before, nothing typed by hand):
  QC/evidence/indep_sens/tol{0.5,1.0,1.5}_{R3H,R2V,R3V}_{velocity,pressure}.json, QC/fix-round/fixround_summary.json,
  RES/b7_subwindow_pair_R2-V_tile512.json, RES/b8_pixel_convergence_R2-V_t1_5.json
OUTPUT FILES: <out>/figure3.{svg,pdf,png}, <out>/figure4.{svg,pdf,png}
ORIGIN: CFD digital-rock manuscript, npj Computational Materials draft rev3, 2026-11-04
"""
import argparse, json
from pathlib import Path
import numpy as np
import figstd as S
from figstd import plt, mpl, COL, MK, MKW, GREY, CASES, FS, ref, letter

QC = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1")
RES = QC / "package/drc_audit_pack/results"
FILE = {"R3-H": "R3H", "R2-V": "R2V", "R3-V": "R3V"}
FR = json.load(open(QC / "fix-round/fixround_summary.json"))


def figure3(out):
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(S.WIDTH, 7.6), gridspec_kw=dict(hspace=0.42))
    tols = [0.5, 1.0, 1.5]
    for cs in CASES:
        r = []
        for t in tols:
            v = json.load(open(QC / f"evidence/indep_sens/tol{t}_{FILE[cs]}_velocity.json"))
            p = json.load(open(QC / f"evidence/indep_sens/tol{t}_{FILE[cs]}_pressure.json"))
            r.append(p["k"]["nodal"]["k_author_def"] / v["k"]["nodal"]["k_author_def"])
        a.plot(tols, r, "-", color=COL[cs], marker=MK[cs], label=cs, **MKW)
    ref(a.axhline(1, color=GREY, lw=0.6, ls=":"))
    a.set_xlabel("Contour tolerance [px]"); a.set_ylabel(r"$k_{\mathrm{pressure}}/k_{\mathrm{velocity}}$")
    a.set_yscale("log"); a.set_yticks([1, 2, 5]); a.set_yticklabels(["1", "2", "5"]); a.set_ylim(0.8, 7)
    a.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    a.set_xticks(tols); a.set_xlim(0.4, 1.6)
    a.legend(frameon=False, loc="upper right", handlelength=1.4, borderaxespad=0.3)
    letter(a, "a")

    for cs in CASES:
        dv = FR["reservoir"][cs]["dev_iface_pct"]; bb = sorted(float(x) for x in dv)
        b.plot([100 * x for x in bb], [dv[f"{x:g}"] for x in bb], "-", color=COL[cs], marker=MK[cs], **MKW)
    ref(b.axhline(0, color=GREY, lw=0.6, ls=":")); ref(b.axvline(10, color=GREY, lw=0.6, ls="--"))
    b.set_xlabel(r"Reservoir width $b/L_s$ [\%]"); b.set_ylabel(r"Deviation [\%]")
    letter(b, "b")

    for cs in CASES:
        dv = FR["tolerance"][cs]["dev_pct"]; tt = sorted(float(x) for x in dv)
        c.plot(tt, [dv[f"{x:g}" if f"{x:g}" in dv else str(x)] for x in tt], "-", color=COL[cs], marker=MK[cs], **MKW)
    ref(c.axhline(0, color=GREY, lw=0.6, ls=":")); ref(c.axvline(1.5, color=GREY, lw=0.6, ls="--"))
    c.set_xlabel("Contour tolerance [px]"); c.set_ylabel(r"Deviation [\%]")
    letter(c, "c")
    fig.align_ylabels([a, b, c])
    S.save(fig, "figure3", out)


def figure4(out):
    fig, (d, e) = plt.subplots(2, 1, figsize=(S.WIDTH, 5.3), gridspec_kw=dict(hspace=0.42))
    for cs in CASES:
        w = FR["window"][cs]; k8 = w["k_iface"][w["band_pct"].index(8.0)]
        d.plot(w["band_pct"], np.array(w["k_iface"]) / k8, "-", color=COL[cs], marker=MK[cs], label=cs, **MKW)
        for bp in w["not_percolating_band_pct"]:
            d.plot([bp], [0.48], marker="x", color=COL[cs], ms=7, mew=1.4, ls="none")
    ref(d.axvline(8, color=GREY, lw=0.6, ls="--")); ref(d.axhline(1, color=GREY, lw=0.6, ls=":"))
    d.text(4.7, 0.48, "not percolating", fontsize=FS["note"], color=COL["R2-V"], ha="left", va="center")
    d.set_xticks([2, 4, 6, 8, 10, 12, 14, 16]); d.set_xlim(1, 17)
    d.set_xlabel(r"Boundary band [\% of the image side]"); d.set_ylabel(r"$k_D/k_D(8\%)$")
    d.set_ylim(0.35, 2.0); d.legend(frameon=False, loc="upper left", ncol=1, handlelength=1.4, borderaxespad=0.3)
    letter(d, "a")

    b7 = json.load(open(RES / "b7_subwindow_pair_R2-V_tile512.json"))
    b8 = json.load(open(RES / "b8_pixel_convergence_R2-V_t1_5.json"))
    sound = [t for t in b7["tiles"] if t.get("status") == "ok" and t.get("sound")]
    kfem = np.array([t["fem"]["5"]["kD_face"] for t in sound]); rd = 100 * np.array([t["rel_diff"] for t in sound])
    e.plot(kfem, rd, "s", color=COL["R2-V"], label="11 sub-window pairs", **MKW)
    t15 = [t for t in sound if t["tag"] == b8["tag"]][0]
    e.plot([t15["fem"]["5"]["kD_face"]], [100 * b8["richardson_vs_fem"]], marker="*", color="black", ms=11, ls="none",
           label=r"one pair, cell size $\rightarrow 0$")
    ref(e.axhline(0, color=GREY, lw=0.6, ls=":"))
    e.set_xscale("log"); e.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    e.set_xlabel(r"$k_D$ of the sub-window, contour FE [px$^2$]")
    e.set_ylabel(r"$k_{\mathrm{pixel}}/k_{\mathrm{contour}}-1$ [\%]")
    e.set_ylim(-36, 6)
    e.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.22), handlelength=1.4, borderaxespad=0.3)
    letter(e, "b")
    fig.align_ylabels([d, e])
    print("sub-windows", len(sound), "median |diff| %.1f%%" % np.median(np.abs(rd)))
    S.save(fig, "figure4", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    figure3(a.out); figure4(a.out)
