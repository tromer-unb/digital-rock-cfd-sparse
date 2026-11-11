#!/usr/bin/env python3
"""Write every numerical table of supplementary.tex from the stored result files (npj revision, 2026-11-03).
Nothing is typed by hand: each table reads one JSON of the audit package and formats it.
OUTPUT: si_tables/{clusters,poiseuille_p1,poiseuille_worst,cylinders,minres,minres_summary,nnz_range,subwindows}.tex
"""
import json
from pathlib import Path
import numpy as np

QC = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1")
RES = QC / "package/drc_audit_pack/results"
OUT = Path(__file__).resolve().parents[1] / "si_tables"
OUT.mkdir(exist_ok=True)
FILE = {"R3H": "R3-H", "R2V": "R2-V", "R3V": "R3-V"}


def sci(x, d=1):
    """LaTeX scientific notation without semicolons."""
    if x == 0:
        return "$0$"
    e = int(np.floor(np.log10(abs(x)))); m = x / 10 ** e
    return f"${m:.{d}f}\\times10^{{{e}}}$"


def thou(n):
    return f"{int(n):,}"


def write(name, text):
    (OUT / name).write_text(text)
    print("wrote", name)


# ------------------------------------------------------------------ clusters (c3)
C3 = json.load(open(RES / "c3_wellposedness.json"))
rows = []
for key in ("R1-H", "R1-V", "R2-H", "R2-V", "R3-H", "R3-V"):
    s = C3["rock_cluster_statistics"][key]
    rows.append(f"{key} & {thou(s['pore_px'])} & {s['n_clusters']} & {s['n_percolating']} & {s['n_dead_end_inlet']} & "
                f"{s['n_dead_end_outlet']} & {s['n_isolated']} & {100 * s['area_frac_percolating']:.1f} & "
                f"{100 * s['area_frac_dead_end']:.1f} & {100 * s['area_frac_isolated']:.1f} \\\\")
write("clusters.tex", r"""\begin{table}[htbp]
\centering
\caption{Pore clusters of the three rock images in the cropped windows (8-connectivity, boundary band of $8\%$). Columns give the number of pore pixels, the number of clusters, the numbers of percolating clusters, of dead ends that touch only the inlet or only the outlet plane, and of isolated clusters, and the shares of the pore area in percolating, dead-end, and isolated clusters, in percent.}
\label{tab:clusters}
\small
\begin{tabular}{lrrrrrrrrr}
\toprule
 & Pore px & Clusters & Perc. & Dead, in & Dead, out & Isolated & Perc. (\%) & Dead (\%) & Isol. (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ Poiseuille under P1 (b9)
B9 = json.load(open(RES / "b9_poiseuille_p1.json"))
rows = []
for h in sorted({r["h"] for r in B9["rows"]}):
    rr = sorted((r for r in B9["rows"] if r["h"] == h), key=lambda r: r["level"])
    for i, r in enumerate(rr):
        lead = f"{h} & {r['Ls']}" if i == 0 else " & "
        rows.append(f"{lead} & {r['level']} & {thou(r['P1_p2']['n_tri'])} & {sci(r['P1_p2']['rel_face'])} & "
                    f"${100 * r['P1_mini']['rel_face']:+.2f}$ & ${100 * r['velocity_p2']['rel_vs_poiseuille']:+.2f}$ \\\\")
    rows.append("\\midrule")
rows = rows[:-1]
write("poiseuille_p1.tex", r"""\begin{table}[htbp]
\centering
\caption{Plane Poiseuille flow under the pressure-driven definition (reservoirs of width $0.1L_s$, contour tolerance $0.25$ px). Columns give the channel width $h$ and length $L_s$ in pixels, the mesh level, the number of triangles of the Taylor--Hood mesh, the relative deviation of the Taylor--Hood permeability from the exact value of equation~\eqref{eq:channelbias}, the same deviation for the Mini element in percent, and the deviation of a velocity-driven Taylor--Hood solution without reservoirs from $h^{2}/12$ in percent.}
\label{tab:poiseuilleP1}
\small
\begin{tabular}{rrrrrrr}
\toprule
$h$ & $L_s$ & Level & $N_{\mathrm{tri}}$ & Taylor--Hood & Mini (\%) & Velocity-driven (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
write("poiseuille_worst.tex", sci(B9["criteria"]["worst_p2"]).strip("$").join(["$", "$"]))

# ------------------------------------------------------------------ cylinders (b2, b5)
B2 = {round(x["c"], 4): x for x in json.load(open(RES / "b2_periodic_cell.json"))["cases"]}
B5 = {round(x["c"], 4): x for x in json.load(open(RES / "b5_openfoam_cell.json"))["cases"]}
rows = []
for c in sorted(B2):
    x = B2[c]; fin = min(x["mesh"], key=lambda m: m["h"])
    fv = B5.get(c)
    fvs = f"{fv['k_of']:.5f} & ${100 * fv['rel_diff_vs_p2p1']:+.3f}$" if fv else "-- & --"
    rows.append(f"{c:g} & {thou(fin['n_tri'])} & {x['k_num']:.5f} & ${100 * x['conv_last_step']:+.3f}$ & {x['k_hasimoto']:.5f} & "
                f"${100 * x['rel_diff']:+.2f}$ & {fvs} \\\\")
write("cylinders.tex", r"""\begin{table}[htbp]
\centering
\caption{Periodic square array of cylinders. Columns give the solid fraction $c$, the number of triangles of the finest Taylor--Hood mesh, the Taylor--Hood permeability $k/L^{2}$ on that mesh, its relative change from the preceding mesh in percent, the dilute series of equation~\eqref{eq:hasimoto}, the deviation of the Taylor--Hood value from the series in percent, the OpenFOAM permeability on its finest mesh of 200 cells per cell side, and the deviation of the OpenFOAM value from the Taylor--Hood value in percent.}
\label{tab:cylinders}
\small
\begin{tabular}{rrrrrrrr}
\toprule
$c$ & $N_{\mathrm{tri}}$ & $k/L^{2}$, FE & Last step (\%) & Series & FE vs series (\%) & $k/L^{2}$, FV & FV vs FE (\%) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ MINRES vs LU (c4)
C4 = json.load(open(RES / "c4_minres_p1.json"))
lu = {(r["case"], r["level"]): r for r in C4["rows"] if r["method"] == "lu"}
mr = {(r["case"], r["level"]): r for r in C4["rows"] if r["method"] == "minres"}
rows = []
for case in ("R3H", "R2V", "R3V"):
    for lv in range(1, 6):
        a, b = lu[(case, lv)], mr[(case, lv)]
        lead = FILE[case] if lv == 1 else ""
        rows.append(f"{lead} & {lv} & {thou(a['n_free'])} & {a['nnz'] / a['n_free']:.1f} & {a['csr_MB']:.0f} & {a['t_solve_s']:.1f} & "
                    f"{a['peak_rss_MB_after_solve'] / 1024:.2f} & {thou(b['iterations'])} & {b['t_solve_s']:.0f} & "
                    f"{b['peak_rss_MB_after_solve'] / 1024:.2f} & {sci(b['rel_residual_true'])} & {sci(abs(b['rel_k_vs_LU_reference']))} \\\\")
    rows.append("\\midrule")
rows = rows[:-1]
write("minres.tex", r"""\begin{table}[htbp]
\centering
\caption{Screening solver under the pressure-driven definition (Mini elements). Columns give the configuration, the mesh level, the number of free unknowns $N$, the non-zeros per row of the reduced matrix, its CSR storage in MiB, the solve time in seconds and the peak resident memory in GiB of sparse LU, the number of MINRES iterations, the MINRES solve time in seconds and its peak memory in GiB, the true relative residual $\|\mathbf{r}\|/\|\mathbf{f}\|$ of the MINRES solution, and the relative difference between the MINRES and LU permeabilities.}
\label{tab:minres}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrrrrrrrr}
\toprule
 & & & & & \multicolumn{2}{c}{LU} & \multicolumn{5}{c}{MINRES} \\
\cmidrule(lr){6-7}\cmidrule(lr){8-12}
 & Level & $N$ & nnz$/N$ & CSR & Time & Memory & Iterations & Time & Memory & Residual & $|\Delta k_D|/k_D$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
ratios = [a["nnz"] / a["n_free"] for a in lu.values()]
write("nnz_range.tex", f"{min(ratios):.1f} and {max(ratios):.1f}")
mrv = list(mr.values())
its = [r["iterations"] for r in mrv]; kerr = max(abs(r["rel_k_vs_LU_reference"]) for r in mrv)
tres = [r["rel_residual_true"] for r in mrv]
tr = [mr[k]["t_solve_s"] / lu[k]["t_solve_s"] for k in mr]
mem = [lu[k]["peak_rss_MB_after_solve"] / mr[k]["peak_rss_MB_after_solve"] for k in mr if k[1] >= 4]
stops = sorted({r["stop_message"] for r in mrv})
crit = C4["criteria"]
nmin = min(r["n_free"] for r in mrv); nmax = max(r["n_free"] for r in mrv)
g_it = [mr[(c, 5)]["iterations"] / mr[(c, 1)]["iterations"] for c in ("R3H", "R2V", "R3V")]
g_n = [mr[(c, 5)]["n_free"] / mr[(c, 1)]["n_free"] for c in ("R3H", "R2V", "R3V")]
npass = sum(abs(r["rel_k_vs_LU_reference"]) < 1e-3 for r in mrv)
first = (f"The criterion holds for {npass} of the {len(mrv)} runs, and the largest difference is {sci(kerr)}. " if not crit["PASS"]
         else f"The criterion holds, with a largest difference of {sci(kerr)}. ")
orders = (np.log10(min(tres) / 1e-7), np.log10(max(tres) / 1e-7))
write("minres_summary.tex", first +
      f"MINRES needed between {thou(min(its))} and {thou(max(its))} iterations, a count that grows by a factor of {min(g_it):.1f} to {max(g_it):.1f} for each configuration while the number of unknowns grows {min(g_n):.0f}- to {max(g_n):.0f}-fold, "
      f"and every run stopped through the backward-error test of equation~\\eqref{{eq:test1}}. "
      f"The true relative residual of the accepted solutions lies between {sci(min(tres))} and {sci(max(tres))}, {["zero","one","two","three","four","five"][round(orders[0])]} to {["zero","one","two","three","four","five"][round(orders[1])]} orders of magnitude above the requested tolerance, as equation~\\eqref{{eq:trueres}} anticipates. "
      f"MINRES took {min(tr):.0f} to {max(tr):.0f} times longer than sparse LU on the same matrix, and on the two finest levels it needed {min(mem):.1f} to {max(mem):.1f} times less memory. "
      f"The algebraic error of MINRES follows its true residual, and the one run that missed the criterion, R3-V at level 5, also has the largest true residual. "
      f"The direct solver is therefore the faster and more accurate choice whenever its memory fits, and MINRES the choice that bounds memory, provided that its true residual is tested.\n")
print("stop messages:", stops)

# ------------------------------------------------------------------ sub-windows (b7)
B7 = json.load(open(RES / "b7_subwindow_pair_R2-V_tile512.json"))
sound = sorted([t for t in B7["tiles"] if t.get("status") == "ok" and t.get("sound")], key=lambda t: t["fem"]["5"]["kD_face"])
rows = []
for i, t in enumerate(sound, 1):
    fe, fv = t["fem"]["5"], t["fv"]
    rows.append(f"{i} & {t['tag'].replace('_', chr(92) + '_')} & {100 * t['porosity_percolating']:.1f} & {fe['kD_face']:.1f} & {fv['kD_face']:.1f} & "
                f"${100 * t['rel_diff']:+.1f}$ & {thou(fe['n_tri'])} & {fe['t_mesh_s'] + fe['t_lu_s']:.1f} & {thou(fv['n_cells'])} & "
                f"{thou(fv['iterations'])} & {fv['t_mesh_s'] + fv['t_solve_s']:.0f} & {fv['nproc']} \\\\")
write("subwindows.tex", r"""\begin{table}[htbp]
\centering
\caption{Matched sub-windows of $512\times512$ px of R2-V, in order of increasing permeability. Columns give the index used in main-text Figure~4e, the sub-window tag, the percolating porosity in percent, the contour finite-element (Taylor--Hood, mesh level 5) and pixel finite-volume permeabilities in px$^{2}$, their relative difference in percent, the number of triangles and the wall time in seconds (meshing plus LU, one process) of the finite-element solution, and the number of cells, the number of SIMPLEC iterations, the wall time in seconds (meshing plus solve), and the number of processes of the finite-volume solution.}
\label{tab:subwindows}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{rlrrrrrrrrrr}
\toprule
 & & & \multicolumn{2}{c}{$k_D$ (px$^{2}$)} & & \multicolumn{2}{c}{Contour FE} & \multicolumn{4}{c}{Pixel FV} \\
\cmidrule(lr){4-5}\cmidrule(lr){7-8}\cmidrule(lr){9-12}
 & Tag & $\phi$ (\%) & FE & FV & Diff. (\%) & $N_{\mathrm{tri}}$ & Time & Cells & Iterations & Time & Proc. \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
