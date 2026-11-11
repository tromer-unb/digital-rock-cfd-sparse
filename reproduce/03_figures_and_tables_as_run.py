#!/usr/bin/env python3
"""
Step 3 of the reproduction: regenerate the remaining figures and the Supplementary Tables of the paper by running the
figure scripts exactly as they were run for the paper, on the data stored in this repository, and compare the results with
the files of the paper in paper/figures and paper/si_tables.

What it regenerates (file names as in paper/figures; figure numbers of the paper in brackets):
  figure1.pdf  [Fig. 1]  percolating domains             reproduce/scripts_as_run/figures/make_fig1_fig2_split.py
  figure3.pdf  [Fig. 2]  definitional sensitivities      reproduce/scripts_as_run/figures/make_fig3_fig4_sens.py
  figure4.pdf  [Fig. 5]  cropping window, representation reproduce/scripts_as_run/figures/make_fig3_fig4_sens.py
  figure5.pdf  [Fig. 3]  verification and validation     reproduce/scripts_as_run/figures/make_fig5_vv.py
  figure6.pdf  [Fig. 6]  solver cost                     reproduce/scripts_as_run/figures/make_fig6_cost.py
  figure7.pdf  [Fig. 7]  accuracy and cost (sub-windows) reproduce/scripts_as_run/figures/make_fig6_cost.py
  figureS1.pdf [Suppl. Fig. 1] Gmsh size field           reproduce/scripts_as_run/figures/make_figures_npj.py S2
  figureS2.pdf [Suppl. Fig. 2] legacy vs reference       reproduce/scripts_as_run/studies/qc_scripts/make_figures_qc.py --fig E3
  si_tables/*.tex  Supplementary Tables 2, 3, 6-8         reproduce/scripts_as_run/figures/make_si_tables.py
(figure2.pdf and figure8-11.pdf are built from the stored results by reproduce/02_figures_from_results.py.)

How: the scripts are kept verbatim and refer to absolute paths of the original project layout
(.../CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1 and .../CFD/PROJETO/convergencia). This runner builds that layout in a
temporary directory from links into data/ (and writes the level-5 domain masks from their compressed archive, checking
the SHA-256 of every reconstructed .npy against the original), copies the scripts there, rewrites ONLY that absolute
path prefix, runs them and compares the outputs with paper/figures and paper/si_tables:
  * tables: byte identity;
  * figures: byte identity of the PDF if the PDF writer is deterministic, otherwise identity of the page rasterized at
    200 dpi (pdftoppm), reported as the number of differing pixels.

Requirements: the PercoFlow environment (numpy, matplotlib 3.11, meshio, pillow), a LaTeX installation with mathptmx
(the figures use Times through usetex) and poppler-utils (pdftoppm).

Usage:  python reproduce/03_figures_and_tables_as_run.py [--keep DIR]
Exit status 0 if every item matches.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
ORIGINAL_PREFIX = "/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO"   # as written in the scripts
PATCHES = {"R1": "patch_y3800_x3800", "R2": "patch_y7600_x19000", "R3": "patch_y7600_x53200"}
LEVELS = ["01_coarse", "02_intermediate", "03_fine", "04_very_fine", "05_extra_fine"]

# (script relative to reproduce/, arguments, outputs to compare: {produced file: file in paper/figures})
JOBS = [
    ("scripts_as_run/figures/make_fig1_fig2_split.py", ["--out", "{FIG}"], {"figure1.pdf": "figure1.pdf"}),
    ("scripts_as_run/figures/make_fig3_fig4_sens.py", ["--out", "{FIG}"], {"figure3.pdf": "figure3.pdf", "figure4.pdf": "figure4.pdf"}),
    ("scripts_as_run/figures/make_fig5_vv.py", ["--out", "{FIG}"], {"figure5.pdf": "figure5.pdf"}),
    ("scripts_as_run/figures/make_fig6_cost.py", ["--out", "{FIG}"], {"figure6.pdf": "figure6.pdf", "figure7.pdf": "figure7.pdf"}),
    ("scripts_as_run/figures/make_figures_npj.py", ["S2"], {"figureS1.pdf": "figureS1.pdf"}),
    ("scripts_as_run/studies/qc_scripts/make_figures_qc.py", ["--fig", "E3", "--out", "{FIG}"], {"QC-E3.pdf": "figureS2.pdf"}),
    ("scripts_as_run/figures/make_si_tables.py", [], {}),
]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def link(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(src, dst)


def build_layout(root: Path) -> Path:
    """Recreate the parts of the original project tree that the scripts read. Returns the PROJETO directory."""
    proj = root / "PROJETO"
    qc = proj / "QUALITY-CHECK" / "QC-20261019-v1"
    # scripts imported by path (copied, not linked: qcdata.py locates the data through Path(__file__).resolve())
    shutil.copytree(REPO / "reproduce/scripts_as_run/studies/qc_scripts", qc / "scripts")
    link(DATA / "results/audit_evidence/indep_sens", qc / "evidence/indep_sens")
    link(DATA / "results/audit_evidence/poiseuille", qc / "evidence/poiseuille")
    for name in ("claims_recomputed.json", "critical_throat.json", "indep_summary.json"):
        link(DATA / "results/audit_evidence" / name, qc / "evidence" / name)
    link(DATA / "results/fix_round/fixround_summary.json", qc / "fix-round/fixround_summary.json")
    link(DATA / "results/fix_round/FR-E1_convergence_tol0.25", qc / "fix-round/FR-E1_convergence_tol0.25")
    link(DATA / "results/benchmarks", qc / "package/drc_audit_pack/results")
    for lev in LEVELS:                      # legacy pipeline on the Poiseuille channels (Supplementary Fig. 2)
        link(DATA / f"results/legacy_poiseuille_rerun/v_{lev}/global_summary.csv",
             qc / f"rerun/poiseuille/v_{lev}/global_summary.csv")
    # the authors' original runs (legacy pipeline): summaries, level-5 percolating domains, images
    conv = proj / "convergencia"
    for lev in LEVELS:
        link(DATA / f"results/legacy_convergence/convergence_mesh_{lev}/global_summary.csv",
             conv / f"convergence_mesh_{lev}/global_summary.csv")
    meta = json.loads((DATA / "results/legacy_convergence/level5_full_percolating_domains.json").read_text())
    arch = np.load(DATA / "results/legacy_convergence/level5_full_percolating_domains.npz")
    for rel, m in meta.items():
        key = "__".join(Path(rel).parts[1:3])
        out = conv / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        np.save(out, arch[key])
        got = sha256(out)
        if got != m["sha256"]:
            raise SystemExit(f"reconstructed {rel}: SHA-256 {got} != original {m['sha256']}")
    for p in PATCHES.values():
        link(DATA / f"images/{p}_c0_mask.png", conv / f"github/structures/{p}_c0_mask.png")
    return proj


def stage_scripts(root: Path, proj: Path) -> Path:
    """Copy the figure scripts next to their original siblings and rewrite only the absolute path prefix."""
    rev3 = root / "rev3"
    shutil.copytree(REPO / "reproduce/scripts_as_run/figures", rev3 / "scripts")
    for py in (rev3 / "scripts").rglob("*.py"):
        txt = py.read_text()
        if ORIGINAL_PREFIX in txt:
            py.write_text(txt.replace(ORIGINAL_PREFIX, str(proj)))
    return rev3


def raster_diff(a: Path, b: Path, work: Path) -> int:
    """Number of differing pixels between the first pages of two PDFs rasterized at 200 dpi."""
    from PIL import Image
    imgs = []
    for p, tag in ((a, "a"), (b, "b")):
        stem = work / f"{p.stem}_{tag}"
        subprocess.run(["pdftoppm", "-r", "200", "-f", "1", "-l", "1", "-png", "-singlefile", str(p), str(stem)], check=True)
        imgs.append(np.asarray(Image.open(f"{stem}.png").convert("RGB")))
    if imgs[0].shape != imgs[1].shape:
        return -1
    return int((imgs[0] != imgs[1]).any(-1).sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--keep", help="build in this directory and keep it (default: a temporary directory)")
    a = ap.parse_args()
    root = Path(a.keep).resolve() if a.keep else Path(tempfile.mkdtemp(prefix="rev3_"))
    if a.keep:
        if root.exists():
            shutil.rmtree(root)
        root.mkdir(parents=True)
    proj = build_layout(root)
    rev3 = stage_scripts(root, proj)
    fig = rev3 / "figures"
    fig.mkdir()
    qc_scripts = proj / "QUALITY-CHECK/QC-20261019-v1/scripts"
    report, ok = [], True
    for script, args, outputs in JOBS:
        path = (rev3 / "scripts" / Path(script).name) if script.startswith("scripts_as_run/figures") else qc_scripts / Path(script).name
        cmd = [sys.executable, str(path)] + [s.replace("{FIG}", str(fig)) for s in args]
        r = subprocess.run(cmd, cwd=path.parent, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-4000:])
            report.append(dict(script=script, status="FAILED"))
            ok = False
            continue
        for produced, target in outputs.items():
            p, t = fig / produced, REPO / "paper/figures" / target
            same = sha256(p) == sha256(t)
            px = 0 if same else raster_diff(p, t, root)
            status = "IDENTICAL" if same else ("SAME RASTER" if px == 0 else f"DIFFERS ({px} px)")
            ok &= px == 0
            report.append(dict(script=script, output=produced, paper=target, status=status))
    for t in sorted((REPO / "paper/si_tables").glob("*.tex")):
        p = rev3 / "si_tables" / t.name
        same = p.exists() and p.read_bytes() == t.read_bytes()
        ok &= same
        report.append(dict(script="scripts_as_run/figures/make_si_tables.py", output=f"si_tables/{t.name}",
                           status="IDENTICAL" if same else "DIFFERS"))
    for row in report:
        print(f"{row['status']:<18} {row.get('output', '-'):<28} <- {row['script']}")
    print("PASS" if ok else "FAIL", f"(work directory: {root})")
    if not a.keep and ok:
        shutil.rmtree(root)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
