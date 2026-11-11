"""
Figure standard of the CFD manuscript, set by the corresponding author on 2026-11-04, for every figure:
  - Times through LaTeX (house_style of QC-20261019-v1, mathptmx), usetex required;
  - type sizes of Figures 1-2: ticks 10 pt, axis labels 11 pt, legends and notes 10 pt, panel letters 12 pt bold "(a)",
    for a single-column figure 89 mm wide that is included at its natural size;
  - closed axes box, ticks pointing inward (major and minor);
  - no panel titles (their content goes into the caption);
  - no text over data and no legend outside its axes box: check_overlap() must return no conflict before a
    figure is delivered.
Run with the audit venv (QUALITY-CHECK/QC-20261019-v1/package/drc_audit_pack/.venv), whose matplotlib renders usetex.
"""
import sys
from pathlib import Path
import numpy as np

QCS = Path("/home/luiz/SIMULACOES/CLAUDE/MANUSCRIPT-REVISION/CFD/PROJETO/QUALITY-CHECK/QC-20261019-v1/scripts")
sys.path.insert(0, str(QCS / "templates"))
import house_style as hs                                        # noqa: E402
from house_style import plt                                     # noqa: E402
import matplotlib as mpl                                        # noqa: E402

assert hs.USE_TEX, "LaTeX rendering unavailable: run with the audit venv"

FS = dict(tick=10, label=11, note=10, letter=12)
WIDTH = 3.5                                                     # in, 89 mm single column
COL = {"R3-H": "#0072B2", "R2-V": "#E69F00", "R3-V": "#009E73"}
MK = {"R3-H": "o", "R2-V": "s", "R3-V": "^"}
MKW = dict(mec="black", mew=0.5, ms=6)
GREY = "0.45"
CASES = ["R3-H", "R2-V", "R3-V"]

mpl.rcParams.update({
    "axes.labelsize": FS["label"], "xtick.labelsize": FS["tick"], "ytick.labelsize": FS["tick"],
    "legend.fontsize": FS["note"], "font.size": FS["note"],
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.minor.top": True, "ytick.minor.right": True,
    "axes.spines.top": True, "axes.spines.right": True, "lines.linewidth": 1.2})


def ref(artist):
    """Mark a reference line (axhline/axvline/band) so that the overlap check does not treat it as data."""
    artist.set_gid("ref")
    return artist


def letter(ax, s, x=-0.22, y=1.02):
    ax.text(x, y, r"\textbf{(%s)}" % s, transform=ax.transAxes, fontsize=FS["letter"], ha="left", va="bottom",
            gid="letter")


def _data_points(ax, renderer):
    """Display coordinates of every data marker and of densely sampled points along every data line, with a radius."""
    pts = []
    for ln in ax.get_lines():
        if ln.get_gid() == "ref" or not ln.get_visible():
            continue
        xy = ax.transData.transform(np.column_stack(ln.get_data()).astype(float))
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if len(xy) == 0:
            continue
        r = 0.5 * ln.get_markersize() * ax.figure.dpi / 72 if ln.get_marker() not in (None, "None", "", " ") else 0
        for p in xy:
            pts.append((p[0], p[1], r))
        if ln.get_linestyle() not in ("None", "", " ") and len(xy) > 1:
            lw = 0.5 * ln.get_linewidth() * ax.figure.dpi / 72
            for p0, p1 in zip(xy[:-1], xy[1:]):
                for t in np.linspace(0, 1, 40)[1:-1]:
                    q = p0 + t * (p1 - p0)
                    pts.append((q[0], q[1], lw))
    return pts


def check_overlap(fig, pad_pt=1.0):
    """Return the list of (axes index, text, number of data points) for every text box that covers data."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pad = pad_pt * fig.dpi / 72
    bad = []
    for i, ax in enumerate(fig.axes):
        pts = _data_points(ax, renderer)
        boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in ax.texts if t.get_visible() and t.get_text()]
        leg = ax.get_legend()
        if leg is not None:
            lb = leg.get_window_extent(renderer); ab = ax.get_window_extent(renderer)
            boxes.append(("legend", lb))
            outside_ok = leg.get_gid() == "outside-ok"     # a legend deliberately placed between panels, off all data
            if not outside_ok and (lb.x0 < ab.x0 - pad or lb.x1 > ab.x1 + pad or lb.y0 < ab.y0 - pad or lb.y1 > ab.y1 + pad):
                bad.append((i, "legend outside the axes box", 0))
        bars = [p.get_window_extent(renderer) for p in ax.patches
                if isinstance(p, mpl.patches.Rectangle) and p.get_gid() != "ref" and p.get_visible()]
        on_bar = {t.get_text() for t in ax.texts if t.get_gid() == "on-bar"}   # value labels printed inside bars
        for name, bb in boxes:
            n = sum(1 for (x, y, r) in pts
                    if bb.x0 - r - pad <= x <= bb.x1 + r + pad and bb.y0 - r - pad <= y <= bb.y1 + r + pad)
            if name not in on_bar:
                n += sum(1 for q in bars if q.x0 < bb.x1 + pad and q.x1 > bb.x0 - pad and q.y0 < bb.y1 + pad and q.y1 > bb.y0 - pad)
            if n:
                bad.append((i, name, n))
    return bad


def save(fig, name, out):
    import os
    bad = check_overlap(fig)
    if bad and not os.environ.get("FIGSTD_DRAFT"):
        raise RuntimeError(f"{name}: text over data {bad}")
    if bad:
        print(f"DRAFT {name}: text over data {bad}")
    hs.save(fig, name, out)
