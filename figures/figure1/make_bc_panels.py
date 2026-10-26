"""Panels b-c of Figure 1: persistence diagrams of the median TYLCV and median mock spot (24 dpi), vector, in the style of the figure."""
import matplotlib
import os
from pathlib import Path as _P
_HERE = _P(__file__).resolve().parent
_B = _P(os.environ.get("BUNDLE_DIR", _HERE.parents[1] / "results" / "bundle"))
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d = pd.read_csv(_B / "04b_h1_panels_data.csv")
dg = d[d.part == "diagram"].copy(); dg["dim"] = dg["dim"].astype(int)
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.major.size": 2.5, "ytick.major.size": 2.5})
H0 = (0.161, 0.502, 0.725); H1 = {"TYLCV": (0.875, 0.361, 0.365), "mock": (0.494, 0.494, 0.494)}
finite = dg[np.isfinite(dg.y)]
TOP = 1.08 * finite.y.max()                       # shared scale for both diagrams
PAGE = {"TYLCV": dict(x0=428.0, x1=570.0, axis_x=454.2), "mock": dict(x0=575.0, x1=706.0, axis_x=591.7)}
Y0, Y1, AX_TOP, AX_BOT, AX_W = 24.5, 206.5, 41.0, 185.0, 108.9
TITLE = {"TYLCV": "Persistence diagram, TYLCV\nmedian spot (24 dpi)", "mock": "Persistence diagram, mock\nmedian spot (24 dpi)"}
for grp, g in PAGE.items():
    W, H = g["x1"] - g["x0"], Y1 - Y0
    fig = plt.figure(figsize=(W / 72, H / 72))
    ax = fig.add_axes([(g["axis_x"] - g["x0"]) / W, (Y1 - AX_BOT) / H, AX_W / W, (AX_BOT - AX_TOP) / H])
    s = dg[dg.group == grp]
    h0 = s[(s.dim == 0) & np.isfinite(s.y)]; h1 = s[s.dim == 1]
    ax.plot([0, TOP], [0, TOP], ls=(0, (3, 2)), lw=0.8, color="#9E9E9E", zorder=1)
    ax.scatter(h0.x, h0.y, s=3.46 ** 2, color=H0, lw=0, zorder=2, label="H$_0$")
    ax.scatter(h1.x, h1.y, s=3.46 ** 2, color=H1[grp], lw=0, zorder=3, label="H$_1$")
    ax.set_xlim(-0.03 * TOP, TOP); ax.set_ylim(-0.03 * TOP, TOP)
    ax.set_xticks(range(0, int(TOP) + 1)); ax.set_yticks(range(0, int(TOP) + 1))
    ax.tick_params(labelsize=5.5, pad=1.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("Birth ε", fontsize=6, labelpad=2); ax.set_ylabel("Death ε", fontsize=6, labelpad=2)
    ax.set_title(TITLE[grp], fontsize=6, pad=3.5, linespacing=1.15)
    ax.legend(frameon=False, loc="lower right", fontsize=5, handletextpad=0.1, borderaxespad=0.2, markerscale=0.9, labelspacing=0.3)
    if grp == "TYLCV":
        ax.text(0.04, 0.985, "Scaled total H$_1$ persistence,\nTYLCV vs mock: −3.5%, P = 0.79\n(8 capture areas)",
                transform=ax.transAxes, ha="left", va="top", fontsize=4.8, style="italic", color="#333333", linespacing=1.2)
    fig.savefig(str(_HERE / "panels" / f"panel_{grp}.pdf"), transparent=True)
    plt.close(fig)
print("panels written; shared axis top =", round(TOP, 2))
