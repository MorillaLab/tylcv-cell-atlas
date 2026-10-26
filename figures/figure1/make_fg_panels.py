"""Panels f-g of Figure 1: UMAP of the PlantAG latent embedding, coloured by training labels (f) and imputed types (g)."""
import matplotlib
import os
from pathlib import Path as _P
_HERE = _P(__file__).resolve().parent
_B = _P(os.environ.get("BUNDLE_DIR", _HERE.parents[1] / "results" / "bundle"))
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

u = pd.read_csv(_B / "08c_plantag_umap.csv")
SCANPY20 = ["#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8", "#ffbb78"]
types = sorted(set(u.pseudo_cell_type) | set(u.predicted_cell_type))
col = dict(zip(types, SCANPY20))
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
RECTS = {"f": (127.6, 118.8, "pseudo_cell_type"), "g": (122.9, 115.7, "predicted_cell_type")}
order = u.sample(frac=1.0, random_state=0).index                      # random draw order: no type hidden under another
for k, (w, h, key) in RECTS.items():
    fig = plt.figure(figsize=(w / 72, h / 72)); ax = fig.add_axes([0.02, 0.02, 0.96, 0.96]); ax.set_axis_off()
    ax.scatter(u.loc[order, "umap1"], u.loc[order, "umap2"], s=1.1, c=u.loc[order, key].map(col), lw=0)
    ax.set_aspect("equal", adjustable="datalim")
    fig.savefig(str(_HERE / "panels" / f"panel_{k}.pdf"), transparent=True); plt.close(fig)
fig = plt.figure(figsize=(78.8 / 72, 131.5 / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
for i, t in enumerate(types):
    y = 1 - (i + 0.8) / (len(types) + 0.6)
    ax.scatter([0.07], [y], s=14, color=col[t], transform=ax.transAxes, clip_on=False)
    ax.text(0.15, y, t, transform=ax.transAxes, va="center", fontsize=6.6)
fig.savefig(str(_HERE / "panels" / "panel_legend.pdf"), transparent=True); plt.close(fig)
print("f, g and legend written for", len(types), "types:", types)
