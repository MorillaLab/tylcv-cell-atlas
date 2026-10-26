# Supplementary Figure 4 | PlantAG-imputed cell types over the brightfield images of the eight capture areas
# Paste into one Colab cell after mounting Google Drive. Writes figures/FigureS4.pdf and figures/FigureS4.png.
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

BASE = Path("/content/drive/MyDrive/spatial_transcriptomics")
ANN = BASE / "processed/adata_combined_squidpy_clean_tda_features_plantag_gnn_embedding.h5ad"
OUT = BASE / "figures"
TYPE_KEY, SAMPLE_KEY = "predicted_cell_type", "sample"
SHOW_IMAGE = True          # False draws the spots alone
SWAP_I7 = False            # True only if the two infected 7-dpi panels show spots off their leaf strips

# Colours of Fig. 1f,g and Fig. 3 (Scanpy default palette, types in alphabetical order)
PALETTE = {"Bundle sheath": "#1f77b4", "Guard": "#ff7f0e", "Hydathodes": "#279e68", "Mesophyll": "#d62728",
           "Mesophyll-Ch": "#aa40fc", "Parenchyma": "#8c564b", "Phloem": "#e377c2", "Vascular bundle": "#b5bd61",
           "Vein": "#17becf", "Xylem": "#aec7e8"}
# Rows: 7 dpi, 24 dpi. Columns: mock 1, mock 2, TYLCV 1, TYLCV 2 (same order as Fig. 3a)
ORDER = [("control_7dpi_1", "M7-1"), ("control_7dpi_2", "M7-2"), ("infected_7dpi_1", "I7-1"), ("infected_7dpi_2", "I7-2"),
         ("control_24dpi_1", "M24-1"), ("control_24dpi_2", "M24-2"), ("infected_24dpi_1", "I24-1"), ("infected_24dpi_2", "I24-2")]

adata = sc.read_h5ad(ANN)
sections = adata.obs[SAMPLE_KEY].astype(str)
types = adata.obs[TYPE_KEY].astype(str)
unknown = sorted(set(types) - set(PALETTE))
if unknown:
    print("types without a colour (drawn grey):", unknown)
spatial = adata.uns.get("spatial", {})
lib_key = next((k for k in ("library_id", "sample_id", "library", SAMPLE_KEY)
                if k in adata.obs and set(adata.obs[k].astype(str)) & set(spatial)), None)


def section_image(sec):
    """Brightfield image, scale from full-resolution pixels to it, spot diameter in full-resolution pixels."""
    lib = adata.obs.loc[sections == sec, lib_key].astype(str).iloc[0] if lib_key else sec
    if SWAP_I7 and lib in ("infected_7dpi_1", "infected_7dpi_2"):
        lib = {"infected_7dpi_1": "infected_7dpi_2", "infected_7dpi_2": "infected_7dpi_1"}[lib]
    d = spatial.get(lib)
    if d is None:
        print(f"{sec}: no image found (libraries: {list(spatial)}); spots drawn alone")
        return None, 1.0, None
    key = "hires" if "hires" in d["images"] else "lowres"
    return np.asarray(d["images"][key]), float(d["scalefactors"][f"tissue_{key}_scalef"]), d["scalefactors"].get("spot_diameter_fullres")


plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 7})
fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.6))
for ax, (sec, label) in zip(axes.ravel(), ORDER):
    m = (sections == sec).to_numpy()
    full = np.asarray(adata.obsm["spatial"], float)[m, :2]
    img, sf, spot_d = section_image(sec) if SHOW_IMAGE else (None, 1.0, None)
    xy = full * sf
    if img is not None:                                       # coordinates are [x, y]; fall back to [y, x] if needed
        inside = lambda p: np.mean((p[:, 0] >= 0) & (p[:, 0] < img.shape[1]) & (p[:, 1] >= 0) & (p[:, 1] < img.shape[0]))
        if inside(xy) < 0.95 and inside(xy[:, ::-1]) >= 0.95:
            xy = xy[:, ::-1]
        ax.imshow(img, cmap="gray" if img.ndim == 2 else None, interpolation="lanczos")
    pitch = float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))   # spot centre-to-centre distance = 100 um
    diam = spot_d * sf if spot_d else 0.55 * pitch                    # spot diameter = 55 um
    cols = types[m].map(PALETTE).fillna("#bbbbbb").to_list()
    ax.add_collection(EllipseCollection(np.full(len(xy), diam), np.full(len(xy), diam), np.zeros(len(xy)), units="xy",
                                        offsets=xy, offset_transform=ax.transData, facecolors=cols, edgecolors="none"))
    pad = 1.5 * pitch
    x0, x1, y0, y1 = xy[:, 0].min() - pad, xy[:, 0].max() + pad, xy[:, 1].min() - pad, xy[:, 1].max() + 2.5 * pad
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title(f"{label} ({int(m.sum())} spots)", fontsize=7, pad=3)
    ax.plot([x0 + pitch, x0 + 11 * pitch], [y1 - 0.9 * pitch] * 2, color="black", lw=1.2, solid_capstyle="butt")  # 1 mm
    ax.text(x0 + 6 * pitch, y1 - 1.4 * pitch, "1 mm", ha="center", va="bottom", fontsize=6)
for row, tp in enumerate(("7 dpi", "24 dpi")):
    fig.text(0.012, 0.73 - row * 0.42, tp, rotation=90, ha="center", va="center", fontsize=8, fontweight="bold")
handles = [Line2D([], [], marker="o", ls="", ms=5, color=c, label=t) for t, c in PALETTE.items()]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=6.5, handletextpad=0.3, columnspacing=1.2)
fig.subplots_adjust(left=0.03, right=0.995, top=0.95, bottom=0.12, wspace=0.05, hspace=0.18)
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "FigureS4.pdf"); fig.savefig(OUT / "FigureS4.png", dpi=300)
plt.show()
print(adata.obs.groupby(SAMPLE_KEY, observed=True)[TYPE_KEY].value_counts().unstack(fill_value=0).loc[[s for s, _ in ORDER]])
