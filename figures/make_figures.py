import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore")
import os
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
B = Path(os.environ.get("BUNDLE_DIR", ROOT / "results" / "bundle"))       # unzipped result bundle (05 for the article)
OUT = Path(os.environ.get("FIG_OUT", HERE / "output")); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300, "pdf.fonttype": 42})
MOCK, INF, GREY = "#0072B2", "#D55E00", "#9E9E9E"
D7, D24 = "#56B4E9", "#E69F00"
ORDER = ["control_7dpi_1", "control_7dpi_2", "infected_7dpi_1", "infected_7dpi_2",
         "control_24dpi_1", "control_24dpi_2", "infected_24dpi_1", "infected_24dpi_2"]
SHORT = {s: ("M" if s.startswith("control") else "I") + s.split("_")[1].replace("dpi", "") + "-" + s[-1] for s in ORDER}
col = lambda s: MOCK if s.startswith("control") else INF
stats = {}


def letter(ax, l, x=-0.16, y=1.06):
    ax.text(x, y, l, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def img(ax, path, crop=(0, 1, 0, 1)):
    im = plt.imread(path)
    h, w = im.shape[:2]
    ax.imshow(im[int(crop[0] * h):int(crop[1] * h), int(crop[2] * w):int(crop[3] * w)], interpolation="lanczos")
    ax.set_axis_off()


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)



# ------------------------------------------------------------------ vector leaf maps (Fig. 2a, Supplementary Fig. 3)
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.colors import to_rgb
LEAF_COL = {"L1": "#0072B2", "L2": "#E69F00", "L3": "#009E73", "fragment": "#BBBBBB"}


def leaf_spots():
    """Analysed spots with array coordinates (um) and the leaf assigned to each spot in the analysis."""
    h = pd.read_csv(B / "04b_h1_persistence_per_spot.csv", index_col=0)
    u = pd.read_csv(B / "08c_plantag_umap.csv", index_col=0)
    assert (h["sample"].astype(str).values == u["sample"].astype(str).values).all()
    lab = pd.DataFrame({"spot": u["spot"].values, "leaf": h["leaf"].astype(str).values})
    m = pd.read_csv(B / "11_viral_per_spot.csv", index_col=0).merge(lab, on="spot", how="left")
    m["leaf_short"] = m["leaf"].str.extract(r"_(L\d+|fragment)$")[0].fillna("fragment")
    m["x"] = m["array_col"] * 50.0
    m["y"] = m["array_row"] * 50.0 * np.sqrt(3)
    return m


def pieces_of(d, radius=101.0):
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    xy = d[["x", "y"]].to_numpy()
    return connected_components(cKDTree(xy).sparse_distance_matrix(cKDTree(xy), radius).tocsr(), directed=False)[1]


def draw_spots(ax, d, colors, title):
    ax.add_collection(EllipseCollection(np.full(len(d), 80.0), np.full(len(d), 80.0), np.zeros(len(d)), units="xy",
                                        offsets=d[["x", "y"]].to_numpy(), offset_transform=ax.transData,
                                        facecolors=colors, edgecolors="none"))
    pad = 140
    x0, x1, y0, y1 = d.x.min() - pad, d.x.max() + pad, d.y.min() - pad, d.y.max() + pad
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_aspect("equal"); ax.set_axis_off()
    ax.plot([x0 + 60, x0 + 1060], [y1 - 60] * 2, color="black", lw=1.0, solid_capstyle="butt")
    ax.text(x0 + 560, y1 - 110, "1 mm", ha="center", va="bottom", fontsize=5)
    ax.set_title(title, fontsize=6.5, pad=2)


def leaf_panels(fig, spec, mode="leaves", letter_tag=None):
    m = leaf_spots(); seg = pd.read_csv(B / "03_leaf_segmentation.csv").set_index("sample")
    sub = GridSpecFromSubplotSpec(2, 4, subplot_spec=spec, hspace=0.18, wspace=0.06)
    for i, s_ in enumerate(ORDER):
        ax = fig.add_subplot(sub[i // 4, i % 4]); d = m[m["sample"] == s_].reset_index(drop=True)
        if mode == "leaves":
            cols = [LEAF_COL[l] for l in d.leaf_short]
        else:   # each tissue piece in its own shade of its leaf's colour
            pc = pieces_of(d); cols = [None] * len(d)
            for lf, g in d.groupby("leaf_short"):
                base = np.array(to_rgb(LEAF_COL[lf]))
                ps = pd.Series(pc[g.index]).value_counts().index          # pieces of this leaf, largest first
                for k, pid in enumerate(ps):
                    shade = base if lf == "fragment" else base + (1 - base) * [0.0, 0.42, 0.62, 0.25, 0.52, 0.72][k % 6]
                    for j in g.index[pc[g.index] == pid]:
                        cols[j] = tuple(np.clip(shade, 0, 1))
        draw_spots(ax, d, cols, f"{SHORT[s_]}  ({int(seg.loc[s_, 'leaves'])} leaves)")
        if i == 0 and letter_tag:
            letter(ax, letter_tag, x=-0.08, y=1.12)
    return m



# Domain of each analysed spot, recovered from 08_domains.png at the exact pixel position of every spot and verified to
# reproduce 08_domain_counts_per_section.csv for all capture areas and domains (no mismatching spot).
DOM_COL = {"0": "#332288", "1": "#DDCC77", "2": "#44AA99", "3": "#CC6677", "4": "#88CCEE", "5": "#AA4499"}   # Tol muted


def domain_panels(fig, spec, letter_tag=None):
    d_all = pd.read_csv(HERE / "data" / "spot_domains.csv", index_col=0, dtype={"domain": str})
    d_all["x"] = d_all["array_col"] * 50.0; d_all["y"] = d_all["array_row"] * 50.0 * np.sqrt(3)
    sub = GridSpecFromSubplotSpec(2, 4, subplot_spec=spec, hspace=0.18, wspace=0.06)
    for i, s_ in enumerate(ORDER):
        ax = fig.add_subplot(sub[i // 4, i % 4]); d = d_all[d_all["sample"] == s_].reset_index(drop=True)
        draw_spots(ax, d, [DOM_COL[k] for k in d["domain"]], SHORT[s_])
        if i == 0 and letter_tag:
            letter(ax, letter_tag, x=-0.08, y=1.12)

# ------------------------------------------------------------------ Figure 2: dataset and data quality
qc = pd.read_csv(B / "01_sample_qc.csv").set_index("sample").loc[ORDER]
fig = plt.figure(figsize=(7.2, 6.6))
gs = GridSpec(2, 3, figure=fig, height_ratios=[1.9, 1], hspace=0.42, wspace=0.55)
leaf_panels(fig, gs[0, :], "leaves", letter_tag="a")
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=l) for l, c in (("leaf 1", LEAF_COL["L1"]), ("leaf 2", LEAF_COL["L2"]),
           ("leaf 3", LEAF_COL["L3"]), ("fragment, not assigned", LEAF_COL["fragment"]))], loc="upper center", bbox_to_anchor=(0.55, 0.955), ncol=4, frameon=False, fontsize=6)
x = np.r_[0:4, 5:9]
def per_section_axis(ax, ylabel, title, tag):
    ax.set_xticks(x); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=60)
    ax.set_ylabel(ylabel); ax.set_title(title); letter(ax, tag)
ax = fig.add_subplot(gs[1, 0])
for xi, s in zip(x, ORDER):
    r = qc.loc[s]
    ax.vlines(xi, r.q25_umis, r.q75_umis, color=col(s), lw=2.2, alpha=0.55)
    ax.plot(xi, r.median_umis, "o", color=col(s), ms=4.5, mec="white", mew=0.5)
ax.set_yscale("log"); ax.set_yticks([20, 50, 100, 200, 500]); ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
per_section_axis(ax, "UMIs per spot (median, IQR)", "Capture depth", "b")
ax = fig.add_subplot(gs[1, 1])
ax.bar(x, qc.median_genes, color=[col(s) for s in ORDER], width=0.75)
per_section_axis(ax, "Genes per spot (median)", "Genes detected", "c")
ax = fig.add_subplot(gs[1, 2])
ax.bar(x, qc.median_pct_mito, color=[col(s) for s in ORDER], width=0.75)
per_section_axis(ax, "Mitochondrial UMIs (%, median)", "Mitochondrial UMIs", "d")
save(fig, "Figure2")

# ------------------------------------------------------------------ Figure 3: PlantAG imputation
fb = pd.read_csv(B / "06_feature_block_share.csv", index_col=0)
ko = pd.read_csv(B / "06_knn_overlap_with_full_graph.csv", index_col=0)
raw = pd.read_csv(B / "07_random_vs_loso_summary.csv", header=None)
hdr0, hdr1 = raw.iloc[0, 1:].tolist(), raw.iloc[1, 1:].tolist()
perf = raw.iloc[3:].set_index(0)
perf.columns = pd.MultiIndex.from_arrays([[h.split(".")[0] for h in hdr0], hdr1])
perf = perf.astype(float)
nul = pd.read_csv(B / "04_label_share_vs_null.csv", index_col=0)
auc = pd.read_csv(B / "04_heldout_marker_auc.csv", index_col=0)
aud = pd.read_csv(B / "04_marker_dictionary_audit.csv").set_index("cell_type")
shp = pd.read_csv(B / "05_predicted_cell_type_share_per_section.csv", index_col=0).loc[ORDER]
pc = pd.read_csv(B / "07_spatial_coherence_of_labels.csv", header=[0, 1], index_col=0)["predicted_cell_type"].astype(float)
TYPES = list(pc.sort_values("n_spots", ascending=False).index)
SCANPY20 = ["#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8",
            "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31"]
TCOL = dict(zip(sorted(TYPES), SCANPY20))        # Scanpy default palette, types in alphabetical order (as in Fig. 1f,g)

fig = plt.figure(figsize=(7.2, 11.0))
gs = GridSpec(4, 3, figure=fig, hspace=0.62, wspace=0.95, height_ratios=[1, 1, 1.1, 0.95])
# a | composition per capture area
ax = fig.add_subplot(gs[0, 0:2])
xpos = np.r_[0:4, 5:9]; bottom = np.zeros(len(ORDER))
for t in TYPES:
    vals = shp[t].to_numpy() if t in shp else np.zeros(len(ORDER))
    ax.bar(xpos, vals, bottom=bottom, color=TCOL[t], width=0.75, lw=0)
    bottom += vals
ax.set_xticks(xpos); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=60)
ax.set_ylabel("Share of spots"); ax.set_ylim(0, 1.0)
ax.text(1.5, -0.3, "7 dpi slide", ha="center", transform=ax.get_xaxis_transform(), fontsize=6.3)
ax.text(6.5, -0.3, "24 dpi slide", ha="center", transform=ax.get_xaxis_transform(), fontsize=6.3)
ax.set_title("Imputed cell types per capture area (colours in b)"); letter(ax, "a", x=-0.09)
# b | spatial coherence per type
ax = fig.add_subplot(gs[0, 2])
pcs = pc.sort_values("ratio_obs_null"); yy = np.arange(len(pcs))
ax.hlines(yy, 1, pcs["ratio_obs_null"], color="#BBBBBB", lw=0.8)
ax.scatter(pcs["ratio_obs_null"], yy, s=22, c=[TCOL[t] for t in pcs.index], zorder=3, edgecolors="black", linewidths=0.3)
ax.axvline(1, ls="--", lw=0.6, color="black")
ax.set_yticks(yy); ax.set_yticklabels([f"{t} ({100 * pc.loc[t, 'n_spots'] / pc['n_spots'].sum():.0f}%)" for t in pcs.index], fontsize=5.6)
ax.set_xlabel("Neighbour purity,\nobserved / shuffled"); ax.set_xlim(0.8, 3.7)
ax.set_title("Spatial coherence"); letter(ax, "b", x=-1.05)
# c | generalisation
ax = fig.add_subplot(gs[1, 0:2])
short_m = {"PlantAG-like GCN | feature graph | original features": "PlantAG (topological graph)",
           "LogReg | TDA-image PCA only (features used to make the labels)": "LogReg, persistence-image PCs",
           "GCN | expression kNN graph | expression": "GCN, expression graph",
           "GCN | spatial graph | expression": "GCN, spatial graph",
           "SGC | spatial graph | expression": "SGC, spatial graph",
           "LogReg | expression PCA": "LogReg, expression PCs"}
order_m = list(short_m); yy = np.arange(len(order_m))[::-1]
ra = perf.loc[order_m, ("random_80_20", "balanced_accuracy")].to_numpy()
lo = perf.loc[order_m, ("leave_one_section_out", "balanced_accuracy")].to_numpy()
ax.barh(yy + 0.2, ra, height=0.38, color="#BBBBBB", label="random 80/20 split of spots")
ax.barh(yy - 0.2, lo, height=0.38, color="#332288", label="whole capture area held out")
for yi, b_ in zip(yy, lo):
    ax.text(b_ + 0.01, yi - 0.2, f"{b_:.2f}", va="center", fontsize=5.5)
NT = len(TYPES)
ax.axvline(1 / NT, ls="--", lw=0.7, color="black"); ax.text(1 / NT + 0.01, -0.95, f"chance (1/{NT})", fontsize=5.5)
ax.set_yticks(yy); ax.set_yticklabels([short_m[m] for m in order_m])
ax.set_xlabel("Balanced accuracy for the imputed types"); ax.set_xlim(0, 1.05); ax.set_ylim(-1.1, len(order_m) - 0.4)
ax.legend(frameon=False, loc="lower right", fontsize=5.8)
ari = json.load(open(B / "07_seed_stability.json"))["seed_ARI_mean"]
ax.set_title(f"Imputation generalises to held-out capture areas (seed ARI {ari:.2f})"); letter(ax, "c", x=-0.3)
# d | what the graph is built from
ax = fig.add_subplot(gs[1, 2])
nd = fb["n_dims"].astype(int)
names = {"expression_PCA": f"Expr.\nPCs\n({nd['expression_PCA']})", "persistence_image_pixels": f"Persist.\nimage\n({nd['persistence_image_pixels']:,})",
         "tda_scalar_features": f"TDA\nscalars\n({nd['tda_scalar_features']})"}
xs = np.arange(3)
ax.bar(xs - 0.2, fb.loc[list(names), "share_after_standardisation"] * 100, width=0.4, color="#CC79A7", label="share of input")
ax.bar(xs + 0.2, ko.loc[list(names), "jaccard_with_full_PlantAG_graph"] * 100, width=0.4, color="#009E73", label="graph overlap (Jaccard ×100)")
ax.set_xticks(xs); ax.set_xticklabels(list(names.values()), fontsize=5.8)
ax.set_ylabel("%"); ax.legend(frameon=False, loc="upper left", fontsize=5.2, handlelength=1.0); ax.set_ylim(0, 140); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_title("PlantAG input and graph"); letter(ax, "d", x=-0.45)
# e | marker coverage
ax = fig.add_subplot(gs[2, 0])
cts = [c for c in aud.index if c in nul.index]
listed = aud.loc[cts, "n_unique"]; indata = aud.loc[cts, "n_in_data"]; det = nul.loc[cts, "n_markers_detected"]
o = np.argsort(listed.to_numpy()); yy = np.arange(len(cts))
ax.barh(yy, listed.to_numpy()[o], color="#DDDDDD", label="listed")
ax.barh(yy, indata.to_numpy()[o], color="#88CCEE", label="in matrix")
ax.barh(yy, det.to_numpy()[o], color="#117733", label="detected")
ax.set_xscale("log"); ax.set_yticks(yy); ax.set_yticklabels(np.array(cts)[o], fontsize=5.6)
ax.set_xlabel("Marker genes (log)"); ax.legend(frameon=False, fontsize=5.2, loc="lower right")
ax.set_title("Literature markers at\nspot resolution"); letter(ax, "e", x=-0.6)
# f | split-half consistency of the marker lists
ax = fig.add_subplot(gs[2, 1])
au = auc.sort_values("auc_heldout_markers_mean")
ax.errorbar(au.auc_heldout_markers_mean, np.arange(len(au)), xerr=au.auc_sd, fmt="o", ms=3.5, color="#882255", ecolor="#AAAAAA", elinewidth=0.6)
ax.axvline(0.5, ls="--", lw=0.7, color="black")
ax.set_yticks(np.arange(len(au))); ax.set_yticklabels(au.index, fontsize=5.6)
ax.set_xlabel("Split-half AUC"); ax.set_xlim(0.4, 0.8); ax.set_title("Internal consistency of\nmarker lists"); letter(ax, "f", x=-0.6)
# g | frequencies vs random-list null
ax = fig.add_subplot(gs[2, 2])
nl = nul.copy()
ax.errorbar(nl.null_share_mean * 100, nl.observed_share * 100, xerr=nl.null_share_sd * 100, fmt="o", ms=3.2,
            color="#332288", ecolor="#AAAAAA", elinewidth=0.6, capsize=0)
lim = max(nl.null_share_mean.max(), nl.observed_share.max()) * 100 * 1.15
ax.plot([0, lim], [0, lim], "--", lw=0.6, color="black")
pr = pearsonr(nl.null_share_mean, nl.observed_share)
ax.text(0.97, 0.04, f"r = {pr[0]:.2f}", transform=ax.transAxes, fontsize=5.8, va="bottom", ha="right")
ax.set_xlabel("Expected with random\nlists of equal size (%)"); ax.set_ylabel("Observed share (%)")
ax.set_title("Type frequencies vs\nrandom-list null"); letter(ax, "g", x=-0.45)
stats["label_share_vs_null_pearson"] = (round(pr[0], 3), float(pr[1]))
stats["label_share_vs_null_spearman"] = round(spearmanr(nl.null_share_mean, nl.observed_share)[0], 3)
# h | infection information in held-out capture areas, capture depth regressed out
inf_s = pd.read_csv(B / "04c_infection_signal_held_out_summary.csv", index_col=0)
FEATS = [("expression PCs (50) | depth removed", "Expression PCs"), ("persistence images | depth removed", "Persistence images"),
         ("TDA scalar descriptors (40) | depth removed", "TDA scalar descriptors"),
         ("H1 loop statistics, scale-free | depth removed", "H$_1$ loop summaries (scaled)"), ("capture depth only", "Capture depth alone")]
ax = fig.add_subplot(gs[3, 0:2])
for k, (tp, c_) in enumerate((("7dpi", "#7570b3"), ("24dpi", "#d95f02"))):
    g_ = inf_s[inf_s.time == tp].set_index("features")
    yy = np.arange(len(FEATS)) + (k - 0.5) * 0.28
    m_ = np.array([g_.loc[f, "mean_auc"] for f, _ in FEATS]); lo_ = np.array([g_.loc[f, "min_auc"] for f, _ in FEATS])
    hi_ = np.array([g_.loc[f, "max_auc"] for f, _ in FEATS])
    ax.errorbar(m_, yy, xerr=[m_ - lo_, hi_ - m_], fmt="o", ms=3.5, color=c_, capsize=1.5, lw=0.8, label=tp.replace("dpi", " dpi"))
ax.axvline(0.5, ls="--", lw=0.6, color="black")
ax.set_yticks(range(len(FEATS))); ax.set_yticklabels([l for _, l in FEATS]); ax.invert_yaxis(); ax.set_xlim(0.3, 1.0)
ax.set_xlabel("ROC AUC for infection in held-out capture areas\n(mean and range over four held-out pairs; capture depth regressed out)")
ax.legend(frameon=False, fontsize=6, loc="lower right")
ax.set_title("Local topology carries infection information"); letter(ax, "h", x=-0.42)
stats["infection_signal"] = inf_s.round(3).to_dict("records")
save(fig, "Figure3")
stats["perf"] = {short_m[m]: {"random_bal": round(float(a), 3), "loso_bal": round(float(b_), 3)} for m, a, b_ in zip(order_m, ra, lo)}
stats["auc_range_non_mesophyll_ch"] = (float(auc.drop("Mesophyll-Ch").auc_heldout_markers_mean.min()), float(auc.drop("Mesophyll-Ch").auc_heldout_markers_mean.max()))
stats["auc_mesophyll_ch"] = float(auc.loc["Mesophyll-Ch", "auc_heldout_markers_mean"])

# ------------------------------------------------------------------ Figure 3: domains, slide structure, 'disruption'
dc = pd.read_csv(B / "08_domain_counts_per_section.csv").set_index("sample").loc[ORDER]
dcoh = pd.read_csv(B / "08_domain_spatial_coherence.csv", index_col=0)
lcoh = pd.read_csv(B / "07_spatial_coherence_of_labels.csv", header=None)
lcoh = lcoh.iloc[2:, :6].set_index(0).astype(float); lcoh.columns = ["n", "purity", "null", "ratio", "z"]
om = pd.read_csv(B / "10_original_metric_vs_within_condition.csv", index_col=0)
dd = pd.read_csv(B / "10_domain_disruption_with_nulls.csv", index_col=0)

fig = plt.figure(figsize=(7.2, 8.4))
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.9, 1, 1], hspace=0.6, wspace=0.55)
domain_panels(fig, gs[0, :], letter_tag="a")
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=f"domain {k}") for k, c in DOM_COL.items()],
           loc="upper center", bbox_to_anchor=(0.55, 0.955), ncol=6, frameon=False, fontsize=6)
ax = fig.add_subplot(gs[1, 0:2])
pal = [DOM_COL[k] for k in sorted(DOM_COL)]   # domains 0-5, same colours as the maps in panel a
fr = dc.div(dc.sum(1), axis=0)
bottom = np.zeros(len(fr))
for k, c in zip(fr.columns, pal):
    ax.bar(np.r_[0:4, 5:9], fr[c].to_numpy() if False else fr[k].to_numpy(), bottom=bottom, color=c, width=0.75, label=f"domain {k}")
    bottom += fr[k].to_numpy()
ax.set_xticks(np.r_[0:4, 5:9]); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=60)
ax.set_ylabel("Share of spots"); ax.set_ylim(0, 1.24); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.legend(frameon=False, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=5.6, columnspacing=0.8, handlelength=1.0, handletextpad=0.3)
ax.text(1.5, -0.42, "7 dpi slide", ha="center", transform=ax.get_xaxis_transform(), fontsize=6.5)
ax.text(6.5, -0.42, "24 dpi slide", ha="center", transform=ax.get_xaxis_transform(), fontsize=6.5)
letter(ax, "b", x=-0.1)
d3 = {t: float(fr.loc[[s for s in ORDER if t in s], "3"].mean()) for t in ("_7dpi", "_24dpi")}
d2 = {t: float(fr.loc[[s for s in ORDER if t in s], "2"].mean()) for t in ("_7dpi", "_24dpi")}
stats["domain3_share"] = d3; stats["domain2_share"] = d2
ax = fig.add_subplot(gs[1, 2])
ax.scatter(np.zeros(len(lcoh)) + np.random.default_rng(1).uniform(-0.12, 0.12, len(lcoh)), lcoh.ratio, s=12, color="#882255", label="PlantAG labels")
ax.scatter(np.ones(len(dcoh)) + np.random.default_rng(2).uniform(-0.12, 0.12, len(dcoh)), dcoh.ratio_obs_null, s=12, color="#117733", label="expression domains")
ax.axhline(1, ls="--", lw=0.6, color="black")
ax.set_xticks([0, 1]); ax.set_xticklabels(["PlantAG\nlabels", "Expression\ndomains"]); ax.set_xlim(-0.5, 1.5)
ax.set_yscale("log"); ax.set_ylabel("Neighbour purity / shuffled"); ax.set_title("Spatial coherence"); letter(ax, "c", x=-0.5)
ax = fig.add_subplot(gs[2, 0:2])
within = om[["cosine_mock1_vs_mock2", "cosine_inf1_vs_inf2"]].max(1)
for t, c, m in (("7dpi", D7, "o"), ("24dpi", D24, "s")):
    sub = om[om.time == t]
    ax.scatter(within[sub.index], sub.cosine_infected_vs_mock, s=16, color=c, marker=m, label=t, edgecolor="white", lw=0.4)
for i, r in om.iterrows():
    if r.cosine_infected_vs_mock > within[i]:
        ax.annotate(f"{r.label} {r.time}", (within[i], r.cosine_infected_vs_mock), fontsize=5, xytext=(2, 2), textcoords="offset points")
lim = max(within.max(), om.cosine_infected_vs_mock.max()) * 1.1
ax.plot([0, lim], [0, lim], "--", lw=0.6, color="black"); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("Largest within-condition distance (mock vs mock, infected vs infected)")
ax.set_ylabel("Infected vs mock distance")
ax.legend(frameon=False, loc="upper left"); ax.set_title("Infected–mock distance of each type vs replicate distance", loc="left")
ax.text(0.98, 0.04, f"{int((om.cosine_infected_vs_mock <= within).sum())} of {int(within.notna().sum())} within replicate range", transform=ax.transAxes, ha="right", fontsize=6)
letter(ax, "d", x=-0.1)
n_above = int((om.cosine_infected_vs_mock > within).sum())
stats["disruption_above_within"] = (n_above, len(om))
ax = fig.add_subplot(gs[2, 2])
dco = dd[dd.metric == "correlation"].dropna(subset=["within_condition_section_mean"])
yy = np.arange(len(dco))[::-1]
for yi, (_, r) in zip(yy, dco.iterrows()):
    ax.plot([r.observed, r.within_condition_section_mean], [yi, yi], color="#BBBBBB", lw=1)
    ax.plot(r.observed, yi, "o", color=INF, ms=4); ax.plot(r.within_condition_section_mean, yi, "o", color=GREY, ms=4)
ax.set_yticks(yy); ax.set_yticklabels([f"domain {int(r.group)}, {r.time}" for _, r in dco.iterrows()], fontsize=5.8)
ax.set_xlabel("1 − Pearson r of mean\nleaf profiles")
ax.plot([], [], "o", color=INF, label="infected vs mock"); ax.plot([], [], "o", color=GREY, label="between sections,\nsame condition")
ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.4, -0.38), ncol=1); ax.set_title("Domain-level distances"); letter(ax, "e", x=-0.55)
save(fig, "Figure4")

# ------------------------------------------------------------------ Figure 4: replicate-aware host response
L7 = pd.read_csv(B / "09b_limma_dupcor_7dpi.csv", index_col=0); L24 = pd.read_csv(B / "09b_limma_dupcor_24dpi.csv", index_col=0)
D7t = pd.read_csv(B / "09_DE_7dpi_all_genes_top300.csv", index_col=0); D24t = pd.read_csv(B / "09_DE_24dpi_all_genes_top300.csv", index_col=0)
ls_ = pd.read_csv(B / "09b_limma_dupcor_summary.csv").set_index("timepoint")

fig = plt.figure(figsize=(7.2, 10.4))
gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.55, height_ratios=[1, 1.35, 1.3])


def volcano(ax, L, title, nlab=()):
    sig = L["adj.P.Val"] < 0.05
    ax.scatter(L.logFC[~sig], -np.log10(L["P.Value"][~sig]), s=6, color="#CCCCCC", lw=0)
    ax.scatter(L.logFC[sig & (L.logFC > 0)], -np.log10(L["P.Value"][sig & (L.logFC > 0)]), s=10, color=INF, lw=0)
    ax.scatter(L.logFC[sig & (L.logFC < 0)], -np.log10(L["P.Value"][sig & (L.logFC < 0)]), s=10, color=MOCK, lw=0)
    for g, (dx, dy) in nlab.items():
        if g not in L.index:
            continue
        r = L.loc[g]
        ax.annotate(g, (r.logFC, -np.log10(r["P.Value"])), fontsize=5.3, xytext=(dx, dy), textcoords="offset points",
                    ha="right" if dx < 0 else "left", arrowprops=dict(arrowstyle="-", lw=0.3, color="#555555"))
    ax.axvline(0, lw=0.5, color="black")
    ax.set_xlabel("log2 fold change (infected/mock)"); ax.set_ylabel("−log10 P")
    ax.set_title(title)


volcano(fig.add_subplot(gs[0, 0]), L7, f"7 dpi: {int((L7['adj.P.Val']<0.05).sum())}/{len(L7)} genes at FDR < 0.05", nlab={"Cab-1D": (8, 2), "Cab-1A": (10, -4), "IPS": (-6, -9), "SlFBA2": (6, 9), "LyesC2p084": (8, 4)})
letter(fig.axes[-1], "a", x=-0.3)
volcano(fig.add_subplot(gs[0, 1]), L24, f"24 dpi: {int((L24['adj.P.Val']<0.05).sum())}/{len(L24)} genes at FDR < 0.05", nlab={"CAB11": (-6, 4), "NR": (-8, 3), "CAB1B": (8, 2), "ltpg1": (8, -5), "SAMDC": (-6, -6), "TLRP": (-4, 9)})
letter(fig.axes[-1], "b", x=-0.3)
ax = fig.add_subplot(gs[0, 2])
rhos = {}
for L, Dt, c, t, m in ((L7, D7t, D7, "7 dpi", "o"), (L24, D24t, D24, "24 dpi", "s")):
    cm = L.index.intersection(Dt.index)
    rhos[t] = spearmanr(L.loc[cm, "logFC"], Dt.loc[cm, "log2FoldChange_leaf"])[0]
    ax.scatter(Dt.loc[cm, "log2FoldChange_leaf"], L.loc[cm, "logFC"], s=8, color=c, marker=m, lw=0, label=f"{t} (ρ = {rhos[t]:.2f}, n = {len(cm)})")
ax.plot([-4, 4], [-4, 4], "--", lw=0.6, color="black"); ax.axhline(0, lw=0.4, color="#888888"); ax.axvline(0, lw=0.4, color="#888888")
ax.set_xlabel("DESeq2, leaf pseudobulks (log2 FC)"); ax.set_ylabel("limma-voom + within-section\ncorrelation (log2 FC)")
ax.legend(frameon=False, fontsize=5.5, loc="upper left"); ax.set_title("Effect sizes agree across models"); letter(ax, "c", x=-0.32)
stats["limma_vs_deseq2_rho"] = {k: round(float(v), 3) for k, v in rhos.items()}
ax = fig.add_subplot(gs[1, 0])
cnt = pd.DataFrame({"DESeq2\nleaf": [int((D7t.padj_leaf < 0.05).sum()), int((D24t.padj_leaf < 0.05).sum())],
                    "DESeq2 leaf\n+ section": [int(D7t.robust.sum()), int(D24t.robust.sum())],
                    "limma\nnested": [int(ls_.loc["7dpi", "fdr05_up"] + ls_.loc["7dpi", "fdr05_down"]), int(ls_.loc["24dpi", "fdr05_up"] + ls_.loc["24dpi", "fdr05_down"])]},
                   index=["7 dpi", "24 dpi"])
xs = np.arange(3)
ax.bar(xs - 0.2, cnt.loc["7 dpi"], width=0.4, color=D7, label="7 dpi"); ax.bar(xs + 0.2, cnt.loc["24 dpi"], width=0.4, color=D24, label="24 dpi")
for xi, a, b_ in zip(xs, cnt.loc["7 dpi"], cnt.loc["24 dpi"]):
    ax.text(xi - 0.2, a + 0.6, str(a), ha="center", fontsize=5.8); ax.text(xi + 0.2, b_ + 0.6, str(b_), ha="center", fontsize=5.8)
ax.set_xticks(xs); ax.set_xticklabels(cnt.columns, fontsize=5.8); ax.set_ylabel("Differentially expressed genes")
ax.legend(frameon=False, fontsize=5.8); ax.set_title("Genes at adjusted P < 0.05"); letter(ax, "d", x=-0.3)
ax = fig.add_subplot(gs[1, 1:])
genes = ["Cab-1A", "Cab-1D", "CAB5", "CAB1B", "CAB11", "CAB13", "psbS", "LyesC2p085", "LyesC2p083", "LyesC2p084", "LyesC2p009",
         "LyesC2t025", "LOC112940383", "IPS", "NR", "TLRP", "dhn", "PAL3", "CHS2", "SAMDC", "ltpg1", "LOC544001", "SlFBA2"]
colsH = [("limma", "7 dpi"), ("DESeq2 leaf", "7 dpi"), ("DESeq2 section", "7 dpi"), ("limma", "24 dpi"), ("DESeq2 leaf", "24 dpi"), ("DESeq2 section", "24 dpi")]
V = np.full((len(genes), 6), np.nan); S = np.zeros((len(genes), 6), bool)
for i, g in enumerate(genes):
    for j, (meth, t) in enumerate(colsH):
        L, Dt = (L7, D7t) if t == "7 dpi" else (L24, D24t)
        if meth == "limma" and g in L.index:
            V[i, j] = L.loc[g, "logFC"]; S[i, j] = L.loc[g, "adj.P.Val"] < 0.05
        elif meth == "DESeq2 leaf" and g in Dt.index:
            V[i, j] = Dt.loc[g, "log2FoldChange_leaf"]; S[i, j] = Dt.loc[g, "padj_leaf"] < 0.05
        elif meth == "DESeq2 section" and g in Dt.index:
            V[i, j] = Dt.loc[g, "log2FoldChange_section"]; S[i, j] = Dt.loc[g, "padj_section"] < 0.05
imh = ax.imshow(np.ma.masked_invalid(V), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
ax.set_facecolor("#F2F2F2")
for i in range(len(genes)):
    for j in range(6):
        if np.isfinite(V[i, j]):
            ax.text(j, i, ("*" if S[i, j] else "") , ha="center", va="center", fontsize=8, color="black")
ax.set_xticks(range(6)); ax.set_xticklabels([f"{m}\n{t}" for m, t in colsH], fontsize=5.6)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=5.8)
ax.axvline(2.5, color="black", lw=0.8)
cb = fig.colorbar(imh, ax=ax, fraction=0.03, pad=0.02); cb.set_label("log2 FC", fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax.set_title("Evidence across models (* adjusted P < 0.05; grey, not tested or not among top 300)"); letter(ax, "e", x=-0.18)
E_t = pd.read_csv(B / "09d_table2_genes_by_type.csv")


def forest(ax, tp, genes_, types_):
    sub = E_t[(E_t.timepoint == tp) & E_t.gene.isin(genes_) & E_t.cell_type.isin(types_)]
    rows = [(g, t) for g in genes_ for t in types_ if ((sub.gene == g) & (sub.cell_type == t)).any()]
    for i, (g, t) in enumerate(rows):
        r = sub[(sub.gene == g) & (sub.cell_type == t)].iloc[0]
        ax.plot([r["CI.L"], r["CI.R"]], [i, i], color=TCOL.get(t, GREY), lw=1.1)
        ax.plot(r["logFC"], i, "o", ms=4, color=TCOL.get(t, GREY), mec="black" if r["adj.P.Val"] < 0.1 else "none", mew=0.7)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{g} · {t}" for g, t in rows], fontsize=5.4); ax.invert_yaxis()
    ax.set_xlabel("log2 fold change, TYLCV vs mock (95% CI)")


ax = fig.add_subplot(gs[2, 0]); forest(ax, "7dpi", ["Cab-1A", "Cab-1D", "IPS", "LyesC2p084"], ["Mesophyll-Ch", "Bundle sheath"])
ax.set_title("7 dpi, within types", loc="left"); letter(ax, "f", x=-1.05)
ax = fig.add_subplot(gs[2, 1:]); forest(ax, "24dpi", ["CAB1B", "ltpg1", "TLRP"], ["Bundle sheath", "Guard", "Hydathodes", "Phloem", "Vein"])
ax.set_title("24 dpi, within types (black outline: adjusted P < 0.1)", loc="left"); letter(ax, "g", x=-0.36, y=1.12)
save(fig, "Figure5")
keyg = {}
for g in genes:
    keyg[g] = {"limma7": (round(float(L7.loc[g, "logFC"]), 2), float(L7.loc[g, "adj.P.Val"])) if g in L7.index else None,
               "limma24": (round(float(L24.loc[g, "logFC"]), 2), float(L24.loc[g, "adj.P.Val"])) if g in L24.index else None,
               "deseq7": (round(float(D7t.loc[g, "log2FoldChange_leaf"]), 2), float(D7t.loc[g, "padj_leaf"]), float(D7t.loc[g, "padj_section"])) if g in D7t.index else None,
               "deseq24": (round(float(D24t.loc[g, "log2FoldChange_leaf"]), 2), float(D24t.loc[g, "padj_leaf"]), float(D24t.loc[g, "padj_section"]), bool(D24t.loc[g, "robust"])) if g in D24t.index else None}
stats["key_genes"] = keyg
stats["limma_summary"] = ls_.reset_index().to_dict("records")

# ------------------------------------------------------------------ Figure 5: virus and velocity
vd = pd.read_csv(B / "11_viral_barcode_diagnostic.csv").set_index("section").loc[ORDER]
pv = pd.read_csv(B / "11_viral_pos_vs_neg_paired_DE_top300.csv", index_col=0)
vin = pd.read_csv(B / "12_velocity_input_qc.csv").set_index("sample").loc[ORDER]
vg = pd.read_csv(B / "12_velocity_genes_by_min_shared_counts.csv", index_col=0)
vconf = pd.read_csv(B / "12_velocity_confidence_per_section.csv", index_col=0).loc[ORDER]
vnull = json.load(open(B / "12_velocity_confidence_vs_null.json"))
uf_l = pd.read_csv(B / "12_unspliced_fraction_test_by_leaf.csv", index_col=0)
uf_s = pd.read_csv(B / "12_unspliced_fraction_test_by_sample.csv", index_col=0)

fig = plt.figure(figsize=(7.2, 6.3))
gs = GridSpec(2, 4, figure=fig, hspace=0.62, wspace=0.95)
x = np.r_[0:4, 5:9]
ax = fig.add_subplot(gs[0, 0])
ax.bar(x, vd.viral_umis_in_csv + 1, color=[col(s) for s in ORDER], width=0.75); ax.set_yscale("log")
for xi, s in zip(x, ORDER):
    ax.text(xi, (vd.loc[s, "viral_umis_in_csv"] + 1) * 1.25, f"{int(vd.loc[s, 'viral_umis_in_csv']):,}", ha="center", va="bottom", fontsize=4.8, rotation=90)
ax.set_xticks(x); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=90, fontsize=5.6); ax.set_ylabel("TYLCV UMIs per library (+1)")
ax.set_ylim(top=2e4); ax.set_title("Viral UMIs per library"); letter(ax, "a", x=-0.45)
ax = fig.add_subplot(gs[0, 1])
sub = vd[vd.viral_umis_in_csv >= 10]
yy = np.arange(len(sub))[::-1]
for yi, (s, r) in zip(yy, sub.iterrows()):
    ax.plot([r.pct_expected_if_random, r.pct_viral_on_loom_tissue], [yi, yi], color="#BBBBBB", lw=1)
    ax.plot(r.pct_expected_if_random, yi, "o", color=GREY, ms=4)
    ax.plot(r.pct_viral_on_loom_tissue, yi, "o", color=col(s), ms=4.5)
ax.set_yticks(yy); ax.set_yticklabels([f"{SHORT[s]} ({int(r.viral_umis_in_csv):,})" for s, r in sub.iterrows()], fontsize=5.8)
ax.set_xlabel("Viral UMIs on tissue spots (%)"); ax.set_xlim(-1, 40)
ax.plot([], [], "o", color=GREY, label="random"); ax.plot([], [], "o", color=INF, label="observed")
ax.legend(frameon=False, fontsize=5.0, loc="lower right", handletextpad=0.2); ax.set_title("On-tissue share"); letter(ax, "b", x=-0.7)
ax = fig.add_subplot(gs[0, 2:])
isv = pv.index.str.startswith("TYLCV")
nlp = -np.log10(pv.pvalue.clip(lower=1e-12))
ax.scatter(pv.log2FoldChange[~isv], nlp[~isv], s=6, color="#CCCCCC", lw=0)
ax.scatter(pv.log2FoldChange[isv], nlp[isv], s=16, color=INF, lw=0)
for (g, r), dy in zip(pv[isv].iterrows(), (4, -9)):
    ax.annotate(f"{g} (adj. P = {r.padj:.3f})", (r.log2FoldChange, nlp[g]), fontsize=5.5, xytext=(-5, dy), textcoords="offset points", ha="right")
ax.text(0.02, 0.62, "300 genes with the smallest P shown\n(P < 10$^{-12}$ drawn at 10$^{-12}$);\nall host genes adj. P > 0.99", transform=ax.transAxes, va="top", fontsize=5.4)
ax.set_ylim(top=nlp.max() * 1.18)
ax.set_xlabel("log2 fold change, viral-positive vs negative spots"); ax.set_ylabel("−log10 P")
ax.set_title("Viral-positive vs viral-negative spots, paired by section")
letter(ax, "c", x=-0.17)
from statsmodels.stats.proportion import proportion_confint
st_ = pd.read_csv(B / "11d_viral_strand_by_capture_area.csv", index_col=0)
stt = pd.read_csv(B / "11d_viral_strand_test.csv", index_col=0).iloc[0]
ist = st_[st_.condition == "infected"].set_index("sample").loc[[s_ for s_ in ORDER if s_.startswith("infected")]]
xx = np.arange(len(ist))
ax = fig.add_subplot(gs[1, 0:2])
for xi, (s_, r) in zip(xx, ist.iterrows()):
    n_ = r.V_all + r.C_all; lo, hi = proportion_confint(r.C_all, n_, method="wilson")
    ax.bar(xi, r.pct_complementary_all, color=col(s_), width=0.65)
    ax.errorbar(xi, r.pct_complementary_all, yerr=[[r.pct_complementary_all - 100 * lo], [100 * hi - r.pct_complementary_all]],
                color="black", capsize=2, lw=0.8)
    ax.text(xi, 2.5, f"{int(n_):,}", ha="center", fontsize=5.3, color="white")
ax.set_xticks(xx); ax.set_xticklabels([SHORT[s_] for s_ in ist.index]); ax.set_ylim(0, 78)
ax.set_ylabel("Complementary-sense share\nof viral UMIs (%, 95% CI)")
ax.text(0.5, 0.97, f"24 vs 7 dpi: odds ratio {stt.odds_ratio:.2f}, P = {stt.p:.3f}", transform=ax.transAxes, ha="center", va="top", fontsize=5.8)
ax.set_title("Viral transcription programme"); letter(ax, "d", x=-0.2)
ax = fig.add_subplot(gs[1, 2:4])
for xi, (s_, r) in zip(xx, ist.iterrows()):
    ax.plot(xi + 0.15, r.pct_complementary_off_tissue, "s", mfc="white", mec=col(s_), ms=4.5, mew=1)
    if np.isfinite(r.pct_complementary_on_tissue):
        ax.plot([xi - 0.15, xi + 0.15], [r.pct_complementary_on_tissue, r.pct_complementary_off_tissue], color="#999999", lw=0.8, zorder=0)
        ax.plot(xi - 0.15, r.pct_complementary_on_tissue, "o", color=col(s_), ms=4.5)
        ax.text(xi, max(r.pct_complementary_on_tissue, r.pct_complementary_off_tissue) + 5, f"P = {r.fisher_p_on_vs_off:.2f}", ha="center", fontsize=5.2)
    else:
        ax.text(xi, r.pct_complementary_off_tissue + 5, "none on\ntissue", ha="center", fontsize=5.0)
ax.plot([], [], "o", color=GREY, label="on tissue"); ax.plot([], [], "s", mfc="white", mec=GREY, label="off tissue")
ax.legend(frameon=False, fontsize=5.5, loc="lower right")
ax.set_xticks(xx); ax.set_xticklabels([SHORT[s_] for s_ in ist.index]); ax.set_xlim(-0.5, len(xx) - 0.5); ax.set_ylim(0, 78)
ax.set_ylabel("Complementary-sense share (%)"); ax.set_title("Strand composition on and off tissue"); letter(ax, "e", x=-0.2)
stats["viral_strand"] = {"by_area": ist[["V_all", "C_all", "pct_complementary_all", "pct_complementary_on_tissue",
                                        "pct_complementary_off_tissue", "fisher_p_on_vs_off"]].round(4).to_dict("index"), "test": stt.to_dict()}
save(fig, "Figure6")
stats["viral"] = vd[["viral_umis_in_csv", "pct_viral_on_loom_tissue", "pct_expected_if_random", "pct_anndata_barcodes_in_loom"]].to_dict("index")
stats["velocity"] = {"unspliced_pct_range": (round(float(vin.unspliced_fraction.min() * 100), 1), round(float(vin.unspliced_fraction.max() * 100), 1)),
                     "genes_by_threshold": vg.n_genes.to_dict(), **vnull,
                     "uf_tests_leaf": int(len(uf_l)), "uf_tests_section": int(len(uf_s)),
                     "uf_min_padj_leaf": float(uf_l.padj_bh.min()), "uf_min_padj_section": float(uf_s.padj_bh.min())}

# ------------------------------------------------------------------ Figure 7: tropism and splicing kinetics by cell type
tr = pd.read_csv(B / "11b_viral_tropism_by_type.csv")
uf_t = pd.read_csv(B / "12b_unspliced_fraction_by_type_leaf.csv", index_col=0)
uf_p = pd.read_csv(B / "12_unspliced_fraction_test_by_leaf.csv", index_col=0)
fig = plt.figure(figsize=(7.2, 5.0))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1], hspace=0.6, wspace=0.75)
# a | observed / expected viral UMIs per imputed type
ax = fig.add_subplot(gs[:, 0])
pooled = tr[tr.scope == "pooled_all"].set_index("cell_type").sort_values("obs_over_exp")
sec_names = {"infected_7dpi_2": ("I7-2", "#E69F00"), "infected_24dpi_1": ("I24-1", "#56B4E9"), "infected_24dpi_2": ("I24-2", "#009E73")}
yy = {t: i for i, t in enumerate(pooled.index)}
for sc, (lab, c) in sec_names.items():
    g = tr[tr.scope == sc]
    g = g[g.cell_type.isin(yy)]
    ax.scatter(np.log2(g.obs_over_exp.clip(lower=1 / 16)), g.cell_type.map(yy) + 0.18, s=12, color=c, alpha=0.9, lw=0, label=lab, zorder=2)
ax.scatter(np.log2(pooled.obs_over_exp), [yy[t] for t in pooled.index], marker="D", s=26, color="black", zorder=3, label="pooled")
for t, r in pooled.iterrows():
    ax.text(3.35, yy[t], f"{r.obs_over_exp:.2f}   {r.p_two_sided:.2f}   {int(r.n_sections_enriched)}/{int(r.n_sections)}", va="center", fontsize=5.3)
ax.axvline(0, color="black", lw=0.6, ls="--")
ax.set_yticks(range(len(pooled))); ax.set_yticklabels(pooled.index)
ax.set_xlim(-4.2, 6.6); ax.set_xticks([-4, -2, 0, 2]); ax.set_xticklabels(["1/16", "1/4", "1", "4"])
ax.axvline(3.15, color="#DDDDDD", lw=0.6)
ax.set_xlabel("Viral UMIs, observed / expected from host UMIs")
ax.text(3.35, len(pooled) - 0.45, "O/E     P     areas", fontsize=5.3, fontweight="bold", va="bottom")
ax.legend(frameon=False, fontsize=5.4, loc="upper left", handletextpad=0.2)
ax.set_title("Viral UMIs by imputed cell type (none at adj. P < 0.05)", loc="left", pad=12); letter(ax, "a", x=-0.38)
# b | all within-type splicing tests
ax = fig.add_subplot(gs[0, 1])
for t_, mk, c in [("7dpi", "o", "#7F7F7F"), ("24dpi", "s", "#222222")]:
    g = uf_t[uf_t.time == t_]
    ax.scatter(g.log_odds_diff, -np.log10(g.p), s=11, marker=mk, color=c, lw=0, label=f"{t_.replace('dpi', ' dpi')} ({len(g)} tests)")
rb = uf_t.loc[uf_t.index == "Rbcs-3B"]
ax.scatter(rb.log_odds_diff, -np.log10(rb.p), s=26, facecolors="none", edgecolors="#D55E00", lw=0.9, label="Rbcs-3B")
for g_, r_ in uf_t[uf_t.padj_all < 0.05].iterrows():
    ax.annotate(f"{g_}, {r_.cell_type} ({r_.time.replace('dpi', ' dpi')})\nFDR = {r_.padj_all:.3f}", (r_.log_odds_diff, -np.log10(r_.p)),
                fontsize=5.0, xytext=(6, -2), textcoords="offset points", ha="left", va="top")
ax.axvline(0, color="black", lw=0.6, ls="--")
ax.set_xlabel("Δ log-odds unspliced, TYLCV vs mock"); ax.set_ylabel("−log10 P (leaves as replicates)")
ax.legend(frameon=False, fontsize=5.0, loc="upper left", bbox_to_anchor=(0.0, 1.0), handletextpad=0.2, borderaxespad=0.2)
ax.set_ylim(top=-np.log10(uf_t.p.min()) * 1.45)
ax.set_title("Splicing kinetics within imputed types", loc="left"); letter(ax, "b", x=-0.3)
# c | Rbcs-3B unspliced fraction per type at 7 dpi
ax = fig.add_subplot(gs[1, 1])
r7 = rb[rb.time == "7dpi"].sort_values("log_odds_diff")
yy7 = np.arange(len(r7))
ax.hlines(yy7, r7.unspliced_frac_infected, r7.unspliced_frac_mock, color="#BBBBBB", lw=1)
ax.scatter(r7.unspliced_frac_mock, yy7, s=16, facecolors="white", edgecolors=MOCK, lw=0.9, label="mock", zorder=3)
ax.scatter(r7.unspliced_frac_infected, yy7, s=16, color=INF, label="TYLCV", zorder=3)
pp = uf_p.loc[(uf_p.index == "Rbcs-3B") & (uf_p.time == "7dpi")].iloc[0]
ax.axvline(pp.unspliced_frac_mock, color=MOCK, lw=0.6, ls=":"); ax.axvline(pp.unspliced_frac_infected, color=INF, lw=0.6, ls=":")
ax.set_yticks(yy7); ax.set_yticklabels(r7.cell_type, fontsize=5.6)
ax.set_xlabel("Rbcs-3B unspliced fraction, 7 dpi"); ax.legend(frameon=False, fontsize=5.2, loc="upper left")
ax.set_title(f"Rbcs-3B: lower in TYLCV in {int((r7.log_odds_diff < 0).sum())} of {len(r7)} types (dotted: pooled)", loc="left"); letter(ax, "c", x=-0.3)
save(fig, "Figure7")

# ------------------------------------------------------------------ Supplementary figures (maps and streams from the audit run)
fig = plt.figure(figsize=(7.2, 3.6)); ax = fig.add_subplot(111); img(ax, B / "11_viral_maps.png"); save(fig, "FigureS1")
fig = plt.figure(figsize=(7.2, 5.2))
gs = GridSpec(2, 3, figure=fig, hspace=0.75, wspace=0.85, width_ratios=[2.1, 1, 1])
x = np.r_[0:4, 5:9]
ax = fig.add_subplot(gs[:, 0]); img(ax, B / "12_velocity_stream_umap.png"); letter(ax, "a", x=0.0, y=1.0)
ax = fig.add_subplot(gs[0, 1])
ax.bar(x, vin.unspliced_fraction * 100, color=[col(s) for s in ORDER], width=0.75)
ax.set_xticks(x); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=90, fontsize=5.6); ax.set_ylabel("Unspliced UMIs (%)")
ax.set_title("Unspliced signal"); letter(ax, "b", x=-0.45)
ax = fig.add_subplot(gs[0, 2])
thr = [int(i.split(">=")[1]) for i in vg.index]
ax.bar(range(len(thr)), vg.n_genes, color="#44AA99"); ax.set_xticks(range(len(thr))); ax.set_xticklabels([f"≥{t}" for t in thr])
ax.axhline(vnull["n_velocity_genes"], ls="--", lw=0.7, color="black")
ax.text(3.4, vnull["n_velocity_genes"] + 8, f"velocity genes: {vnull['n_velocity_genes']}", ha="right", fontsize=5.5)
ax.set_xlabel("Shared spliced+unspliced counts"); ax.set_ylabel("Genes"); ax.set_title("Genes usable for velocity"); letter(ax, "c", x=-0.45)
ax = fig.add_subplot(gs[1, 1])
ax.bar(x, vconf.pct_spots_above_null_q95, color=[col(s) for s in ORDER], width=0.75)
ax.axhline(5, ls="--", lw=0.7, color="black"); ax.text(8.6, 7, "null: 5%", fontsize=5.5, ha="right")
ax.set_xticks(x); ax.set_xticklabels([SHORT[s] for s in ORDER], rotation=90, fontsize=5.6)
ax.set_ylabel("Spots above shuffled-null 95th pct (%)")
ax.set_title("Velocity confidence")
letter(ax, "d", x=-0.45)
ax = fig.add_subplot(gs[1, 2])
for d_, c, lab in ((uf_l, "#332288", "leaf"), (uf_s, "#CC6677", "section")):
    p = np.sort(d_.p.dropna().to_numpy()); n = len(p)
    ax.plot(-np.log10((np.arange(1, n + 1) - 0.5) / n), -np.log10(p), "o", ms=2.5, color=c, label=f"{lab} (n={n})")
m_ = 2.6; ax.plot([0, m_], [0, m_], "--", lw=0.6, color="black")
ax.set_xlabel("Expected −log10 P"); ax.set_ylabel("Observed −log10 P"); ax.legend(frameon=False, fontsize=5.5)
ax.set_title("Unspliced-fraction tests"); letter(ax, "e", x=-0.45)
save(fig, "FigureS2")
fig = plt.figure(figsize=(7.2, 4.2)); leaf_panels(fig, GridSpec(1, 1, figure=fig, top=0.93, bottom=0.02)[0, 0], "pieces")
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=l) for l, c in (("pieces of leaf 1", LEAF_COL["L1"]), ("pieces of leaf 2", LEAF_COL["L2"]),
           ("pieces of leaf 3", LEAF_COL["L3"]), ("fragments", LEAF_COL["fragment"]))], loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=6)
save(fig, "FigureS3")

stats["spots_total_kept"] = int(pd.read_csv(B / "03_qc_filter.csv").kept.sum())
stats["spots_total_raw"] = int(qc.n_spots.sum())
stats["leaf_seg"] = pd.read_csv(B / "03_leaf_segmentation.csv")[["sample", "pieces", "leaves", "leaf_sizes"]].to_dict("records")
stats["seed"] = json.load(open(B / "07_seed_stability.json"))
json.dump(stats, open(OUT / "stats.json", "w"), indent=1, default=str)
print(json.dumps(stats, indent=1, default=str)[:9000])


# ------------------------------------------------------------------ placeholders (authors' figures, pending results)
def placeholder(name, title, lines, size=(7.2, 4.2)):
    fig = plt.figure(figsize=size); ax = fig.add_axes([0.02, 0.02, 0.96, 0.96]); ax.set_axis_off()
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="#F4F4F4", edgecolor="#999999", lw=1, ls="--"))
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=11, fontweight="bold", color="#444444")
    ax.text(0.5, 0.42, lines, ha="center", va="center", fontsize=8, color="#555555", linespacing=1.6)
    fig.savefig(OUT / f"{name}.png", dpi=200, facecolor="white"); plt.close(fig)
placeholder("FigureS4_placeholder", "Supplementary Figure 4 | PlantAG-imputed cell types in space (authors' figure)",
            "Insert the spatial maps of the ten imputed types in the eight capture areas (colours as in Fig. 1g).")
