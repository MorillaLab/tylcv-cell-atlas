"""
plantag_audit.py — reliability audit and replicate-aware re-analysis of the TYLCV Visium dataset.

Written for Google Colab with the Drive layout of `test1_tylcv_7_24dpi.ipynb`, but runs anywhere.
Required: numpy, scipy, pandas, scikit-learn, anndata, scanpy (+ leidenalg/igraph), statsmodels,
          matplotlib, h5py. Optional: pydeseq2 (pseudobulk DE), harmonypy (integration),
          gudhi + persim (re-computing persistence images).
No torch / torch_geometric is needed: the GCN used to audit PlantAG is implemented in numpy/scipy and
follows PyG GCNConv propagation (symmetric normalisation, self-loops, edge weights, bias after
aggregation) with the PlantAG architecture (3 x GCN(256) -> Linear, dropout 0.5, Adam lr 5e-3, wd 1e-4).
"""
from __future__ import annotations

import itertools
import json
import os
import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------
# 0. Small utilities
# ----------------------------------------------------------------------------------------------

class Bundle:
    """Collects small, shareable outputs (CSV / PNG / JSON). Heavy objects stay on Drive."""

    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def table(self, df, name, index=True):
        path = self.root / f"{name}.csv"
        df.to_csv(path, index=index)
        return path

    def fig(self, fig, name, dpi=130, show=True):
        path = self.root / f"{name}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        if show:                                       # also show it inline when running in Jupyter/Colab
            try:
                from IPython import get_ipython
                from IPython.display import display
                if get_ipython() is not None:
                    display(fig)
            except Exception:
                pass
        plt.close(fig)
        return path

    def json(self, obj, name):
        path = self.root / f"{name}.json"
        path.write_text(json.dumps(obj, indent=2, default=str))
        return path

    def zip(self, max_mb=29.0):
        archive = shutil.make_archive(str(self.root), "zip", self.root)
        size_mb = os.path.getsize(archive) / 1e6
        if size_mb > max_mb:
            warnings.warn(f"Bundle is {size_mb:.1f} MB; split it (split -b 25M) before uploading.")
        return archive, round(size_mb, 2)


def _as_str(series):
    """Labels as an object array of str; missing values become "nan" (pandas 3 keeps NaN after astype(str))."""
    s = pd.Series(series, copy=False)
    return s.astype(object).where(s.notna(), "nan").astype(str).to_numpy(dtype=object)


def get_counts(adata, layer=None):
    """Raw counts as CSR; refuses normalised / scaled matrices."""
    X = adata.layers[layer] if layer else adata.X
    X = sp.csr_matrix(X)
    if X.nnz and not np.allclose(X.data[: min(X.nnz, 200_000)], np.round(X.data[: min(X.nnz, 200_000)])):
        raise ValueError("This matrix is not integer counts. Load the raw object or pass the counts layer.")
    return X


def log_normalise(counts, target_sum=None):
    """Library-size normalisation + log1p (same as sc.pp.normalize_total(target_sum=None) + log1p)."""
    lib = np.asarray(counts.sum(1)).ravel()
    target = np.median(lib[lib > 0]) if target_sum is None else target_sum
    scale = target / np.maximum(lib, 1)
    out = sp.diags(scale) @ counts
    out = sp.csr_matrix(out)
    out.data = np.log1p(out.data)
    return out


def barcode16(names):
    """16-nt Visium barcode from obs_names or velocyto CellIDs ('prefix:BARCODEx', 'BARCODE-1_sample', ...)."""
    pat = re.compile(r"[ACGT]{16}")
    out = []
    for name in names:
        tail = str(name).split(":")[-1]
        m = pat.search(tail)
        out.append(m.group(0) if m else None)
    return np.array(out, dtype=object)


def physical_coords(obs, row_key="array_row", col_key="array_col", pitch_um=100.0):
    """Visium array indices -> micrometres (adjacent spots are 100 µm apart, rows 86.6 µm)."""
    x = obs[col_key].to_numpy(float) * pitch_um / 2.0
    y = obs[row_key].to_numpy(float) * pitch_um * np.sqrt(3) / 2.0
    return np.c_[x, y]


def _match_any(index, patterns):
    mask = np.zeros(len(index), bool)
    for p in patterns:
        mask |= np.asarray(pd.Index(index).astype(str).str.contains(p, regex=True), dtype=bool)
    return mask


# ----------------------------------------------------------------------------------------------
# 1. Provenance and QC
# ----------------------------------------------------------------------------------------------

def organelle_genes_from_gtf(gtf_path,
                             mito_regex=r"(?i)(mito|^MT$|^chrM$|NC_035963)",
                             plastid_regex=r"(?i)(chloro|plastid|^Pt$|^Pltd$|^chrC$|NC_007898)"):
    """Gene names/ids on organelle sequences of the GTF used by Space Ranger (most reliable source)."""
    import gzip
    mito, plastid = set(), set()
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            is_m, is_p = re.search(mito_regex, f[0]), re.search(plastid_regex, f[0])
            if not (is_m or is_p):
                continue
            names = set(re.findall(r'(?:gene_id|gene_name|gene)\s+"([^"]+)"', f[8]))
            (mito if is_m else plastid).update(names)
    return mito, plastid


def flag_organelles(adata, mito_genes=(), plastid_genes=(), mito_patterns=(r"^CQW56_",), plastid_patterns=()):
    """adata.var['organelle'] in {'mito','plastid','nuclear'}. CQW56_ = tomato mitochondrial locus tags."""
    v = pd.Index(adata.var_names.astype(str))
    m = v.isin(list(mito_genes)) | _match_any(v, mito_patterns)
    p = (v.isin(list(plastid_genes)) | _match_any(v, plastid_patterns)) & ~m
    adata.var["organelle"] = np.select([m, p], ["mito", "plastid"], "nuclear")
    return adata.var["organelle"].value_counts()


def spot_qc(adata, counts, sample_key="sample"):
    tot = np.asarray(counts.sum(1)).ravel()
    org = adata.var["organelle"].to_numpy() if "organelle" in adata.var else np.full(adata.n_vars, "nuclear")
    def frac(kind):
        cols = np.where(org == kind)[0]
        if len(cols) == 0:
            return np.zeros(len(tot))
        return np.asarray(counts[:, cols].sum(1)).ravel() / np.maximum(tot, 1) * 100
    df = pd.DataFrame({"sample": _as_str(adata.obs[sample_key]), "total_counts": tot,
                       "n_genes": np.diff(counts.indptr), "pct_mito": frac("mito"),
                       "pct_plastid": frac("plastid")}, index=adata.obs_names)
    if "in_tissue" in adata.obs:
        df["in_tissue"] = adata.obs["in_tissue"].to_numpy()
    return df


def sample_qc_summary(qc):
    g = qc.groupby("sample")
    out = pd.DataFrame({
        "n_spots": g.size(),
        "n_in_tissue_0": g["in_tissue"].apply(lambda s: int((s == 0).sum())) if "in_tissue" in qc else np.nan,
        "median_umis": g["total_counts"].median(),
        "q25_umis": g["total_counts"].quantile(0.25),
        "q75_umis": g["total_counts"].quantile(0.75),
        "pct_spots_lt_50_umis": g["total_counts"].apply(lambda s: 100 * (s < 50).mean()),
        "median_genes": g["n_genes"].median(),
        "median_pct_mito": g["pct_mito"].median(),
        "median_pct_plastid": g["pct_plastid"].median(),
    })
    return out.round(2)


def plot_qc(qc):
    samples = sorted(qc["sample"].unique())
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, col, lab in zip(axes, ["total_counts", "pct_mito", "pct_plastid"],
                            ["log10 UMIs per spot", "% mitochondrial UMIs", "% plastid UMIs"]):
        data = [np.log10(qc.loc[qc["sample"] == s, col] + 1) if col == "total_counts"
                else qc.loc[qc["sample"] == s, col] for s in samples]
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(samples) + 1))
        ax.set_xticklabels(samples, rotation=60, ha="right", fontsize=8)
        ax.set_title(lab)
    fig.tight_layout()
    return fig


def qc_keep_mask(qc, mad_k=3.0, min_counts=None, min_genes=None, require_in_tissue=True, max_pct_organelle=None):
    """Removes out-of-tissue spots and per-section low-library outliers (median - k*MAD on log counts)."""
    from scipy.stats import median_abs_deviation
    keep = np.ones(len(qc), bool)
    if require_in_tissue and "in_tissue" in qc:
        keep &= qc["in_tissue"].to_numpy() == 1
    logc = np.log1p(qc["total_counts"])
    lo = logc.groupby(qc["sample"]).transform(lambda s: s.median() - mad_k * median_abs_deviation(s, scale="normal"))
    keep &= (logc >= lo).to_numpy()
    if min_counts:
        keep &= qc["total_counts"].to_numpy() >= min_counts
    if min_genes:
        keep &= qc["n_genes"].to_numpy() >= min_genes
    if max_pct_organelle is not None:
        keep &= (qc["pct_mito"] + qc["pct_plastid"]).to_numpy() <= max_pct_organelle
    return keep


# ----------------------------------------------------------------------------------------------
# 2. Velocyto / loom provenance (the original integration matched loom rows to spots by position)
# ----------------------------------------------------------------------------------------------

def _read_loom(path, genes=None):
    import h5py
    with h5py.File(path, "r") as f:
        ids = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["col_attrs"]["CellID"][:]])
        gnames = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["row_attrs"]["Gene"][:]])
        layers = {}
        for name in ("spliced", "unspliced", "ambiguous"):
            if name in f["layers"]:
                layers[name] = f["layers"][name]
        total = np.zeros(len(ids))
        step = 2000
        for name, ds in layers.items():
            if name == "ambiguous":
                continue
            for i in range(0, ds.shape[0], step):
                total += ds[i:i + step, :].sum(0)
        out = {"ids": ids, "genes": gnames, "total": total}
        if genes is not None:
            sel = np.where(np.isin(gnames, genes))[0]
            out["genes_sel"] = gnames[sel]
            for name, ds in layers.items():
                if name == "ambiguous":
                    continue
                out[name] = sp.csr_matrix(ds[np.sort(sel), :].T)
            out["genes_sel"] = gnames[np.sort(sel)]
    return out


def audit_loom_mapping(adata, counts, loom_map, loom_dir, sample_key="sample"):
    """Per sample: does the loom hold the same barcodes, in the same order, with matching counts?
    Visium barcodes encode array position and repeat across capture areas, so identity is tested with
    the per-barcode count correlation between each loom and each AnnData section (not barcode overlap)."""
    from scipy.stats import spearmanr
    samples = list(loom_map)
    tot = np.asarray(counts.sum(1)).ravel()
    ad_bc, ad_tot = {}, {}
    for s in samples:
        m = _as_str(adata.obs[sample_key]) == s
        ad_bc[s] = barcode16(adata.obs_names[m])
        ad_tot[s] = pd.Series(tot[m], index=ad_bc[s]).groupby(level=0).sum()
    looms, rows = {}, []
    for s, fn in loom_map.items():
        path = Path(loom_dir) / fn
        if not path.exists():
            rows.append({"sample": s, "loom": fn, "status": "missing"})
            continue
        lo = _read_loom(path)
        bc = barcode16(lo["ids"])
        looms[s] = pd.Series(lo["total"], index=bc).groupby(level=0).sum()
        a = ad_bc[s]
        n = min(len(a), len(bc))
        rows.append({"sample": s, "loom": fn, "status": "ok", "n_spots_adata": len(a), "n_cells_loom": len(bc),
                     "barcode_overlap": round(float(np.isin(a, bc).mean()), 3),
                     "positional_match_used_by_original_code": round(float(np.mean(a[:n] == bc[:n])), 3)})
    corr = pd.DataFrame(index=[f"adata:{s}" for s in samples], columns=[f"loom:{s}" for s in looms], dtype=float)
    for s in samples:
        for t, lt in looms.items():
            common = ad_tot[s].index.intersection(lt.index)
            if len(common) >= 20:
                corr.loc[f"adata:{s}", f"loom:{t}"] = spearmanr(ad_tot[s][common], lt[common])[0]
    summary = pd.DataFrame(rows).set_index("sample")
    best = corr.astype(float).idxmax(axis=1).str.replace("loom:", "", regex=False)
    summary["loom_that_best_matches_counts"] = [best.get(f"adata:{s}") for s in summary.index]
    return summary, corr.round(3)


def integrate_velocity_by_barcode(adata, loom_map, loom_dir, sample_key="sample"):
    """Correct alternative to the original integration: matches each spot to its own barcode in the
    loom that actually belongs to its section. Returns spliced/unspliced CSR matrices aligned to adata."""
    genes = adata.var_names.to_numpy()
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    S = sp.lil_matrix((adata.n_obs, len(genes)), dtype=np.float32)
    U = sp.lil_matrix((adata.n_obs, len(genes)), dtype=np.float32)
    report = {}
    for s, fn in loom_map.items():
        rows = np.where(_as_str(adata.obs[sample_key]) == s)[0]
        lo = _read_loom(Path(loom_dir) / fn, genes=genes)
        lb = pd.Series(np.arange(len(lo["ids"])), index=barcode16(lo["ids"]))
        lb = lb[~lb.index.duplicated()]
        bc = barcode16(adata.obs_names[rows])
        ok = np.isin(bc, lb.index)
        cols = gpos[lo["genes_sel"]].to_numpy()
        src = lb[bc[ok]].to_numpy()
        S[rows[ok][:, None], cols] = lo["spliced"][src].toarray()
        U[rows[ok][:, None], cols] = lo["unspliced"][src].toarray()
        report[s] = {"spots": len(rows), "matched_by_barcode": int(ok.sum())}
    return S.tocsr(), U.tocsr(), pd.DataFrame(report).T


# ----------------------------------------------------------------------------------------------
# 3. Replicate structure: leaves inside each capture area
# ----------------------------------------------------------------------------------------------

def spatial_adjacency(adata, sample_key="sample", radius_um=101.0):
    """Hexagonal neighbour graph (within sections only). 101 µm = first ring; 175 µm bridges 1-spot gaps."""
    from sklearn.neighbors import radius_neighbors_graph
    blocks, order = [], []
    samples = _as_str(adata.obs[sample_key])
    for s in pd.unique(samples):
        idx = np.where(samples == s)[0]
        xy = physical_coords(adata.obs.iloc[idx])
        blocks.append(radius_neighbors_graph(xy, radius=radius_um, mode="connectivity", include_self=False))
        order.append(idx)
    A = sp.block_diag(blocks).tocsr()
    perm = np.concatenate(order)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    return A[inv][:, inv].tocsr()


def _weakest_cut(A, min_spots):
    """Fiedler-vector sweep on a connected piece: the cut with the lowest conductance that leaves >= min_spots
    on each side. Returns (conductance, boolean side) or None if the piece is too small to cut."""
    from scipy.sparse.csgraph import laplacian
    A = sp.csr_matrix(A, dtype=float)
    n = A.shape[0]
    if n < 2 * min_spots:
        return None
    L = laplacian(A, normed=True)
    if n <= 3000:
        f = np.linalg.eigh(L.toarray())[1][:, 1]
    else:
        from scipy.sparse.linalg import eigsh
        w, V = eigsh(L, k=2, sigma=-1e-4, which="LM")
        f = V[:, np.argsort(w)[1]]
    order = np.argsort(f, kind="stable")
    deg = np.asarray(A.sum(1)).ravel()
    vol_tot = deg.sum()
    inside = np.zeros(n, bool)
    cut = vol = 0.0
    best_phi, best_k = np.inf, None
    for k, v in enumerate(order[:-1], start=1):
        lo, hi = A.indptr[v], A.indptr[v + 1]
        cut += deg[v] - 2 * A.data[lo:hi][inside[A.indices[lo:hi]]].sum()
        vol += deg[v]
        inside[v] = True
        if min_spots <= k <= n - min_spots:
            phi = cut / max(min(vol, vol_tot - vol), 1e-12)
            if phi < best_phi:
                best_phi, best_k = phi, k
    if best_k is None:
        return None
    side = np.zeros(n, bool)
    side[order[:best_k]] = True
    return float(best_phi), side


def segment_leaves(adata, sample_key="sample", min_spots=15, radius_um=101.0, n_expected=None,
                   split=(), max_conductance=None, group_pieces=True, attach_fragments_um=200.0, leaf_groups=None):
    """Leaves = connected tissue pieces in each capture area, largest first (<section>_L1, _L2, ...).
    Touching leaves merge into one piece. Nothing is cut unless you ask:
      * every piece >= 2*min_spots in a section with fewer than `n_expected` leaves gets a SUGGESTED cut: the
        weakest spectral cut (lowest conductance) leaving >= min_spots per side, shown in adata.obs['suggested_cut']
        and in the piece table;
      * `split`: piece names from that table (e.g. 'infected_7dpi_2_P1') to cut; list a name twice to cut twice;
      * `max_conductance`: also cut any piece whose weakest cut is <= this value. In tests touching or slightly
        overlapping leaves scored 0.006-0.014 and single leaves 0.022-0.066, but large, narrow single leaves can
        score lower, so check the plot before relying on it.
    Broken leaves: if a section has MORE pieces than `n_expected` and group_pieces=True, pieces are merged into
    n_expected leaves by single linkage on the smallest gap between them (pieces of one leaf lie closer together
    than different leaves); the gaps bridged are reported. Fragments (< min_spots) within attach_fragments_um of a
    leaf join it. `leaf_groups` overrides this per section: {"infected_24dpi_1": [["P1", "P4"], ["P2", "P3"], ["P5"]]}
    using names from the piece table (cut parts are named P1.1, P1.2).
    Returns (section_table, piece_table)."""
    from collections import Counter
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    want = Counter(split) if not isinstance(split, dict) else Counter(split)
    leaf_groups = leaf_groups or {}
    A = spatial_adjacency(adata, sample_key, radius_um)
    samples = _as_str(adata.obs[sample_key])
    leaf = np.empty(adata.n_obs, dtype=object)
    cutlab = np.full(adata.n_obs, "-", dtype=object)
    sec_rows, piece_rows = [], []
    for s in pd.unique(samples):
        idx = np.where(samples == s)[0]
        xy = physical_coords(adata.obs.iloc[idx])
        As = A[idx][:, idx].tocsr()
        ncomp, lab = connected_components(As, directed=False)
        sizes = np.bincount(lab)
        order = np.argsort(-sizes, kind="stable")
        parts, frags = [], [np.where(lab == c)[0] for c in order if sizes[c] < min_spots]
        big = [c for c in order if sizes[c] >= min_spots]
        suggest = bool(n_expected) and len(big) < n_expected
        for rank, c in enumerate(big, start=1):
            name = f"{s}_P{rank}"
            members = np.where(lab == c)[0]
            n_cuts = want.pop(name, 0)
            res = _weakest_cut(As[members][:, members], min_spots) if (suggest or n_cuts or max_conductance) else None
            auto = res is not None and max_conductance is not None and res[0] <= max_conductance
            n_cuts = max(n_cuts, 1 if auto else 0)
            if res is not None:
                cutlab[idx[members[res[1]]]] = f"{name}: side A"
                cutlab[idx[members[~res[1]]]] = f"{name}: side B"
            pieces = [members]
            for _ in range(n_cuts if res else 0):
                pieces.sort(key=len, reverse=True)
                target = pieces.pop(0)
                r = _weakest_cut(As[target][:, target], min_spots)
                if r is None:
                    pieces.append(target)
                    break
                for side in (target[r[1]], target[~r[1]]):
                    nc, sl = connected_components(As[side][:, side], directed=False)
                    for k in range(nc):
                        (pieces if (sl == k).sum() >= min_spots else frags).append(side[sl == k])
            pieces.sort(key=len, reverse=True)
            pnames = [name] if len(pieces) == 1 else [f"{name}.{k}" for k in range(1, len(pieces) + 1)]
            parts += list(zip(pnames, pieces))
            piece_rows.append({"section": s, "piece": name, "spots": len(members),
                               "weakest_cut_conductance": round(res[0], 4) if res else np.nan,
                               "cut_sizes": f"{res[1].sum()}/{(~res[1]).sum()}" if res else "",
                               "cuts_applied": n_cuts if res else 0})
        # ---- group parts into leaves
        short = lambda n: n[len(s) + 1:]
        how, gaps = "pieces = leaves", []
        if s in leaf_groups:
            lookup = {short(n): m for n, m in parts}
            groups, used = [], set()
            for g in leaf_groups[s]:
                missing = [x for x in g if x not in lookup]
                if missing:
                    raise ValueError(f"{s}: unknown piece names {missing}; valid: {sorted(lookup)}")
                groups.append(([x for x in g], np.concatenate([lookup[x] for x in g])))
                used.update(g)
            groups += [([k], m) for k, m in lookup.items() if k not in used]
            how = "manual leaf_groups"
        elif group_pieces and n_expected and len(parts) > n_expected:
            D = np.full((len(parts), len(parts)), np.inf)
            trees = [cKDTree(xy[m]) for _, m in parts]
            for a in range(len(parts)):
                for b in range(a + 1, len(parts)):
                    D[a, b] = D[b, a] = trees[a].query(xy[parts[b][1]])[0].min()
            clusters = [[k] for k in range(len(parts))]
            while len(clusters) > n_expected:
                best = min(((D[np.ix_(clusters[a], clusters[b])].min(), a, b)
                            for a in range(len(clusters)) for b in range(a + 1, len(clusters))))
                gaps.append(round(float(best[0])))
                clusters[best[1]] += clusters[best[2]]
                del clusters[best[2]]
            groups = [([short(parts[k][0]) for k in cl], np.concatenate([parts[k][1] for k in cl])) for cl in clusters]
            how = f"{len(parts)} pieces grouped into {n_expected} by proximity: verify in 03_leaves"
        else:
            groups = [([short(n)], m) for n, m in parts]
        if attach_fragments_um and groups and frags:
            gtrees = [cKDTree(xy[m]) for _, m in groups]
            still = []
            for f in frags:
                d = [t.query(xy[f])[0].min() for t in gtrees]
                k = int(np.argmin(d))
                if d[k] <= attach_fragments_um:
                    groups[k] = (groups[k][0], np.concatenate([groups[k][1], f]))
                else:
                    still.append(f)
            frags = still
        groups.sort(key=lambda g: len(g[1]), reverse=True)
        part_to_leaf = {}
        for rank, (names_g, g) in enumerate(groups, start=1):
            leaf[idx[g]] = f"{s}_L{rank}"
            for nm in names_g:
                part_to_leaf[nm.split(".")[0]] = part_to_leaf.get(nm.split(".")[0], []) + [f"L{rank}"]
        for f in frags:
            leaf[idx[f]] = f"{s}_fragment"
        for row in piece_rows:
            if row["section"] == s:
                row["leaf"] = ",".join(sorted(set(part_to_leaf.get(short(row["piece"]), [])))) or "fragment"
        nl = len(groups)
        if nl == 0:
            status = f"no piece with >= {min_spots} spots: check in_tissue and array_row/array_col for this section"
        elif n_expected and nl < n_expected:
            status = f"{nl} leaves, expected {n_expected}: check '03_suggested_cuts'"
        elif n_expected and nl > n_expected:
            status = f"{nl} leaves, expected {n_expected}: check the plot"
        else:
            status = "ok" if how == "pieces = leaves" else how
        if gaps and max(gaps) > 1000:
            status += f" | largest gap bridged {max(gaps)} um: check"
        sec_rows.append({"sample": s, "pieces": ncomp, "leaves": nl, "grouping": how,
                         "gaps_bridged_um": gaps, "fragment_spots": int(sum(len(f) for f in frags)),
                         "leaf_sizes": sorted([len(g) for _, g in groups], reverse=True),
                         "leaf_members": " | ".join("+".join(n) for n, _ in groups), "status": status})
    if want:
        warnings.warn(f"split names not found (use the 'piece' column): {sorted(want)}")
    adata.obs["leaf"] = pd.Categorical(leaf)
    adata.obs["suggested_cut"] = pd.Categorical(cutlab)
    return pd.DataFrame(sec_rows).set_index("sample"), pd.DataFrame(piece_rows)


def plot_by_section(adata, color_key, sample_key="sample", ncols=4, s=6, cmap=None, title=None):
    samples = list(pd.unique(_as_str(adata.obs[sample_key])))
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.8 * nrows), squeeze=False)
    vals = adata.obs[color_key]
    is_cat = not pd.api.types.is_numeric_dtype(vals)
    if is_cat:
        cats = list(pd.unique(_as_str(vals)))
        pal = plt.get_cmap("tab20")
        cmap_d = {c: pal(i % 20) for i, c in enumerate([c for c in cats if c != "-"])}
        cmap_d["-"] = (0.85, 0.85, 0.85, 1.0)          # light grey for 'nothing to show here
    for ax, smp in zip(axes.ravel(), samples):
        m = _as_str(adata.obs[sample_key]) == smp
        xy = physical_coords(adata.obs[m])
        if is_cat:
            ax.scatter(xy[:, 0], -xy[:, 1], c=[cmap_d[v] for v in _as_str(vals[m])], s=s, linewidths=0)
        else:
            sca = ax.scatter(xy[:, 0], -xy[:, 1], c=vals[m].to_numpy(float), s=s, cmap=cmap or "viridis", linewidths=0)
            fig.colorbar(sca, ax=ax, fraction=0.04)
        ax.set_title(smp, fontsize=9)
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    if is_cat:
        handles = [plt.Line2D([], [], marker="o", ls="", color=cmap_d[c], label=c) for c in cats[:30]]
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7)
    fig.suptitle(title or color_key)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# 4. Marker dictionary audit and the pseudo-labelling rule
# ----------------------------------------------------------------------------------------------

STRESS_REGEX = (r"(?i)^(PR-?1|PR1|NIA|NR$|ACX1|OPR|LOX|AOS|AOC|JAZ|MYC2|COI1|WRKY|PAL|PDF1|NPR1|EDS1|"
                r"SID2|ICS|PAD4|PBS3|FMO1|ALD1|SARD|HSP|TCTP|MT[0-9])")


def audit_marker_dict(markers, var_names):
    var = set(map(str, var_names))
    usage = {}
    for ct, genes in markers.items():
        for g in genes:
            usage.setdefault(g, set()).add(ct)
    rows = []
    for ct, genes in markers.items():
        uniq = list(dict.fromkeys(genes))
        rows.append({
            "cell_type": ct, "n_listed": len(genes), "n_unique": len(uniq),
            "n_in_data": sum(g in var for g in uniq),
            "likely_missing_comma": [g for g in uniq if re.match(r"^LOC\d+[A-Za-z]", g)],
            "not_in_data": [g for g in uniq if g not in var][:25],
            "duplicates": [g for g in uniq if genes.count(g) > 1],
            "shared_with": {g: sorted(usage[g] - {ct}) for g in uniq if len(usage[g]) > 1},
            "stress_or_defence_named": [g for g in uniq if re.match(STRESS_REGEX, g)],
        })
    return pd.DataFrame(rows).set_index("cell_type")


def zmean_scores(Xnorm, var_names, markers):
    """Exactly the original rule: per-gene z-score over all spots (NaN->0), then mean over the type's genes."""
    vi = pd.Series(np.arange(len(var_names)), index=pd.Index(map(str, var_names)))
    vi = vi[~vi.index.duplicated()]
    out = {}
    for ct, genes in markers.items():
        cols = [vi[g] for g in dict.fromkeys(genes) if g in vi.index]
        if not cols:
            continue
        E = Xnorm[:, cols]
        E = E.toarray() if sp.issparse(E) else np.asarray(E)
        sd = E.std(0)
        sd[sd == 0] = np.nan
        out[ct] = np.nan_to_num((E - E.mean(0)) / sd).mean(1)
    return pd.DataFrame(out)


def cluster_argmax(scores, clusters):
    clusters = np.asarray(clusters).astype(str)
    lab = {c: scores[clusters == c].mean(0).idxmax() for c in pd.unique(clusters)}
    return np.array([lab[c] for c in clusters], dtype=object)


def null_label_composition(Xnorm, var_names, markers, clusters, n_iter=100, seed=0, n_bins=20):
    """Re-runs the cluster-argmax rule with random, expression-matched gene sets of the *same sizes*.
    Undetected markers (all-zero in every spot) are matched with other undetected genes, so their dilution of the
    mean z-score is reproduced too. If the observed composition resembles this null, labels mostly reflect the
    marker sets (their size and detectability), not the tissue."""
    rng = np.random.default_rng(seed)
    names = np.array(list(map(str, var_names)))
    vi = pd.Series(np.arange(len(names)), index=names)
    vi = vi[~vi.index.duplicated()]
    mean_expr = np.asarray(Xnorm.mean(0)).ravel()
    used = {g for gs in markers.values() for g in gs}
    is_marker = np.isin(names, list(used))
    detected = mean_expr > 0
    bins = np.full(len(names), -1)                       # -1 = never detected
    nb = int(min(n_bins, max(detected.sum(), 1)))
    if detected.any():
        bins[detected] = pd.qcut(pd.Series(mean_expr[detected]).rank(method="first"), nb, labels=False).to_numpy()
    pool = {b: np.where((bins == b) & ~is_marker)[0] for b in np.unique(bins)}
    nonempty = np.array(sorted(b for b, p in pool.items() if len(p)))
    if not len(nonempty):
        raise ValueError("no non-marker genes available to build the null")

    def pick(b):
        if len(pool.get(b, ())) == 0:                    # borrow the nearest bin that has candidates
            b = nonempty[np.argmin(np.abs(nonempty - b))]
        return rng.choice(pool[b])

    found = {ct: [vi[g] for g in dict.fromkeys(gs) if g in vi.index] for ct, gs in markers.items()}
    found = {ct: c for ct, c in found.items() if c}
    shares = []
    for _ in range(n_iter):
        fake = {ct: list(names[[pick(bins[c]) for c in cols]]) for ct, cols in found.items()}
        lab = cluster_argmax(zmean_scores(Xnorm, names, fake), clusters)
        shares.append(pd.Series(lab).value_counts(normalize=True))
    null = pd.concat(shares, axis=1).fillna(0)
    return pd.DataFrame({"n_markers_in_data": {ct: len(c) for ct, c in found.items()},
                         "n_markers_detected": {ct: int(detected[c].sum()) for ct, c in found.items()},
                         "null_share_mean": null.mean(1), "null_share_sd": null.std(1)}).fillna(0)


def marker_holdout_auc(Xnorm, var_names, markers, clusters, n_splits=10, seed=0, min_markers=4):
    """Label clusters with half A of every marker set, then ask whether the held-out half B scores
    those spots higher (one-vs-rest AUC). ~0.5 = labels not supported by independent markers."""
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    vset = set(map(str, var_names))
    present = {ct: [g for g in dict.fromkeys(gs) if g in vset] for ct, gs in markers.items()}
    present = {ct: gs for ct, gs in present.items() if len(gs) >= 2}
    res = {ct: [] for ct in present if len(present[ct]) >= min_markers}
    for _ in range(n_splits):
        A, B = {}, {}
        for ct, gs in present.items():
            perm = list(rng.permutation(gs))
            h = len(perm) // 2
            A[ct], B[ct] = perm[:max(h, 1)], perm[max(h, 1):]
        labA = cluster_argmax(zmean_scores(Xnorm, var_names, A), clusters)
        scB = zmean_scores(Xnorm, var_names, {k: v for k, v in B.items() if v})
        for ct in res:
            y = labA == ct
            if 0 < y.sum() < len(y) and ct in scB:
                res[ct].append(roc_auc_score(y, scB[ct]))
    return pd.DataFrame({ct: {"auc_heldout_markers_mean": np.mean(v) if v else np.nan,
                              "auc_sd": np.std(v) if v else np.nan, "n_splits_labelled": len(v),
                              "n_markers": len(present[ct])} for ct, v in res.items()}).T


# ----------------------------------------------------------------------------------------------
# 5. Are labels confounded with infection?
# ----------------------------------------------------------------------------------------------

def label_condition_test(obs, label_key, sample_key="sample", cond_key="condition", time_key="timepoint",
                         infected="infected"):
    """Section-level quasi-binomial GLM per label: logit(p) ~ condition + time (8 sections = 8 observations)."""
    import statsmodels.api as sm
    from scipy.stats import t as tdist
    tab = pd.crosstab(_as_str(obs[sample_key]), _as_str(obs[label_key]))
    meta = obs.groupby(_as_str(obs[sample_key]))[[cond_key, time_key]].first().loc[tab.index]
    n = tab.sum(1)
    design = pd.DataFrame({"const": 1.0, "infected": (meta[cond_key].astype(str) == infected).astype(float),
                           "time": (meta[time_key].astype(str) == sorted(meta[time_key].astype(str).unique())[-1]).astype(float)})
    rows = []
    for lab in tab.columns:
        k = tab[lab].to_numpy(float)
        y = np.c_[k + 0.5, (n - k).to_numpy(float) + 0.5]
        try:
            fit = sm.GLM(y, design, family=sm.families.Binomial()).fit()
            ntr = y.sum(1)
            mu = np.asarray(fit.fittedvalues)
            x2 = np.sum((y[:, 0] - ntr * mu) ** 2 / (ntr * mu * (1 - mu)))
            df = max(len(ntr) - design.shape[1], 1)
            disp = max(x2 / df, 1.0)                       # quasi-binomial: inflate SE by overdispersion
            lor = fit.params["infected"]
            tval = lor / (fit.bse["infected"] * np.sqrt(disp))
            p = 2 * tdist.sf(abs(tval), df)
        except Exception:
            lor, p, disp = np.nan, np.nan, np.nan
        inf = design["infected"] == 1
        rows.append({"label": lab, "share_mock": (k[~inf] / n[~inf]).mean(), "share_infected": (k[inf] / n[inf]).mean(),
                     "log_odds_infected_vs_mock": lor, "overdispersion": disp, "p_quasibinomial": p})
    out = pd.DataFrame(rows).set_index("label").sort_values("p_quasibinomial")
    return out, (tab.T / n).T.round(3)


# ----------------------------------------------------------------------------------------------
# 6. PlantAG diagnostics and evaluation (numpy GCN = PlantAG architecture)
# ----------------------------------------------------------------------------------------------

def block_dominance(blocks):
    dims = {k: int(v.shape[1]) for k, v in blocks.items()}
    tot = sum(dims.values())
    return pd.DataFrame({"n_dims": dims, "share_after_standardisation": {k: round(d / tot, 4) for k, d in dims.items()}})


def knn_sets(X, k=30, metric="cosine"):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric).fit(X)
    return nn.kneighbors(X, return_distance=False)[:, 1:]


def knn_jaccard(nA, nB):
    return float(np.mean([len(set(a) & set(b)) / len(set(a) | set(b)) for a, b in zip(nA, nB)]))


def knn_graph(X, k=30, metric="cosine", symmetric=True):
    """Weighted kNN graph like PlantAG (weight = cosine similarity), symmetrised."""
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k, metric=metric).fit(X)
    dist, ind = nn.kneighbors(X)
    w = np.clip(1.0 - dist, 1e-6, None) if metric == "cosine" else np.exp(-dist / (dist.mean() + 1e-12))
    n = X.shape[0]
    A = sp.csr_matrix((w.ravel(), (np.repeat(np.arange(n), k), ind.ravel())), shape=(n, n))
    return A.maximum(A.T).tocsr() if symmetric else A


def gcn_norm(A, self_loop_weight=1.0):
    """D^-1/2 (A + I) D^-1/2 as in PyG GCNConv (add_remaining_self_loops with fill value 1)."""
    A = sp.csr_matrix(A, dtype=np.float64).tolil()
    diag = A.diagonal()
    A.setdiag(np.where(diag > 0, diag, self_loop_weight))
    A = A.tocsr()
    d = np.asarray(A.sum(1)).ravel()
    dinv = 1.0 / np.sqrt(np.maximum(d, 1e-12))
    return (sp.diags(dinv) @ A @ sp.diags(dinv)).tocsr()


class NumpyGCN:
    """Full-batch GCN with PlantAG's architecture, trained with Adam (L2 weight decay as torch.optim.Adam)."""

    def __init__(self, n_in, n_classes, hidden=256, n_layers=3, dropout=0.5, lr=0.005, weight_decay=1e-4, seed=0):
        self.rng = np.random.default_rng(seed)
        dims = [n_in] + [hidden] * n_layers
        glorot = lambda a, b: self.rng.uniform(-np.sqrt(6 / (a + b)), np.sqrt(6 / (a + b)), (a, b))
        self.params = [glorot(dims[i], dims[i + 1]) for i in range(n_layers)] + [np.zeros(hidden) for _ in range(n_layers)]
        bound = 1 / np.sqrt(hidden)
        self.params += [self.rng.uniform(-bound, bound, (hidden, n_classes)), self.rng.uniform(-bound, bound, n_classes)]
        self.L, self.p, self.lr, self.wd = n_layers, dropout, lr, weight_decay
        self.m = [np.zeros_like(q) for q in self.params]
        self.v = [np.zeros_like(q) for q in self.params]
        self.t = 0

    def _forward(self, S, X, train):
        L = self.L
        Ws, bs, Wc, bc = self.params[:L], self.params[L:2 * L], self.params[2 * L], self.params[2 * L + 1]
        H, cache = X, []
        for l in range(L):
            HW = H @ Ws[l]
            Z = S @ HW + bs[l]
            A = np.maximum(Z, 0)
            mask = None
            if train and self.p > 0 and l < L - 1:          # PlantAG: dropout after conv1 and conv2 only
                mask = (self.rng.random(A.shape) >= self.p) / (1 - self.p)
                A = A * mask
            cache.append((H, Z, mask))
            H = A
        return H @ Wc + bc, H, cache

    def embed(self, S, X):
        return self._forward(S, X, train=False)[1]

    def predict_proba(self, S, X):
        logits = self._forward(S, X, train=False)[0]
        e = np.exp(logits - logits.max(1, keepdims=True))
        return e / e.sum(1, keepdims=True)

    def _loss_and_grads(self, S, X, y, idx, class_w):
        L = self.L
        logits, H, cache = self._forward(S, X, train=True)
        e = np.exp(logits - logits.max(1, keepdims=True))
        P = e / e.sum(1, keepdims=True)
        w = class_w[y[idx]]
        loss = -(w * np.log(P[idx, y[idx]] + 1e-12)).sum() / w.sum()
        dlog = np.zeros_like(P)
        dlog[idx] = P[idx]
        dlog[idx, y[idx]] -= 1
        dlog[idx] *= (w / w.sum())[:, None]
        Ws, Wc = self.params[:L], self.params[2 * L]
        gW, gb = [None] * L, [None] * L
        gWc, gbc = H.T @ dlog, dlog.sum(0)
        dH = dlog @ Wc.T
        for l in reversed(range(L)):
            Hin, Z, mask = cache[l]
            dA = dH * mask if mask is not None else dH
            dZ = dA * (Z > 0)
            gb[l] = dZ.sum(0)
            StdZ = S.T @ dZ
            gW[l] = Hin.T @ StdZ
            if l > 0:
                dH = StdZ @ Ws[l].T
        grads = gW + gb + [gWc, gbc]
        return loss, [g + self.wd * q for g, q in zip(grads, self.params)]

    def step(self, S, X, y, idx, class_w, b1=0.9, b2=0.999, eps=1e-8):
        loss, grads = self._loss_and_grads(S, X, y, idx, class_w)
        self.t += 1
        for i, g in enumerate(grads):
            self.m[i] = b1 * self.m[i] + (1 - b1) * g
            self.v[i] = b2 * self.v[i] + (1 - b2) * g * g
            mh, vh = self.m[i] / (1 - b1 ** self.t), self.v[i] / (1 - b2 ** self.t)
            self.params[i] -= self.lr * mh / (np.sqrt(vh) + eps)
        return loss


def fit_gcn(S, X, y, train_idx, val_idx=None, n_classes=None, epochs=200, patience=30, seed=0, eval_every=5, **kw):
    """Early stopping on validation loss (checked every `eval_every` epochs; patience in epochs)."""
    n_classes = n_classes or int(y.max()) + 1
    counts = np.bincount(y[train_idx], minlength=n_classes).astype(float)
    class_w = np.where(counts > 0, len(train_idx) / (n_classes * np.maximum(counts, 1)), 0.0)  # sklearn 'balanced'
    model = NumpyGCN(X.shape[1], n_classes, seed=seed, **kw)
    best, best_params, since = np.inf, None, 0
    for epoch in range(1, epochs + 1):
        model.step(S, X, y, train_idx, class_w)
        if val_idx is not None and len(val_idx) and epoch % eval_every == 0:
            P = model.predict_proba(S, X)
            vloss = -np.mean(np.log(P[val_idx, y[val_idx]] + 1e-12))
            if vloss < best - 1e-4:
                best, best_params, since = vloss, [q.copy() for q in model.params], 0
            else:
                since += eval_every
                if since >= patience:
                    break
    if best_params is not None:
        model.params = best_params
    model.epochs_run = epoch
    return model


def sgc_features(S, X, K=2):
    H = X
    for _ in range(K):
        H = S @ H
    return np.asarray(H)


def _logreg():
    from sklearn.linear_model import LogisticRegression
    return LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0)


def _metrics(y_true, y_pred, n_classes):
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
    return {"accuracy": accuracy_score(y_true, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
            "macro_f1": f1_score(y_true, y_pred, average="macro", labels=np.arange(n_classes), zero_division=0)}


def evaluate_models(specs, y, groups, scheme="loso", seed=0, epochs=200, patience=30, test_size=0.2):
    """specs: list of dicts {name, method in {'gcn','sgc','logreg'}, X, S (normalised adjacency or None)}.
    scheme='random' reproduces PlantAG's protocol (stratified 80/20 split of spots, transductive);
    scheme='loso' trains on 7 sections and predicts the held-out one."""
    from sklearn.model_selection import train_test_split
    n_classes = int(y.max()) + 1
    rng = np.random.default_rng(seed)
    if scheme == "random":
        tr, te = train_test_split(np.arange(len(y)), test_size=test_size, stratify=y, random_state=seed)
        folds = [("random_80_20", tr, te)]
    else:
        folds = [(g, np.where(groups != g)[0], np.where(groups == g)[0]) for g in pd.unique(groups)]
    rows, preds = [], {}
    for spec in specs:
        pred_all = np.full(len(y), -1)
        for fold, tr, te in folds:
            if spec["method"] == "gcn":
                val = rng.choice(tr, size=max(1, int(0.1 * len(tr))), replace=False)
                trn = np.setdiff1d(tr, val)
                model = fit_gcn(spec["S"], spec["X"], y, trn, val, n_classes, epochs=epochs, patience=patience, seed=seed)
                p = model.predict_proba(spec["S"], spec["X"]).argmax(1)
            else:
                F = spec["X"] if spec["method"] == "logreg" else sgc_features(spec["S"], spec["X"])
                clf = _logreg().fit(F[tr], y[tr])
                p = np.full(len(y), -1)
                p[te] = clf.predict(F[te])
            pred_all[te] = p[te]
            m = _metrics(y[te], p[te], n_classes)
            rows.append({"model": spec["name"], "scheme": scheme, "fold": fold, "n_test": len(te), **m})
        preds[spec["name"]] = pred_all
    df = pd.DataFrame(rows)
    summary = df.groupby("model")[["accuracy", "balanced_accuracy", "macro_f1"]].agg(["mean", "min", "max"]).round(3)
    return df, summary, preds


def spatial_coherence(labels, A_spatial, sections, n_perm=200, seed=0):
    """Neighbour purity per label vs within-section permutation. Real tissue domains -> ratio >> 1."""
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels).astype(str)
    cats = list(pd.unique(labels))
    A = sp.csr_matrix(A_spatial)
    deg = np.asarray(A.sum(1)).ravel()

    def purity(lab):
        L = sp.csr_matrix((np.ones(len(lab)), (np.arange(len(lab)), pd.Categorical(lab, categories=cats).codes)),
                          shape=(len(lab), len(cats)))
        same = np.asarray((A @ L).multiply(L).sum(1)).ravel()
        return np.array([same[lab == c].sum() / max(deg[lab == c].sum(), 1) for c in cats])

    obs = purity(labels)
    sec = np.asarray(sections).astype(str)
    null = []
    for _ in range(n_perm):
        perm = labels.copy()
        for s in pd.unique(sec):
            m = np.where(sec == s)[0]
            perm[m] = rng.permutation(perm[m])
        null.append(purity(perm))
    null = np.array(null)
    return pd.DataFrame({"n_spots": [int((labels == c).sum()) for c in cats], "neighbour_purity": obs,
                         "null_mean": null.mean(0), "ratio_obs_null": obs / np.maximum(null.mean(0), 1e-9),
                         "z": (obs - null.mean(0)) / np.maximum(null.std(0), 1e-9)}, index=cats).round(3)


def seed_stability(spec, y, train_idx, seeds=(0, 1, 2, 3, 4), epochs=200, patience=30):
    from sklearn.metrics import adjusted_rand_score
    n_classes = int(y.max()) + 1
    preds = []
    for s in seeds:
        rng = np.random.default_rng(s)
        val = rng.choice(train_idx, size=max(1, int(0.1 * len(train_idx))), replace=False)
        trn = np.setdiff1d(train_idx, val)
        model = fit_gcn(spec["S"], spec["X"], y, trn, val, n_classes, epochs=epochs, patience=patience, seed=s)
        preds.append(model.predict_proba(spec["S"], spec["X"]).argmax(1))
    aris = [adjusted_rand_score(a, b) for a, b in itertools.combinations(preds, 2)]
    return float(np.mean(aris)), float(np.min(aris)), preds


# ----------------------------------------------------------------------------------------------
# 7. Robust alternative: spatially-aware tissue domains
# ----------------------------------------------------------------------------------------------

def tissue_domains(adata, counts, sample_key="sample", A_spatial=None, n_hvg=2000, n_pcs=30, lam=0.5,
                   resolution=0.4, use_harmony=True, exclude_organelle=True, seed=0):
    """BANKSY-like domains: [own expression PCs, neighbour-mean PCs] -> kNN -> Leiden."""
    import anndata as ad
    import scanpy as sc
    keep_genes = np.ones(adata.n_vars, bool)
    if exclude_organelle and "organelle" in adata.var:
        keep_genes = adata.var["organelle"].to_numpy() == "nuclear"
    a = ad.AnnData(X=counts[:, keep_genes].astype(np.float32), obs=adata.obs[[sample_key]].copy())
    a.var_names = adata.var_names[keep_genes]
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars), batch_key=sample_key, flavor="seurat")
    sc.pp.pca(a, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)
    rep = "X_pca"
    if use_harmony:
        try:
            import harmonypy
            ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, sample_key, random_state=seed, verbose=False)
            Z = np.asarray(ho.Z_corr)
            a.obsm["X_harmony"] = Z.T if Z.shape[0] != a.n_obs else Z
            rep = "X_harmony"
        except Exception as err:
            warnings.warn(f"Harmony skipped ({err}); domains may follow sections.")
    E = a.obsm[rep]
    if A_spatial is None:
        A_spatial = spatial_adjacency(adata, sample_key)
    deg = np.maximum(np.asarray(A_spatial.sum(1)).ravel(), 1)
    nbr = sp.diags(1 / deg) @ A_spatial @ E
    a.obsm["X_banksy"] = np.hstack([np.sqrt(1 - lam) * E, np.sqrt(lam) * nbr])
    sc.pp.neighbors(a, use_rep="X_banksy", n_neighbors=15, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, key_added="domain", random_state=seed, flavor="igraph",
                 n_iterations=2, directed=False)
    return a.obs["domain"].astype(str).to_numpy(), a.obsm["X_banksy"]


def domain_marker_profile(Xnorm, var_names, markers, domains):
    sc_ = zmean_scores(Xnorm, var_names, markers)
    sc_["domain"] = np.asarray(domains).astype(str)
    return sc_.groupby("domain").mean().round(3)


# ----------------------------------------------------------------------------------------------
# 8. Pseudobulk and replicate-aware differential expression
# ----------------------------------------------------------------------------------------------

def pseudobulk(counts, obs, group_cols, var_names, min_spots=10, meta_cols=()):
    keys = np.array(["|".join(t) for t in zip(*[_as_str(obs[c]) for c in group_cols])], dtype=object)
    codes, uniq = pd.factorize(keys)
    M = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    pb = pd.DataFrame(np.asarray((M @ counts).todense()), index=uniq, columns=list(map(str, var_names)))
    meta = obs.assign(_k=keys).groupby("_k")[list(group_cols) + list(meta_cols)].first().loc[uniq]
    meta["n_spots"] = np.bincount(codes)
    keep = meta["n_spots"] >= min_spots
    return pb.loc[keep], meta.loc[keep]


def run_deseq2(pb, meta, cond_col="condition", ref="control", alt="infected", covariates=(), min_total=10,
               min_samples=2, n_cpus=2):
    """pydeseq2 wrapper (works with the 0.4 and >=0.5 APIs)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    counts = pb.round().astype(int)
    keep = (counts.sum(0) >= min_total) & ((counts > 0).sum(0) >= min_samples)
    counts = counts.loc[:, keep]
    md = meta[[cond_col, *covariates]].astype(str).copy()
    md.index = counts.index
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            dds = DeseqDataSet(counts=counts, metadata=md, design="~" + " + ".join([*covariates, cond_col]),
                               refit_cooks=True, n_cpus=n_cpus, quiet=True)
        except TypeError:
            dds = DeseqDataSet(counts=counts, metadata=md, design_factors=[*covariates, cond_col],
                               ref_level=[cond_col, ref], refit_cooks=True, n_cpus=n_cpus, quiet=True)
        dds.deseq2()
        try:
            ds = DeseqStats(dds, contrast=[cond_col, alt, ref], n_cpus=n_cpus, quiet=True)
        except TypeError:
            ds = DeseqStats(dds, contrast=[cond_col, alt, ref], quiet=True)
        ds.summary()
    return ds.results_df.sort_values("pvalue")


def two_level_de(counts, obs, var_names, time_value, time_col="timepoint", sample_col="sample", leaf_col="leaf",
                 cond_col="condition", domain_col=None, domain_value=None, exclude=None, min_spots=10, **kw):
    """Leaf-level DE (6 vs 6, anti-conservative: leaves share their capture area) and capture-area-level
    DE (2 vs 2, conservative). 'robust' = padj_leaf < 0.05 and p_section < 0.05 with the same sign."""
    m = obs[time_col].astype(str).to_numpy() == str(time_value)
    if domain_col is not None:
        m &= obs[domain_col].astype(str).to_numpy() == str(domain_value)
    if leaf_col in obs:
        m &= ~obs[leaf_col].astype(str).str.endswith("_fragment").to_numpy()
    cols = np.ones(len(var_names), bool) if exclude is None else ~np.asarray(exclude)
    C, o, v = counts[m][:, cols], obs.loc[m], np.asarray(var_names)[cols]
    pb_l, md_l = pseudobulk(C, o, [sample_col, leaf_col], v, min_spots, meta_cols=[cond_col])
    pb_s, md_s = pseudobulk(C, o, [sample_col], v, min_spots, meta_cols=[cond_col])
    info = {"n_leaf_pseudobulks": md_l[cond_col].value_counts().to_dict(),
            "n_section_pseudobulks": md_s[cond_col].value_counts().to_dict(),
            "leaf_pseudobulks_per_section": md_l[sample_col].value_counts().sort_index().to_dict()}
    usable = lambda md: md[cond_col].astype(str).nunique() == 2 and md[cond_col].astype(str).value_counts().min() >= 2
    if not usable(md_l):
        info["skipped"] = ("needs >= 2 leaf pseudobulks (>= %d spots each) in BOTH conditions; check 03_leaf_segmentation "
                           "and the in_tissue flag for the sections listed above" % min_spots)
        return pd.DataFrame(), info
    res_l = run_deseq2(pb_l, md_l, cond_col, **kw)
    out = res_l[["baseMean", "log2FoldChange", "pvalue", "padj"]].add_suffix("_leaf")
    if usable(md_s):
        res_s = run_deseq2(pb_s, md_s, cond_col, **kw)
        out = out.join(res_s[["log2FoldChange", "pvalue", "padj"]].add_suffix("_section"), how="left")
        out["robust"] = ((out["padj_leaf"] < 0.05) & (out["pvalue_section"] < 0.05)
                         & (np.sign(out["log2FoldChange_leaf"]) == np.sign(out["log2FoldChange_section"])))
    else:
        for c in ("log2FoldChange_section", "pvalue_section", "padj_section"):
            out[c] = np.nan
        out["robust"] = False
        info["section_level"] = "skipped: fewer than 2 sections per condition"
    return out.sort_values("pvalue_leaf"), info


LIMMA_R = r'''
# 11b | Optional: nested-design DE with limma-voom + duplicateCorrelation, one fit per time point.
# Leaves are nested in capture areas (sections): leaf-level DESeq2 treats them as independent (anti-conservative),
# section-level DESeq2 has 2 vs 2 (underpowered). duplicateCorrelation estimates how similar leaves of the same
# section are and lmFit uses that correlation, which sits between the two. Inputs: pb_counts.csv, pb_meta.csv (cell 11).
OUT    <- "."     # run from inside the plantag_audit folder:  Rscript limma_voom_dupcor.R
BUNDLE <- file.path(OUT, "bundle")
EXCLUDE_REGEX <- "^TYLCV|^CQW56_"     # viral features + mitochondrial genes; add plastid IDs once confirmed, e.g. "|^LyesC2"

for (p in c("limma", "edgeR")) if (!requireNamespace(p, quietly = TRUE)) {            # first run per session: a few minutes
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(p, update = FALSE, ask = FALSE)
}
suppressMessages({library(limma); library(edgeR)})

cnt_all  <- t(as.matrix(read.csv(file.path(OUT, "pb_counts.csv"), row.names = 1, check.names = FALSE)))
meta_all <- read.csv(file.path(OUT, "pb_meta.csv"), row.names = 1, check.names = FALSE)
use      <- !grepl("_fragment$|^out_of_tissue$", meta_all$leaf)                      # whole leaves only
cnt_all  <- cnt_all[!grepl(EXCLUDE_REGEX, rownames(cnt_all)), use, drop = FALSE]
meta_all <- meta_all[use, , drop = FALSE]

summ <- list()
for (tp in sort(unique(meta_all$timepoint))) {
  m <- meta_all[meta_all$timepoint == tp, , drop = FALSE]
  m$condition <- factor(m$condition, levels = c("control", "infected"))
  n_pb  <- table(m$condition)
  n_sec <- tapply(m$sample, m$condition, function(x) length(unique(x)))
  if (any(is.na(n_sec) | n_sec < 2 | n_pb < 2)) { cat(tp, ": skipped, needs >= 2 sections and >= 2 leaves per condition\n"); next }
  y <- DGEList(cnt_all[, rownames(m), drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, group = m$condition), , keep.lib.sizes = FALSE])
  design <- model.matrix(~ condition, m)
  v   <- voom(y, design)                                                            # round 1
  rho <- duplicateCorrelation(v, design, block = m$sample)$consensus
  v   <- voom(y, design, block = m$sample, correlation = rho)                       # round 2 (standard)
  rho <- duplicateCorrelation(v, design, block = m$sample)$consensus
  rho_used <- if (is.na(rho)) 0 else max(rho, 0)                                    # negative/NA -> leaves treated as independent
  fit <- eBayes(lmFit(v, design, block = m$sample, correlation = rho_used))
  res <- topTable(fit, coef = "conditioninfected", number = Inf, sort.by = "P")
  write.csv(res, file.path(BUNDLE, sprintf("09b_limma_dupcor_%s.csv", tp)))
  summ[[tp]] <- data.frame(timepoint = tp, leaves_control = n_pb[["control"]], leaves_infected = n_pb[["infected"]],
                           sections_control = n_sec[["control"]], sections_infected = n_sec[["infected"]],
                           genes_tested = nrow(res), within_section_correlation = round(rho, 3), correlation_used = round(rho_used, 3),
                           fdr05_up = sum(res$adj.P.Val < 0.05 & res$logFC > 0), fdr05_down = sum(res$adj.P.Val < 0.05 & res$logFC < 0))
  print(summ[[tp]], row.names = FALSE)
}
if (length(summ)) write.csv(do.call(rbind, summ), file.path(BUNDLE, "09b_limma_dupcor_summary.csv"), row.names = FALSE)
'''


# ----------------------------------------------------------------------------------------------
# 9. "Disruption" with proper baselines and nulls
# ----------------------------------------------------------------------------------------------

def _logcpm(pb):
    lib = pb.sum(1).to_numpy(float)
    return np.log1p(pb.to_numpy(float) / lib[:, None] * 1e6)


def _dist(a, b, metric, centre=None):
    if metric == "euclidean":
        return float(np.linalg.norm(a - b))
    if metric == "correlation":                      # 1 - Pearson r between the two mean profiles
        a, b = a - a.mean(), b - b.mean()
    return float(1 - a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


def disruption_with_nulls(pb, meta, group_col, time_col="timepoint", cond_col="condition", section_col="sample",
                          infected="infected", metrics=("cosine_raw", "correlation", "euclidean"),
                          max_perm=5000, seed=0, min_cpm_genes=1.0):
    """For each group (domain or label) x time: distance between infected and mock mean profiles
    (pseudobulk log-CPM per leaf), compared with (i) mock-vs-mock and infected-vs-infected section
    distances and (ii) permutation nulls at leaf level and at section level (2 vs 2 -> 3 splits only)."""
    rng = np.random.default_rng(seed)
    rows = []
    for (g, t), sub in meta.groupby([group_col, time_col], observed=True):
        cond = sub[cond_col].astype(str).to_numpy() == infected
        if cond.sum() < 2 or (~cond).sum() < 2:
            continue
        X = _logcpm(pb.loc[sub.index])
        expressed = np.exp(X).mean(0) - 1 >= min_cpm_genes
        X = X[:, expressed]
        centre = X.mean(0)
        secs = sub[section_col].astype(str).to_numpy()
        for metric in metrics:
            d = lambda i, j: _dist(X[i].mean(0), X[j].mean(0), metric, centre)
            obs_d = d(cond, ~cond)
            n, k = len(cond), int(cond.sum())
            combos = list(itertools.combinations(range(n), k)) if n <= 16 else None
            if combos is not None and len(combos) <= max_perm:
                idxs = combos
            else:
                idxs = [tuple(rng.choice(n, k, replace=False)) for _ in range(max_perm)]
            null = []
            for c in idxs:
                mask = np.zeros(n, bool)
                mask[list(c)] = True
                null.append(d(mask, ~mask))
            null = np.array(null)
            us = pd.unique(secs)
            sec_cond = {s: cond[secs == s][0] for s in us}
            inf_s = [s for s in us if sec_cond[s]]
            mock_s = [s for s in us if not sec_cond[s]]
            within = []
            for pair in (mock_s, inf_s):
                if len(pair) == 2:
                    within.append(d(secs == pair[0], secs == pair[1]))
            sec_null = []
            for comb in itertools.combinations(us, len(inf_s)):
                if 2 * len(inf_s) == len(us) and us[0] not in comb:      # a split and its complement are identical
                    continue
                mask = np.isin(secs, comb)
                if 0 < mask.sum() < n:
                    sec_null.append(d(mask, ~mask))
            rows.append({"group": g, "time": t, "metric": metric, "n_inf_units": k, "n_mock_units": n - k,
                         "observed": obs_d, "within_condition_section_mean": np.mean(within) if within else np.nan,
                         "leaf_perm_null_mean": null.mean(), "leaf_perm_p": (1 + (null >= obs_d - 1e-12).sum()) / (1 + len(null)),
                         "section_level_distinct_splits": len(sec_null),
                         "section_level_p": (int((np.array(sec_null) >= obs_d - 1e-12).sum()) / len(sec_null)) if sec_null else np.nan})
    return pd.DataFrame(rows)


def original_metric_with_baselines(Xmat, obs, label_key, time_col="timepoint", cond_col="condition",
                                   section_col="sample", infected="infected", n_shuffle=200, seed=0):
    """The original spot-level metric (cosine of mean profiles of infected vs mock spots per label),
    next to the same metric for mock-section-1 vs mock-section-2, infected-1 vs infected-2, and a
    spot-level condition shuffle. If within-condition distances are as large, the 'hierarchy' is noise."""
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(Xmat) if sp.issparse(Xmat) else np.asarray(Xmat)
    mean = lambda m: np.asarray(X[np.where(m)[0]].mean(0)).ravel()
    cos = lambda a, b: float(1 - a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
    lab, tim = obs[label_key].astype(str).to_numpy(), obs[time_col].astype(str).to_numpy()
    con, sec = obs[cond_col].astype(str).to_numpy() == infected, obs[section_col].astype(str).to_numpy()
    rows = []
    for l in pd.unique(lab):
        for t in pd.unique(tim):
            m = (lab == l) & (tim == t)
            if m[con].sum() < 3 or m[~con].sum() < 3:
                continue
            obs_d = cos(mean(m & con), mean(m & ~con))
            within = []
            for flag in (False, True):
                ss = [s for s in pd.unique(sec[m & (con == flag)])]
                if len(ss) == 2:
                    a, b = m & (sec == ss[0]), m & (sec == ss[1])
                    if a.sum() >= 3 and b.sum() >= 3:
                        within.append(cos(mean(a), mean(b)))
            idx = np.where(m)[0]
            k = int((m & con).sum())
            shuf = []
            for _ in range(n_shuffle):
                pick = np.zeros(len(con), bool)
                pick[rng.choice(idx, k, replace=False)] = True
                shuf.append(cos(mean(pick), mean(m & ~pick)))
            rows.append({"label": l, "time": t, "n_infected_spots": k, "n_mock_spots": int((m & ~con).sum()),
                         "cosine_infected_vs_mock": obs_d,
                         "cosine_mock1_vs_mock2": within[0] if len(within) > 0 else np.nan,
                         "cosine_inf1_vs_inf2": within[1] if len(within) > 1 else np.nan,
                         "spot_shuffle_null_mean": float(np.mean(shuf))})
    return pd.DataFrame(rows).round(4)


# ----------------------------------------------------------------------------------------------
# 10. Viral reads (after re-running Space Ranger with TYLCV in the reference)
# ----------------------------------------------------------------------------------------------

def viral_summary(adata, counts, viral_features, sample_key="sample", cond_key="condition", infected="infected",
                  bg_quantile=0.999):
    """Viral UMIs per spot; mock sections define the background (index hopping / ambient)."""
    from scipy.spatial import cKDTree
    vi = [list(map(str, adata.var_names)).index(f) for f in viral_features]
    v = np.asarray(counts[:, vi].sum(1)).ravel()
    tot = np.asarray(counts.sum(1)).ravel()
    adata.obs["viral_umis"] = v
    adata.obs["viral_per_10k"] = v / np.maximum(tot, 1) * 1e4
    mock = adata.obs[cond_key].astype(str).to_numpy() != infected
    thr = max(1.0, float(np.quantile(v[mock], bg_quantile))) if mock.any() else 1.0
    adata.obs["viral_positive"] = v > thr
    xy = physical_coords(adata.obs)
    dist = np.full(adata.n_obs, np.nan)
    samples = _as_str(adata.obs[sample_key])
    for s in pd.unique(samples):
        m = samples == s
        pos = m & adata.obs["viral_positive"].to_numpy()
        if pos.any():
            dist[m] = cKDTree(xy[pos]).query(xy[m])[0]
    adata.obs["dist_to_viral_um"] = dist
    tab = adata.obs.groupby(samples).agg(n_spots=("viral_umis", "size"), viral_umis=("viral_umis", "sum"),
                                         pct_positive=("viral_positive", lambda s: 100 * s.mean()),
                                         max_umis=("viral_umis", "max"))
    tab["background_threshold"] = thr
    return tab.round(3)


def paired_viral_de(counts, obs, var_names, leaf_col="leaf", status_col="viral_positive", min_spots=8, **kw):
    """Within infected leaves: viral-positive vs viral-negative spots, paired by leaf (design ~ leaf + status)."""
    o = obs.copy()
    o["status"] = np.where(o[status_col].to_numpy(), "pos", "neg")
    pb, md = pseudobulk(counts, o, [leaf_col, "status"], var_names, min_spots)
    both = md.groupby(leaf_col, observed=True)["status"].nunique()
    md = md[md[leaf_col].isin(both[both == 2].index)]
    if md[leaf_col].nunique() < 2:
        raise ValueError("Need >= 2 infected leaves with both viral-positive and -negative spots.")
    return run_deseq2(pb.loc[md.index], md, "status", ref="neg", alt="pos", covariates=(leaf_col,), **kw)


# ----------------------------------------------------------------------------------------------
# 11. RNA velocity, fixed (barcode matching, correct loom per section, within-section smoothing,
#     replicate-aware condition tests and a permutation null)
# ----------------------------------------------------------------------------------------------

PAPER_VELOCITY_GENES = ["Rbcs-3B", "rbcS1", "PSBR", "CAB7", "ca2", "PsbQ", "SlFBA1", "cTPT",
                        "LOC101263969", "LOC101249777", "LOC543974"]


def assign_looms(loom_summary, loom_corr, loom_map, min_corr=0.9):
    """Pick, for every section, the loom whose per-barcode totals correlate with that section's counts
    (decided by the data, not the file name). Sections without a convincing match are excluded."""
    rows, mapping = [], {}
    for s in loom_summary.index:
        best = loom_summary.loc[s].get("loom_that_best_matches_counts")
        key = (f"adata:{s}", f"loom:{best}")
        r = float(loom_corr.loc[key]) if best is not None and key[1] in loom_corr.columns else np.nan
        ok = best is not None and np.isfinite(r) and r >= min_corr
        mapping[s] = loom_map[best] if ok else None
        rows.append({"sample": s, "loom_file_used": mapping[s], "count_correlation": round(r, 4) if np.isfinite(r) else r,
                     "reassigned_from_filename": bool(ok and best != s), "excluded": not ok})
    return mapping, pd.DataFrame(rows).set_index("sample")


def integrate_velocity_by_barcode(adata, loom_for_sample, loom_dir, sample_key="sample"):
    """Spliced/unspliced counts aligned to adata spot-by-spot via the 16-nt barcode, from the loom that
    belongs to each section (see assign_looms). Replaces the original positional mapping."""
    genes = pd.Index(list(map(str, adata.var_names)))
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    gpos = gpos[~gpos.index.duplicated()]
    parts = {"spliced": ([], [], []), "unspliced": ([], [], [])}
    report, samples = {}, _as_str(adata.obs[sample_key])
    for s, fn in loom_for_sample.items():
        rows = np.where(samples == s)[0]
        if fn is None or len(rows) == 0:
            report[s] = {"loom": fn, "spots": len(rows), "matched_by_barcode": 0}
            continue
        lo = _read_loom(Path(loom_dir) / fn, genes=genes.to_numpy())
        gsel = pd.Index(lo["genes_sel"])
        keepg = np.where(~gsel.duplicated())[0]
        cols = gpos[gsel[keepg]].to_numpy()
        lb = pd.Series(np.arange(len(lo["ids"])), index=barcode16(lo["ids"]))
        lb = lb[~lb.index.duplicated()]
        bc = barcode16(adata.obs_names[rows])
        ok = pd.Index(bc).isin(lb.index)
        src = lb.loc[bc[ok]].to_numpy()
        for layer in parts:
            M = lo[layer][src][:, keepg].tocoo()
            parts[layer][0].append(rows[ok][M.row])
            parts[layer][1].append(cols[M.col])
            parts[layer][2].append(M.data)
        report[s] = {"loom": fn, "spots": len(rows), "matched_by_barcode": int(ok.sum())}
    shape = (adata.n_obs, len(genes))
    out = []
    for layer in ("spliced", "unspliced"):
        r, c, d = parts[layer]
        out.append(sp.csr_matrix((np.concatenate(d) if d else [], (np.concatenate(r) if r else [], np.concatenate(c) if c else [])),
                                 shape=shape, dtype=np.float32))
    return out[0], out[1], pd.DataFrame(report).T


def velocity_input_qc(adata, S, U, counts, sample_key="sample", thresholds=(5, 10, 20, 30)):
    """How much unspliced signal is there? Also checks that integrated totals track the spot counts."""
    from scipy.stats import spearmanr
    s_tot, u_tot = np.asarray(S.sum(1)).ravel(), np.asarray(U.sum(1)).ravel()
    c_tot = np.asarray(counts.sum(1)).ravel()
    samples = _as_str(adata.obs[sample_key])
    rows = []
    for smp in pd.unique(samples):
        m = samples == smp
        has = (s_tot[m] + u_tot[m]) > 0
        rows.append({"sample": smp, "spots": int(m.sum()), "spots_with_loom_data": int(has.sum()),
                     "median_spliced_umis": float(np.median(s_tot[m])), "median_unspliced_umis": float(np.median(u_tot[m])),
                     "pct_spots_with_unspliced": 100 * float((u_tot[m] > 0).mean()),
                     "unspliced_fraction": float(u_tot[m].sum() / max(s_tot[m].sum() + u_tot[m].sum(), 1)),
                     "spearman_loom_total_vs_counts": float(spearmanr(s_tot[m] + u_tot[m], c_tot[m])[0]) if has.sum() > 10 else np.nan})
    nz = (S > 0).multiply(U > 0)
    shared = np.asarray(nz.multiply(S + U).sum(0)).ravel()       # scVelo's definition of 'shared counts'
    genes_ok = {f"genes_with_shared_counts>={t}": int((shared >= t).sum()) for t in thresholds}
    return pd.DataFrame(rows).set_index("sample").round(3), pd.Series(genes_ok)


def within_section_neighbors(v, emb, sample_key="sample", n_neighbors=30, seed=0):
    """kNN graph computed separately in every section, so moments never average infected with mock spots."""
    import anndata as ad
    import scanpy as sc
    samples = _as_str(v.obs[sample_key])
    Ds, Cs, order = [], [], []
    for s in pd.unique(samples):
        idx = np.where(samples == s)[0]
        sub = ad.AnnData(X=np.zeros((len(idx), 1), np.float32))
        sub.obsm["emb"] = np.asarray(emb[idx], dtype=np.float32)
        sc.pp.neighbors(sub, use_rep="emb", n_neighbors=int(min(n_neighbors, len(idx) - 1)), random_state=seed)
        Ds.append(sub.obsp["distances"])
        Cs.append(sub.obsp["connectivities"])
        order.append(idx)
    perm = np.concatenate(order)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    v.obsp["distances"] = sp.block_diag(Ds).tocsr()[inv][:, inv].tocsr()
    v.obsp["connectivities"] = sp.block_diag(Cs).tocsr()[inv][:, inv].tocsr()
    v.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances",
                          "params": {"n_neighbors": int(n_neighbors), "method": "umap", "metric": "euclidean"}}


def run_velocity(adata, S, U, sample_key="sample", n_top_genes=2000, min_shared_counts=10, n_neighbors=30,
                 n_pcs=30, mode="stochastic", nuclear_only=True, seed=0):
    """scVelo on barcode-matched layers. Neighbours: within-section kNN on PCA of the velocity genes'
    expression (the original used the PlantAG embedding across all sections, i.e. circular labels and
    infected/mock mixing). Organellar genes are excluded (they are not spliced and dominate the counts)."""
    import anndata as ad
    import scanpy as sc
    import scvelo as scv
    scv.settings.verbosity = 1
    has = (np.asarray(S.sum(1)).ravel() + np.asarray(U.sum(1)).ravel()) > 0
    genes = np.ones(adata.n_vars, bool)
    if nuclear_only and "organelle" in adata.var:
        genes = adata.var["organelle"].to_numpy() == "nuclear"
    v = ad.AnnData(X=S[has][:, genes].copy(), obs=adata.obs.loc[has].copy())
    v.var_names = adata.var_names[genes]
    v.layers["spliced"] = S[has][:, genes].copy()
    v.layers["unspliced"] = U[has][:, genes].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import inspect
        kw = {"min_shared_counts": min_shared_counts}
        if "log" in inspect.signature(scv.pp.filter_and_normalize).parameters:   # scVelo < 0.3.3 also logs X
            kw["log"] = False
        scv.pp.filter_and_normalize(v, **kw)                                     # filters + normalises layers
        sc.pp.log1p(v)
        if n_top_genes and v.n_vars > n_top_genes:
            sc.pp.highly_variable_genes(v, n_top_genes=n_top_genes, flavor="seurat")
            v = v[:, v.var["highly_variable"].to_numpy()].copy()
        sc.pp.pca(v, n_comps=int(min(n_pcs, v.n_vars - 1, v.n_obs - 1)), random_state=seed)
        g = sc.pp.neighbors(v, n_neighbors=15, random_state=seed, copy=True)      # global graph: display only
        sc.tl.umap(g, random_state=seed)
        v.obsm["X_umap"] = g.obsm["X_umap"]
        within_section_neighbors(v, v.obsm["X_pca"], sample_key, n_neighbors, seed)
        scv.pp.moments(v, n_neighbors=None, n_pcs=None)
        scv.tl.velocity(v, mode=mode)
        scv.tl.velocity_graph(v, n_jobs=1, show_progress_bar=False)
        scv.tl.velocity_confidence(v)
    return v


def velocity_gene_table(v, cond_key="condition", infected="infected"):
    vg = v.var["velocity_genes"].to_numpy(bool) if "velocity_genes" in v.var else np.ones(v.n_vars, bool)
    V = np.asarray(v.layers["velocity"])
    inf = v.obs[cond_key].astype(str).to_numpy() == infected
    tab = pd.DataFrame({"velocity_gene": vg,
                        "gamma": v.var.get("velocity_gamma", pd.Series(np.nan, index=v.var_names)).to_numpy(),
                        "fit_r2": v.var.get("velocity_r2", pd.Series(np.nan, index=v.var_names)).to_numpy(),
                        "mean_velocity_all": np.nanmean(V, 0),
                        "mean_velocity_infected": np.nanmean(V[inf], 0) if inf.any() else np.nan,
                        "mean_velocity_mock": np.nanmean(V[~inf], 0) if (~inf).any() else np.nan,
                        "mean_Ms": np.asarray(v.layers["Ms"]).mean(0), "mean_Mu": np.asarray(v.layers["Mu"]).mean(0)},
                       index=v.var_names)
    return tab.sort_values("fit_r2", ascending=False)


def unspliced_fraction_test(S, U, obs, var_names, unit_key="leaf", cond_key="condition", time_key="timepoint",
                            infected="infected", genes=None, min_unspliced=20, min_spliced=50):
    """Calibrated, expression-independent test of altered splicing kinetics (replaces condition tests on velocity).
    Per gene and time point: raw unspliced (u) and spliced (s) counts summed per unit (leaf or section), then a
    quasi-binomial model of the unspliced fraction u/(u+s) ~ condition with units as replicates.
    If kinetics are unchanged the fraction does not move when expression goes up or down, whereas raw velocity
    (u - gamma*s) scales with expression and per-cell normalisation of the sparse unspliced layer distorts ratios.
    log_odds_diff > 0: higher unspliced fraction in infected (e.g. induction faster than steady state)."""
    from scipy.stats import t as tdist
    from statsmodels.stats.multitest import multipletests
    names = np.asarray(list(map(str, var_names)))
    cols = np.arange(len(names)) if genes is None else np.where(np.isin(names, list(map(str, genes))))[0]
    units = _as_str(obs[unit_key])
    valid = ~(pd.Series(units).str.endswith("_fragment") | pd.Series(units).isin(["out_of_tissue", "nan"])).to_numpy()
    codes, uniq = pd.factorize(units)
    rows = np.where(valid)[0]
    M = sp.csr_matrix((np.ones(len(rows)), (codes[rows], rows)), shape=(len(uniq), len(units)))
    s_sum = np.asarray((M @ sp.csr_matrix(S)[:, cols]).todense(), float)
    u_sum = np.asarray((M @ sp.csr_matrix(U)[:, cols]).todense(), float)
    meta = obs.assign(_u=units).groupby("_u")[[cond_key, time_key]].first().loc[uniq]
    used = np.bincount(codes[rows], minlength=len(uniq)) > 0
    out = []
    for t in sorted(meta[time_key].astype(str).unique()):
        sel = (meta[time_key].astype(str).to_numpy() == t) & used
        inf = meta[cond_key].astype(str).to_numpy()[sel] == infected
        if inf.sum() < 2 or (~inf).sum() < 2:
            continue
        u, s = u_sum[sel] + 0.5, s_sum[sel] + 0.5
        n = u + s
        ok = (u_sum[sel].sum(0) >= min_unspliced) & (s_sum[sel].sum(0) >= min_spliced)
        Ui, Si, Um, Sm = u[inf].sum(0), s[inf].sum(0), u[~inf].sum(0), s[~inf].sum(0)
        lor = np.log(Ui / Si) - np.log(Um / Sm)
        se = np.sqrt(1 / Ui + 1 / Si + 1 / Um + 1 / Sm)
        fi, fm = Ui / (Ui + Si), Um / (Um + Sm)
        phat = np.where(inf[:, None], fi[None, :], fm[None, :])
        x2 = ((u - n * phat) ** 2 / (n * phat * (1 - phat))).sum(0)
        df = int(sel.sum()) - 2
        disp = np.maximum(x2 / max(df, 1), 1.0)
        p = 2 * tdist.sf(np.abs(lor / (se * np.sqrt(disp))), max(df, 1))
        padj = np.full(len(p), np.nan)
        if ok.any():
            padj[ok] = multipletests(p[ok], method="fdr_bh")[1]
        out.append(pd.DataFrame({"time": t, "unit": unit_key, "n_infected_units": int(inf.sum()),
                                 "n_mock_units": int((~inf).sum()), "unspliced_frac_infected": fi,
                                 "unspliced_frac_mock": fm, "log_odds_diff": lor, "dispersion": disp,
                                 "tested": ok, "p": np.where(ok, p, np.nan), "padj_bh": padj},
                                index=pd.Index(names[cols], name="gene")))
    if not out:
        return pd.DataFrame()
    res = pd.concat(out)
    return res[res["tested"]].drop(columns="tested").sort_values("p")


def velocity_confidence_null(v, sample_key="sample", n_perm=3, seed=0, mode="stochastic"):
    """Same genes, graph and spliced data, but unspliced values shuffled among spots within each section
    (independently per gene). If observed confidence ~ null, the velocity field is not supported."""
    import scvelo as scv
    rng = np.random.default_rng(seed)
    samples = _as_str(v.obs[sample_key])
    obs_conf = v.obs["velocity_confidence"].to_numpy(float)
    null_meds, null_all = [], []
    for _ in range(n_perm):
        w = v.copy()
        Uu = w.layers["unspliced"]
        Uu = Uu.toarray() if sp.issparse(Uu) else np.array(Uu)
        for s in pd.unique(samples):
            idx = np.where(samples == s)[0]
            order = np.argsort(rng.random((len(idx), Uu.shape[1])), axis=0)
            Uu[idx] = np.take_along_axis(Uu[idx], order, axis=0)
        w.layers["unspliced"] = sp.csr_matrix(Uu)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scv.pp.moments(w, n_neighbors=None, n_pcs=None)
            scv.tl.velocity(w, mode=mode)
            scv.tl.velocity_confidence(w)
        c = w.obs["velocity_confidence"].to_numpy(float)
        null_meds.append(np.nanmedian(c))
        null_all.append(c)
    null_all = np.concatenate(null_all)
    q95 = np.nanquantile(null_all, 0.95)
    per_sec = pd.DataFrame({"observed_median_confidence": pd.Series(obs_conf).groupby(samples).median(),
                            "pct_spots_above_null_q95": pd.Series(obs_conf > q95).groupby(samples).mean() * 100})
    summary = {"observed_median_confidence": float(np.nanmedian(obs_conf)),
               "null_median_confidence_mean": float(np.mean(null_meds)), "null_q95": float(q95),
               "pct_spots_above_null_q95": float(100 * np.nanmean(obs_conf > q95)),
               "n_velocity_genes": int(v.var["velocity_genes"].sum()) if "velocity_genes" in v.var else None}
    return summary, per_sec.round(3)


# ----------------------------------------------------------------------------------------------
# 12. Viral counts from the quick route (count_viral_umis.py) -> extra features in the matrix
# ----------------------------------------------------------------------------------------------

def attach_viral_counts(adata, counts, csv_for_sample, sample_key="sample", features=("TYLCV_V", "TYLCV_C")):
    """Adds per-spot viral UMI columns (from count_viral_umis.py CSVs: barcode,TYLCV_V,TYLCV_C) as new
    features. Sections without a CSV get zeros and are reported."""
    import anndata as ad
    extra = np.zeros((adata.n_obs, len(features)), np.float32)
    samples = _as_str(adata.obs[sample_key])
    tissue = {t: pd.Index(barcode16(adata.obs_names[samples == t])) for t in pd.unique(samples)}
    report = {}
    for s in pd.unique(samples):
        rows = np.where(samples == s)[0]
        path = csv_for_sample.get(s)
        if path is None or not Path(path).exists():
            report[s] = {"csv": None, "matched_spots": 0}
            continue
        tab = pd.read_csv(path)
        tab["bc16"] = barcode16(tab["barcode"])
        tab = tab.dropna(subset=["bc16"]).groupby("bc16")[list(features)].sum()
        bc = pd.Index(barcode16(adata.obs_names[rows]))
        hit = bc.isin(tab.index)
        extra[rows[hit]] = tab.loc[bc[hit], list(features)].to_numpy(np.float32)
        # identity check: whose tissue do this file's viral UMIs fall on? (barcodes encode array position)
        tot = tab[list(features)].sum(1)
        frac = {t: float(tot[tot.index.isin(ix)].sum() / tot.sum()) for t, ix in tissue.items()} if tot.sum() > 0 else {}
        best = max(frac, key=frac.get) if frac else None
        status = "n/a (no viral UMIs)" if not frac else ("ok" if frac[best] <= frac[s] + 0.05 else f"MISMATCH: looks like {best}")
        if status.startswith("MISMATCH"):
            warnings.warn(f"{Path(path).name}: viral UMIs fall mostly on the tissue of {best} ({100 * frac[best]:.0f}%) "
                          f"rather than {s} ({100 * frac[s]:.0f}%). Check which BAM/run belongs to which section.")
        report[s] = {"csv": str(path), "matched_spots": int(hit.sum()), "viral_umis": float(extra[rows].sum()),
                     "pct_csv_umis_on_this_tissue": round(100 * frac.get(s, np.nan), 1) if frac else np.nan,
                     "identity_check": status}
    new_counts = sp.hstack([counts, sp.csr_matrix(extra)]).tocsr()
    var = pd.concat([adata.var, pd.DataFrame(index=list(features))])
    if "organelle" in var:
        var.loc[list(features), "organelle"] = "viral"
    new = ad.AnnData(X=new_counts, obs=adata.obs.copy(), var=var)
    new.obsm = adata.obsm.copy()
    return new, new_counts, pd.DataFrame(report).T


# ----------------------------------------------------------------------------------------------
# 13. Viral counts from the merged-reference Space Ranger re-runs -> per-section CSVs for cell 14
# ----------------------------------------------------------------------------------------------

def viral_csvs_from_spaceranger(sr_dir, raw, counts, out_dir, sample_key="sample", viral_prefix="TYLCV",
                                min_corr=0.95, h5_name="outs/filtered_feature_bc_matrix.h5"):
    """Reads every <sr_dir>/<run>/outs/filtered_feature_bc_matrix.h5, takes the TYLCV_V / TYLCV_C columns and
    writes <out_dir>/<section>.csv (barcode,TYLCV_V,TYLCV_C,total), the format cell 14 attaches.
    Each run is assigned to the existing section whose per-barcode tomato totals it reproduces (Spearman),
    so folder names do not matter and crossed labels (e.g. the 7-dpi infected pair) are caught.
    A run is skipped if its best correlation is < min_corr or if two runs claim the same section."""
    import scanpy as sc
    from scipy.stats import spearmanr
    sr_dir, out_dir = Path(sr_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tot = np.asarray(counts.sum(1)).ravel()
    samples = _as_str(raw.obs[sample_key])
    old = {s: pd.Series(tot[samples == s], index=barcode16(raw.obs_names[samples == s])).groupby(level=0).sum()
           for s in pd.unique(samples)}
    rows, tables = [], {}
    for h5 in sorted(sr_dir.glob(f"*/{h5_name}")):
        run = h5.relative_to(sr_dir).parts[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = sc.read_10x_h5(str(h5))
        names = pd.Index(a.var_names.astype(str))
        vmask = names.str.startswith(viral_prefix)
        X = sp.csr_matrix(a.X)
        bc = barcode16(a.obs_names)
        host = pd.Series(np.asarray(X[:, np.where(~vmask)[0]].sum(1)).ravel(), index=bc).groupby(level=0).sum()
        corr = {}
        for s, o in old.items():
            common = o.index.intersection(host.index)
            corr[s] = float(spearmanr(o[common], host[common])[0]) if len(common) >= 20 else np.nan
        ranked = sorted([(v, s) for s, v in corr.items() if np.isfinite(v)], reverse=True)
        best_r, best = ranked[0] if ranked else (np.nan, None)
        V = pd.DataFrame(X[:, np.where(vmask)[0]].toarray(), columns=list(names[vmask]))
        for f in ("TYLCV_V", "TYLCV_C"):
            if f not in V:
                V[f] = 0
        V.insert(0, "barcode", list(a.obs_names.astype(str)))
        V["total"] = V["TYLCV_V"] + V["TYLCV_C"]
        matched = int(pd.Index(bc).isin(old[best].index).sum()) if best else 0
        rows.append({"run": run, "section": best, "host_count_spearman": round(best_r, 4),
                     "second_best": f"{ranked[1][1]} ({ranked[1][0]:.3f})" if len(ranked) > 1 else "",
                     "spots_in_run": a.n_obs, "spots_matched": matched, "viral_umis": int(V["total"].sum()),
                     "pct_spots_with_virus": round(100 * float((V["total"] > 0).mean()), 2),
                     "viral_features": ",".join(names[vmask])})
        tables[run] = V[["barcode", "TYLCV_V", "TYLCV_C", "total"]]
    tab = pd.DataFrame(rows)
    if tab.empty:
        raise FileNotFoundError(f"no {h5_name} found under {sr_dir}")
    dup = tab["section"].duplicated(keep=False) & tab["section"].notna()
    tab["status"] = np.select([tab["host_count_spearman"] < min_corr, dup], ["low_correlation", "duplicate_section"], "ok")
    for _, r in tab[tab["status"] == "ok"].iterrows():
        tables[r["run"]].to_csv(out_dir / f"{r['section']}.csv", index=False)
    missing = sorted(set(old) - set(tab.loc[tab["status"] == "ok", "section"]))
    if missing:
        warnings.warn(f"sections without a matching re-run: {missing}")
    return tab.set_index("run")


# ----------------------------------------------------------------------------------------------
# 14. Where do the viral UMIs land? (AnnData spots vs the BAM's own tissue spots)
# ----------------------------------------------------------------------------------------------

SECTION_CODES = {"control_7dpi_1": "N0524", "control_7dpi_2": "N0525", "control_24dpi_1": "N0036",
                 "control_24dpi_2": "N0037", "infected_7dpi_1": "N0523", "infected_7dpi_2": "N0522",
                 "infected_24dpi_1": "N0034", "infected_24dpi_2": "N0035"}   # same pairing as make_bams_tsv.sh


def viral_barcode_diagnostic(raw, loom_dir, csv_dir, section_codes=None, sample_key="sample", n_array_spots=4992):
    """For each section, compares the barcodes in viral_counts/<section>.csv with
      (a) that section's spots in the AnnData (all spots, before any filtering), and
      (b) the tissue spots Space Ranger reported for the SAME BAM: the barcodes velocyto used, read from the loom
          whose name carries the same facility code (N0522...).
    Reading the table:
      pct_anndata_barcodes_in_loom ~100  -> the AnnData section really is that BAM's capture area.
      pct_viral_on_loom_tissue >> pct_expected_if_random -> viral reads come from tissue: route B is fine.
      pct_viral_on_loom_tissue >> pct_viral_on_anndata_spots -> the AnnData barcodes do not match that BAM.
      pct_viral_on_loom_tissue ~ pct_expected_if_random -> viral reads are not tissue-associated: check route B."""
    import h5py
    section_codes = section_codes or SECTION_CODES
    samples = _as_str(raw.obs[sample_key])
    looms = sorted(Path(loom_dir).glob("*.loom"))
    rows = []
    for s, code in section_codes.items():
        ad_bc = pd.Index(barcode16(raw.obs_names[samples == s])).dropna()
        lf = [p for p in looms if re.search(rf"(^|[^0-9]){code}([^0-9]|$)", p.name)]
        lo_bc = pd.Index([])
        if len(lf) == 1:
            with h5py.File(lf[0], "r") as f:
                ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["col_attrs"]["CellID"][:]]
            lo_bc = pd.Index(barcode16(ids)).dropna()
        v = pd.Series(dtype=float)
        csv = Path(csv_dir) / f"{s}.csv"
        if csv.exists():
            t = pd.read_csv(csv)
            t["bc"] = barcode16(t["barcode"])
            v = t.dropna(subset=["bc"]).groupby("bc")["total"].sum()
        tot = float(v.sum())
        pct = lambda ix: round(100 * float(v[v.index.isin(ix)].sum()) / tot, 1) if tot and len(ix) else np.nan
        rows.append({"section": s, "code": code, "loom": lf[0].name if len(lf) == 1 else f"{len(lf)} looms match",
                     "anndata_spots": len(ad_bc), "loom_tissue_spots": len(lo_bc),
                     "pct_anndata_barcodes_in_loom": round(100 * float(ad_bc.isin(lo_bc).mean()), 1)
                     if len(ad_bc) and len(lo_bc) else np.nan,
                     "viral_umis_in_csv": int(tot), "pct_viral_on_anndata_spots": pct(ad_bc),
                     "pct_viral_on_loom_tissue": pct(lo_bc),
                     "pct_expected_if_random": round(100 * len(lo_bc) / n_array_spots, 1) if len(lo_bc) else np.nan})
    return pd.DataFrame(rows).set_index("section")


# ----------------------------------------------------------------------------------------------
# 15. Cell-type-resolved analyses with the PlantAG labels: viral tropism and splicing kinetics
# ----------------------------------------------------------------------------------------------

def transfer_labels(target_obs, source_obs, cols=("predicted_cell_type",), sample_key="sample"):
    """Copies label columns from the annotated object to the analysis object: by obs name, else by
    (section, 16-nt barcode). Unlabelled spots get "nan" (skipped downstream).
    Returns (DataFrame aligned to target_obs.index, fraction of spots that received a label)."""
    cols = list(cols)
    out = pd.DataFrame({c: pd.Series(_as_str(source_obs[c]), index=source_obs.index).reindex(target_obs.index)
                        for c in cols}, index=target_obs.index)
    labelled = lambda df: float((df[cols[0]].notna() & (df[cols[0]] != "nan")).mean())
    if labelled(out) < 0.9:
        bt = pd.Series(barcode16(target_obs.index), index=target_obs.index)
        bs = pd.Series(barcode16(source_obs.index), index=source_obs.index)
        kt = pd.Index(_as_str(target_obs[sample_key])) + "|" + pd.Index(_as_str(bt))
        ks = pd.Index(_as_str(source_obs[sample_key])) + "|" + pd.Index(_as_str(bs))
        m = pd.DataFrame({c: _as_str(source_obs[c]) for c in cols}, index=ks)[bs.notna().to_numpy()]
        m = m[~m.index.duplicated()]
        alt = pd.DataFrame({c: m[c].reindex(kt).to_numpy() for c in cols}, index=target_obs.index)
        alt.loc[bt.isna().to_numpy(), cols] = np.nan            # no barcode: never matched by the fallback
        if labelled(alt) > labelled(out):
            out = alt
    out = pd.DataFrame({c: _as_str(out[c]) for c in cols}, index=target_obs.index)
    return out, labelled(out)


def viral_tropism_by_type(obs, viral, host_total, type_key="predicted_cell_type", sample_key="sample",
                          cond_key="condition", infected="infected", time_key="timepoint", leaf_key="leaf",
                          n_perm=5000, seed=0, min_spots=5):
    """Do some cell types carry more on-tissue viral UMIs than their sequencing depth predicts?
    In each infected section: O = viral UMIs on spots of type k; E = viral UMIs of the section x (host UMIs on type k /
    host UMIs of the section); O/E > 1 means enrichment beyond depth. Viral RNA spread evenly over spots (ambient RNA,
    index hopping) gives O/E ~ 1 for every type: it dilutes a real tropism but cannot create one.
    Significance: labels permuted among spots within each section (each spot keeps its viral and host UMIs), statistic
    recomputed; sections pooled by summing O and E (stratified permutation), per time point and over all infected
    sections. Spots are not independent (labels and viral foci are both spatially clustered), so permutation P values
    are optimistic: read them with n_sections_enriched and n_leaves_enriched, which count replicate-level agreement."""
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(seed)
    inf = _as_str(obs[cond_key]) == infected
    o = obs.loc[inf]
    v = np.asarray(viral, float)[inf]
    t = np.asarray(host_total, float)[inf]
    sec = _as_str(o[sample_key])
    tp = _as_str(o[time_key])
    lf = _as_str(o[leaf_key]) if leaf_key in o else np.array(["-"] * len(o))
    typ = _as_str(o[type_key])
    keep = typ != "nan"
    v, t, sec, tp, lf, typ = v[keep], t[keep], sec[keep], tp[keep], lf[keep], typ[keep]
    types = np.array(sorted(pd.unique(typ), key=str), dtype=object)
    K = len(types)
    code = pd.Categorical(typ, categories=types).codes
    sections = [s for s in pd.unique(sec) if v[sec == s].sum() > 0]
    idx = {s: np.where(sec == s)[0] for s in sections}

    def stats(codes):
        O = np.zeros((len(sections), K)); E = np.zeros((len(sections), K))
        for j, s in enumerate(sections):
            ii = idx[s]
            c = codes[ii]
            O[j] = np.bincount(c, weights=v[ii], minlength=K)
            E[j] = v[ii].sum() * np.bincount(c, weights=t[ii], minlength=K) / max(t[ii].sum(), 1.0)
        return O, E

    O, E = stats(code)
    scopes = {s: [j] for j, s in enumerate(sections)}
    for time in sorted(pd.unique(tp)):
        js = [j for j, s in enumerate(sections) if tp[idx[s][0]] == time]
        if len(js) > 1:
            scopes[f"pooled_{time}"] = js
    if len(sections) > 1:
        scopes["pooled_all"] = list(range(len(sections)))
    lor = lambda O_, E_: np.log((O_ + 0.5) / (E_ + 0.5))
    obs_stat = {sc: lor(O[js].sum(0), E[js].sum(0)) for sc, js in scopes.items()}
    ge = {sc: np.zeros(K) for sc in scopes}
    le = {sc: np.zeros(K) for sc in scopes}
    for _ in range(n_perm):
        perm = code.copy()
        for s in sections:
            ii = idx[s]
            perm[ii] = code[ii][rng.permutation(len(ii))]
        Op, Ep = stats(perm)
        for sc, js in scopes.items():
            st = lor(Op[js].sum(0), Ep[js].sum(0))
            ge[sc] += st >= obs_stat[sc] - 1e-12
            le[sc] += st <= obs_stat[sc] + 1e-12
    # replicate-level agreement: sections and leaves with O/E > 1 for the type
    sec_enr = (O > E) & (O + E > 0)
    leaf_rows = []
    for (s, l), ii in pd.Series(np.arange(len(v))).groupby([sec, lf]):
        ii = ii.to_numpy()
        if v[ii].sum() <= 0 or l.endswith("_fragment") or l in ("-", "nan", "out_of_tissue"):
            continue
        c = code[ii]
        Ol = np.bincount(c, weights=v[ii], minlength=K)
        El = v[ii].sum() * np.bincount(c, weights=t[ii], minlength=K) / max(t[ii].sum(), 1.0)
        n_k = np.bincount(c, minlength=K)
        leaf_rows.append((s, tp[ii[0]], Ol > El, n_k >= min_spots))
    rows = []
    for sc, js in scopes.items():
        Os, Es = O[js].sum(0), E[js].sum(0)
        p_up = (ge[sc] + 1) / (n_perm + 1)
        p_dn = (le[sc] + 1) / (n_perm + 1)
        p2 = np.minimum(1.0, 2 * np.minimum(p_up, p_dn))
        in_scope = [sections[j] for j in js]
        n_spots = np.bincount(code[np.isin(sec, in_scope)], minlength=K)
        lr = [r for r in leaf_rows if r[0] in in_scope]
        for k in range(K):
            if n_spots[k] < min_spots:
                continue
            rows.append({"scope": sc, "cell_type": types[k], "n_spots": int(n_spots[k]),
                         "viral_umis_observed": float(Os[k]), "viral_umis_expected": round(float(Es[k]), 2),
                         "obs_over_exp": round(float((Os[k] + 0.5) / (Es[k] + 0.5)), 3),
                         "p_enriched": float(p_up[k]), "p_depleted": float(p_dn[k]), "p_two_sided": float(p2[k]),
                         "n_sections": len(js), "n_sections_enriched": int(sec_enr[js, k].sum()),
                         "n_leaves_evaluable": int(sum(r[3][k] for r in lr)),
                         "n_leaves_enriched": int(sum(r[2][k] and r[3][k] for r in lr))})
    out = pd.DataFrame(rows)
    if len(out):
        out["padj_bh"] = np.nan
        for sc, g in out.groupby("scope"):
            out.loc[g.index, "padj_bh"] = multipletests(g["p_two_sided"], method="fdr_bh")[1]
    return out


def unspliced_fraction_by_type(S, U, obs, var_names, type_key="predicted_cell_type", unit_key="leaf",
                               min_spots_per_unit=5, **kw):
    """Splicing kinetics per PlantAG cell type: the calibrated unspliced-fraction test (quasi-binomial, replicate units,
    never spots) run separately on the spots of each type. Units with fewer than min_spots_per_unit spots of the type
    are dropped; a type is tested at a time point only with >= 2 infected and >= 2 mock units. padj_bh corrects within
    each type and time point; padj_all corrects across every type x gene x time test."""
    from statsmodels.stats.multitest import multipletests
    S, U = sp.csr_matrix(S), sp.csr_matrix(U)
    types = _as_str(obs[type_key])
    res = []
    for k in sorted(pd.unique(types), key=str):
        if k == "nan":
            continue
        rows = np.where(types == k)[0]
        if len(rows) < 4 * min_spots_per_unit:
            continue
        o = obs.iloc[rows].copy()
        u = pd.Series(_as_str(o[unit_key]))
        small = u.map(u.value_counts()).to_numpy() < min_spots_per_unit
        o[unit_key] = np.where(small, "nan", u.to_numpy())
        r = unspliced_fraction_test(S[rows], U[rows], o, var_names, unit_key, **kw)
        if len(r):
            res.append(r.assign(cell_type=k, n_spots_type=len(rows)))
    if not res:
        return pd.DataFrame()
    out = pd.concat(res)
    out["padj_all"] = multipletests(out["p"].to_numpy(), method="fdr_bh")[1]
    return out.sort_values("p")


# ----------------------------------------------------------------------------------------------
# 16. Topology of local expression neighbourhoods: H1 persistence, infected vs mock (capture areas as replicates)
# ----------------------------------------------------------------------------------------------

def h1_persistence_by_spot(X, k=30, method="ripser"):
    """Vietoris-Rips persistence of each spot's local neighbourhood (the spot and its k nearest neighbours in X), as in
    PlantAG. Returns (DataFrame with one row per spot, list of (H0, H1) diagrams). Columns: n_h1 (number of loops),
    total_h1 and max_h1 (summed and largest persistence), local_scale (mean distance to the k neighbours) and
    total_h1_scaled = total_h1 / local_scale, which removes the dependence of persistence on neighbourhood size.
    method="plantag" reproduces PlantAG's compute_local_diagrams + clean_diagrams exactly: k points including the spot
    itself, a Rips complex truncated at the distance to the farthest of them (GUDHI), and infinite deaths replaced by
    birth + the largest finite death of that spot's diagram. method="ripser" (default) uses the spot and its k neighbours
    without truncation."""
    if method == "plantag":
        return _h1_plantag(np.asarray(X, float), k)
    from sklearn.neighbors import NearestNeighbors
    try:
        from ripser import ripser

        def dgm(P):
            d = ripser(P, maxdim=1)["dgms"]
            return d[0], d[1]
    except ImportError:
        import gudhi

        def dgm(P):
            st = gudhi.RipsComplex(points=P).create_simplex_tree(max_dimension=2)
            st.compute_persistence()
            h0 = np.asarray(st.persistence_intervals_in_dimension(0)).reshape(-1, 2)
            h1 = np.asarray(st.persistence_intervals_in_dimension(1)).reshape(-1, 2)
            return h0, h1
    X = np.asarray(X, float)
    dist, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    rows, diagrams = [], []
    for i in range(len(X)):
        h0, h1 = dgm(X[idx[i]])
        pers = h1[:, 1] - h1[:, 0] if len(h1) else np.zeros(0)
        scale = float(dist[i, 1:].mean())
        rows.append({"n_h1": int(len(pers)), "total_h1": float(pers.sum()), "max_h1": float(pers.max()) if len(pers) else 0.0,
                     "local_scale": scale, "total_h1_scaled": float(pers.sum() / scale) if scale > 0 else np.nan})
        diagrams.append((h0, h1))
    return pd.DataFrame(rows), diagrams


def _h1_plantag(X, k=30):
    import gudhi as gd
    from sklearn.neighbors import NearestNeighbors
    dist, idx = NearestNeighbors(n_neighbors=k, algorithm="ball_tree").fit(X).kneighbors(X)
    rows, diagrams = [], []
    for i in range(len(X)):
        st = gd.RipsComplex(points=X[idx[i]], max_edge_length=np.max(dist[i])).create_simplex_tree(max_dimension=2)
        pers = st.persistence()
        finite = [d for _, (b, d) in pers if np.isfinite(d)]
        top = max(finite) if finite else 1.0
        dg = {0: [], 1: []}
        for dim, (b, d) in pers:
            d = b + top if not np.isfinite(d) else d
            if d > b and dim in dg:
                dg[dim].append((b, d))
        h0, h1 = (np.asarray(dg[q], float).reshape(-1, 2) for q in (0, 1))
        p1 = h1[:, 1] - h1[:, 0] if len(h1) else np.zeros(0)
        scale = float(dist[i, 1:].mean())
        rows.append({"n_h1": int(len(p1)), "total_h1": float(p1.sum()), "max_h1": float(p1.max()) if len(p1) else 0.0,
                     "local_scale": scale, "total_h1_scaled": float(p1.sum() / scale) if scale > 0 else np.nan})
        diagrams.append((h0, h1))
    return pd.DataFrame(rows), diagrams


def h1_persistence_test(spot_stats, obs, sample_key="sample", cond_key="condition", time_key="timepoint",
                        infected="infected", leaf_key=None, depth=None,
                        stats=("total_h1_scaled", "n_h1", "max_h1", "total_h1")):
    """Infected vs mock with capture areas as replicates. Spot statistics are averaged per capture area and compared by
    (i) OLS on capture-area means, value ~ infected + time point (8 areas, 5 residual df); (ii) an exact permutation of
    condition labels within time points (36 relabellings for 2 + 2 areas per time point, so P >= 2/36); and (iii) if
    leaves are given, a leaf-level mixed model with a random intercept per capture area. Each statistic is tested raw
    and, if depth is given, after regressing log10 UMIs per spot out across all spots. The primary statistic is the first
    in `stats`; the others are secondary."""
    import itertools
    import statsmodels.formula.api as smf
    d = spot_stats.reset_index(drop=True).copy()
    d["area"] = _as_str(obs[sample_key]); d["time"] = _as_str(obs[time_key])
    d["is_inf"] = (_as_str(obs[cond_key]) == infected).astype(int)
    if leaf_key:
        d["leaf"] = _as_str(obs[leaf_key])
    variants = [("raw", None)]
    if depth is not None:
        variants.append(("depth_adjusted", np.log10(np.asarray(depth, float) + 1)))
    rows, area_tables = [], []
    for st in stats:
        for vname, ld in variants:
            y = d[st].astype(float).to_numpy()
            if ld is not None:
                ok = np.isfinite(y) & np.isfinite(ld)
                coef = np.polyfit(ld[ok], y[ok], 1)
                y = y - np.polyval(coef, ld) + np.nanmean(y)
            d["_y"] = y
            A = d.groupby("area").agg(y=("_y", "mean"), is_inf=("is_inf", "first"), time=("time", "first"),
                                      n_spots=("_y", "size")).reset_index()
            fit = smf.ols("y ~ is_inf + C(time)", data=A).fit()

            def contrast(lab):
                diffs = [g.y[lab[g.index] == 1].mean() - g.y[lab[g.index] == 0].mean()
                         for _, g in A.groupby("time") if len(set(lab[g.index])) == 2]
                return float(np.mean(diffs))
            lab0 = A["is_inf"].to_numpy()
            obs_c = contrast(lab0)
            groups = [list(itertools.combinations(g.index, int(g["is_inf"].sum()))) for _, g in A.groupby("time")]
            perm = []
            for combo in itertools.product(*groups):
                lab = np.zeros(len(A), int)
                lab[[i for c in combo for i in c]] = 1
                perm.append(contrast(lab))
            perm = np.asarray(perm)
            p_mixed = np.nan
            if leaf_key:
                L = d[~d["leaf"].str.endswith("_fragment") & ~d["leaf"].isin(["nan", "out_of_tissue"])]
                L = L.groupby(["area", "leaf"]).agg(y=("_y", "mean"), is_inf=("is_inf", "first"), time=("time", "first")).reset_index()
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        p_mixed = float(smf.mixedlm("y ~ is_inf + C(time)", L, groups=L["area"]).fit(reml=True).pvalues["is_inf"])
                except Exception:
                    pass
            mock_mean = float(A.y[A.is_inf == 0].mean())
            rows.append({"statistic": st, "adjustment": vname, "primary": st == stats[0],
                         "mock_mean": mock_mean, "infected_minus_mock": float(fit.params["is_inf"]),
                         "relative_change_pct": 100 * float(fit.params["is_inf"]) / mock_mean if mock_mean else np.nan,
                         "p_capture_area_ols": float(fit.pvalues["is_inf"]), "within_time_contrast": obs_c,
                         "p_exact_permutation": float(np.mean(np.abs(perm) >= abs(obs_c) - 1e-12)),
                         "n_permutations": int(len(perm)), "p_leaf_mixed_model": p_mixed})
            area_tables.append(A.assign(statistic=st, adjustment=vname))
    return pd.DataFrame(rows), pd.concat(area_tables, ignore_index=True)


def plot_h1_panels(spot_stats, diagrams, obs, area_table, test_table, cond_key="condition", time_key="timepoint",
                   infected="infected", stat="total_h1_scaled", time_value=None):
    """Representative persistence diagrams (the infected and the mock spot closest to their group's median `stat` at one
    time point) and capture-area means of `stat` with the replicate-level P values. Returns the figure."""
    cond = _as_str(obs[cond_key]); tim = _as_str(obs[time_key])
    tsort = lambda ts: sorted(ts, key=lambda t: (int(re.sub(r"\D", "", t) or 0), t))      # 7dpi before 24dpi
    time_value = time_value or tsort(pd.unique(tim))[-1]
    v = spot_stats[stat].to_numpy()
    picks = {}
    for lab, is_inf in (("TYLCV-infected", True), ("mock", False)):
        m = np.where(((cond == infected) == is_inf) & (tim == time_value) & np.isfinite(v))[0]
        picks[lab] = m[np.argmin(np.abs(v[m] - np.median(v[m])))]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.0))
    finite_max = lambda a: float(a[np.isfinite(a)].max()) if np.isfinite(a).any() else 1.0
    top = 1.08 * max(max(finite_max(diagrams[i][0][:, 1]), finite_max(diagrams[i][1][:, 1]) if len(diagrams[i][1]) else 0) for i in picks.values())
    for ax, (lab, i) in zip(axes[:2], picks.items()):
        h0, h1 = diagrams[i]
        h0 = h0[np.isfinite(h0[:, 1])]
        ax.plot([0, top], [0, top], ls="--", lw=0.7, color="grey")
        ax.scatter(h0[:, 0], h0[:, 1], s=14, color="#1F77B4", label="H$_0$")
        ax.scatter(h1[:, 0], h1[:, 1], s=14, color="#D62728" if lab.startswith("TYLCV") else "#7F7F7F", label="H$_1$")
        ax.set_xlim(-0.03 * top, top); ax.set_ylim(-0.03 * top, top)
        ax.set_xlabel("Birth ε"); ax.set_ylabel("Death ε"); ax.legend(frameon=False, loc="lower right", fontsize=7)
        ax.set_title(f"{lab} spot ({time_value}; median {stat})", fontsize=8)
    ax = axes[2]
    A = area_table[(area_table.statistic == stat) & (area_table.adjustment == "raw")]
    for (t, inf_), g in A.groupby(["time", "is_inf"]):
        x = tsort(pd.unique(A.time)).index(t) + (0.15 if inf_ else -0.15)
        ax.scatter(np.full(len(g), x), g.y, s=22, color="#D62728" if inf_ else "#1F77B4", label=("TYLCV" if inf_ else "mock"))
    h, l = ax.get_legend_handles_labels(); keep = dict(zip(l, h))
    ax.legend(keep.values(), keep.keys(), frameon=False, fontsize=7)
    ax.set_xticks(range(len(pd.unique(A.time)))); ax.set_xticklabels(tsort(pd.unique(A.time)))
    r = test_table[(test_table.statistic == stat) & (test_table.adjustment == "raw")].iloc[0]
    ax.set_ylabel(f"{stat}, capture-area mean")
    ax.set_title(f"TYLCV vs mock: {r.relative_change_pct:+.1f}%\nP = {r.p_capture_area_ols:.3f} (areas), {r.p_exact_permutation:.3f} (exact)", fontsize=8)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# 17. Excluding spots that are not leaf tissue (e.g. displaced material under a section)
# ----------------------------------------------------------------------------------------------

def keep_listed_spots(adata, csv_path, sample_key="sample"):
    """Keeps, for each section listed in csv_path (columns: section, spot), only the listed spots; sections that are not
    listed are left whole. Spots are matched by obs name, else by (section, 16-nt barcode).
    Returns (filtered AnnData, summary table, figure of kept and removed spots in the affected sections)."""
    tab = pd.read_csv(csv_path)
    assert {"section", "spot"} <= set(tab.columns), "the CSV needs the columns: section, spot"
    samples = _as_str(adata.obs[sample_key])
    names = pd.Index(adata.obs_names.astype(str))
    bc = pd.Series(barcode16(adata.obs_names), index=adata.obs_names).astype(str)
    drop = np.zeros(adata.n_obs, bool)
    rows = []
    for sec, g in tab.groupby("section"):
        in_sec = samples == str(sec)
        if not in_sec.any():
            raise ValueError(f"section {sec!r} from {csv_path} is not in the data")
        listed = pd.Index(g["spot"].astype(str))
        keep = in_sec & names.isin(listed)
        if keep.sum() < len(listed):                          # fall back to barcodes for spots not matched by name
            want = set(pd.Series(barcode16(listed)).dropna().astype(str))
            keep |= in_sec & bc.isin(want).to_numpy()
        if keep.sum() == 0:
            raise ValueError(f"none of the {len(listed)} spots listed for {sec!r} matched the data (names or barcodes)")
        drop |= in_sec & ~keep
        rows.append({"section": sec, "spots_before": int(in_sec.sum()), "spots_listed": len(listed),
                     "spots_kept": int(keep.sum()), "spots_removed": int((in_sec & ~keep).sum())})
    summary = pd.DataFrame(rows).set_index("section")
    fig = _plot_kept_removed(adata, samples, drop, list(summary.index))
    return adata[~drop].copy(), summary, fig


def propose_block_exclusion(adata, counts, section, sample_key="sample", k=2, min_size=200, depth_quantile=0.10):
    """Fallback when no list of kept spots exists. In `section`, finds regions thicker than about 2k + 1 spots
    (morphological opening of the spot lattice: k erosions, k dilations) with at least min_size spots: a filled block,
    unlike the narrow strips of transverse leaf sections. Spots inside the block are proposed for removal unless their
    UMI count reaches the depth_quantile of the section's strip spots, i.e. they look like leaf tissue lying over the
    block. Returns (DataFrame with columns section, spot to keep, figure). Check the figure against the tissue image."""
    from scipy.sparse.csgraph import connected_components
    samples = _as_str(adata.obs[sample_key])
    idx = np.where(samples == section)[0]
    sub = adata[idx]
    A = spatial_adjacency(sub, sample_key, 101.0).tocsr()
    full = np.ones(len(idx), bool)
    cur = full.copy()
    for _ in range(k):                                         # erosion: keep spots whose 6 lattice neighbours are all present
        deg = np.asarray(A[:, cur].sum(1)).ravel()
        cur = cur & (deg >= 6)
    for _ in range(k):                                         # dilation within the section
        cur = cur | (np.asarray(A[:, cur].sum(1)).ravel() > 0) & full
    n, lab = connected_components(A[cur][:, cur], directed=False)
    block = np.zeros(len(idx), bool)
    sizes = np.bincount(lab) if n else np.zeros(0)
    pos = np.where(cur)[0]
    for c in np.where(sizes >= min_size)[0]:
        block[pos[lab == c]] = True
    umis = np.asarray(counts[idx].sum(1)).ravel()
    thr = np.quantile(umis[~block], depth_quantile) if (~block).any() else np.inf
    remove = block & (umis < thr)
    keep_names = sub.obs_names[~remove]
    drop = np.zeros(adata.n_obs, bool); drop[idx[remove]] = True
    fig = _plot_kept_removed(adata, samples, drop, [section],
                             title=f"{section}: {int(remove.sum())} spots proposed for removal (block, < {thr:.0f} UMIs)")
    return pd.DataFrame({"section": section, "spot": keep_names.astype(str)}), fig


def _plot_kept_removed(adata, samples, drop, sections, title=None):
    fig, axes = plt.subplots(1, len(sections), figsize=(4.2 * len(sections), 4.0), squeeze=False)
    for ax, sec in zip(axes[0], sections):
        m = samples == sec
        xy = physical_coords(adata.obs[m])
        d = drop[m]
        ax.scatter(xy[~d, 0], xy[~d, 1], s=6, color="#1f77b4", label=f"kept ({int((~d).sum())})")
        ax.scatter(xy[d, 0], xy[d, 1], s=6, color="#d9d9d9", label=f"removed ({int(d.sum())})")
        ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title or sec, fontsize=8); ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# 18. Interactive spot editor (Colab / Jupyter): choose by hand which spots of a section are leaf tissue
# ----------------------------------------------------------------------------------------------

def _section_image(adata, section, sample_key="sample", library_id=None):
    """Tissue image of `section` from adata.uns['spatial'] (Space Ranger / Squidpy layout).
    Returns (image as uint8 array, scale from full-resolution pixels to that image, spot diameter in full-resolution
    pixels or None, library id) or None, with a message saying why."""
    sp_all = adata.uns.get("spatial", {}) if hasattr(adata, "uns") else {}
    if not sp_all:
        return None, "no adata.uns['spatial'] (no tissue images stored)"
    samples = _as_str(adata.obs[sample_key])
    obs = adata.obs[samples == section]
    lib = library_id
    if lib is None:
        for col in ("library_id", "library", "batch", sample_key):
            if col in obs and str(obs[col].iloc[0]) in sp_all:
                lib = str(obs[col].iloc[0]); break
    if lib is None:
        hits = [k for k in sp_all if section in str(k)]
        lib = hits[0] if len(hits) == 1 else (next(iter(sp_all)) if len(sp_all) == 1 else None)
    if lib is None or lib not in sp_all:
        return None, f"no image matched to {section!r}; pass library_id=<one of {list(sp_all)}>"
    d = sp_all[lib]
    imgs, sf = d.get("images", {}), d.get("scalefactors", {})
    key = "hires" if "hires" in imgs else ("lowres" if "lowres" in imgs else None)
    if key is None:
        return None, f"library {lib!r} has no hires/lowres image"
    img = np.asarray(imgs[key])
    if img.dtype != np.uint8:
        img = np.clip(img * (255 if img.max() <= 1.0 else 1), 0, 255).astype(np.uint8)
    if img.ndim == 2:
        img = np.stack([img] * 3, -1)
    return (img[..., :3], float(sf.get(f"tissue_{key}_scalef", 1.0)), sf.get("spot_diameter_fullres"), lib), f"image of library {lib!r}"


def spot_editor(adata, counts, section, out_csv, sample_key="sample", library_id=None, start="proposal",
                k=2, min_size=200, max_px=1000, return_html=False):
    """Interactive editor for the spots of one section. Spots are drawn over the section's tissue image when
    adata.uns['spatial'] holds it (else on the spot grid), starting from the automatic proposal with the thick block
    removed entirely (start='proposal') or with every spot kept (start='all').
    Click a spot to toggle it (blue = kept, red ring = removed); drag a loop to apply the chosen mode (remove / keep)
    to every spot inside it; Undo reverts the last action. 'Save to Drive' writes out_csv (columns: section, spot) in
    Colab; 'Download CSV' saves the same file through the browser anywhere."""
    import base64, io, json, uuid
    from PIL import Image
    samples = _as_str(adata.obs[sample_key])
    idx = np.where(samples == section)[0]
    if not len(idx):
        raise ValueError(f"section {section!r} not found")
    names = adata.obs_names[idx].astype(str).tolist()
    prop, fig = propose_block_exclusion(adata, counts, section, sample_key, k=k, min_size=min_size, depth_quantile=1.0)
    plt.close(fig)
    keep0 = np.isin(names, prop["spot"].astype(str).to_numpy()) if start == "proposal" else np.ones(len(idx), bool)
    res, note = _section_image(adata, section, sample_key, library_id)
    img_b64, xy, r = None, None, None
    if res is not None and "spatial" in adata.obsm:
        img, scalef, spot_d, lib = res
        full = np.asarray(adata.obsm["spatial"])[idx][:, :2].astype(float)
        for cand in (full, full[:, ::-1]):                       # [x, y] expected; try [y, x] if spots fall off the image
            p = cand * scalef
            inside = (p[:, 0] >= 0) & (p[:, 0] < img.shape[1]) & (p[:, 1] >= 0) & (p[:, 1] < img.shape[0])
            if inside.mean() > 0.95:
                xy = p; break
        if xy is not None:
            nn = np.sort(np.linalg.norm(xy[:, None] - xy[None, :min(len(xy), 400)], axis=-1), axis=0)[1]
            r = (spot_d * scalef / 2) if spot_d else 0.42 * float(np.median(nn))
            m = 6 * r
            x0, y0 = int(max(0, xy[:, 0].min() - m)), int(max(0, xy[:, 1].min() - m))
            x1, y1 = int(min(img.shape[1], xy[:, 0].max() + m)), int(min(img.shape[0], xy[:, 1].max() + m))
            crop = Image.fromarray(img[y0:y1, x0:x1])
            f = min(1.0, max_px / max(crop.size))
            if f < 1.0:
                crop = crop.resize((int(crop.size[0] * f), int(crop.size[1] * f)), Image.LANCZOS)
            buf = io.BytesIO(); crop.save(buf, format="PNG"); img_b64 = base64.b64encode(buf.getvalue()).decode()
            xy = (xy - [x0, y0]) * f; r = r * f
            W, H = crop.size
        else:
            note = f"{note}, but spot coordinates did not fall on it; drawing the spot grid instead"
    if xy is None:                                               # spot grid (array coordinates, 100-um pitch)
        pc = physical_coords(adata.obs.iloc[idx])
        pc = pc - pc.min(0)
        f = max_px / max(pc.max(0).max(), 1.0)
        xy = pc * f + 40; r = 27.5 * f
        W, H = int(xy[:, 0].max() + 40), int(xy[:, 1].max() + 40)
    uid = uuid.uuid4().hex[:8]
    data = json.dumps({"x": np.round(xy[:, 0], 2).tolist(), "y": np.round(xy[:, 1], 2).tolist(), "r": round(float(r), 2),
                       "names": names, "keep": keep0.astype(bool).tolist(), "section": section})
    image_tag = (f'<image id="im-{uid}" href="data:image/png;base64,{img_b64}" x="0" y="0" width="{W}" height="{H}" '
                 f'style="pointer-events:none"/>') if img_b64 else ""
    html = f"""
<div id="ed-{uid}" style="font-family:Arial,sans-serif;font-size:13px;line-height:1.5">
 <div><b>{section}</b>: click a spot to switch it between <span style="color:#1f77b4">kept</span> and
  <span style="color:#d62728">removed</span>; drag a loop around spots to apply the mode below to all of them. ({note})</div>
 <div style="margin:6px 0">Loop mode:
  <label><input type="radio" name="mode-{uid}" value="remove" checked> remove</label>
  <label><input type="radio" name="mode-{uid}" value="keep"> keep</label>
  &nbsp;<button id="undo-{uid}">Undo</button> <button id="prop-{uid}">Start from proposal</button>
  <button id="all-{uid}">Keep all</button>
  {'<label><input type="checkbox" id="imgbox-' + uid + '" checked> tissue image</label>' if img_b64 else ''}
  &nbsp;<button id="save-{uid}" style="font-weight:bold">Save to Drive</button> <button id="dl-{uid}">Download CSV</button>
  &nbsp;<span id="cnt-{uid}"></span> <span id="msg-{uid}" style="color:#2a7f2a"></span></div>
 <svg id="svg-{uid}" viewBox="0 0 {W} {H}" style="width:100%;max-width:{max_px}px;border:1px solid #bbb;cursor:crosshair;user-select:none;background:#fff">
  {image_tag}<g id="g-{uid}" style="pointer-events:none"></g>
  <polyline id="lasso-{uid}" points="" fill="rgba(214,39,40,0.10)" stroke="#d62728" stroke-width="{max(1.0, r / 4):.2f}" style="pointer-events:none"/>
 </svg>
</div>
<script>
(function () {{
  const D = {data};
  const uid = "{uid}", NS = "http://www.w3.org/2000/svg", n = D.x.length;
  const svg = document.getElementById("svg-" + uid), g = document.getElementById("g-" + uid);
  const lasso = document.getElementById("lasso-" + uid), cnt = document.getElementById("cnt-" + uid), msg = document.getElementById("msg-" + uid);
  let keep = D.keep.slice(); const hist = [];
  const circles = D.x.map(function (x, i) {{
    const c = document.createElementNS(NS, "circle");
    c.setAttribute("cx", x); c.setAttribute("cy", D.y[i]); c.setAttribute("r", D.r); g.appendChild(c); return c; }});
  function paint() {{
    let k = 0;
    circles.forEach(function (c, i) {{
      if (keep[i]) {{ k++; c.setAttribute("fill", "#1f77b4"); c.setAttribute("fill-opacity", "0.85"); c.setAttribute("stroke", "none"); }}
      else {{ c.setAttribute("fill", "none"); c.setAttribute("stroke", "#d62728"); c.setAttribute("stroke-width", Math.max(0.6, D.r * 0.28)); }} }});
    cnt.textContent = "kept " + k + " / removed " + (n - k);
  }}
  function toSvg(e) {{ const p = svg.createSVGPoint(); p.x = e.clientX; p.y = e.clientY; return p.matrixTransform(svg.getScreenCTM().inverse()); }}
  function inPoly(px, py, poly) {{
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {{
      const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
      if (((yi > py) !== (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) inside = !inside; }}
    return inside;
  }}
  let drawing = false, path = [], down = null, travel = 0;
  svg.addEventListener("mousedown", function (e) {{ e.preventDefault(); drawing = true; travel = 0; const p = toSvg(e); path = [[p.x, p.y]]; down = [e.clientX, e.clientY]; }});
  svg.addEventListener("mousemove", function (e) {{
    if (!drawing) return; const p = toSvg(e); path.push([p.x, p.y]);
    travel = Math.max(travel, Math.hypot(e.clientX - down[0], e.clientY - down[1]));
    lasso.setAttribute("points", path.map(function (q) {{ return q[0] + "," + q[1]; }}).join(" ")); }});
  window.addEventListener("mouseup", function (e) {{
    if (!drawing) return; drawing = false; lasso.setAttribute("points", "");
    const before = keep.slice();
    if (travel < 4) {{                                                          // a click (a loop ends where it began)
      let best = -1, bd = Infinity;
      for (let i = 0; i < n; i++) {{ const d = Math.hypot(D.x[i] - path[0][0], D.y[i] - path[0][1]); if (d < bd) {{ bd = d; best = i; }} }}
      if (best >= 0 && bd <= 1.6 * D.r) keep[best] = !keep[best];
    }} else if (path.length > 2) {{                                             // a loop: apply the mode inside it
      const k = document.querySelector('input[name="mode-' + uid + '"]:checked').value === "keep";
      for (let i = 0; i < n; i++) if (inPoly(D.x[i], D.y[i], path)) keep[i] = k;
    }}
    if (keep.some(function (v, i) {{ return v !== before[i]; }})) hist.push(before);
    paint();
  }});
  document.getElementById("undo-" + uid).onclick = function () {{ if (hist.length) {{ keep = hist.pop(); paint(); }} }};
  document.getElementById("prop-" + uid).onclick = function () {{ hist.push(keep.slice()); keep = D.keep.slice(); paint(); }};
  document.getElementById("all-" + uid).onclick = function () {{ hist.push(keep.slice()); keep = keep.map(function () {{ return true; }}); paint(); }};
  const box = document.getElementById("imgbox-" + uid), im = document.getElementById("im-" + uid);
  if (box && im) box.onchange = function () {{ im.style.display = box.checked ? "" : "none"; }};
  function keptNames() {{ return D.names.filter(function (_, i) {{ return keep[i]; }}); }}
  function csvText() {{ return "section,spot\\n" + keptNames().map(function (s) {{ return D.section + "," + s; }}).join("\\n") + "\\n"; }}
  document.getElementById("dl-" + uid).onclick = function () {{
    const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([csvText()], {{type: "text/csv"}}));
    a.download = D.section + "_kept_spots.csv"; document.body.appendChild(a); a.click(); a.remove(); msg.textContent = "CSV downloaded"; }};
  document.getElementById("save-" + uid).onclick = async function () {{
    if (!(window.google && google.colab && google.colab.kernel)) {{ msg.textContent = "not in Colab: use Download CSV"; return; }}
    msg.textContent = "saving...";
    try {{
      const res = await google.colab.kernel.invokeFunction("plantag.save_kept_spots", [D.section, keptNames()], {{}});
      const j = res && res.data && res.data["application/json"]; msg.textContent = (j && j.msg) || "saved";
    }} catch (err) {{ msg.textContent = "save failed (" + err + "): use Download CSV"; }}
  }};
  paint();
}})();
</script>"""
    if return_html:
        return html
    from IPython.display import HTML, JSON, display
    try:
        from google.colab import output

        def _save(sec, kept):
            pd.DataFrame({"section": sec, "spot": list(kept)}).to_csv(out_csv, index=False)
            return JSON({"msg": f"saved {len(kept)} kept spots to {out_csv}"})
        output.register_callback("plantag.save_kept_spots", _save)
    except ImportError:
        pass
    display(HTML(html))


def h1_panels_data(spot_stats, diagrams, obs, area_table, cond_key="condition", time_key="timepoint", infected="infected",
                   stat="total_h1_scaled", time_value=None):
    """Data behind Fig. 1b,c: persistence diagrams of the infected and the mock spot closest to their group's median `stat`
    at one time point (the same rule as plot_h1_panels), plus the capture-area means of `stat`."""
    cond = _as_str(obs[cond_key]); tim = _as_str(obs[time_key])
    tsort = lambda ts: sorted(ts, key=lambda t: (int(re.sub(r"\D", "", t) or 0), t))
    time_value = time_value or tsort(pd.unique(tim))[-1]
    v = spot_stats[stat].to_numpy()
    rows = []
    for grp, is_inf in (("TYLCV", True), ("mock", False)):
        m = np.where(((cond == infected) == is_inf) & (tim == time_value) & np.isfinite(v))[0]
        i = m[np.argmin(np.abs(v[m] - np.median(v[m])))]
        for dim, D in enumerate(diagrams[i]):
            rows += [{"part": "diagram", "group": grp, "spot": str(obs.index[i]), "dim": dim, "x": float(b), "y": float(d)} for b, d in D]
    A = area_table[(area_table.statistic == stat) & (area_table.adjustment == "raw")]
    rows += [{"part": "area_mean", "group": "TYLCV" if r.is_inf else "mock", "spot": r.area, "dim": r.time, "x": np.nan, "y": r.y}
             for r in A.itertuples()]
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------
# 19. Infection information in a feature set, tested on capture areas held out of training
# ----------------------------------------------------------------------------------------------

def infection_signal_held_out(feature_sets, obs, sample_key="sample", cond_key="condition", time_key="timepoint",
                              infected="infected", depth=None, n_pcs=50, seed=0):
    """Within each time point, a logistic regression is trained on one infected and one mock capture area and scored
    (ROC AUC) on the other infected-mock pair, for every pairing, so that the score measures infection information that
    generalises across biological replicates. If `depth` (log10 UMIs per spot) is given, it is scored as a feature set
    on its own and regressed out of every other feature set inside each training fold ('... | depth removed').
    Returns per-pair results and a summary (mean, min, max AUC, pairs above 0.5)."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression, LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    sec, cond, tim = _as_str(obs[sample_key]), _as_str(obs[cond_key]), _as_str(obs[time_key])
    y = (cond == infected).astype(int)
    sets = {k: np.asarray(v, float).reshape(len(sec), -1) for k, v in feature_sets.items()}
    d = None if depth is None else np.asarray(depth, float).reshape(-1, 1)
    if d is not None:
        sets["capture depth only"] = d
    tsort = lambda ts: sorted(ts, key=lambda t: (int(re.sub(r"\D", "", t) or 0), t))
    rows = []
    for tp in tsort(pd.unique(tim)):
        I = sorted(pd.unique(sec[(tim == tp) & (y == 1)])); M = sorted(pd.unique(sec[(tim == tp) & (y == 0)]))
        if len(I) < 2 or len(M) < 2:
            continue
        for a in I:
            for b in M:
                tr = np.isin(sec, [a, b]); te = (tim == tp) & ~tr
                for name, A in sets.items():
                    for resid in ((False, True) if d is not None and name != "capture depth only" else (False,)):
                        Xtr, Xte = A[tr], A[te]
                        if resid:
                            lr = LinearRegression().fit(d[tr], Xtr)
                            Xtr, Xte = Xtr - lr.predict(d[tr]), Xte - lr.predict(d[te])
                        sc = StandardScaler().fit(Xtr); Xtr, Xte = np.nan_to_num(sc.transform(Xtr)), np.nan_to_num(sc.transform(Xte))
                        if Xtr.shape[1] > n_pcs:
                            p = PCA(n_pcs, random_state=seed).fit(Xtr); Xtr, Xte = p.transform(Xtr), p.transform(Xte)
                        clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(Xtr, y[tr])
                        rows.append({"time": tp, "features": name + (" | depth removed" if resid else ""),
                                     "trained_on": f"{a} + {b}", "auc": roc_auc_score(y[te], clf.predict_proba(Xte)[:, 1])})
    res = pd.DataFrame(rows)
    summ = (res.groupby(["time", "features"], sort=False)
               .agg(mean_auc=("auc", "mean"), min_auc=("auc", "min"), max_auc=("auc", "max"), n_pairs=("auc", "size"),
                    pairs_above_half=("auc", lambda s: int((s > 0.5).sum()))).reset_index())
    return res, summ


def plot_infection_signal(summ):
    feats = list(dict.fromkeys(summ["features"]))
    times = list(dict.fromkeys(summ["time"]))
    fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(feats) + 1.0))
    for k, tp in enumerate(times):
        g = summ[summ["time"] == tp].set_index("features").reindex(feats)
        yy = np.arange(len(feats)) + (k - (len(times) - 1) / 2) * 0.25
        ax.errorbar(g["mean_auc"], yy, xerr=[g["mean_auc"] - g["min_auc"], g["max_auc"] - g["mean_auc"]], fmt="o", ms=4,
                    capsize=2, lw=0.8, label=tp)
    ax.axvline(0.5, color="grey", ls="--", lw=0.7)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels(feats, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("ROC AUC for infection in held-out capture areas (mean, range over pairs)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# 20. Viral UMIs on tissue: one shared scale for all capture areas
# ----------------------------------------------------------------------------------------------

def _short_label(name):
    m = re.match(r"(control|mock|infected)_(\d+)dpi_(\d+)", str(name))
    return f"{'I' if m.group(1) == 'infected' else 'M'}{m.group(2)}-{m.group(3)}" if m else str(name)


def plot_viral_maps(obs, viral, sample_key="sample", order=None):
    """Spots of every capture area on the array grid: spots without viral UMIs in light grey, spots with viral UMIs
    coloured on one logarithmic scale shared by all capture areas. Panels in the order of Fig. 3a."""
    from matplotlib.colors import LogNorm
    sec = _as_str(obs[sample_key]); v = np.asarray(viral, float)
    if order is None:
        key = lambda s: (int((re.search(r"_(\d+)dpi", s) or [0, 0])[1]), "infected" in s, s)
        order = sorted(pd.unique(sec), key=key)
    pc = physical_coords(obs)
    vmax = max(1.0, float(np.nanmax(v))) if np.isfinite(v).any() else 1.0
    norm = LogNorm(vmin=1, vmax=vmax)
    ncol = 4; nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.6 * nrow + 0.6), squeeze=False)
    sm = None
    for ax, s in zip(axes.ravel(), order):
        m = sec == s; p = pc[m]; vv = v[m]; pos = vv > 0
        ax.scatter(p[~pos, 0], p[~pos, 1], s=3, color="#d9d9d9", lw=0)
        sm = ax.scatter(p[pos, 0], p[pos, 1], s=9, c=vv[pos], cmap="magma_r", norm=norm, lw=0)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_axis_off()
        ax.set_title(f"{_short_label(s)}\n{int(vv.sum())} viral UMIs on {int(pos.sum())} of {int(m.sum())} spots", fontsize=7)
    for ax in axes.ravel()[len(order):]:
        ax.set_axis_off()
    fig.subplots_adjust(left=0.01, right=0.9, top=0.92, bottom=0.03, wspace=0.06, hspace=0.3)
    if sm is not None:
        cb = fig.colorbar(sm, cax=fig.add_axes([0.925, 0.25, 0.012, 0.5]))
        cb.set_label("TYLCV UMIs per spot (grey: none)", fontsize=7); cb.ax.tick_params(labelsize=6)
    return fig


# ----------------------------------------------------------------------------------------------
# 21. Pre-specified analyses added before the final run (gene sets, thresholds and groupings fixed here)
# ----------------------------------------------------------------------------------------------

# Regular expressions on NCBI product descriptions (case-insensitive); nuclear genes only (organelle loci excluded).
PATHWAY_RULES = {
    "Photosynthesis, light reactions": r"chlorophyll a-b binding|photosystem|cytochrome b6-f|plastocyanin|oxygen-evolving enhancer|ferredoxin--NADP",
    "Calvin-Benson cycle": (r"ribulose bisphosphate carboxylase|ribulose-1,5-bisphosphate carboxylase|rubisco activase|phosphoribulokinase"
                            r"|sedoheptulose-1,7-bisphosphatase|fructose-1,6-bisphosphatase, chloroplastic"
                            r"|glyceraldehyde-3-phosphate dehydrogenase [AB], chloroplastic|phosphoglycerate kinase, chloroplastic"
                            r"|transketolase, chloroplastic"),
    "DNA replication and S phase": (r"proliferating cell nuclear antigen|DNA polymerase|DNA replication licensing factor"
                                    r"|minichromosome maintenance|ribonucleoside-diphosphate reductase|replication factor C"
                                    r"|replication protein A|\bhistone H[34]\b|\bhistone H2[AB]\b|\bcyclin-[ABD]|E2F"),
    "RNA silencing": (r"argonaute|endoribonuclease Dicer|dicer-like|RNA-dependent RNA polymerase|suppressor of gene silencing"
                      r"|small RNA 2'-O-methyltransferase|double-stranded RNA-binding protein"),
    "Salicylic acid pathway": r"pathogenesis-related protein 1|regulatory protein NPR|isochorismate synthase|protein EDS1|lipase-like PAD4",
    "Jasmonate pathway": (r"lipoxygenase|allene oxide synthase|allene oxide cyclase|12-oxophytodienoate reductase|jasmonate"
                          r"|protein TIFY|coronatine-insensitive|proteinase inhibitor"),
    "Phenylpropanoids and flavonoids": (r"phenylalanine ammonia-lyase|trans-cinnamate 4-monooxygenase|4-coumarate--CoA ligase"
                                        r"|chalcone synthase|chalcone--flavon[oa]ne isomerase|chalcone isomerase|flavanone 3-hydroxylase"
                                        r"|naringenin,2-oxoglutarate 3-dioxygenase|cinnamoyl-CoA reductase|cinnamyl alcohol dehydrogenase"
                                        r"|caffeoyl-CoA O-methyltransferase"),
    "Heat-shock proteins": r"heat shock|heat-shock|\bHSP\d",
}


def gene_products_from_gtf(gtf_path):
    """Product description of each gene from the GTF of the Space Ranger reference (NCBI attributes 'product')."""
    import gzip
    prod = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = re.search(r'product "([^"]+)"', f[8])
            if not m:
                continue
            for key in ("gene_id", "gene_name", "gene"):
                g = re.search(key + r' "([^"]+)"', f[8])
                if g and g.group(1) not in prod:
                    prod[g.group(1)] = m.group(1)
    return pd.Series(prod, name="product", dtype=object)


def pathway_gene_sets(products, genes, rules=None, exclude_regex=r"^(CQW56_|LyesC2|TYLCV)", min_size=5):
    """Gene sets from product descriptions with the pre-specified PATHWAY_RULES; returns (long table, set sizes)."""
    rules = PATHWAY_RULES if rules is None else rules
    genes = pd.Index([str(g) for g in genes])
    prod = products.reindex(genes).fillna("").astype(str)
    keep = ~genes.str.match(exclude_regex)
    rows = []
    for name, rx in rules.items():
        hit = prod.str.contains(rx, case=False, regex=True).to_numpy() & keep
        rows += [(name, g, prod[g]) for g in genes[hit]]
    df = pd.DataFrame(rows, columns=["gene_set", "gene", "product"])
    sizes = df["gene_set"].value_counts().reindex(list(rules)).fillna(0).astype(int).rename("n_genes_in_matrix")
    return df[df["gene_set"].isin(sizes[sizes >= min_size].index)].reset_index(drop=True), sizes


def thin_counts(counts, obs, sample_key="sample", time_key="timepoint", target=None, seed=0):
    """Binomial thinning. Default: within each time point (slide), every capture area is thinned to the median UMIs
    per spot of the shallowest capture area of that slide; with `target`, every capture area to that median.
    Returns thinned counts (float32 CSR) and the factors used."""
    rng = np.random.default_rng(seed)
    C = sp.csr_matrix(counts)
    sec, tim = _as_str(obs[sample_key]), _as_str(obs[time_key])
    tot = np.asarray(C.sum(1)).ravel()
    med = pd.Series(tot).groupby(sec).median()
    fac = np.ones(C.shape[0]); rows = []
    for tp in pd.unique(tim):
        secs = list(pd.unique(sec[tim == tp]))
        tgt = float(target) if target is not None else float(min(med[s] for s in secs))
        for s in secs:
            f = min(1.0, tgt / med[s]) if med[s] > 0 else 1.0
            fac[sec == s] = f
            rows.append({"timepoint": tp, "sample": s, "median_umis_before": float(med[s]), "target_median": tgt, "factor": f})
    Cc = C.tocoo()
    data = rng.binomial(np.rint(Cc.data).astype(np.int64), fac[Cc.row])
    T = sp.csr_matrix((data.astype(np.float32), (Cc.row, Cc.col)), shape=C.shape)
    T.eliminate_zeros()
    out = pd.DataFrame(rows)
    tot2 = np.asarray(T.sum(1)).ravel()
    out["median_umis_after"] = [float(np.median(tot2[sec == s])) for s in out["sample"]]
    return T, out


def compare_de(orig, new, genes=(), alpha=0.05):
    """Agreement between two two_level_de results (original vs re-analysis)."""
    common = orig.index.intersection(new.index)
    so, sn = set(orig.index[orig["padj_leaf"] < alpha]), set(new.index[new["padj_leaf"] < alpha])
    ro, rn = set(orig.index[orig["robust"]]), set(new.index[new["robust"]])
    from scipy.stats import spearmanr
    top = orig.index[:300].intersection(new.index)
    summary = {"padj_leaf_lt_0.05_original": len(so), "padj_leaf_lt_0.05_reanalysis": len(sn), "shared": len(so & sn),
               "robust_original": len(ro), "robust_reanalysis": len(rn), "robust_shared": len(ro & rn),
               "logFC_spearman_top300_original": float(spearmanr(orig.loc[top, "log2FoldChange_leaf"], new.loc[top, "log2FoldChange_leaf"])[0]) if len(top) > 2 else np.nan,
               "genes_compared": len(common)}
    g = [x for x in genes if x in orig.index or x in new.index]
    detail = pd.DataFrame({"log2FC_original": orig["log2FoldChange_leaf"].reindex(g), "padj_original": orig["padj_leaf"].reindex(g),
                           "log2FC_reanalysis": new["log2FoldChange_leaf"].reindex(g), "padj_reanalysis": new["padj_leaf"].reindex(g)})
    return summary, detail


def leaf_jackknife(counts, obs, var_names, genes, time_value, time_col="timepoint", sample_col="sample", leaf_col="leaf",
                   cond_col="condition", min_spots=10):
    """Leaf-level DESeq2 refitted with each leaf left out in turn; effect sizes of `genes` in every refit."""
    m = (_as_str(obs[time_col]) == str(time_value)) & ~pd.Series(_as_str(obs[leaf_col])).str.endswith("_fragment").to_numpy()
    pb, md = pseudobulk(counts[m], obs.loc[m], [sample_col, leaf_col], var_names, min_spots, meta_cols=[cond_col])
    full = run_deseq2(pb, md, cond_col)
    rows = [{"timepoint": time_value, "left_out": "none", "gene": g, "log2FC": full.loc[g, "log2FoldChange"],
             "padj": full.loc[g, "padj"]} for g in genes if g in full.index]
    for leaf in md.index:
        keep = md.index != leaf
        if md.loc[keep, cond_col].astype(str).value_counts().min() < 2:
            continue
        r = run_deseq2(pb.loc[keep], md.loc[keep], cond_col)
        rows += [{"timepoint": time_value, "left_out": leaf, "gene": g, "log2FC": r.loc[g, "log2FoldChange"],
                  "padj": r.loc[g, "padj"]} for g in genes if g in r.index]
    res = pd.DataFrame(rows)
    if res.empty:
        return res, res
    f = res[res["left_out"] == "none"].set_index("gene")
    j = res[res["left_out"] != "none"]
    summ = j.groupby("gene").agg(n_refits=("log2FC", "size"), log2FC_min=("log2FC", "min"), log2FC_max=("log2FC", "max"),
                                 padj_max=("padj", "max"), refits_padj_lt_0_05=("padj", lambda s: int((s < 0.05).sum())))
    summ.insert(0, "log2FC_all_leaves", f["log2FC"].reindex(summ.index))
    summ.insert(1, "padj_all_leaves", f["padj"].reindex(summ.index))
    summ["same_sign_in_all_refits"] = [bool(np.all(np.sign(j.loc[j.gene == g, "log2FC"]) == np.sign(summ.loc[g, "log2FC_all_leaves"])))
                                       for g in summ.index]
    summ.insert(0, "timepoint", time_value)
    return res, summ.reset_index()


def markers_by_type(counts, obs, var_names, type_key="predicted_cell_type", sample_key="sample", min_spots=10,
                    lfc=1.0, fdr=0.05, exclude=None, min_capture_areas=3):
    """De novo markers with capture areas as replicates: in every capture area with >= min_spots spots of the type and
    of the rest, one pseudobulk of the type and one of all other spots; DESeq2 ~ capture area + type (paired), type vs
    rest. Marker: padj < fdr and log2 fold change > lfc."""
    lab, sec = _as_str(obs[type_key]), _as_str(obs[sample_key])
    ok = ~pd.Series(lab).isin(["nan", "None", ""]).to_numpy()
    cols = np.ones(len(var_names), bool) if exclude is None else ~np.asarray(exclude, bool)
    C, v = sp.csr_matrix(counts)[:, cols], np.asarray(var_names)[cols]
    out, info = [], []
    for t in sorted(pd.unique(lab[ok])):
        o = pd.DataFrame({"capture_area": sec, "group": np.where(lab == t, "type", "rest")}, index=obs.index)
        pb, md = pseudobulk(C[ok], o.loc[ok], ["capture_area", "group"], v, min_spots)
        both = md.groupby("capture_area")["group"].nunique()
        areas = both[both == 2].index
        info.append({"cell_type": t, "capture_areas_used": len(areas)})
        if len(areas) < min_capture_areas:
            continue
        sel = md["capture_area"].isin(areas).to_numpy()
        r = run_deseq2(pb.loc[sel], md.loc[sel], cond_col="group", ref="rest", alt="type", covariates=("capture_area",))
        r = r[["baseMean", "log2FoldChange", "pvalue", "padj"]].copy()
        r.insert(0, "cell_type", t); r["n_capture_areas"] = len(areas)
        out.append(r)
    res = pd.concat(out) if out else pd.DataFrame()
    if len(res):
        res["marker"] = (res["padj"] < fdr) & (res["log2FoldChange"] > lfc)
    return res, pd.DataFrame(info)


def marker_overlap(markers, marker_dict):
    """For each type: de novo markers among the genes tested, literature markers among them, and the overlap
    (hypergeometric P that the overlap is at least as large by chance)."""
    from scipy.stats import hypergeom
    rows = []
    for t, g in markers.groupby("cell_type"):
        universe = set(map(str, g.index)); de_novo = set(map(str, g.index[g["marker"]]))
        lit = set(map(str, marker_dict.get(t, ()))) & universe
        k, N, K, n = len(de_novo & lit), len(universe), len(lit), len(de_novo)
        rows.append({"cell_type": t, "genes_tested": N, "de_novo_markers": n, "literature_markers_tested": K, "overlap": k,
                     "fold_enrichment": (k / n) / (K / N) if n and K else np.nan,
                     "p_hypergeometric": float(hypergeom.sf(k - 1, N, K, n)) if n and K else np.nan,
                     "overlapping_genes": ", ".join(sorted(de_novo & lit))})
    return pd.DataFrame(rows)


def viral_strand_by_library(csv_for_sample, obs, sample_key="sample", cond_key="condition", time_key="timepoint",
                            infected="infected"):
    """Virion-sense (TYLCV_V) and complementary-sense (TYLCV_C) UMIs per library, for all barcodes and split into
    barcodes of analysed spots (on tissue) and all others; plus a quasi-binomial test of time point (infected libraries,
    capture areas as units) and, per library, Fisher's test of on- vs off-tissue strand composition."""
    import statsmodels.api as sm
    from scipy.stats import fisher_exact
    sec = _as_str(obs[sample_key]); rows = []
    for s in pd.unique(sec):
        path = csv_for_sample.get(s)
        if path is None or not Path(path).exists():
            continue
        tab = pd.read_csv(path); tab["bc16"] = barcode16(tab["barcode"])
        tab = tab.dropna(subset=["bc16"]).groupby("bc16")[["TYLCV_V", "TYLCV_C"]].sum()
        on = tab.index.isin(pd.Index(barcode16(obs.index[sec == s])))
        o1 = obs.loc[sec == s].iloc[0]
        r = {"sample": s, "condition": str(o1[cond_key]), "timepoint": str(o1[time_key])}
        for part, msk in (("all", np.ones(len(tab), bool)), ("on_tissue", on), ("off_tissue", ~on)):
            v_, c_ = float(tab.loc[msk, "TYLCV_V"].sum()), float(tab.loc[msk, "TYLCV_C"].sum())
            r.update({f"V_{part}": v_, f"C_{part}": c_, f"pct_complementary_{part}": 100 * c_ / (v_ + c_) if v_ + c_ else np.nan})
        r["fisher_p_on_vs_off"] = (float(fisher_exact([[r["V_on_tissue"], r["C_on_tissue"]], [r["V_off_tissue"], r["C_off_tissue"]]])[1])
                                   if min(r["V_on_tissue"] + r["C_on_tissue"], r["V_off_tissue"] + r["C_off_tissue"]) > 0 else np.nan)
        rows.append(r)
    tab = pd.DataFrame(rows)
    tests = []
    inf = tab[(tab["condition"] == infected) & ((tab["V_all"] + tab["C_all"]) > 0)]
    if inf["timepoint"].nunique() == 2 and len(inf) >= 3:
        tsort = sorted(inf["timepoint"].unique(), key=lambda t: int(re.sub(r"\D", "", t) or 0))
        X = sm.add_constant((inf["timepoint"] == tsort[1]).astype(float).to_numpy())
        y = np.c_[inf["C_all"].to_numpy(), inf["V_all"].to_numpy()]
        from scipy.stats import t as t_dist
        fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
        phi = max(1.0, float(fit.pearson_chi2 / fit.df_resid)) if fit.df_resid > 0 else np.nan   # never below binomial
        tval = float(fit.params[1] / (fit.bse[1] * np.sqrt(phi)))
        tests.append({"test": f"complementary-sense share, {tsort[1]} vs {tsort[0]} (quasi-binomial, capture areas as units)",
                      "odds_ratio": float(np.exp(fit.params[1])), "dispersion": phi, "df": float(fit.df_resid),
                      "p": float(2 * t_dist.sf(abs(tval), fit.df_resid)), "n_capture_areas": len(inf),
                      f"pct_complementary_{tsort[0]}": 100 * inf.loc[inf.timepoint == tsort[0], "C_all"].sum() / inf.loc[inf.timepoint == tsort[0], ["V_all", "C_all"]].values.sum(),
                      f"pct_complementary_{tsort[1]}": 100 * inf.loc[inf.timepoint == tsort[1], "C_all"].sum() / inf.loc[inf.timepoint == tsort[1], ["V_all", "C_all"]].values.sum()})
    return tab, pd.DataFrame(tests)


def plot_viral_strand(tab, infected="infected"):
    from statsmodels.stats.proportion import proportion_confint
    t = tab[(tab["condition"] == infected) & ((tab["V_all"] + tab["C_all"]) > 0)].copy()
    t["order"] = [int(re.sub(r"\D", "", x) or 0) for x in t["timepoint"]]
    t = t.sort_values(["order", "sample"])
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    for i, r in enumerate(t.itertuples()):
        n = r.V_all + r.C_all
        lo, hi = proportion_confint(r.C_all, n, method="wilson")
        ax.bar(i, r.pct_complementary_all, color="#7570b3" if r.order == t["order"].min() else "#d95f02", alpha=0.8)
        ax.errorbar(i, r.pct_complementary_all, yerr=[[r.pct_complementary_all - 100 * lo], [100 * hi - r.pct_complementary_all]],
                    color="black", capsize=2, lw=0.8)
        for part, mk in (("on_tissue", "o"), ("off_tissue", "s")):
            val = getattr(r, f"pct_complementary_{part}")
            if np.isfinite(val):
                ax.scatter(i + (0.18 if part == "off_tissue" else -0.18), val, marker=mk, s=14, color="white", edgecolor="black", zorder=3)
    ax.set_xticks(range(len(t))); ax.set_xticklabels([_short_label(s) for s in t["sample"]])
    ax.set_ylabel("Complementary-sense share of\nviral UMIs (%)"); ax.set_ylim(0, 100)
    ax.set_title("Bars: all barcodes (95% CI); circles: on tissue; squares: off tissue", fontsize=7)
    fig.tight_layout()
    return fig
