"""Reproduce Figure 4a-b — cell-type cosine distances."""
import scanpy as sc, pandas as pd, numpy as np
from scipy.spatial.distance import cosine
import scipy.sparse as sp

REPORTED_24DPI = {
    "Parenchyma": 0.62, "Vein": 0.58, "Mesophyll": 0.41,
    "Xylem": 0.38, "Epidermal": 0.35, "Bundle sheath": 0.33,
    "Subsidiary": 0.18, "Phloem": 0.15, "Stomatal": 0.12, "Hydathodes": 0.09,
}

def compute(adata, tp="24dpi"):
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    rows = []
    for ct in adata.obs["planttag_celltype"].unique():
        m = adata.obs["planttag_celltype"] == ct
        t = adata.obs["timepoint"] == tp
        inf  = X[(m) & (t) & (adata.obs["condition"]=="TYLCV")]
        ctrl = X[(m) & (t) & (adata.obs["condition"]=="mock")]
        if len(inf) and len(ctrl):
            rows.append({"cell_type": ct, "d": cosine(inf.mean(0), ctrl.mean(0))})
    return pd.DataFrame(rows).sort_values("d", ascending=False)

if __name__ == "__main__":
    adata = sc.read_h5ad("../data/tylcv_spatial_1740cells_planttag_annotated.h5ad")
    df = compute(adata)
    print(df.to_string(index=False))
