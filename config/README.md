# Settings used for the article

Notebook 02, cell 1:

| Setting | Value |
|---|---|
| `RAW_H5AD` | `processed/adata_combined_squidpy.h5ad` |
| `ANN_H5AD` | `processed/adata_combined_squidpy_clean_tda_features_plantag_gnn_embedding.h5ad` |
| `KEEP_SPOTS_CSV` | `M7-1_kept_spots.csv` (this folder) |
| `GTF` | `genes/genes.gtf.gz` of the Space Ranger reference |
| `LEAF_SPLITS`, `LEAF_GROUPS` | empty (automatic leaf grouping) |
| `SEED`, `FAST` | 0, False |

Notebook 01, cell 1: `SEED = 0`, `N_EPOCHS = 100`, `DEVICE = "cpu"`, `THIN_TO_MEDIAN_UMIS = None` (100 for the
sensitivity run).

**To add:** `M7-1_kept_spots.csv`, the 253 spots of capture area M7-1 (control_7dpi_1) kept after excluding the
displaced material (columns `section,spot`), as saved by the editor in cell 3a.
