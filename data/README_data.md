# Data access

## Processed AnnData (GitHub releases)
`tylcv_spatial_4259cells_planttag_annotated.h5ad`

Key slots:
- `adata.obs["planttag_celltype"]`  — 11 cell types
- `adata.obs["condition"]`          — TYLCV or mock
- `adata.obs["timepoint"]`          — 7dpi or 24dpi
- `adata.obsm["X_planttag"]`        — 64-dim latent embeddings

## Raw FASTQ
GEO accession: **pending**
Samples: A1 A2 (TYLCV) + A3 A4 (mock) at 7 and 24 dpi
