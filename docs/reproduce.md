# Reproducing the article

All notebooks run on Google Colab (CPU runtime is sufficient). Paths are set in the first cells.

## 1. Count matrices (Space Ranger)

Space Ranger v2.0.0 (7-dpi slide) and v2.0.1 (24-dpi slide) against SL3.1 / NCBI Annotation Release 103.
`reference/make_tylcv_reference.py` writes the TYLCV FASTA and GTF; `reference/merge_references.sh` and
`reference/spaceranger_tylcv.sh` build and run a combined tomato + TYLCV reference if you want viral features in the
matrix (route A).

## 2. Viral UMIs (route B, used in the article)

`viral/tylcv_viral_reads.sh` extracts the reads Space Ranger left unmapped, aligns them to TYLCV with minimap2 and calls
`viral/count_viral_umis.py`, which writes one CSV per library with virion-sense (`TYLCV_V`) and complementary-sense
(`TYLCV_C`) UMIs per barcode. `viral/bams.tsv` maps BAM files to capture areas (the two infected 7-dpi libraries are
crossed in the facility file names; `make_bams_tsv.sh` applies the correction). `viral/check_laptop.sh` checks a local
machine (Apple Silicon supported through conda).

## 3. Spot exclusion and PlantAG training

1. In notebook 02, run cells 1-2 and the interactive editor in cell 3a, or use `config/M7-1_kept_spots.csv`.
   Set `KEEP_SPOTS_CSV` in cell 1, then run cell 3c to write `adata_combined_squidpy_clean.h5ad`.
2. Run `notebooks/01_PlantAG_training.ipynb` on the clean object (about 15-30 min on CPU; seeds fixed).
   It writes `..._clean_tda_features_plantag_gnn_embedding.h5ad`.
3. Depth-matched sensitivity run (Supplementary Data 15): set `THIN_TO_MEDIAN_UMIS = 100` in cell 1 of notebook 01 and
   run it again; outputs get a `_depth100` suffix.

## 4. Analyses

Notebook 02 with the settings in `config/README.md`. Run all cells, including the optional R cells at the end (limma,
camera). Everything is written to a result bundle (`bundle.zip`), about 95 files.

## 5. Figures

```bash
unzip bundle.zip -d results/bundle
python figures/make_figures.py                       # Figures 2-7, Supplementary Figures 1-3 -> figures/output
cd figures/figure1 && python make_bc_panels.py && python make_fg_panels.py && python build_figure1_final.py
```

Supplementary Figure 4 needs the brightfield images and is drawn in Colab with `figures/make_supplementary_figure4.py`.
`BUNDLE_DIR`, `FIG_OUT` and `FIG1_SOURCE_PDF` override the default paths.
