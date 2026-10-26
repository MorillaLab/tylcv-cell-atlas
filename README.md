# 🍅 Spatial transcriptomics of TYLCV-infected tomato leaves

Code to reproduce every result and figure of:

> Luna, A. M., Aberbache, M. et al. *Spatial transcriptomics with PlantAG reveals a tissue-wide response of tomato
> leaves to Tomato yellow leaf curl virus* (submitted).

Cell types are imputed with **PlantAG** (Plant Tissue Annotation via Graph Topology), developed in our group and
maintained separately: [MorillaLab/PlantAG](https://github.com/MorillaLab/PlantAG).

## Data

- **Visium data** (8 capture areas; tomato cv. Moneymaker, TYLCV-infected and mock, 7 and 24 dpi): GEO [accession to be added].
- **Result bundles, trained PlantAG model and processed objects:** Zenodo, https://doi.org/10.5281/zenodo.23664403.
- **Reference:** tomato genome SL3.1 (GCF_000188115.5) with NCBI *Solanum lycopersicum* Annotation Release 103, plus TYLCV.

## Workflow

| Step | What | Where |
|---|---|---|
| 1 | Space Ranger with the tomato reference (and, optionally, TYLCV) | `reference/` |
| 2 | Viral UMIs per barcode from unmapped reads, strand-resolved | `viral/` |
| 3 | Exclusion of spots outside leaf tissue (M7-1) and PlantAG training | `config/`, `notebooks/01_PlantAG_training.ipynb` |
| 4 | Analyses: QC, leaves, PlantAG evaluation, topology, domains, differential expression, virus, splicing kinetics | `notebooks/02_PlantAG_reliability_audit.ipynb`, `src/plantag_audit.py` |
| 5 | Figures 1-7 and Supplementary Figures 1-4 from the result bundle | `figures/` |

Step by step, with settings and run times: [docs/reproduce.md](docs/reproduce.md). What each file of the result bundle
contains and which manuscript item it supports: [docs/bundle_files.md](docs/bundle_files.md).

## Layout

```
reference/    TYLCV reference, merged tomato+TYLCV reference, Space Ranger runs
viral/        unmapped reads -> TYLCV alignment -> UMIs per barcode (TYLCV_V, TYLCV_C)
config/       settings used for the article and the list of retained M7-1 spots
notebooks/    01 PlantAG training (Colab) | 02 analysis notebook (Colab) | optional limma notebook
src/          plantag_audit.py, the analysis module used by notebook 02
figures/      scripts that draw the article figures from the result bundle
results/      place the unzipped result bundle here (results/bundle)
docs/         reproduction guide and bundle file dictionary
```

## Citation

See [CITATION.cff](CITATION.cff). Please also cite PlantAG.

## Licence

Code: MIT ([LICENSE](LICENSE)). Data: see GEO and Zenodo.
