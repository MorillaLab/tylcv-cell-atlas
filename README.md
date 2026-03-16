<div align="center">

# 🍅 TYLCV Cell Atlas

### The first spatially resolved, cell-type-annotated transcriptomic atlas<br>of a plant–geminivirus interaction

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-green.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.9](https://img.shields.io/badge/Python-3.9-3776AB?logo=python&logoColor=white)](https://python.org)
[![PlantAG](https://img.shields.io/badge/Annotated%20with-PlantAG-1D9E75)](https://github.com/MorillaLab/PlantAG)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](https://doi.org)
[![GEO](https://img.shields.io/badge/GEO-pending-blue)](https://ncbi.nlm.nih.gov/geo)

**Luna\*, Aberbache\*, Albendea, Montano, Castillo, Bejarano, Doblas & Morilla**
<br>
*IHSM-UMA-CSIC, University of Malaga · Université Sorbonne Paris Nord*
<br>
*\* Equal contributors · Correspondence: ian.morilla@ihsm.uma-csic.es*

<br>

> *"TYLCV silences immunity in the cells it travels through,<br>
> and devastates the cells it colonises — a spatial logic<br>
> only visible at single-cell resolution."*

</div>

---

## The atlas

This repository is the **data and analysis home** for the first cell-type-resolved spatial transcriptomics atlas of *Tomato yellow leaf curl virus* (TYLCV) infection in *Solanum lycopersicum*. It contains the complete dataset, the experimental protocol that produced it, and every analysis notebook needed to reproduce the paper's figures.

The computational engine that annotates cell types — **PlantAG** — lives in its own repository: [github.com/MorillaLab/PlantAG](https://github.com/MorillaLab/PlantAG).

---

## Why this atlas matters

TYLCV causes near-total yield losses in susceptible tomato worldwide. Despite decades of molecular virology, we never knew **which cells** the virus targets, which defend themselves, which it exploits as a silent highway, and which it leaves alone. Bulk transcriptomics averages that information away. This atlas puts it back.

| Discovery | What it means |
|-----------|--------------|
| **Parenchyma is the primary target** (d = 0.62) | SA + JA immune co-activation + photosynthetic shutdown in the same cells — the molecular basis of chlorosis |
| **Vein is the secondary target** (d = 0.58) | Vascular stress + defence gene induction; cell wall remodelling at the infection front |
| **Phloem is paradoxically stable** (d = 0.15) | Active viral conduit, yet transcriptionally silent — direct spatial evidence for effector-mediated immune suppression |
| **A directed pseudotime gradient** | Vein → Bundle sheath → Parenchyma mirrors the physical path of viral cell-to-cell spread |
| **SA + JA co-activation** | Challenges canonical SA–JA antagonism in a geminiviral infection context |
| **Chlorosis is cell-autonomous** | Individual Parenchyma cells simultaneously activate immunity and suppress SlRBCS, SlLHCB1, SlPSBA |

---

## Dataset

| Sample | Condition | Time point | Cells |
|--------|-----------|------------|------:|
| A1 | TYLCV-infected | 7 dpi | ~435 |
| A2 | TYLCV-infected | 24 dpi | ~435 |
| A3 | Mock control | 7 dpi | ~435 |
| A4 | Mock control | 24 dpi | ~435 |
| **Total** | | | **1,740** |

**Platform:** 10x Genomics Visium Spatial Gene Expression v1
**Sequencing:** Illumina NovaSeq 6000 · 150 bp paired-end
**Reference:** *S. lycopersicum* SL4.0 / ITAG4.0

**Raw data + processed AnnData:** GEO accession pending → `data/README_data.md`

---

## Ten annotated cell types

<div align="center">

| Cell type | Colour | Viral impact (24 dpi) | Key markers |
|-----------|--------|----------------------|-------------|
| Parenchyma | 🔴 `#D62728` | **HIGH** d = 0.62 | SlWRKY33, SlPR1b, SlNPR1 |
| Vein | 🟠 `#E07B39` | **HIGH** d = 0.58 | SlSUC2, SlSWEET11, SlVND7 |
| Mesophyll | 🔵 `#1F77B4` | Moderate d = 0.41 | SlRBCS, SlLHCB1, SlCAB3 |
| Xylem | 🟡 `#BCBD22` | Moderate d = 0.38 | SlXCP2, SlIRX3, SlLAC4 |
| Epidermal | 🟣 `#9467BD` | Moderate d = 0.35 | SlCER6, SlKCS20, SlFDH |
| Bundle sheath | 🩵 `#17BECF` | Moderate d = 0.33 | SlPDLP1, SlCALS3, SlSEOR1 |
| Subsidiary | 🟤 `#8C564B` | Low d = 0.18 | SlML1, SlBRX, SlWUS2 |
| Phloem | 🟢 `#2CA02C` | **Stable** d = 0.15 | SlSUC1, SlAPL, SlPTF1 |
| Stomatal | 🩷 `#E377C2` | Stable d = 0.12 | SlEPF1, SlFAMA, SlSPCH |
| Hydathodes | ⚫ `#7F7F7F` | Stable d = 0.09 | SlHAP2, SlNAC55, SlDREB1 |

</div>

---

## Repository structure

```
tylcv-cell-atlas/
│
├── data/
│   ├── README_data.md          # GEO accession, download instructions
│   └── sample_metadata.csv     # Sample table (condition, dpi, n_cells)
│
├── protocol/
│   ├── tissue_preparation.md   # Full wet-lab protocol (step-by-step)
│   ├── permeabilisation.md     # PVP40/DTT optimisation — the key innovation
│   └── library_prep.md         # Visium library prep parameters
│
├── analysis/
│   ├── 01_preprocessing.py     # QC, normalisation, HVG, clustering
│   ├── 02_annotation.py        # PlantAG annotation (calls PlantAG package)
│   ├── 03_validation.py        # Marker gene validation, spatial coherence
│   ├── 04_cosine_disruption.py # Cell-type cosine distances, impact tiers
│   ├── 05_trajectory.py        # Diffusion maps, DPT, RNA velocity
│   ├── 06_deg_enrichment.py    # Wilcoxon DE, GO enrichment (g:Profiler)
│   ├── cosine_validation/      # Supplementary Note 1 — robustness analyses
│   └── trajectory/             # scVelo dynamical model
│
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_planttag_annotation.ipynb
│   ├── 03_validation.ipynb
│   ├── 04_trajectory_analysis.ipynb
│   ├── 05_viral_tropism.ipynb
│   └── 06_supplementary_figures.ipynb
│
├── environment.yml
├── requirements.txt
└── README.md
```

---

## Quickstart

```bash
git clone https://github.com/MorillaLab/tylcv-cell-atlas.git
cd tylcv-cell-atlas
conda env create -f environment.yml
conda activate tylcv-atlas

# Download processed data (GEO, pending)
# Place tylcv_spatial_1740cells_planttag_annotated.h5ad in data/

# Run the full analysis
jupyter nbconvert --to notebook --execute notebooks/01_preprocessing.ipynb
jupyter nbconvert --to notebook --execute notebooks/02_planttag_annotation.ipynb
# ... continue through 06
```

---

## Reproducing key figures

```python
import scanpy as sc
import pandas as pd
from analysis.cosine_disruption import compute_cosine_disruption, REPORTED_DISTANCES

# Load atlas
adata = sc.read_h5ad("data/tylcv_spatial_1740cells_planttag_annotated.h5ad")

# Replicate Figure 4a — cell-type viral disruption
df = compute_cosine_disruption(
    adata,
    celltype_key="planttag_celltype",
    condition_key="condition",
    infected_label="TYLCV",
    control_label="mock",
)
print(df.to_string(index=False))

# Compare with reported values (Table 2)
for _, row in df.iterrows():
    reported = REPORTED_DISTANCES["24dpi"].get(row["cell_type"], "n/a")
    print(f"{row['cell_type']:<16} computed={row['cosine_distance']:.3f}  reported={reported}")
```

---

## The wet-lab breakthrough

Tomato leaves are hostile to standard Visium protocols. We solved three barriers:

```
Problem 1: Polyphenols degrade RNA and inhibit reverse transcription
Solution:  2% PVP40 sequesters phenolics before enzyme treatment

Problem 2: Cellulosic cell wall blocks mRNA access
Solution:  Extended 40 min Visium enzyme permeabilisation
           after PVP40/DTT pre-treatment

Problem 3: Waxy cuticle disrupts cryosectioning
Solution:  3:1 ethanol:glacial acetic acid fixation preserves integrity
```

Full protocol: `protocol/tissue_preparation.md`
This is the **first validated Visium workflow for tomato leaf** and is
directly transferable to other polyphenol-rich crop species.

---

## Citation

```bibtex
@article{luna_aberbache_2026_tylcv_atlas,
  title   = {Cell-type-resolved spatial transcriptomics reveals stratified
             viral tropism and emergent defence coordination in
             {TYLCV}-infected tomato leaves},
  author  = {Luna, Ana M\*. and Aberbache\*, Melissa and Albendea, Mar{\'{i}}a
             and Montano, Jos{\'{e}} A. and Castillo, Araceli and Bejarano, Eduardo R. and Doblas, Ver{\'{o}}nica G.
             and Morilla, Ian},
  journal = {in preparation},
  year    = {2026},
  url     = {https://github.com/MorillaLab/tylcv-cell-atlas}
}
```

If you use PlantAG for annotation, please also cite:

```bibtex
@software{morilla_2026_planttag,
  title   = {{PlantAG}: Plant Tissue Annotation via Graph Topology},
  author  = {Aberbache, Melissa and Morilla, Ian},
  year    = {2026},
  doi     = {https://doi.org/10.5281/zenodo.19048604},
  url     = {https://github.com/MorillaLab/PlantAG}
}
```

---

## Licence

Data and protocol: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
Analysis code: [GPL-3.0](LICENSE)

---

<div align="center">
<b>MorillaLab · IHSM-UMA-CSIC · University of Malaga</b><br>
<a href="mailto:ian.morilla@ihsm.uma-csic.es">ian.morilla@ihsm.uma-csic.es</a>
</div>
