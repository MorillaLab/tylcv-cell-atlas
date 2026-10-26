# Figures

- `make_figures.py`: Figures 2-7 and Supplementary Figures 1-3 from the result bundle (`BUNDLE_DIR`, default `../results/bundle`).
- `figure1/`: Figure 1 = the authors' original figure (`Fig1_authors_original.pdf`) with panels b, c, f and g redrawn
  from the bundle and the specifications of panel e updated (`make_bc_panels.py`, `make_fg_panels.py`, `build_figure1_final.py`).
- `make_supplementary_figure4.py`: Supplementary Figure 4 (run in Colab, it needs the brightfield images).
- `data/spot_domains.csv`: expression domain of each analysed spot (Figure 4a). It was recovered from the domain map of the
  bundle (`08_domains.png`) at the pixel position of every spot and reproduces `08_domain_counts_per_section.csv` exactly
  for all capture areas and domains.
