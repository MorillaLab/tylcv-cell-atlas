# Result bundle: file dictionary

Files of bundle 05 (the article), grouped by prefix.

| Prefix | Analysis | Manuscript | Files |
|---|---|---|---|
| 00 | Spot exclusion in capture area M7-1 | Methods (quality control); Supplementary Data 3 | `00_spot_exclusion.csv`<br>`00_spot_exclusion.png`<br>`00_spot_exclusion_proposal.png` |
| 01 | Quality control per capture area; most abundant features | Fig. 2b-d; Results (workflow) | `01_qc_boxplots.png`<br>`01_sample_qc.csv`<br>`01_top30_genes_by_umis.csv` |
| 02 | Velocyto loom provenance and barcode matching | Methods (RNA velocity) | `02_loom_mapping.csv`<br>`02_loom_vs_section_count_correlation.csv` |
| 03 | Tissue pieces, leaf grouping, spot filtering | Fig. 2a; Supplementary Fig. 3; Supplementary Data 3 | `03_leaf_pieces.csv`<br>`03_leaf_segmentation.csv`<br>`03_leaves.png`<br>`03_qc_filter.csv`<br>`03_suggested_cuts.png` |
| 04 | Marker lists: dictionary audit, coverage, random-list null, split-half AUC | Fig. 3e-g; Table 1; Supplementary Data 4 | `04_heldout_marker_auc.csv`<br>`04_label_share_vs_null.csv`<br>`04_label_share_vs_null.png`<br>`04_marker_dictionary_audit.csv` |
| 04b | Local H1 persistence and the infected-versus-mock test | Fig. 1b,c; Supplementary Data 14 | `04b_h1_panels_data.csv`<br>`04b_h1_persistence.png`<br>`04b_h1_persistence_by_capture_area.csv`<br>`04b_h1_persistence_per_spot.csv`<br>`04b_h1_persistence_test.csv` |
| 04c | Infection information in held-out capture areas | Fig. 3h; Supplementary Data 16 | `04c_infection_signal_held_out.png`<br>`04c_infection_signal_held_out_pairs.csv`<br>`04c_infection_signal_held_out_summary.csv` |
| 05 | Imputed types per capture area and association with infection | Fig. 3a; Table 1; Supplementary Data 5 | `05_predicted_cell_type_share_per_section.csv`<br>`05_predicted_cell_type_vs_infection.csv`<br>`05_pseudo_cell_type_share_per_section.csv`<br>`05_pseudo_cell_type_vs_infection.csv` |
| 06 | PlantAG input blocks and graph overlap | Fig. 3d | `06_feature_block_share.csv`<br>`06_knn_overlap_with_full_graph.csv` |
| 07 | PlantAG evaluation: held-out capture areas, baselines, seeds, spatial coherence | Fig. 3b,c; Table 1 | `07_loso_per_section.csv`<br>`07_random_vs_loso_summary.csv`<br>`07_seed_stability.json`<br>`07_spatial_coherence_of_labels.csv` |
| 07c | De novo markers and overlap with the literature lists | Supplementary Data 17 | `07c_denovo_markers.csv`<br>`07c_denovo_markers_capture_areas.csv`<br>`07c_denovo_markers_vs_literature.csv` |
| 08 | Expression domains: maps, composition, coherence, infection | Fig. 4a-c | `08_domain_counts_per_section.csv`<br>`08_domain_marker_profile.csv`<br>`08_domain_spatial_coherence.csv`<br>`08_domain_vs_infection.csv`<br>`08_domains.png` |
| 08c | UMAP of the PlantAG embedding | Fig. 1f,g | `08c_plantag_umap.csv` |
| 09 | DESeq2 at leaf and capture-area level | Fig. 5; Table 2 | `09_DE_24dpi_all_genes_top300.csv`<br>`09_DE_24dpi_no_organelle_top300.csv`<br>`09_DE_7dpi_all_genes_top300.csv`<br>`09_DE_7dpi_no_organelle_top300.csv` |
| 09b | limma-voom with duplicateCorrelation | Fig. 5; Table 2 | `09b_limma_dupcor_24dpi.csv`<br>`09b_limma_dupcor_7dpi.csv`<br>`09b_limma_dupcor_summary.csv`<br>`09b_limma_vs_deseq2.csv` |
| 09c | Depth-matched differential expression | Supplementary Data 18 | `09c_DE_depth_matched_24dpi_top300.csv`<br>`09c_DE_depth_matched_7dpi_top300.csv`<br>`09c_depth_matched_summary.csv`<br>`09c_depth_matched_table2_genes.csv`<br>`09c_depth_matching_factors.csv`<br>`09c_limma_depth_matched_24dpi.csv`<br>`09c_limma_depth_matched_7dpi.csv`<br>`09c_limma_depth_matched_summary.csv` |
| 09d | Table 2 genes within imputed types; heterogeneity | Fig. 5f,g; Supplementary Data 19 | `09d_table2_genes_by_type.csv`<br>`09d_table2_genes_type_heterogeneity.csv` |
| 09e | Pathway gene sets and camera tests | Supplementary Data 19 | `09e_camera_gene_sets.csv`<br>`09e_gene_set_sizes.csv`<br>`09e_gene_sets.csv` |
| 09f | Leave-one-leaf-out refits | Supplementary Data 18 | `09f_leaf_jackknife_refits.csv`<br>`09f_leaf_jackknife_summary.csv` |
| 10 | Infected-versus-mock distances against replicate variation | Fig. 4d,e | `10_domain_disruption_with_nulls.csv`<br>`10_original_metric_vs_within_condition.csv` |
| 11 | Viral UMIs: per library, on tissue, maps, viral-positive spots | Fig. 6a-c; Supplementary Fig. 1 | `11_viral_barcode_diagnostic.csv`<br>`11_viral_csv_attach.csv`<br>`11_viral_maps.png`<br>`11_viral_per_section.csv`<br>`11_viral_per_spot.csv`<br>`11_viral_pos_vs_neg_paired_DE_top300.csv`<br>`11_viral_positive_share_by_domain.csv` |
| 11b | Viral tropism by imputed type | Fig. 7a | `11b_viral_tropism_by_type.csv`<br>`11b_viral_tropism_by_type.png` |
| 11c | Vascular-group tropism test | Supplementary Data 20 | `11c_viral_tropism_vascular_group.csv` |
| 11d | Viral strand programme | Fig. 6d,e; Supplementary Data 20 | `11d_viral_strand.png`<br>`11d_viral_strand_by_capture_area.csv`<br>`11d_viral_strand_test.csv` |
| 12 | RNA velocity and pooled splicing-kinetics tests | Supplementary Fig. 2 | `12_unspliced_fraction_test_by_leaf.csv`<br>`12_unspliced_fraction_test_by_sample.csv`<br>`12_velocity_barcode_matching.csv`<br>`12_velocity_confidence_maps.png`<br>`12_velocity_confidence_per_section.csv`<br>`12_velocity_confidence_vs_null.json`<br>`12_velocity_gene_table.csv`<br>`12_velocity_genes_by_min_shared_counts.csv`<br>`12_velocity_input_qc.csv`<br>`12_velocity_loom_assignment.csv`<br>`12_velocity_stream_umap.png` |
| 12b | Splicing kinetics within imputed types | Fig. 7b,c | `12b_unspliced_fraction_by_type_leaf.csv`<br>`12b_unspliced_fraction_by_type_sample.csv` |

`00_spot_exclusion_proposal.png` is a leftover of an earlier run in the same Drive folder; notebook 02 does not write it.
