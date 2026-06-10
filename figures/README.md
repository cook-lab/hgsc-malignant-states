# Figure Panel Map

Maps every manuscript panel to the repository script that generates it. This table is the source of truth for the panel→script mapping; each script's header docstring lists its exact inputs, outputs, and the panel(s) it supports. Paths are relative to the repo root (`hgsc-malignant-states/`).

**Legend:** SecA = Secretory A (progenitor-like, NMF Factor 3); SecB = Secretory B (differentiated/adaptive, NMF Factor 2); Int = Intermediate. Per-figure scripts read cached intermediates produced by the `atlas/` and `spatial/` backends; the `figures/_prep/` helpers build shared caches (`meta.parquet`, atlas obs/UMAP parquets). External/known-gap panels are flagged in **Notes**.

## Output location

Every figure script writes its panel file(s) to a **per-figure subdirectory** under `OUTPUT_ROOT/figures/`:
`figure1/` … `figure7/`, `supplementary/` (SF1–SF14), and the placeholder dirs `figure_icon7_bevacizumab/` /
`figure_ccle_smac_mimetics/`. Nothing is written to the `figures/` root. Each panel is exported as both `.svg`
(vector, for assembly) and `.png` (raster preview); spatial/ROI scripts may also emit `.pdf`.

## Main figures

| Panel | Script | Notes |
|---|---|---|
| Fig 1A | `figures/figure1/01_atlas_level1_umap.py` | |
| Fig 1B,1C,1D,1E,1I | `figures/figure1/02_atlas_study_contribution.py` | study contribution + epithelial composition |
| Fig 1F | `figures/figure1/03_atlas_epithelial_secretory_ciliated_umap.py` | |
| Fig 1G | `figures/figure1/04_atlas_SecA_nmf_factor_umap.py` | Factor_3 / SecA |
| Fig 1H | `figures/figure1/05_atlas_SecB_nmf_factor_umap.py` | Factor_2 / SecB |
| Fig 1J | `figures/figure1/06_atlas_cnv_alluvial.py` | consumes `atlas/05_cnv` within-clone coexistence cache |
| Fig 2A | `figures/figure2/01_atlas_radar_characterization.py` | PROGENy/Hallmark/DoRothEA/Flux radar |
| Fig 2B | `figures/figure2/02_atlas_g2m_barplot.py` | cell-cycle cache from NMF-characterization prep helper |
| Fig 2C | `figures/figure2/03_atlas_volcano_secA_secB.py` | SecA/SecB markers from `shared/signatures.yml` |
| Fig 2D,2E | `figures/figure2/04_atlas_site_secb_shift.py` | |
| Fig 2F,2G | `figures/figure2/05_atlas_treatment_secb_shift.py` | |
| Fig 3A | `figures/figure3/01_organoids_models_umap.R` | **EXTERNAL** organoid data (`organoids_root`); TRUST-EXISTING per author decision |
| Fig 3B | `figures/figure3/02_organoids_secB_ucell_by_model.py` | **EXTERNAL** organoid dependency |
| Fig 3C | *(no in-repo generator)* | **EXTERNAL / KNOWN GAP**: flow cytometry, no in-repo script |
| Fig 3D | `figures/figure3/03_organoids_secB_timecourse.py` | **EXTERNAL** organoid dependency |
| Fig 3E | `figures/figure3/04_organoids_secB_perturbations.py` | **EXTERNAL** organoid dependency |
| Fig 3F | `figures/figure3/05_organoids_g2m_barplot.py` | **EXTERNAL** organoid dependency |
| Fig 3G | `figures/figure3/06_organoids_radar_characterization.py` | **EXTERNAL** organoid dependency |
| Fig 3H | `figures/figure3/07_atlas_tcga_subtype_epitype.py` | atlas/ConsensusOV, in-repo (n=96) |
| Fig 4A,4B | `figures/figure4/01_xenium_composition_by_tissue_dendro.py` | |
| Fig 4C | `figures/figure4/02_xenium_whole_tissue_snapshot.R` | also drives SF10A via `all` arg |
| Fig 4D | `figures/figure4/03_fig2d_roi_celltype.R` | **FLAG**: colours by in-object SingleR/pre-06f labels, not 06f reclassification |
| Fig 4E | `figures/figure4/04_lee_bivariate_segregation.py` | |
| Fig 4F | `figures/figure4/05_bilisa_regime_by_label.py` | |
| Fig 4G | `figures/figure4/06_gam_microenvironment_polarization.R` | |
| Fig 4H | `figures/figure4/07_xenium_vascular_distance_paired.py` | |
| Fig 4I | `figures/figure4/08_xenium_roi_OTB_2384_vascular_distance.R` | |
| Fig 4J | `figures/figure4/09_xenium_epi_density_paired.py` | |
| Fig 4K | `figures/figure4/10_xenium_roi_OTB_2384_pathway_panels.R` | |
| Fig 5A | `figures/figure5/01_gam_epithelial_pathways_polarization.R` | |
| Fig 5B,5C | `figures/figure5/02_roi_OTB_2384_roi06_svgs.R` | |
| Fig 5D,5E | `figures/figure5/03_morphometrics_paired.py` | |
| Fig 5F | `figures/figure5/04_atlas_seca_secb_autocrine_shift.py` | **KNOWN ISSUE**: generator reproduces the panel exactly, but thresholds on expression magnitude (`lrscore>0.5`), not significance (see REPRODUCIBILITY.md → Fig 5F) |
| Fig 5G | `figures/figure5/05_gam_epithelial_genes_polarization.R` | |
| Fig 5H | `figures/figure5/06_roi_SP24_24824_zoom_genes.R` | cell-type map over the full SP24_24824 ROI (`SP24_24824_roi_celltype_full.*`); same script renders 5I |
| Fig 5I | `figures/figure5/06_roi_SP24_24824_zoom_genes.R` | 4-gene full-ROI maps (CTNNB1/ITGB5/MMP7/ICAM1); shares script with 5H |
| Fig 6A,6C | `figures/figure6/01_tma_hypoxia_gradient_cores.R` | writes `gradient_metrics.csv` consumed by scripts 03/04 |
| Fig 6B | `figures/figure6/02_niche_immune_composition_gam.R` | |
| Fig 6D | `figures/figure6/03_tma_gradient_cores_immune_highlight_svgs.R` (immune) + `04_tma_gradient_cores_secb_highlight.R` (SecB) | |
| Fig 6E,6F | `figures/figure6/05_fig_exhaustion_ucell.R` | |
| Fig 6G | `figures/figure6/06_gam_neighborhood_immunomod_genes_polarization.R` | |
| Fig 6H | `figures/figure6/07_gam_macrophage_focal_polarization.R` | |
| Fig 6I,6J | `figures/figure6/08_fig_macrophage_apoptosis_prolif_morphology.R` | **KNOWN ISSUE (6J)**: reproduces 6/8 p=0.042 vs published 8/8 p=0.008 (see REPRODUCIBILITY.md P1-b); migrated faithfully, not changed |
| Fig 6K,6L | `figures/figure6/09_xenium_roi_whole_tissue_investigation.R` | **FLAG**: requires `<sample>` CLI arg; renders nothing without it |
| Fig 7A,7B | `figures/figure7/01_xenium_forest_cox.R` | Cox forest OS/PFS |
| Fig 7C,7D | `figures/figure7/02_xenium_protein_correlation_row.py` | |
| Fig 7E,7F,7G | `figures/figure7/03_tcga_km_forest.R` | TCGA KM + stepwise Cox; BayesPrism covariate-set caveat (REPRODUCIBILITY.md P2-b) |

## Pending / placeholder panels (not yet assigned a manuscript number)

| Panel | Script | Notes |
|---|---|---|
| *TBD (expanded Fig 7)* | `figures/figure_icon7_bevacizumab/01_icon7_bev_slope_reversal.R` | **PLACEHOLDER** — ICON7/GSE140082 bevacizumab external validation (SecB prognostic under chemo, abolished by bev). Reads `cfg_obj("icon7_cohort")`; in-script per-arm Cox. Rename dir + assign panel once decided. Source module + full report: `…/2026_final_xenium_analysis/davids side quests/ICON7/`. |
| *TBD (supplementary / mechanistic)* | `figures/figure_ccle_smac_mimetics/02_prism_smac_mimetics.R` (main) + `03_crispr_targeted_nfkb.R` (support); prereq `01_score_hgsc_lines.R` | **PLACEHOLDER** — CCLE/DepMap 24Q4 + PRISM 25Q2. SecB-polarized HGSC lines are selectively sensitive to SMAC mimetics (LCL-161 ρ=−1.00 p=0.017; GDC-0152 ρ=−0.71 p=0.088; birinapant ρ=−0.50), but NOT to chemo controls; single-gene BIRC3 Chronos is null/opposite (ρ=+0.22) — IAP family redundancy interpretation — while upstream TRAF2 trends (ρ=−0.45). Reads `cfg_obj("ccle_line_scores")` + raw Chronos/PRISM. Source module + report: `…/davids side quests/ccle_depmap/`. |

## Supplementary figures

| Panel | Script | Notes |
|---|---|---|
| SF1A-C | `figures/supplementary/01_SF1_atlas_qc_metrics_by_study.py` | cache via `_prep/00_extract_atlas_obs.py` |
| SF2A | `figures/supplementary/02_SF2A_atlas_study_umap.py` | |
| SF2B | `figures/supplementary/03_SF2B_atlas_integration_study_composition.py` | cache via `_prep/00b_extract_integration_umaps.py` |
| SF2C | `figures/supplementary/04_SF2C_atlas_integration_per_study_umaps.py` | |
| SF3A-C | `figures/supplementary/05_SF3A_atlas_clinical_umaps.py` | |
| SF3D-F | `figures/supplementary/06_SF3D_atlas_genomic_umaps.py` | |
| SF4A | `figures/supplementary/07_SF4A_atlas_annotation_validation.py` | |
| SF4B | `figures/supplementary/08_SF4B_atlas_canonical_markers_dotplot.py` | reads `09f` dotplot_stats + lilra4_stats |
| SF4C | `figures/supplementary/09_SF4C_atlas_cnv_aneuploid_by_celltype.py` | n=251 |
| SF5 | `figures/supplementary/10_SF5_atlas_nmf_factor_umaps.py` | |
| SF6A-B | `figures/supplementary/11_SF6_atlas_secAB_expression_umaps.py` | SecA panels now use FBXO21 (signatures.yml) in place of source WT1 |
| SF7 | `figures/supplementary/12_SF7_atlas_cnv_validation.py` | **KNOWN ISSUE**: caption "3735 clones" is a typo for **375** (see REPRODUCIBILITY.md P2-a) |
| SF8 | `figures/supplementary/13_SF8_atlas_functional_heatmap.py` | scFEA flux label-map fixed |
| SF9A-E | `figures/supplementary/14_SF9_atlas_epitype_by_metadata.py` | **FLAG**: source emits 6 panels; published PDF shows 5 (met-site dropped) |
| SF10A | `figures/figure4/02_xenium_whole_tissue_snapshot.R` (`all` arg) | shared with Fig 4C |
| SF10B | `figures/supplementary/15_SF10B_tma_composition_spatial.R` | |
| SF10C | `figures/supplementary/16_SF10C_xenium_roi_and_core_celltype.R` | **FLAG**: low-confidence; no dedicated FTE-core gallery (cores mode defaults to TMA cores) |
| SF11 | `figures/supplementary/17_SF11_atlas_xenium_wt_polarization_violins.py` | **FLAG**: WT-only strip; combined 3-cohort assembler missing (per LINEAGE) |
| SF12,SF13,SF14 | `figures/supplementary/18_SF12_SF13_SF14_suppl_per_sample_gams.R` | SF12A-D / SF13A-S / SF14A-H |

## Supplementary data tables

| Item | Script | Notes |
|---|---|---|
| Supp Data 1 | `tables/05_supp_data_1_cell_metadata.py` | newly authored stub; run after `atlas/02_annotation/12a` treatment_status fix |
| Supp Data 2 | `tables/06_supp_data_2_level1_de.py` | newly authored stub (level-1 Wilcoxon DE) |
| Supp Data 3 | `tables/07_supp_data_3_ciliated_vs_secretory_de.py` | newly authored stub (logic from excluded 11v4 sandbox) |
| Supp Data 4 | `tables/01_supp_data_4_nmf_loadings.py` | NMF gene loadings (tidy reshape of `11d`) |
| Supp Data 5 | `tables/02_supp_data_5_xenium_gene_panel.R` | Xenium gene panel + probe QC |
| Supp Data 6 | `tables/03_supp_data_6_pathway_gene_sets.R` | 37 UCell pathway modules |
| Supp Data 7 | `tables/04_supp_data_7_autocrine_lr.py` | autocrine SecA/SecB L-R pairs; curated merge (Fig 5F caveat applies) |

## Prep helpers (shared caches, not panels)

- `figures/_prep/00_extract_atlas_obs.py` — atlas obs QC caches (SF1)
- `figures/_prep/00b_extract_integration_umaps.py` — integration UMAP caches (SF2B/C)
- `figures/_prep/01_export_tma_barcode_patient_map.py` — TMA barcode->patient map for Fig 4A/B
- `figures/_prep/fig_secretory_polarization_00_prepare_data.py` — shared `schema_nmf` `meta.parquet` (Fig 1F/G/H/I, 2C/D/E/G, SF5/6/9)

## Panel-map completeness

All numbered manuscript panels Fig 1A..7G, SF1..SF14, and Supp Data 1-7 are accounted for. The only panels **without an in-repo generator** are the documented external/known-gap panels: **Fig 3A,3B,3D,3E,3F,3G** (external organoid data, TRUST-EXISTING) and **Fig 3C** (external flow cytometry). **Fig 5F** has a generator and reproduces exactly, but carries a statistical-threshold caveat (see `docs/REPRODUCIBILITY.md`). See `docs/REPRODUCIBILITY.md` for the full audit verdict and known-issues list.
