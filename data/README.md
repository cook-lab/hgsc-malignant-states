# Data

No data are stored in this repository.

Data deposits:

- Processed atlas objects, Xenium, Visium HD, PDO snRNA-seq and KPCA.B scRNA-seq: **[accession to be added]**.
- The 13 public scRNA-seq studies of the atlas: their original accessions (Table S1).
- TCGA-OV: GDC.
- ICON7: GSE140082.

## How the code finds the data

Set `DATA_ROOT` to the folder that holds the deposit. All code outside the copies in `as_run/`, `models/` and
`icon7/` reads its inputs under `DATA_ROOT`, through `config/config.yml` or a path relative to `DATA_ROOT` written
in the script, and writes to `OUTPUT_ROOT`. The copies read the lab project folders they ran in; their READMEs list
those paths.

The code expects the deposit to keep the layout of the lab project folders. Paths below are relative to
`DATA_ROOT`:

| Folder | Contents |
|---|---|
| `2026_final_atlas/` | the atlas project (lab folder `202601_final_atlas`; Table S7 uses that name, see below) |
| `2026_final_xenium_analysis/` | the Xenium analysis project: probe QC, the TMA A clinical table, the TMA A labels object, the immunofluorescence workbook, the ICON7 folder |
| `202609_final_hgsc_spatial/` | raw Xenium runs, TMA B de-array tables, clinical workbooks, the TMA B clinical extract |
| `202603_organoids/` | the organoid project; `ORGANOIDS_ROOT` must point here (it does not follow `DATA_ROOT`) |
| `202608_kpca_secb_mouse/` | the KPCA.B mouse project |
| `202605_epitype_manuscript/` | the Xenium gene-panel table |
| `20260606_murine_secb_ihc/` | the KRT7 IHC QuPath project (an optional check of Fig 3l) |

Every file and folder named after a sample must carry its study code (WT01–WT08, TMA-P###, FTE01–FTE15), because the
code looks samples up by code. On the laboratory's data these names carry source identifiers.

**[not yet in the deposit plan]** marks an input that the code reads but the deposit plan does not list yet.

## Entry objects (`config.yml`, `objects:`)

| Key | Path | Role |
|---|---|---|
| `atlas_scanvi` | `2026_final_atlas/hgsc_atlas_scanvi.h5ad` | integrated atlas (scANVI) after the second doublet pass and the atlas UMAP; input of the level-1 annotation (`atlas/02_annotation/as_run/06_celltype_level1_leiden.py`); no figure reads it |
| `atlas_final` | `2026_final_atlas/hgsc_atlas_final.h5ad` | final annotated atlas (1,980,703 cells): every atlas figure, Fig 2f, Figure S4b, Tables S1–S3 |
| `xenium_ref` | `2026_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad` | SingleR reference (`spatial/03_annotation_polarization/01_annotation_singler.R`, `05_clean_split_rctd.R`); Figure S14; Table S7 |
| `sfe_dir` | `2026_final_xenium_analysis/output/sfe` | Xenium SFEs `sfe_WT01`–`sfe_WT08`, `sfe_tma` (TMA A), `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` (TMA B); see "Xenium" below |
| `icon7_cohort` | `2026_final_xenium_analysis/davids side quests/ICON7/data/processed/cohort_filtered.tsv` | ICON7 per-patient table (`icon7/scripts/04_descriptive.R`): Fig 7k–n, Table S12 |
| `icon7_results_dir` | `2026_final_xenium_analysis/davids side quests/ICON7/results/tables` | `17_chemo_arm_continuous_polarization.tsv`, `17b_continuous_polarization_per_arm.tsv` (`icon7/scripts/15_chemo_arm_analysis.R`): the model A check of Table S12. **[not yet in the deposit plan]** |
| `icon7_expr` | `2026_final_xenium_analysis/davids side quests/ICON7/data/processed/expr_gene_symbol.rds` | ICON7 gene × array matrix (`icon7/scripts/02_build_gene_matrix.R`): gene presence in Table S6. **[not yet in the deposit plan]** |

## What the deposit must hold, by data type

**Atlas** (`2026_final_atlas/`):
- lineage subsets `celltype_h5ad/hgsc_atlas_final_<lineage>.h5ad`. The epithelial subset is read by the Fig 1g prep,
  Fig 2a (scFEA), 2b, 2g, Fig 6h (`figures/figure6/10_fig6h_atlas_epithelial_pct.py`), Figure S6, Table S6 and the
  NMF re-implementation; the macrophage subset by Fig 6f,h (see "Inputs placed in `OUTPUT_ROOT`").
- NMF: `output/11d_epithelial_nmf/{11d_nmf_usage,11d_nmf_loadings,11d_factor_scores}.csv`; gene classes
  `output/11e_nmf_characterization/11e_gene_classification.csv`; consensus programs `output/11f_consensus_nmf/` and
  `output/11g_secAB_consensus_nmf/` (Table S5)
- epithelial UMAP and states: `output/fig_secretory_polarization/data/meta.parquet` (`UMAP1`/`UMAP2` for every
  epithelial UMAP; `schema_nmf` for Table S6 and Fig 6h)
- per-cell PROGENy / Hallmark / DoRothEA scores: `output/21_epitype_functional_characterization/` (Fig 2a, Figure S8)
- atlas UCell export: `output/18_ucell_atlas/`. It gives the Xenium polarization thresholds
  (`spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R`, unless `ATLAS_UCELL_DIR` is set), the atlas
  reference of Fig 3b–e and the gene space of Table S6.
- CopyKAT results: `output/19_cnv/per_sample/` (Fig 1g, 1l, 1m; Figures S4c and S7)
- TCGA-OV scores: `output/22_tcga_deconvolution/22d_signature_scores.csv` (Fig 7g–j)
- scFEA model files: `tools/scFEA/data/` (Fig 2a, Figure S8)
- data preparation, in `pre-processing and integration/`:
  - `processed/atlas_concat_counts_only_X.h5ad` and `processed/atlas_concatenated_filtered.h5ad` (Figure S1)
  - `processed/nath_2021.h5ad`, `atlas_metadata.xlsx` and `raw/nath_2021/GSE158722.cell_annotations.txt.gz`
    (`figures/_prep/nath2021_fix.py`)
  - `raw/<study>/`: Table S1 reads the GEO sample accessions from the raw file names.
- marker tables: `output/06_celltype_level1/markers/`, `output/08_celltype_level2/markers/` and
  `output/supplemental/{T1_atlas_metadata.csv, T2_celltypelevel1_markers.csv.csv, T3_epithelial_markeres.csv}`
  (Tables S1–S4). **[not yet in the deposit plan]**
- published per-cell labels: `output/spectrum_epithelial_crosstab/{GSE180661_GEO_cells.tsv.gz,
  geistlinger_cell_annotations.csv}` (Figure S9f, g). **[not yet in the deposit plan]**
- TCGA-OV expression: `data/cibersort_data_prev/tcga_ecotyper.txt` (Table S6, gene presence; see "TCGA-OV")
- `03_preprocess_hvg.py`, `04_cellassign.py` and `07_finalize.py` (`atlas/01_preprocess_qc/`) read
  `processed/atlas_concatenated_filtered.h5ad` directly under `2026_final_atlas/`, not under
  `pre-processing and integration/`; link the file there before running them.

**Xenium:**
- raw Xenium output folders, one per run, in `202609_final_hgsc_spatial/data/HGSC xenium /` (the folder name ends
  with a space): `WT02`–`WT04`, `WT06`–`WT08`, `WT05_WT01` (one slide, two tissues), `CHUM_TMA_1` and `CHUM_TMA_2`
  (TMA A), `COEUR_TMA_G` and `COEUR_TMA_H` (TMA B). Read by `spatial/01_build_sfe/` (without the morphology images;
  `RAW_XENIUM_DIR` overrides the folder there), Table S7 (`gene_panel.json`, `cell_feature_matrix.h5`), and the Fig 6
  ROI picker and Figure 6i (`cell_boundaries.parquet` of `WT05_WT01` and `WT02`).
- TMA B de-array tables: `202609_final_hgsc_spatial/data/dearray/{G,H}_cell_core_assignment.parquet` (`DEARRAY_DIR`)
- TMA A core and patient labels: the earlier processed object
  `2026_final_xenium_analysis/data/previously processed/xenium_tma_final_sfe_hdf5/` (`se.rds`, read by
  `spatial/01_build_sfe/05_rebuild_existing_cohort_rerun.R`). **[not yet in the deposit plan]**
- probe-QC outputs: `2026_final_xenium_analysis/output/05_probe_qc/{probe_qc_full.csv, flag_qc.csv,
  gene_exclusion_decisions.csv, genes_exclude_singler.txt}`, written by `spatial/02_qc/as_run/explore_probe_qc.Rmd`.
  Read by Table S7 and `spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R`; the exclusion list is copied into
  `OUTPUT_ROOT` (see below).
- the Xenium gene-panel table `202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv`
  (Tables S6 and S9; written by `tables/as_run/generate_xenium_gene_panel_table.R`)
- Table S7 also reads `2026_final_xenium_analysis/metadata/panel_markers.csv`, the panel design sheet
  `2026_final_xenium_analysis/data/Updated Xenium_2024_10.xlsx`, the SingleR reference under the lab folder name
  (`202601_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad`) and
  `202609_final_hgsc_spatial/hgsc-malignant-states/shared/signatures.yml` (the same file as this repository's
  `shared/signatures.yml`). **[not yet in the deposit plan]**
- SFEs (config `sfe_dir`). `load_sfe()` reads `$OUTPUT_ROOT/sfe` first, then `sfe_dir` (`RERUN_SFE_DIR` when it is
  set). The backend builds the SFEs from the raw runs into `$OUTPUT_ROOT/sfe`, so `sfe_dir` is needed only to start
  the Xenium figures without the backend. It must then hold the SFEs of the final run (with the SingleR, UCell and
  `pathway_*` columns), named as above. On the lab share the folder at this path holds SFEs from before the final
  re-run. Set `RERUN_SFE_DIR=$OUTPUT_ROOT/sfe` when running the backend: `spatial/02_qc/04_qc_filter_normalize.R`
  rewrites the SFEs in `sfe_dir` in place.

**Visium HD:** Space Ranger outputs of the WT01 section at 16 µm and 8 µm bins, and of a WT03 section at 16 µm.
`figures/figure6/01_fig6abcd_visium_hd_bins.py` and `03_fig6h_visium_hd_candidates.py` read both sections; the WT03
section enters only the per-gene concordance calls of `figures/figure6/11_fig6h_concordance.py`. The manuscript
reports the WT01 section only. The code reads the Space Ranger outputs from `OUTPUT_ROOT` (see below).

**Models:**
- PDOs and OPTO98 (`202603_organoids/`, through `ORGANOIDS_ROOT`):
  `output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds` (Fig 3a),
  `output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv` (Fig 3b) and
  `output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv` (Fig 3c–e). Fig 3f–h is drawn by
  `models/organoids/scripts/12h_panels_fgh_layout.py` from `output/12_GSE284507_spheroid_organoid/12_scores_all_sources.csv`.
- KPCA.B (`202608_kpca_secb_mouse/output/`): the handoff object `33_handoff/kpca_final_annotated.h5ad` (Fig 3i–k,
  Figure S10) and the NMF and alignment tables `31_fulltumour_nmf/` and `32_secb_align/` (Figure S10);
  `02_orthologs/{ortholog_map.csv, 02_coverage_summary.csv}` (Fig 3j, Figure S10) **[not yet in the deposit plan]**
- the KRT7 IHC snapshot of Fig 3l, given to `figures/figure3/06_fig3l_krt7_ihc_scale_bar.py` on the command line. The
  same script reads the scan sizes in `20260606_murine_secb_ihc/output/qupath_k7/data/*/server.json` for an optional
  check. **[not yet in the deposit plan: the QuPath project]**

**Clinical** (Figure 7a–f, Tables S8 and S11). **[Access route to be added.]** These tables carry the laboratory's
patient identifiers.
- TMA A clinical table: `2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv`
  (`spatial/06_clinical/as_run/`)
- TMA B clinical extract: `202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv`
- TMA clinical workbooks in `202609_final_hgsc_spatial/data/clinical data/`: `chum_tma_clinical_data_updated.xlsx`
  (TMA A) and `coeur_clinical_tmagh.xlsx` (TMA B). `spatial/02_qc/06b_final_tma_cores.py` reads them from
  `data/clinical data/` in the folder that contains the repository, not from `DATA_ROOT`.
- whole-tissue workbooks in the same folder, read by Table S8: the verified whole-tissue cohort workbook and
  `cooklab_hgsc_tissue_expansion.xlsx` (file names in `tables/08_tableS8_xenium_cohort.py`).
  **[not yet in the deposit plan]**
- immunofluorescence MFI per TMA A core: `2026_final_xenium_analysis/data/Copy of MFI all tissues TMA valid all markers_.xlsx`,
  sheet `Pour SPSS` (Figure 7e,f)

**TCGA-OV:** Fig 7g–j reads the scored table `2026_final_atlas/output/22_tcga_deconvolution/22d_signature_scores.csv`.
The bulk expression and clinical tables, `2026_final_atlas/data/cibersort_data_prev/{tcga_ecotyper.txt,
tcga_hla_clinical.csv}`, are the inputs of `atlas/07_deconvolution_survival/as_run/`; Table S6 also reads
`tcga_ecotyper.txt`.

**ICON7:** the folder `2026_final_xenium_analysis/davids side quests/ICON7/` (config keys above). The raw GEO files
`data/raw/GSE140082_{series_matrix.txt.gz, geo.normdata.gz}` are the inputs of `icon7/scripts/`.

**Scripts whose gene lists Table S6 parses** (under `DATA_ROOT`): `2026_final_atlas/scripts/{18b_ucell_atlas_scoring.R,
11d_epithelial_nmf.py, 18b_v2_ucell_atlas_scoring.R}`, `202603_organoids/scripts/{00_setup.R,
04v5_secB_only_comparison.py}` and `2026_final_xenium_analysis/davids side quests/ICON7/scripts/03_score_signatures.R`.
The copies in `atlas/03_epithelial_nmf/as_run/`, `models/organoids/scripts/` and `icon7/scripts/` carry the same
lists; `18b_v2_ucell_atlas_scoring.R` has no copy in this repository. **[not yet in the deposit plan]**

## Inputs placed in `OUTPUT_ROOT`

The code reads these from `OUTPUT_ROOT`. They are copies, files made by hand or interactively, or files that no
repository script writes in full. The folder READMEs give the commands.

| File or folder under `$OUTPUT_ROOT` | Source | Read by |
|---|---|---|
| `05_probe_qc/genes_exclude_singler.txt` | copy of the deposited probe-QC list | `spatial/03_annotation_polarization/01_annotation_singler.R`, `05_clean_split_rctd.R`; Figure S14; Table S7 |
| `_logs/singler_parallel.log` | console log of `spatial/03_annotation_polarization/01_annotation_singler.R` (`tee`) | Table S7 |
| `fig6_singler/inputs/pathway_gene_sets_v2.csv` | copy of `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv` | `figures/figure6/05_fig6fh_xenium_reciprocal_effects.R` |
| `fig6_singler/visium_local/{WT01B1, WT03B1}/` (16 µm: `filtered_feature_bc_matrix.h5`, `tissue_positions.parquet`) and `WT01B1_8um/` (`filtered_feature_bc_matrix.h5`, `spatial/tissue_positions.parquet`) | the deposited Visium HD Space Ranger outputs, staged by hand | `figures/figure6/01`–`03` |
| `atlas_local/hgsc_atlas_final_macrophage.h5ad` | a copy of the deposited macrophage subset (the default path; the subset's path can be given as the first argument instead) | `figures/figure6/08_fig6fh_atlas_macrophage_pseudobulk.py` |
| `fig6_singler/fig6_rois.json` | the Fig 6 ROI picker (`figures/figure6/fig6_roi_picker/`), interactive; `figures/figure6/README.md` gives the two entries used. **[not yet in the deposit plan]** | `figures/figure6/13_fig6_draw_panels.py` (the Visium window of 6a–d and the 6i ROI) |
| `fig6_singler/atlas_subtype_dotmatrix.csv` | written by `figures/figure6/08b_fig6f_atlas_subtype_dotmatrix.py` from the `atlas_pb/` sums of step 08; the final run's file. **[not yet in the deposit plan]** | `figures/figure6/13_fig6_draw_panels.py` (Fig 6f dot plot) |
| `fig6_singler/visium_n2_bins_WT01.csv.gz` | the final run's file: `figures/figure6/01_fig6abcd_visium_hd_bins.py`, then `01b_fig6c_visium_hd_bins_extra_genes.py` (CA9, CTSL, TIMP1, ITGB6 and ITGAV columns) rewrite it. **[not yet in the deposit plan]** | `figures/figure6/13_fig6_draw_panels.py` (Fig 6c) |
| `final_panels/figure7/Fig7_ABCD_{cox_2core_results, forest_fdr, km_logrank}.csv` | copies of Figure 7a–d outputs | Tables S8, S11 |
| `final_panels/supplementary/figures/SuppFig*` | copies of the supplementary renders (the S10 files are written there directly) | `figures/supplementary/23_figS1-18_composites.py` |
| `private_keys/SuppTable8_code_key_PRIVATE.csv` | written by Table S8 on the laboratory's data; private | scripts that meet source identifiers (not needed on coded data) |

Table S6 also reads the Fig 2g DE table from the repository's own `output/figure2/` folder and asserts that it equals
`$OUTPUT_ROOT/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv`; with `OUTPUT_ROOT` set elsewhere, copy it there.

## Intermediate outputs

Backend stages and figure-specific analyses write their results to `OUTPUT_ROOT`. Later scripts read them from
there, so run each figure folder in the order of its README. The largest are:
- re-run Xenium objects: `$OUTPUT_ROOT/sfe/`
- SingleR label tables: `$OUTPUT_ROOT/singler_recompute/`
- final TMA cores: `$OUTPUT_ROOT/07_core_qc/`
- Xenium raw-count export: `$OUTPUT_ROOT/fig6_viewer/degdata/`
- Figure 6 analysis tables: `$OUTPUT_ROOT/fig6_singler/`
- atlas loader caches: `$OUTPUT_ROOT/_prep_caches/`
- supplementary tables and the Figure 6 and S10 panel files: `$OUTPUT_ROOT/final_panels/`; other renders:
  `$OUTPUT_ROOT/figures/`
