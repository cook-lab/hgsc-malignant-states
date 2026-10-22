# tables: Supplementary Tables S1–S12

One script per table (one script writes S2–S4). Each table is an `.xlsx` workbook with a README sheet first, plus one
CSV per data sheet: `SuppTable<N>_<name>.xlsx` and `SuppTable<N>_<name>_<sheet>.csv`, written to
`$OUTPUT_ROOT/final_panels/supplementary/tables/` (`OUT` in `supp_xlsx.py`; `01` has its own writer and the same
folder). The scripts read finished upstream outputs. Most of them also check their values against the outputs of the
figure scripts and stop if they differ.

## Running

```bash
bash tables/build_all_tables.sh        # S1 → S12, in order; PY chooses the Python interpreter
```

- **Environment.** Export `DATA_ROOT` (the deposit root; the default is the placeholder `/path/to/deposit`),
  `OUTPUT_ROOT`, `LANG=en_US.UTF-8` and `LC_ALL=en_US.UTF-8` first. Give `OUTPUT_ROOT` as an absolute path, the one the
  figure scripts wrote to: `build_all_tables.sh` runs from `tables/`, and the default (`output`) is relative to the
  working directory. S6, S11 and S12 call `Rscript` from the `PATH`.
- **Paths.** Every script builds its paths from `DATA_ROOT`, `OUTPUT_ROOT` and the repository root. The deposit
  folders read are `$DATA_ROOT/2026_final_atlas/` (written `ATLAS` below; S7 reads the same atlas folder under its lab
  name, `202601_final_atlas/`), `2026_final_xenium_analysis/`, `202609_final_hgsc_spatial/`, `202605_epitype_manuscript/`
  and `202603_organoids/scripts/`. The "Reads" column gives the files.
- **Final TMA cores.** S8 and S11 read `$OUTPUT_ROOT/07_core_qc/final_tma_{cores,patients}.csv`
  (`spatial/02_qc/06b_final_tma_cores.py`). `FINAL_TMA_CORES` overrides the cores file.
- **Study codes.** S8 writes the code → source-identifier key to `$OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv`,
  never to the repository or the table folder. S10 and the figure scripts read the key from there, or from
  `$STUDY_CODE_KEY`, only when their inputs still carry source identifiers.
- **Order.** Run the upstream steps in the "Reads" column first. In particular:
  - S1 reads the Nath 2021 correction `$OUTPUT_ROOT/_prep_caches/nath2021_cell_correction.parquet`: run
    `figures/_prep/nath2021_fix.py` first.
  - S8 and S11 compare against the Figure 7a–d CSVs. After `figures/figure7/01`–`02`, copy `cox_2core_results.csv`,
    `xenium_tma_forest_fdr.csv` and `xenium_tma_km_logrank.csv` from `$OUTPUT_ROOT/figures/figure7/` to
    `$OUTPUT_ROOT/final_panels/figure7/` as `Fig7_ABCD_cox_2core_results.csv`, `Fig7_ABCD_forest_fdr.csv` and
    `Fig7_ABCD_km_logrank.csv`. S8 reads the first one, S11 all three.
  - S12 compares against `figures/figure7/06_fig7klmn_icon7_survival_row.R` (its CSVs in
    `$OUTPUT_ROOT/figures/figure7/`) and `icon7/scripts/15_chemo_arm_analysis.R` (its tables in the ICON7 folder,
    config `objects.icon7_results_dir`).
  - S10 needs `figures/figure6/01`, `02` (16 and 8 µm) and `04`; S9 needs `spatial/04_pathway_scoring/04_pathway_scoring.R`,
    `04c_pathway_scoring_allcells.R` and the count export `spatial/05_count_export/as_run/02_deg_precompute.R`.
  - Downstream: Figures S11 and S16–S18 read the S8 sheet `SuppTable8_xenium_clinical_cohort_Whole_tissues.csv`, and
    Figure S10 reads two S6 CSVs, from the table folder.
- **Runtime.** S6 reads all 575,366 epithelial cells and calls `Rscript`; S7 streams the raw Xenium count matrices
  (3–5 min); S11 and S12 run embedded R code through `Rscript`.

## Run order

| Step | Script | Table (sheets) | Reads | Upstream writer |
|---|---|---|---|---|
| 1 | `01_tableS1_atlas_studies.py` | S1, study characteristics (Studies, Samples, Cells_by_lineage) | `ATLAS/output/supplemental/T1_atlas_metadata.csv` (per-cell metadata); `ATLAS/hgsc_atlas_final.h5ad` (obs, cross-check); `ATLAS/pre-processing and integration/atlas_metadata.xlsx`; raw file names under `ATLAS/pre-processing and integration/raw/<study>/` (GSM accessions); `$OUTPUT_ROOT/_prep_caches/nath2021_cell_correction.parquet` | T1 export **[code missing]**, corrected by `atlas/02_annotation/as_run/fix_metastatic_site_docs.py`; atlas object: `atlas/02_annotation/as_run/09d_publication_figures.py`, `atlas/02_annotation/as_run/12a_fix_zhang2022_treatment_status.py`, `atlas/02_annotation/as_run/fix_metastatic_site_vazquez.py`; cache: `figures/_prep/nath2021_fix.py` |
| 2 | `02_tableS2-4_annotation_markers.py` | S2, level-1 markers (Top50, All_markers); S3, level-2 markers (Top50, All_markers, Cell_types); S4, ciliated vs secretory (All_DEGs, Canonical_markers) | S2: `ATLAS/output/06_celltype_level1/markers/level1_<type>.csv` (13 files), checked against `ATLAS/output/supplemental/T2_celltypelevel1_markers.csv.csv`. S3: `ATLAS/output/08_celltype_level2/markers/*.csv`, restricted to the 51 final level-2 types of `ATLAS/hgsc_atlas_final.h5ad`. S4: `ATLAS/output/supplemental/T3_epithelial_markeres.csv` | S2: **[code missing]** (level-1 Wilcoxon export; level-1 annotation in `atlas/02_annotation/as_run/06_celltype_level1_leiden.py`, `06b_remap_level1.py`). S3: `atlas/02_annotation/as_run/08_finalize_celltype_level2.py` writes the marker files; they are renamed by the atlas project's `08b_rename_level2_labels.py` (not in this repository), `atlas/02_annotation/as_run/08c_rename_level2_labels_v2.py` and `10c_rename_inflammatory_macrophage.py`. S4: **[code missing]** |
| 3 | `05_tableS5_epithelial_nmf_programs.py` | S5, NMF programs (Factor_annotation, Loadings_long, Top50_wide, Top50_specific_wide, SecA_SecB_classification, Consensus_metaprograms, Consensus_MP_genes, SecAB_restricted_MPs, SecAB_restricted_MP_genes) | `ATLAS/output/11d_epithelial_nmf/{11d_nmf_loadings,11d_factor_scores}.csv`; `ATLAS/output/11e_nmf_characterization/11e_gene_classification.csv` (recomputed from 11d and asserted equal); `ATLAS/output/11f_consensus_nmf/`, `ATLAS/output/11g_secAB_consensus_nmf/`; factor labels parsed (not run) from `figures/supplementary/11_figS5_nmf_factor_umaps.py` and `figures/figure1/07_fig1j_program_factor_umaps.py` | `atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py` (re-implemented by `atlas/03_epithelial_nmf/01_epithelial_nmf.py` and verified against the deposit by `01b_verify_against_deposited.py`), `11e_nmf_characterization.py`, `11f_consensus_nmf.py`, `11g_secAB_consensus_nmf.py` |
| 4 | `06_tableS6_secA_secB_signatures.py` | S6, SecA/SecB signatures (Scoring_signatures, Signature_usage, Extended_lists, Other_named_markers) | `shared/signatures.yml`; 11d loadings and 11e classes (177 + 177); the Fig 2g DESeq2 table `$OUTPUT_ROOT/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv`, asserted identical to the copy in the repository's `output/figure2/`; `ATLAS/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` counts and `ATLAS/output/fig_secretory_polarization/data/meta.parquet` (`schema_nmf`; % expressing); the Xenium panel table `$DATA_ROOT/202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv`; `$OUTPUT_ROOT/06d_annotation_noBCAM/shared_gene_list_xenium_atlas.txt` (457 genes); `ATLAS/output/18_ucell_atlas/shared_gene_list.txt` (14,271 genes); TCGA `ATLAS/data/cibersort_data_prev/tcga_ecotyper.txt` and ICON7 `expr_gene_symbol.rds` (path of config `objects.icon7_expr`; gene presence only); the gene lists of the scripts listed under S6 below | `figures/figure2/04_fig2g_pseudobulk_volcano.py`; `atlas/03_epithelial_nmf/as_run/fig_secretory_polarization_00_prepare_data.py`; `tables/as_run/generate_xenium_gene_panel_table.R`; `spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R`; `atlas/03_epithelial_nmf/as_run/18b_ucell_atlas_scoring.R`; `icon7/scripts/02_build_gene_matrix.R` |
| 5 | `07_tableS7_xenium_gene_panel.py` | S7, Xenium panel and probe QC (Gene_panel, Summary) | `gene_panel.json` and `cell_feature_matrix.h5` of the raw Xenium runs in `$DATA_ROOT/202609_final_hgsc_spatial/data/HGSC xenium /`; `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz`; `$DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/{probe_qc_full,flag_qc,gene_exclusion_decisions}.csv`; `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (copy of the notebook's list); the SingleR reference `$DATA_ROOT/202601_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad` (gene names); the SingleR console log `$OUTPUT_ROOT/_logs/singler_parallel.log`; `$DATA_ROOT/2026_final_xenium_analysis/metadata/panel_markers.csv` and `$DATA_ROOT/2026_final_xenium_analysis/data/Updated Xenium_2024_10.xlsx` (design sheet); `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv`; `signatures.yml`, read from `$DATA_ROOT/202609_final_hgsc_spatial/hgsc-malignant-states/shared/` (the same file as this repository's `shared/signatures.yml`) | `spatial/02_qc/as_run/explore_probe_qc.Rmd`; `spatial/03_annotation_polarization/01_annotation_singler.R` (log), `06_extract_singler.R`; `atlas/08_xenium_reference/as_run/xenium_celltype_markers_and_report.py`; `spatial/04_pathway_scoring/04_pathway_scoring.R` |
| 6 | `08_tableS8_xenium_cohort.py` | S8, Xenium cohort (Whole_tissues, TMA_patients, TMA_cores, TMA_core_selection, FTE_cores, Summary) | SingleR cell tables of WT01–WT08, `tma`, `COEUR_TMA_G`, `COEUR_TMA_H`; final TMA core tables (through `figures/_prep/final_tma_cores.py`); in `$DATA_ROOT/202609_final_hgsc_spatial/data/clinical data/`: the verified whole-tissue cohort workbook, `cooklab_hgsc_tissue_expansion.xlsx` (histology check), `chum_tma_clinical_data_updated.xlsx`, `coeur_clinical_tmagh.xlsx`; the TMA A and TMA B clinical tables used by Figure 7 (`$DATA_ROOT/2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv`, `$DATA_ROOT/202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv`); `$OUTPUT_ROOT/final_panels/figure7/Fig7_ABCD_cox_2core_results.csv` | `spatial/03_annotation_polarization/06_extract_singler.R`; `spatial/02_qc/06b_final_tma_cores.py`; `spatial/06_clinical/as_run/explore_clinical_associations.Rmd`; TMA B table **[code missing]**; `figures/figure7/01`–`02` + hand copy |
| 7 | `09_tableS9_ucell_gene_sets.py` | S9, gene sets (Module_summary, Gene_sets_long) | `$OUTPUT_ROOT/9b_scoring/{pathway_gene_sets_v2,pathway_coverage_allcells,scoring_summary_allcells}.csv`; gene lists parsed (not run) from `spatial/04_pathway_scoring/04b_pathway_gene_sets.R` (asserted identical), `figures/figure5/01_fig5a_epithelial_programs_build.R` (DDR, EPIS), `figures/figure5/08_fig5i_macrophage_programs_build.R` (macrophage sets), `figures/figure6/06_fig6g_macrophage_secB_neighbours.py` (GLYCO, RESID); the Xenium panel table; `$OUTPUT_ROOT/fig6_viewer/degdata/WT01_genes.txt` (gene-list check) | `spatial/04_pathway_scoring/04_pathway_scoring.R` (gene sets), `04c_pathway_scoring_allcells.R` (coverage, summary); `tables/as_run/generate_xenium_gene_panel_table.R`; `spatial/05_count_export/as_run/02_deg_precompute.R` |
| 8 | `10_tableS10_secb_niche_macrophage_degs.py` | S10, SecB-niche macrophage DE (Visium_16um, Visium_8um, Visium_top_macrophage_genes, Xenium_niche_vs_control, Xenium_macrophage_counts) | `$OUTPUT_ROOT/fig6_singler/visiumWT01_niche_vs_distal_volcano.csv` and `_8um.csv`; `visium_n2_bins_WT01.csv.gz` (bin-class counts); `mac_signature_deg_WT.csv`, `mac_counts_WT.csv` | `figures/figure6/02_fig6e_visium_hd_niche_volcano.py` (16, 8), `01_fig6abcd_visium_hd_bins.py`, `04_fig6fhi_xenium_niche_macrophages.R` |
| 9 | `11_tableS11_xenium_tma_cox.py` | S11, TMA multivariable Cox (Multivariable_Cox, KM_tertiles, Full_models, PH_tests, Cohort) | The inputs of `figures/figure7/01_fig7abcd_tma_cox_km_tables.R` (paths in the embedded R code); `spatial/00_setup/final_tma_cores.R`; `$OUTPUT_ROOT/final_panels/figure7/Fig7_ABCD_{cox_2core_results,forest_fdr,km_logrank}.csv` | `figures/figure7/01`–`02` + hand copies |
| 10 | `12_tableS12_icon7_cox.py` | S12, ICON7 Cox (Models, Full_coefficients, Group_sizes, KM_medians) | `cohort_filtered.tsv` (config `objects.icon7_cohort`); `17_chemo_arm_continuous_polarization.tsv`, `17b_continuous_polarization_per_arm.tsv` (config `objects.icon7_results_dir`); `icon7_km_pfs_within_group_hr.csv`, `icon7_timevarying_forest_15mo_data.csv` in `$OUTPUT_ROOT/figures/figure7/` | `icon7/scripts/04_descriptive.R`, `15_chemo_arm_analysis.R`; `figures/figure7/06_fig7klmn_icon7_survival_row.R` |

`supp_xlsx.py` is the shared writer (`write_table()`) of S2–S12.

## What each table checks or adds

- **S1.** Publication fields (accession, platform, chemistry, tissue processing) are written in the script
  (`STUDY_INFO`), each with the source it was verified against. Totals: 13 studies, 158 patients, 362 samples,
  1,980,703 cells.
- **S2–S4.** These are the marker tables made during atlas annotation, not a re-test on the final atlas: 51,592
  level-1 pairs; 197,160 level-2 pairs over 49 types (pericytes and neutrophils have a single level-2 type); 10,029
  ciliated vs secretory DEGs. The level-2 Wilcoxon code is in `atlas/02_annotation/as_run/08_finalize_celltype_level2.py`.
- **S5.** k = 10 on 3,007 genes; 177 SecA-specific, 177 SecB-specific and 124 shared genes; consensus NMF k = 5–10,
  30 runs per k. Renaming the two parsed figure scripts breaks S5.
- **S6.** Parses the SecA/SecB lists of the scripts that scored them and asserts each equals `signatures.yml`
  (`Signature_usage` sheet): `figures/figure6/01_fig6abcd_visium_hd_bins.py`, and the copies in the lab project folders
  under `$DATA_ROOT`: `2026_final_atlas/scripts/18b_ucell_atlas_scoring.R`, `202603_organoids/scripts/00_setup.R` and
  `04v5_secB_only_comparison.py`, and `2026_final_xenium_analysis/davids side quests/ICON7/scripts/03_score_signatures.R`.
  It also reads, without asserting, the lists of `2026_final_atlas/scripts/11d_epithelial_nmf.py` and
  `18b_v2_ucell_atlas_scoring.R`. The copies of these scripts in this repository (`atlas/03_epithelial_nmf/as_run/`,
  `models/organoids/scripts/`, `icon7/scripts/`) carry the same lists; `18b_v2_ucell_atlas_scoring.R` has no copy here.
- **S7.** 477 genes; 31 fail probe QC (≥ 3 flags); 28 are excluded from SingleR (429 genes used); 20 are Xenium-only.
  It asserts that every run folder under `$DATA_ROOT/202609_final_hgsc_spatial/data/HGSC xenium /` carries the same 477
  genes, and that the gene counts in the SingleR log match.
- **S8.** Builds the study codes. On the laboratory's data it assigns codes with a fixed seed (20260928) and writes the
  key as described under "Study codes"; when the key exists it asserts the codes are unchanged. Age is released in
  10-year bands and OS/PFS as months. It re-implements the Figure 7a–d cohort rules and asserts OS n = 104 and
  PFS n = 103. 378 patients / 503 cores; TMA A 101 and TMA B 317 patients profiled. The neoadjuvant (n = 5) and
  non-HGSC (n = 4) exclusions of the Methods are recorded in sheet `TMA_patients`.
- **S9.** 37 UCell modules (325 gene entries) and the custom sets of Figures 5a, 5i and 6g.
- **S10.** Applies the same macrophage-content adjustment, carry-over rule and stromal class as the 6e volcano. The
  Xenium sheets are the pseudobulk limma-voom model of `04` (7 whole tissues with ≥ 50 macrophages per group; WT07, with
  39 SecB-niche macrophages, is left out). Sections are labelled through `config/study_codes.py`.
- **S11.** Runs a copy of the model code of `01_fig7abcd_tma_cox_km_tables.R` (embedded R; `ST11_WORK` and
  `FINAL_TMA_CORES_R` are set by the script), then exports every coefficient, the feature SDs, Schoenfeld-residual
  tests (`cox.zph`) and the tertile Kaplan–Meier comparison. Hazard ratios, q values and KM rows must match the three
  `Fig7_ABCD_*` CSVs to 1e-8. Intermediates go to `$OUTPUT_ROOT/supp_tables/ST11/`.
- **S12.** Embedded R (`ST12_WORK`, `ST12_COHORT` set by the script), on the n = 191 FIGO III/IV complete cases; every
  model is adjusted for FIGO stage, debulking status and age:
  - model A: the per-arm continuous-polarization Cox (OS and PFS), rebuilt from section B.2 of
    `icon7/scripts/15_chemo_arm_analysis.R` and checked against its two output tables;
  - model B: the bevacizumab effect on PFS within each polarization group, with KM medians and log-rank p (Figure
    7k,l), checked against `icon7_km_pfs_within_group_hr.csv`;
  - model C: the 15-month time-varying bevacizumab effect with Benjamini–Hochberg q and the period interaction
    (Figure 7m,n), checked against `icon7_timevarying_forest_15mo_data.csv`;
  - model D: the treatment × polarization-group interaction for PFS with its likelihood-ratio test.

  Intermediates go to `$OUTPUT_ROOT/supp_tables/ST12/`.

## as_run/

`as_run/generate_xenium_gene_panel_table.R` is a copy of the script that ran in the lab project folder
`202605_epitype_manuscript` (`supplemental_material/scripts/`). It is not trimmed: all of it builds, summarises or
writes the panel table (its summary prints the 477 genes, the 28 excluded and the 20 Xenium-only). It reads the
probe-QC outputs in the Xenium project folder `2026_final_xenium_analysis`
(`output/05_probe_qc/{gene_exclusion_decisions,probe_qc_full}.csv`) and writes
`supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv` in `202605_epitype_manuscript`: the 477 panel genes with
QC status and `in_analysis_panel`. S6 and S9 read that table; S7 rebuilds its own panel table from the same probe-QC
outputs. Its paths are those of the original project.

## Missing code and manual steps

- **[code missing]** The interactive export that wrote `ATLAS/output/supplemental/T1_atlas_metadata.csv`, the per-cell
  metadata behind S1 (later corrected in place by `atlas/02_annotation/as_run/fix_metastatic_site_docs.py`).
- The level-1 markers `ATLAS/output/06_celltype_level1/markers/level1_<type>.csv` (S2) are written by
  `atlas/02_annotation/as_run/xenium_coverage_report_level1.py`; `T2_celltypelevel1_markers.csv.csv`, which S2 checks
  against, is the same files stacked (the stacking step is not recorded).
- **[code missing]** The export that wrote `ATLAS/output/supplemental/T3_epithelial_markeres.csv` (S4). Its
  Ciliated-up rows equal the ciliated-cluster markers of the atlas project's `11v2_deg.py` / `11v4_deg.py`; the export
  that added the Secretory-up rows and the percentages is not recorded (see `atlas/02_annotation/README.md`).
- The TMA B (COEUR) clinical table read by S8 and Figure 7a–d is written by
  `spatial/06_clinical/as_run/coeur_clinical_extract.py`. S8 re-derives the same fields from the COEUR workbook and
  asserts that they are equal.
- **[code missing]** Download and preparation of the TCGA-OV expression table `tcga_ecotyper.txt` (S6 gene presence;
  Figure 7g–j).
- Hand copies:
  - the three `Fig7_ABCD_*.csv` files in `$OUTPUT_ROOT/final_panels/figure7/` (S8, S11; see "Order");
  - the Fig 2g DESeq2 table in the repository's `output/figure2/` (S6), unless `OUTPUT_ROOT` is that `output/` folder;
  - `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt`, a copy of the probe-QC notebook's list (S7);
  - `$OUTPUT_ROOT/_logs/singler_parallel.log`, the console output of `spatial/03_annotation_polarization/01_annotation_singler.R`
    (S7). No script writes it.
