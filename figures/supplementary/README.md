# Supplementary Figures S1–S18

This folder holds the generators of all 18 supplementary figures and the script that assembles the lettered pages.
The figures fall into three groups:
- **Figures S1–S9** are atlas figures. They use the same loaders as Figures 1–2 (`figures/_prep/`).
- **Figure S10** compares the KPCA.B mouse tumour programs with the human atlas programs. It depends on the mouse
  pipeline in [`../../models/kpca_mouse/`](../../models/kpca_mouse/).
- **Figures S11–S18** are Xenium figures. They depend on the Xenium backend in [`../../spatial/`](../../spatial/) and,
  for S16–S18, on the per-cell tables of Figures 4 and 5.

| Figure | Content |
|---|---|
| S1 | atlas QC metrics by study |
| S2 | study contribution |
| S3 | clinical and genomic metadata on the atlas UMAP |
| S4 | annotation validation |
| S5 | NMF factor usage |
| S6 | SecA and SecB genes on the epithelial UMAP |
| S7 | states within CopyKAT clones |
| S8 | functional heatmaps |
| S9 | state composition by metadata and published labels |
| S10 | KPCA.B programs vs human programs |
| S11–S13 | Xenium cell-type maps (whole tissues; TMA A and FTE cores; TMA B) |
| S14 | Xenium annotation dot plot |
| S15 | polarization distribution across the Xenium cohort |
| S16–S18 | per-tissue GAMs (environment, epithelium, macrophage) |

## Run order

**Requires.**
- **Figures S1–S9.**
  - The deposited atlas objects ([`../../atlas/README.md`](../../atlas/README.md)).
  - The loaders `figures/_prep/nath2021_fix.py` → `figures/_prep/atlas_final_tables.py`
    ([`../_prep/README.md`](../_prep/README.md)).
  - For S8, the Fig 2a prep (`figures/_prep/fig2a_scfea_all_cells.py` → `figures/_prep/fig2a_functional_state_summaries.py`).
    S8 rebuilds the per-state z-scores itself; when the Fig 2a prep tables exist, it checks that they agree and draws
    those. The per-cell scFEA flux always comes from `fig2a_scfea_all_cells.py`.
- **Figure S10.**
  - The KPCA.B pipeline in `models/kpca_mouse/scripts/`:
    - `00_setup.py`
    - QC and annotation: `01`, `15b`, `18`, `24`, `27`
    - the hand-off object: `33_consolidate_handoff.py` → `33_handoff/kpca_final_annotated.h5ad`
    - orthologues: `02_orthologs.py` → `02_orthologs/ortholog_map.csv`
    - the full-tumour NMF: `31_fulltumour_nmf.py` with `04_mouse_nmf.py`
    - the alignment to the human factors: `32_identify_secb_align.py` with `05_crossspecies_alignment.py`
  - The human `11d_nmf_loadings.csv` (deposited; `atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py`).
  - Table S6 (`tables/06_tableS6_secA_secB_signatures.py`), which writes the two CSVs that S10 reads.
- **Figures S11–S18.** The Xenium backend, up to the SingleR label tables (see `spatial/`):
  - SFE build: `spatial/01_build_sfe/`
  - cell QC: `spatial/02_qc/04_qc_filter_normalize.R`
  - probe QC: `spatial/02_qc/as_run/explore_probe_qc.Rmd`. Its exclusion list is copied by hand to
    `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (see "Manual steps").
  - SingleR and UCell polarization: `spatial/03_annotation_polarization/01`–`04`
  - label tables: `spatial/03_annotation_polarization/06_extract_singler.R`, run once per SFE →
    `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz`
  - the final TMA core set: `spatial/02_qc/05b_core_qc.R`, run once per array → `spatial/02_qc/06b_final_tma_cores.py`
    → `$OUTPUT_ROOT/07_core_qc/final_tma_{cores,patients}.csv`
  - for S16–S18: the pathway scores of `spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R`, and the Fig 4f and
    Fig 5 per-cell builders listed in step 22
  - for S11 and S16–S18: Table S8 (`tables/08_tableS8_xenium_cohort.py`), which writes the whole-tissue sheet
    `$OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable8_xenium_clinical_cohort_Whole_tissues.csv` that these
    scripts read. Table S8 needs the Figure 7a–d results first (see [`../../tables/README.md`](../../tables/README.md)).

Unless a different folder is given, renders go to `$OUTPUT_ROOT/figures/supplementary/` as `.svg` and `.png`. Paths
without a root are under `$DATA_ROOT/2026_final_atlas/`.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `00_figS1_extract_atlas_obs.py` | `pre-processing and integration/processed/atlas_concat_counts_only_X.h5ad` (before QC) and `atlas_concatenated_filtered.h5ad` (after the UMI/gene filter, before the 0.25 Scrublet cut), read in chunks | `$OUTPUT_ROOT/_prep_caches/atlas_obs_{prefilter,postfilter}.parquet` | S1 (numbers) |
| 2 | `01_figS1_atlas_qc_metrics_by_study.py` | step 1 caches | `SF1A_atlas_qc_total_counts_by_study`, `SF1B_atlas_qc_genes_by_study`, `SF1C_atlas_qc_doublet_score_by_study` | S1a–c |
| 3 | `02_figS2a_atlas_study_umap.py` | `load_obs()` | `SF2A_atlas_study_umap` | S2a |
| 4 | `03_figS2b_study_composition.py` | `load_obs()` | `SF2B_atlas_integration_study_composition` | S2b |
| 5 | `04_figS2c_per_study_umaps.py` | `load_obs()` | `SF2C_atlas_integration_per_study_umaps` | S2c |
| 6 | `05_figS3abc_clinical_umaps.py` | `load_obs()` | `SF3_treatment_status`, `SF3_anatomic_site`, `SF3_metastatic_site` | S3a, b, c |
| 7 | `06_figS3def_genomic_umaps.py` | `load_obs()` | `SF3_tp53_status`, `SF3_hrd_status`, `SF3_brca_status` | S3d, e, f |
| 8 | `07_figS4a_annotation_validation.py` | `load_obs()`: `celltype_pred`, `celltype_probability`, `celltype_confident` | `SF4A_atlas_annotation_validation` | S4a |
| 9 | `08_figS4b_marker_stats.py` | `hgsc_atlas_final.h5ad` X (streamed in chunks); `load_obs()` | `$OUTPUT_ROOT/02_annotation/09f_dotplot_canonical_markers/{dotplot_stats,lilra4_stats}.csv` | S4b (numbers) |
| 10 | `09_figS4b_canonical_markers_dotplot.py` | step 9 | `SF4B_atlas_canonical_markers_dotplot` | S4b |
| 11 | `10_figS4c_cnv_aneuploid_by_state.py` | `output/19_cnv/per_sample/*/copykat_prediction.csv` (251 samples); `load_obs()`, `load_epithelial()` | `SuppFig4c_atlas_cnv_aneuploid_epithelial_vs_reference`; `$OUTPUT_ROOT/figures/supplementary/tables/SuppFig4c_cnv_aneuploid_by_epithelial_state.csv` | S4c |
| 12 | `11_figS5_nmf_factor_umaps.py` | `load_epithelial()` | `SF5_atlas_nmf_factor_umaps` | S5 |
| 13 | `12_figS6_secAB_expression_umaps.py` | `load_epithelial()`; epithelial expression from `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad`; `shared/signatures.yml` | `SF6a_atlas_SecA_signature_umaps`, `SF6b_atlas_SecB_signature_umaps` | S6a, b |
| 14 | `13_figS7_cnv_clone_validation.py` | CopyKAT `copykat_subclones.csv`, `barcodes.csv`, `copykat_subclone_qc.csv`, `DONE.txt`; `load_obs()`, `load_epithelial()` | `SF7_atlas_cnv_validation`; `$OUTPUT_ROOT/_prep_caches/sf7_cnv/{within_clone_coexistence, per_sample_verdict, sample_manifest}.csv` | S7. Run it before `figures/_prep/text_state_coexistence.py`. |
| 15 | `14_figS8_functional_heatmap.py` | `output/21_epitype_functional_characterization/21_{progeny,hallmark,dorothea}_per_cell.parquet`; `$OUTPUT_ROOT/figure2/fig2a_scfea_all/21_flux_per_cell.parquet` and `$OUTPUT_ROOT/figure2/fig2a_functional/21_*_zscored.csv` (Fig 2a prep); `load_epithelial()`; `tools/scFEA/data/Human_M168_information.symbols.csv` | `$OUTPUT_ROOT/_prep_caches/sf8_functional/21_*_{means,zscored}.csv`; `SF8_atlas_functional_heatmap` | S8 |
| 16 | `15_figS9_composition_by_metadata.py` | `load_epithelial()`; `output/spectrum_epithelial_crosstab/GSE180661_GEO_cells.tsv.gz` (f); `output/spectrum_epithelial_crosstab/geistlinger_cell_annotations.csv` (g) | `SF9{a..g}_atlas_epitype_by_<metadata>`; `SF9_atlas_epitype_by_metadata_counts.csv` (n and % per bar) | S9a–g |
| 17 | `16_figS10_kpca_mouse_vs_human_programs.py` | mouse project outputs under `$DATA_ROOT/202608_kpca_secb_mouse/output/` (`33_handoff/kpca_final_annotated.h5ad`, `31_fulltumour_nmf/`, `32_secb_align/`, `02_orthologs/`); `output/11d_epithelial_nmf/11d_nmf_loadings.csv`; `shared/signatures.yml`; `$OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable6_secA_secB_scoring_signatures_{Scoring_signatures,Extended_lists}.csv` (Table S6) | in `$OUTPUT_ROOT/final_panels/supplementary/figures/`: `SuppFig10{a,b,c,d}_*`; in `$OUTPUT_ROOT/final_panels/supplementary/tables/`: `SuppFig10_kpca_{ortholog_expression, program_alignment}.csv` (source data of a and c). Both go to `$SF10_OUT` when it is set. | S10a–d |
| 18 | `17_figS11-13_xenium_celltype_maps.R [all\|wt\|tmaA\|fte\|tmaB\|legend]` | `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz`; final cores (`spatial/00_setup/final_tma_cores.R`); the Table S8 whole-tissue CSV (`wt`) | in `$OUTPUT_ROOT/figures/supplementary/SF10_SF11_singler/`: `SuppFig11_xenium_whole_tissue_celltype` (`wt`), `SuppFig12a_xenium_TMA_A_paired_cores_celltype` (`tmaA`), `SuppFig12b_xenium_FTE_cores_celltype` (`fte`), `SuppFig13_xenium_TMA_B_cores_celltype` (`tmaB`), `SuppFig1{1,2,3}_celltype_legend` (`legend`), `SF10_SF11_SF12_counts.csv` | S11, S12a, S12b, S13 |
| 19 | `18_figS14_annotation_validation_extract.R` | SFEs (`colData` `singler_pruned`; `logcounts`); label tables | `$OUTPUT_ROOT/_prep_caches/sf12/<sample>_{agg_groups, agg_sumlog, agg_nexpr}` | S14 (numbers) |
| 20 | `19_figS14_annotation_validation_plot.py` | step 19 caches; `obj("xenium_ref")` (SingleR reference, for the Xenium-only genes); `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt`; final cores (`figures/_prep/final_tma_cores.py`) | `SF12B_xenium_canonical_markers_dotplot` (**S14**); `$OUTPUT_ROOT/_prep_caches/sf12/sf12b_dotplot_stats.csv` | S14 |
| 21 | `20_figS15_polarization_distribution.py` | label tables; final cores (`figures/_prep/final_tma_cores.py`) | in `$OUTPUT_ROOT/figures/supplementary/SF13_singler/`: `SuppFig15_xenium_polarization_distribution` (both panels), `SF13_units.csv`, `SF13_arms.csv`, `SF13_tests.csv` | S15a, b |
| 22 | `21_figS16-18_per_sample_gams.R` | `$OUTPUT_ROOT/singler_recompute/singler_gam_percell.csv` (`figures/figure4/08_fig4f_gam_features.R`); `$OUTPUT_ROOT/split_recompute/fig5a_epi_program_scores_presplit.csv` (`figures/figure5/01_fig5a_epithelial_programs_build.R`) and `fig5e_macrophage_percell_presplit.csv` (`figures/figure5/08_fig5i_macrophage_programs_build.R`); SFE `pathway_*` scores and counts; for a pooled-vs-main check, the main-figure curves in `$OUTPUT_ROOT/figures/figure4/gam_curves.csv` (`figures/figure4/09_fig4f_gam_fit.R`) and `$OUTPUT_ROOT/figures/figure5/{gam_curves_5a,fig5e_mac_gam_curves}.csv` (`figures/figure5/02_fig5a_epithelial_programs_fit.R`, `figures/figure5/09_fig5i_macrophage_programs_fit.R`); the Table S8 whole-tissue CSV | `$OUTPUT_ROOT/supplementary/sf14_16_per_sample_gams/{curves_pooled, curves_per_sample, per_sample_consistency, pooled_vs_main, sample_meta}.csv`; `extra_{epi,mac}_percell.csv.gz` (cache) | S16–S18 (numbers) |
| 23 | `22_figS16-18_per_sample_gams_plot.py [out_dir]` | step 22 | in `out_dir` (default `$OUTPUT_ROOT/supplementary/sf14_16_per_sample_gams/figures/`): `SuppFig16_per_sample_gams_environment`, `SuppFig17_per_sample_gams_epithelial`, `SuppFig18_per_sample_gams_macrophage` (each with its tissue legend); and `$OUTPUT_ROOT/supplementary/sf14_16_per_sample_gams/SuppFig16-18_per_sample_consistency.csv` | S16, S17, S18 |
| 24 | manual: copy the renders into `$OUTPUT_ROOT/final_panels/supplementary/figures/` under the names in "Manual steps" | – | – | – |
| 25 | `23_figS1-18_composites.py [N …]` | the `SuppFig<N>_*.png` files in `$OUTPUT_ROOT/final_panels/supplementary/figures/` | `$OUTPUT_ROOT/final_panels/supplementary/figures/composites/SuppFig<N>_composite.png` (7.5-inch page, 300 dpi; it adds the bold panel letters). A figure whose files are missing is skipped. | S1–S18 pages |

**Order.**
- Steps 1–16 (atlas) can run in any order, with three exceptions: 1 comes before 2, 9 comes before 10, and the
  Fig 2a prep comes before 15.
- Step 17 needs Table S6 (`tables/06_…`).
- Steps 18–23 need the Xenium backend first. Step 19 comes before 20, and 22 before 23. Step 18 (`wt`) and step 22
  need Table S8 (`tables/08_…`).
- Step 25 needs the copies of step 24; the S10 files are already in place.

**Script names and output names.** The script names and header comments carry the final figure numbers. Some
internal labels and output names still carry the earlier numbers; each script header gives the mapping.

| Script | Output names it writes | Final figure |
|---|---|---|
| `17_…` | folder `SF10_SF11_singler/`, `SF10_SF11_SF12_counts.csv` (its `panel` column: SF10, SF11a, SF11b, SF12) | S11–S13 |
| `18_…`, `19_…` | `SF12B_xenium_canonical_markers_dotplot`, cache folder `sf12/` | S14 |
| `20_…` | folder `SF13_singler/`, tables `SF13_*.csv` | S15 |
| `21_…`, `22_…` | internal ids SF14 / SF15 / SF16, folder `sf14_16_per_sample_gams/` | S16 / S17 / S18 (the plot script writes the final numbers) |

## Environment

| Variable | Needed by | Value |
|---|---|---|
| `DATA_ROOT` | all | the deposit root, which holds `2026_final_atlas/` and `202608_kpca_secb_mouse/` (default `/path/to/deposit`) |
| `OUTPUT_ROOT` | all | where outputs are written (default `./output`, relative to the working directory). Use one root for the atlas and Xenium figures and the tables: S8 reads the Fig 2a prep from it, S10 the Table S6 CSVs, and S11 and S16–S18 the Table S8 CSV. `00_figS1_extract_atlas_obs.py` stops if it is not set. |
| `SF10_OUT` | S10 | optional: one folder for the S10 figures and tables, instead of `$OUTPUT_ROOT/final_panels/supplementary/{figures,tables}/` |
| `STUDY_CODE_KEY` | S11–S13, S15, S16–S18 | optional: the private study-code key (default `$OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv`, written by Table S8). It is needed only when the inputs still carry the lab's source identifiers; study codes in the deposited data map to themselves. The key is never written to an output. |
| `FINAL_TMA_CORES` | S12–S15 | optional: the final core table (default `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv`) |
| `RERUN_SFE_DIR` | S14 (step 19), S16–S18 (step 22) | optional: see `load_sfe()` below |
| `SF12_SAMPLES` | S14 (step 19) | optional: comma-separated subset of the 11 SFEs to aggregate |
| `LANG`, `LC_ALL` | every R script | `en_US.UTF-8`, so that `config.yml` is read correctly |

**Script-specific settings.**
- `00_figS1_extract_atlas_obs.py` reads the two concatenated objects from
  `$DATA_ROOT/2026_final_atlas/pre-processing and integration/processed/`.
- `16_figS10_…` builds its input paths from `DATA_ROOT` (`MOUSE`, `H_LOAD`) and the repository root (`REPO`), not
  through `config/`. It uses seed 0, the seed of the Fig 3i export.
- `load_sfe()` (`spatial/00_setup/00_setup.R`) looks for an SFE in `$OUTPUT_ROOT/sfe`, then in `$RERUN_SFE_DIR` if it
  is set, otherwise in `$DATA_ROOT/2026_final_xenium_analysis/output/sfe` (`config.yml` `objects.sfe_dir`, the SFEs
  from before the re-run). The re-run SFEs are written to `$OUTPUT_ROOT/sfe` by the backend.
- `23_figS1-18_composites.py` draws the panel letters with the macOS system font Arial Bold (path on its `FONT` line).
- Python steps ran in E2 and R steps in E1.

## Panel → script

| Panel | What it shows | Drawn by | Numbers computed by |
|---|---|---|---|
| **S1a** | total UMIs per cell by study, cut at 500 | `01_figS1_…` | `00_figS1_extract_atlas_obs.py` (all cells of the pre-QC object); QC in `atlas/01_preprocess_qc/as_run/02_concat_qc_doublets.ipynb` |
| **S1b** | genes per cell by study, cut at 300 | `01_figS1_…` | as S1a |
| **S1c** | Scrublet score by study, cut at 0.25 | `01_figS1_…` | `00_figS1_…` (scores in `atlas_concatenated_filtered.h5ad`); the 0.25 cut is applied in `atlas/01_preprocess_qc/as_run/03b_refilter_umap.py` |
| **S2a** | atlas UMAP (1,980,703 cells) by study | `02_figS2a_…` | `figures/_prep/atlas_final_tables.py` |
| **S2b** | study contribution to each of the 13 level-1 types | `03_figS2b_…` | loader |
| **S2c** | one UMAP per study | `04_figS2c_…` | loader |
| **S3a** | atlas UMAP by treatment status | `05_figS3abc_…` | loader (with the Nath 2021 fix) |
| **S3b** | atlas UMAP by anatomic site | `05_figS3abc_…` | loader |
| **S3c** | atlas UMAP by broad metastatic site | `05_figS3abc_…` | loader |
| **S3d** | atlas UMAP by TP53 status | `06_figS3def_…` | loader |
| **S3e** | atlas UMAP by HRD status | `06_figS3def_…` | loader |
| **S3f** | atlas UMAP by BRCA status | `06_figS3def_…` | loader |
| **S4a** | per-cell predicted class and its probability, by class, with cells per class | `07_figS4a_…` | `celltype_pred` / `celltype_probability`, written by `atlas/01_preprocess_qc/07_finalize.py` from the CellAssign run of `atlas/01_preprocess_qc/04_cellassign.py` |
| **S4b** | canonical markers across the 13 level-1 types | `09_figS4b_…` | `08_figS4b_marker_stats.py` (all 1,980,703 cells) |
| **S4c** | % CopyKAT-aneuploid cells per sample for the epithelial compartment and the three reference lineages (fibroblast, endothelial, T/NK); pooled % | `10_figS4c_…` | the same script (deposited CopyKAT calls); the epithelial states are in the CSV only |
| **S5** | usage of the 10 NMF factors on the epithelial UMAP (575,366 cells) | `11_figS5_…` | deposited NMF (`atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py`). `tables/05_tableS5_epithelial_nmf_programs.py` reads the factor labels (`FACTOR_NAMES`) from this script by path, so keep its name. |
| **S6a** | Factor 3 (SecA) usage and the 7 SecA genes on the epithelial UMAP | `12_figS6_…` | loader and epithelial X |
| **S6b** | Factor 2 (SecB) usage and the 7 SecB genes | `12_figS6_…` | as S6a |
| **S7** | state composition of every CopyKAT subclone with ≥ 20 cells (378 subclones in 251 samples), grouped by sample, with a sample clonality banner (143 monoclonal, 54 clonally driven, 54 mixed) | `13_figS7_…` | the same script. Its verdicts are in `per_sample_verdict.csv`: monoclonal if one clone passes; clonally driven if one clone holds ≥ 70% of the sample's SecA cells and one holds ≥ 70% of its SecB cells; otherwise mixed. |
| **S8** | heatmaps of SecA, Intermediate, SecB and ciliated: PROGENy (14), Hallmark (50), scFEA (165), DoRothEA (294) | `14_figS8_…` | Fig 2a prep (`figures/_prep/fig2a_*`), checked against the script's own rebuild; constant scFEA modules are dropped |
| **S9a** | state composition by treatment status | `15_figS9_…` | loader |
| **S9b** | by anatomic site | `15_figS9_…` | loader |
| **S9c** | by TP53 status | `15_figS9_…` | loader |
| **S9d** | by HRD status | `15_figS9_…` | loader |
| **S9e** | by BRCA status | `15_figS9_…` | loader |
| **S9f** | by Vázquez-García (SPECTRUM) epithelial cluster | `15_figS9_…` | the published GEO cell table GSE180661, joined on the atlas cell id |
| **S9g** | by Geistlinger per-cell consensusOV subtype | `15_figS9_…` | `geistlinger_cell_annotations.csv` (published per-cell labels; see "Missing code" below) |
| **S10a** | human SecA / SecB markers and their mouse orthologues in KPCA.B tumour epithelium (12,429 cells) | `16_figS10_…` | the same script; % expressing in the atlas comes from Table S6 |
| **S10b** | UMAP of the orthologues expressed in ≥ 20% of tumour cells | `16_figS10_…` | the same script (embedding from the mouse hand-off object) |
| **S10c** | Spearman correlation of the four KPCA.B NMF programs with the human factors and polarization; ns = FDR > 0.05 | `16_figS10_…` | `models/kpca_mouse/scripts/32_identify_secb_align.py` (gene-label permutation null, BH over 44 tests, via `05_crossspecies_alignment.py`). `16_…` recomputes the matrix and asserts it matches. |
| **S10d** | usage of each KPCA.B NMF program on the tumour UMAP | `16_figS10_…` | `models/kpca_mouse/scripts/31_fulltumour_nmf.py` |
| **S11** | SingleR cell types (18) in the 8 whole tissues (WT01–WT04 treatment-naive, WT05–WT08 post-chemotherapy; 1,898,772 cells) | `17_figS11-13_… wt` | cell count: `figures/figure4/01_fig4a_composition_table.py`. Treatment status: Table S8. |
| **S12a** | TMA A, 95 patients / 177 cores, duplicate cores side by side | `17_… tmaA` | `spatial/02_qc/06b_final_tma_cores.py` |
| **S12b** | 15 fallopian-tube epithelium cores (34,281 cells) | `17_… fte` | cell count: `figures/figure4/01_fig4a_composition_table.py` |
| **S13** | TMA B, 283 patients / 326 cores (43 two-core pairs) | `17_… tmaB` | `spatial/02_qc/06b_final_tma_cores.py` |
| **S14** | canonical-marker dot plot across the 18 Xenium types. † marks genes excluded from SingleR by probe QC (`genes_exclude_singler.txt`); ‡ marks Xenium-only genes. | `19_figS14_…` (`SF12B_…`) | `18_figS14_annotation_validation_extract.R` |
| **S15a** | polarization score (SecB − SecA UCell) per unit: 15 FTE cores, 8 whole tissues, 378 TMA patients (503 cores); dashed lines at the SecA / SecB thresholds (0.109429 / 0.210501) | `20_figS15_…` | the same script; `SF13_units.csv` |
| **S15b** | all secretory cells pooled per arm | `20_figS15_…` | `SF13_arms.csv`. The HGSC-vs-FTE test (two-sided Mann–Whitney on per-unit median polarization) is in `SF13_tests.csv`. |
| **S16a–d** | per-tissue GAMs against focal polarization: a cell density, b distance to the nearest vascular cell, c matrix remodelling, d neighbourhood hypoxia; pooled curve over 1,235,350 secretory cells | `22_figS16-18_…_plot.py` | `21_figS16-18_per_sample_gams.R`: all cells, tissue-weighted `bam(y ~ s(x, k = 10))`, per-tissue curves, the number of tissues moving in the pooled direction, and the paired Wilcoxon of the top vs bottom 20% of polarization (n = 8) |
| **S17a–n** | the same for epithelium: a proliferation, b Myc, c PI3K, d Wnt, e DNA damage repair, f NF-κB, g glycolysis (UCell); h IL32, i CD55, j LDHA, k CTNNB1, l ICAM1, m MMP7, n TGM2 (expression) | `22_…` | `21_…`. Proliferation, DNA damage repair, NF-κB, glycolysis and LDHA come from the Fig 5a per-cell table (`figures/figure5/01_…`); Myc, PI3K and Wnt are SFE `pathway_*` scores; the other genes are log1p CP10k from the SFE counts. |
| **S18a–h** | the same for macrophages against local epithelial polarization within 50 µm: a glycolysis, b NF-κB activation (UCell); c A2M, d C1QC, e CXCL10, f CXCL11, g CTSL, h ICAM1 | `22_…` | `21_…`; the programs and the x axis come from `figures/figure5/08_…` |

**Numbers in the text from these figures.**

| Text | Value | Where it is computed |
|---|---|---|
| Discussion | 93.5% of patients | `13_figS7_…` → `figures/_prep/text_state_coexistence.py` |
| Results | p = 0.003, HGSC above FTE | `SF13_tests.csv` from `20_figS15_…` |
| S16 and S17 legends | 1,235,350 secretory cells | the rows of `singler_gam_percell.csv` (`figures/figure4/08_fig4f_gam_features.R`); also `n_cells` in the S16–S18 statistics table |

## Not in this folder

**Upstream pipelines.**
- The loaders and prep scripts are in `figures/_prep/`.
- The atlas backend is in `atlas/`, the mouse pipeline in `models/kpca_mouse/scripts/`, and the Xenium backend in
  `spatial/`.
- The Fig 4 and Fig 5 per-cell builders are in `figures/figure4/` and `figures/figure5/`.
- Table S6 and Table S8 are built in `tables/`.

**Manual steps.**
- **Panel files for the composites.** Copy each render (`.png`, and `.svg` for the deposit) into
  `$OUTPUT_ROOT/final_panels/supplementary/figures/`, under the name that `23_figS1-18_composites.py` expects. The S10
  files are written there directly.

  | Render | Copy as |
  |---|---|
  | `SF1A_…`, `SF1B_…`, `SF1C_…` | `SuppFig1a_atlas_qc_total_counts_by_study`, `SuppFig1b_atlas_qc_genes_by_study`, `SuppFig1c_atlas_qc_doublet_score_by_study` |
  | `SF2A_atlas_study_umap` | `SuppFig2a_atlas_umap_by_study` |
  | `SF2B_atlas_integration_study_composition` | `SuppFig2b_atlas_study_contribution_level1` |
  | `SF2C_atlas_integration_per_study_umaps` | `SuppFig2c_atlas_per_study_umaps` |
  | `SF3_treatment_status`, `SF3_anatomic_site`, `SF3_metastatic_site` | `SuppFig3a_atlas_umap_treatment_status`, `SuppFig3b_atlas_umap_anatomic_site`, `SuppFig3c_atlas_umap_metastatic_site` |
  | `SF3_tp53_status`, `SF3_hrd_status`, `SF3_brca_status` | `SuppFig3d_atlas_umap_TP53_status`, `SuppFig3e_atlas_umap_HRD_status`, `SuppFig3f_atlas_umap_BRCA_status` |
  | `SF4A_atlas_annotation_validation` | `SuppFig4a_atlas_scanvi_prediction_probability` |
  | `SF4B_atlas_canonical_markers_dotplot` | `SuppFig4b_atlas_canonical_markers_dotplot` |
  | `SuppFig4c_atlas_cnv_aneuploid_epithelial_vs_reference` | same name |
  | `SF5_atlas_nmf_factor_umaps` | `SuppFig5_atlas_nmf_factor_umaps` |
  | `SF6a_…`, `SF6b_…` | `SuppFig6a_atlas_SecA_signature_umaps`, `SuppFig6b_atlas_SecB_signature_umaps` |
  | `SF7_atlas_cnv_validation` | `SuppFig7_atlas_clone_epitype_composition` |
  | `SF8_atlas_functional_heatmap` | `SuppFig8_atlas_functional_heatmap` |
  | `SF9{a..g}_atlas_epitype_by_<metadata>` | `SuppFig9{a..g}_atlas_epitype_by_<metadata>` |
  | `SuppFig11*`–`SuppFig13*` (in `SF10_SF11_singler/`) | same names |
  | `SF12B_xenium_canonical_markers_dotplot` | `SuppFig14b_xenium_canonical_markers_dotplot` |
  | `SuppFig15_xenium_polarization_distribution` (in `SF13_singler/`) | same name |
  | `SuppFig16*`–`SuppFig18*` (in the step-23 `out_dir`) | same names; or give that folder as `out_dir` |

  The statistics tables that go with the figures are placed in `$OUTPUT_ROOT/final_panels/supplementary/tables/`
  the same way: `SuppFig4c_cnv_aneuploid_by_epithelial_state.csv`, `SF13_tests.csv` (as
  `SuppFig15_fte_vs_hgsc_polarization_tests.csv`) and `SuppFig16-18_per_sample_consistency.csv`. The
  `SuppFig10_kpca_*.csv` tables are written there directly. No script reads these copies.
- **Probe QC.** `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (Figure S14 †) is a copy of
  `$DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/genes_exclude_singler.txt`, written by the probe-QC
  notebook `spatial/02_qc/as_run/explore_probe_qc.Rmd`. No script copies it into `$OUTPUT_ROOT`.

**Missing code.**
- `as_run/geistlinger_cell_annotations.R` writes `output/spectrum_epithelial_crosstab/geistlinger_cell_annotations.csv`,
  the Geistlinger per-cell consensusOV subtypes behind Figure S9g, from the per-cell `colData` of the authors' published
  SingleCellExperiment objects. It ran inline in the atlas project on 2026-06-12 and is recovered from the analysis
  session log; re-run on those objects it writes a file byte-identical to the deposited table (checked 2026-10-22).
- **Download (no script).** The GEO cell table of GSE180661 (Figure S9f), read as
  `output/spectrum_epithelial_crosstab/GSE180661_GEO_cells.tsv.gz`.
- The KPCA.B tumour-epithelium UMAP (`epi_umap_1/2`, used by S10b and S10d) is computed by
  `models/kpca_mouse/scripts/27b_epithelial_umap.py`.
  The coordinates exist only in the deposited mouse hand-off object (see `models/kpca_mouse/`).
- The TMA de-array code behind the core and patient assignments of S12–S15 is covered in `spatial/`.
