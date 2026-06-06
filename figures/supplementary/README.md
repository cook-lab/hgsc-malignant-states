# Supplementary figures (SF1–SF14)

Cleaned reproduction scripts for the manuscript supplementary figures. Each
script imports the central config (Python: `from config.config import ...`;
R: `source("../../config/config.R")`), sets the global seed from config at any
stochastic step, uses the epithelial label **Intermediate** (formerly
"Transitioning"), and loads SecA/SecB signatures from `shared/signatures.yml`
where a marker list is needed.

All figure outputs are written to
`output_root/figures/supplementary/` and the shared parquet caches to
`output_root/figures/data/` (and `.../figures/data_fig1/` for the full-atlas meta
cache produced by the Fig 1 extraction step).

## Ordered scripts

| # | Script | Panel(s) | Inputs → Outputs |
|---|--------|----------|------------------|
| 00 | `00_extract_atlas_obs.py` | (cache for SF1) | atlas_scanvi h5ad + pre-QC concat → `figures/data/atlas_obs_{pre,post}filter.parquet` |
| 00b | `00b_extract_integration_umaps.py` | (cache for SF2B/C) | integration h5ads + atlas_scanvi → `figures/data/*_umap.parquet` |
| 01 | `01_SF1_atlas_qc_metrics_by_study.py` | SF1A-C | obs parquets → `SF1_atlas_qc_metrics_by_study.{svg,png}` |
| 02 | `02_SF2A_atlas_study_umap.py` | SF2A | `figures/data_fig1/meta.parquet` → `SF2A_atlas_study_umap.{svg,png}` |
| 03 | `03_SF2B_atlas_integration_study_composition.py` | SF2B | `figures/data/atlas_final_umap.parquet` → `SF2B_...{svg,png}` |
| 04 | `04_SF2C_atlas_integration_per_study_umaps.py` | SF2C | `figures/data/atlas_final_umap.parquet` → `SF2C_...{svg,png}` |
| 05 | `05_SF3A_atlas_clinical_umaps.py` | SF3A-C | `atlas_final` h5ad → `SF3_{anatomic_site,metastatic_site,treatment_celltype}.{svg,png}` |
| 06 | `06_SF3D_atlas_genomic_umaps.py` | SF3D-F | `atlas_final` h5ad → `SF3_{tp53,hrd,brca}_status.{svg,png}` |
| 07 | `07_SF4A_atlas_annotation_validation.py` | SF4A | `atlas_final` h5ad → `SF4A_atlas_annotation_validation.{svg,png}` |
| 08 | `08_SF4B_atlas_canonical_markers_dotplot.py` | SF4B | `09f_dotplot_canonical_markers/{dotplot_stats,lilra4_stats}.csv` → `SF4B_...{svg,png}` |
| 09 | `09_SF4C_atlas_cnv_aneuploid_by_celltype.py` | SF4C | `_data_cnv_aneuploid_by_celltype.csv` (committed) → `SF4C_...{svg,png}` |
| 10 | `10_SF5_atlas_nmf_factor_umaps.py` | SF5 | `fig_secretory_polarization/.../meta.parquet` + `11d_nmf_usage.csv` → `SF5_...{svg,png}` |
| 11 | `11_SF6_atlas_secAB_expression_umaps.py` | SF6A-B | meta.parquet + `atlas_epithelial` h5ad + `11d_nmf_usage.csv` + `shared/signatures.yml` → `SF6_...{svg,png}` |
| 12 | `12_SF7_atlas_cnv_validation.py` | SF7 | `19_cnv/tables/{within_clone_coexistence,per_sample_verdict,sample_manifest}.csv` → `SF7_...{svg,png}` |
| 13 | `13_SF8_atlas_functional_heatmap.py` | SF8 | `21_epitype_functional_characterization/21_*_zscored.csv` + scFEA M168 info → `SF8_...{svg,png}` |
| 14 | `14_SF9_atlas_epitype_by_metadata.py` | SF9A-E | meta.parquet + `atlas_epithelial` h5ad → `SF9_...{svg,png}` |
| 15 | `15_SF10B_tma_composition_spatial.R` | SF10B | `sfe_tma_filtered` → `SF10B_tma_composition_spatial.{pdf,png,svg}` |
| 16 | `16_SF10C_xenium_roi_and_core_celltype.R` | SF10C | `sfe_tma_filtered` + 06f → `SF10C_xenium_by_core_celltype/core_*.png` (mode `cores`) |
| 17 | `17_SF11_atlas_xenium_wt_polarization_violins.py` | SF11 | `18_ucell_atlas/...fulllabels.csv` + 06f reclassification → `SF11_...{png,svg}` |
| 18 | `18_SF12_SF13_SF14_suppl_per_sample_gams.R` | SF12, SF13, SF14 | 28/22/16b_v2/19d/19e caches + `9b_scoring/pathway_gene_sets_v2.csv` → `SF1{2,3,4}_suppl_gam_{env,epi,mac}.{svg,pdf}` |

## Notes / cross-area dependencies

- **SF10A** is rendered by the shared Fig 4C script
  (`xenium_whole_tissue_snapshot`, mode `all`) and lives in the main figures
  area — it is not duplicated here.
- The R scripts (15, 16, 18) `source()` the migrated spatial setup
  (`spatial/00_setup/00_setup.R` for `load_sfe`, `ref_palette`, `out_dir`,
  `theme_lab`) and, for script 18, `spatial/00_setup/36_helpers.R` for
  `compute_signature_score`. The spatial shard owns those files; the
  `ref_palette` key must use **"Intermediate epithelium"** (renamed there).
- SF6 loads the canonical 7-gene **noBCAM** SecA/SecB lists from
  `shared/signatures.yml`. This substitutes the upstream display set's WT1 for
  FBXO21 in the SecA marker panel (the inline divergent list was removed per
  convention #6).
- SF8 reads the scFEA M168 module-info table to make metabolic-flux module
  labels human-readable (resolves the int-`Module_id` vs `M_x` quirk).
- SF11 reads the canonical `18_ucell_atlas` scoring (NOT `18b_v2`), matching the
  xenium noBCAM signature (Q5).
- The shared per-figure parquet caches (`figures/data/`, `figures/data_fig1/`)
  are produced by the `00*` helpers here and by the Fig 1 data-extraction step
  in the main figures area.
