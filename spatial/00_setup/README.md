# `spatial/00_setup/`: shared R setup and helpers

Three files that other scripts `source()`; none of them is run on its own. `00_setup.R` is sourced at the top of the
spatial R stage scripts and of most Xenium figure R scripts. The other two give figure and table scripts the final TMA
core set and the SingleR cell labels.

## Files

| Step | Script | Reads | Writes / defines | Used by |
|---|---|---|---|---|
| helper | `00_setup.R` | `config/config.R` + `config.yml`; `shared/signatures.yml` | creates the output folders under `$OUTPUT_ROOT` (`figures/`, `sfe/`, `03_04_qc/`, `05_probe_qc/`, `06d_annotation_noBCAM/`, `06f_reclassification_polarization/`, `06g_clean_split/`, `07_core_qc/`); defines `load_sfe()`, `save_sfe()`, `log_session()`, `sfe_dir`, `out_dir`, `data_dir`, `sfe_names_wt` / `sfe_names_all`, `SECA_GENES` / `SECB_GENES`, `ref_palette`, `theme_lab()`, `svg_inline_fill()`, `spatial_point_layer()` | the R stage scripts in `spatial/` (except `spatial/04_pathway_scoring/04b_pathway_gene_sets.R`) and `spatial/05_count_export/as_run/02_deg_precompute.R`; the R scripts of `figures/figure4/`, `figures/figure5/`, `figures/figure7/` (except `01`) and `figures/supplementary/` (`17`, `18`, `21`). `tables/06_tableS6_secA_secB_signatures.py` checks that it reads `shared/signatures.yml` |
| helper | `final_tma_cores.R` | `$FINAL_TMA_CORES` (default `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv`) and `final_tma_patients.csv` next to it, written by `spatial/02_qc/06b_final_tma_cores.py` | `final_core_ids(array)`, `final_tma_cores()`, `final_tma_patients()` | `figures/figure4/04_fig4d_bivariate_moran.R`; `figures/figure5/03_fig5b_nuclear_area.R`, `figures/figure5/10_fig5j_macrophage_area.R`; `figures/figure7/01_fig7abcd_tma_cox_km_tables.R`, `figures/figure7/03_fig7ef_protein_correlation.R`; `figures/supplementary/17_figS11-13_xenium_celltype_maps.R`; `tables/11_tableS11_xenium_tma_cox.py` (embedded R) |
| helper | `singler_labels.R` | `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` (from `spatial/03_annotation_polarization/06_extract_singler.R`) | `singler_labels(sfe, name)` → one row per SFE cell: `barcode`, `lab` (= `cell_label_singler`), `polarization`; alias `presplit_labels()` | `figures/figure4/04_fig4d_bivariate_moran.R`, `figures/figure4/06_fig4e_bilisa.R`; `figures/figure5/01_fig5a_epithelial_programs_build.R`, `03_fig5b_nuclear_area.R`, `05_fig5cdegh_roi_export.R`, `06_fig5f_immune_niche_build.R`, `08_fig5i_macrophage_programs_build.R`, `10_fig5j_macrophage_area.R` |

The Python counterpart of `final_tma_cores.R` is `figures/_prep/final_tma_cores.py` (same file, same rules).

## Notes

- **Sample lists.** `sfe_names_wt` is the 8 whole tissues of `config.yml` (`cohort.whole_tissue`, WT01–WT08).
  `sfe_names_all` is `sfe_tma` (TMA A), the arrays of `cohort.tma_new` (`sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H`) and the
  8 whole tissues. The backend steps loop over `sfe_names_all`.
- **Where SFEs are read and written.** `load_sfe(name)` looks in `$OUTPUT_ROOT/sfe`, then in `sfe_dir`. `sfe_dir` is
  `RERUN_SFE_DIR` if set, otherwise `$DATA_ROOT/2026_final_xenium_analysis/output/sfe` (config `objects.sfe_dir`),
  which on the lab's data holds SFEs from before the final re-run. The final run set `RERUN_SFE_DIR=$OUTPUT_ROOT/sfe`.
  `save_sfe()` writes only under `$OUTPUT_ROOT/sfe`. `data_dir` is `$DATA_ROOT/2026_final_xenium_analysis/data`.
- **Final core set.** `final_core_ids("sfe_tma" | "COEUR_TMA_G" | "COEUR_TMA_H")` returns the raw core ids of that
  array that are `final == TRUE`.
- **SingleR labels.** `singler_labels()` renames "Transitioning epithelium" to "Intermediate epithelium" and stops if
  the sample's table is missing. Its header explains why `colData(sfe)$cell_label` is not used: that column holds the
  RCTD/SPLIT labels.
- `config/config.R` needs `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8` to read `config.yml`.
