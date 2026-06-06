# Figure 4 — Targeted spatial transcriptomics; spatial gradient

One reproduction script per panel. R figure scripts source the shared spatial
setup (`spatial/00_setup/00_setup.R`), which loads the central config and
provides `load_sfe`, `ref_palette` (with **Intermediate**), `theme_lab`, and the
seed. Python scripts import `config.config`. Epithelial label standardized to
**Intermediate** (was "Transitioning"). SFE/output paths are config-resolved;
06f overrides read from `output_root/06f_reclassification_polarization/`.

| # | Script | Panel | Inputs -> Output |
|---|---|---|---|
| 01 | `01_xenium_composition_by_tissue_dendro.py` | 4A,4B | `figures/cell_type_counts_by_sample.csv` + `06f_.../reclassified_xenium_scores.csv` + core/barcode maps -> `xenium_composition_by_tissue_dendro.{png,svg}` |
| 02 | `02_xenium_whole_tissue_snapshot.R` | 4C (also SF10A) | `sfe_<sample>` + 06f -> `xenium_whole_tissue_snapshot_<sample>.{png,svg}` |
| 03 | `03_fig2d_roi_celltype.R` | 4D | `sfe_OTB_2384` (ROI_C) -> `fig2d_ROI_C_wide_{celltype,macrophage}.{svg,png}` |
| 04 | `04_lee_bivariate_segregation.py` | 4E | `44_spatial_autocorrelation/{interpretation_summary,tma_patient_level_lee}.csv` -> `lee_bivariate_segregation.{png,svg}` |
| 05 | `05_bilisa_regime_by_label.py` | 4F | `44_spatial_autocorrelation/bilisa_vs_label_crosstab.csv` -> `bilisa_regime_by_label.{png,svg}` |
| 06 | `06_gam_microenvironment_polarization.R` | 4G | `16b_niche_succession_gams_v2/neighborhood_features.rds` + SFEs -> `gam_microenvironment_polarization.{pdf,png,svg}` |
| 07 | `07_xenium_vascular_distance_paired.py` | 4H | `22_vascular_proximity/vascular_distance_summary.csv` -> `xenium_vascular_distance_paired.{png,svg}` |
| 08 | `08_xenium_roi_OTB_2384_vascular_distance.R` | 4I | `sfe_OTB_2384` + 06f (ROI_C) -> `xenium_roi_OTB_2384_vascular_distance.{png,svg}` |
| 09 | `09_xenium_epi_density_paired.py` | 4J | `10_clinical_v2/per_patient_features_v2.csv` -> `xenium_epi_density_paired.{png,svg}` |
| 10 | `10_xenium_roi_OTB_2384_pathway_panels.R` | 4K | `sfe_OTB_2384` + 06f (ROI_C; UCell pathway cols) -> `xenium_roi_OTB_2384_{matrix_remodeling,hypoxia,oxidative_stress}.{png,svg}` |

Flags (LINEAGE.md):
- `03_fig2d_roi_celltype.R` colours by the in-object `cell_label` (pre-06f), not
  the 06f reclassification used by the rest of Fig 4 — verify against the
  published 4D. The original ran from the xenium-root cwd; here it is rebased on
  the shared spatial setup so SFE paths resolve from config.
- `02_xenium_whole_tissue_snapshot.R` also drives SF10A via the `all` CLI arg.

Outputs land in `figures_dir`.
