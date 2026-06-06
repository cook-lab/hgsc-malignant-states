# spatial/05_gradients_gams

GAM-based spatial gradient modeling along the SecA→SecB polarization axis, plus
vascular-proximity and glandular-architecture gradient metrics. `01/04/05` source
`../00_setup/00_setup.R`; `02/03` are self-contained and source only `../../config/config.R`
(they route bare `output/`/`sfe` paths through `cfg_path`/`cfg_obj`). Seed from `CFG$seed`.
Epithelial label standardized to **Intermediate**.

## Ordered scripts

| # | Script | Inputs → Outputs | Panels |
|---|--------|------------------|--------|
| 01 | `01_niche_succession_gams.R` | SFEs (cell_label, pathway_*, polarization_UCell) → `16b_niche_succession_gams/neighborhood_features.rds` + GAM tables/figures | **Fig 4G, 6G, SF12** |
| 02 | `02_gene_polarization_gams.R` | `16b_niche_succession_gams/neighborhood_features.rds` + SFEs → `19d_gene_polarization_gams/{epithelial_expression_polarization,gene_gam_results}.rds`, summaries | **Fig 5A, 5G, SF13** |
| 03 | `03_gene_polarization_gams_all_celltypes.R` | SFEs (RANN nearest-epi polarization) → `19e_gene_gams_all_celltypes/tme_expression_polarization.rds`, per-celltype GAMs | **Fig 6H, SF14** |
| 04 | `04_vascular_proximity.R` | SFEs (coords) → `22_vascular_proximity/{vascular_distance_summary,vascular_distance_all_cells}.csv` | **Fig 4H, SF12** |
| 05 | `05_glandular_architecture.R` | sfe_tma_filtered + WT SFEs → `28_glandular_architecture/per_cell_architecture_wt.rds`, per-sample/per-patient medians | **SF12** |

## Note
Per the migration conventions the `_v2` output suffix is dropped: `02` reads the
`neighborhood_features.rds` written by `01` from `16b_niche_succession_gams/` (the
base script supersedes v1/v2). The original read from `16b_niche_succession_gams_v2/`.
