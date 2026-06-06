# spatial/04_neighborhood

Xenium spatial neighborhood / niche backend. All scripts source `../../config/config.R`
(paths, `CFG$seed`) and `../00_setup/00_setup.R` (`load_sfe`/`save_sfe`, palettes, `nb_names`).
Epithelial label standardized to **Intermediate** (was "Transitioning").

## Ordered scripts

| # | Script | Inputs → Outputs | Role / panels |
|---|--------|------------------|---------------|
| 01 | `01_colocalization.R` | SFEs (cell_label) → `output/08_colocalization/` tables + figures | Cell-type colocalization (Fig 4 context) |
| 02 | `02_neighborhood_k10_production.R` | `09_neighborhood/neighborhood_feature_matrix.rds` + SFEs → `neighborhood_assignments_k10.csv`, centers/composition/mapping; writes `neighborhood`/`neighborhood_name` to SFEs | **Canonical k=10 assignment** (feeds Fig 4G/6G) |
| 03 | `03_neighborhood_finalize_names.R` | SFEs (nb_1..nb_10) → SFEs with refreshed `neighborhood_name` | Name-only writeback |
| 04 | `04_pathway_scoring.R` | SFEs → SFE `pathway_*` cols + `9b_scoring/pathway_gene_sets_v2.csv` | **Supp Data 6**; feeds GAMs |
| 05 | `05_macrophage_niche_analysis.R` | SFEs (cell_label + neighborhood) → `13_macrophage_niche/` | Niche-conditioned macrophage DEG/correlation |
| 06 | `06_tcell_niche_analysis.R` | SFEs → `13_macrophage_niche/` | T-cell niche characterization |
| 07 | `07_bcell_niche_analysis.R` | SFEs → `13_macrophage_niche/` | B-cell niche characterization |
| 08 | `08_nkcell_niche_analysis.R` | SFEs → `13_macrophage_niche/` | NK-cell niche characterization |
| 09 | `09_epithelial_neighborhood_characterization.R` | SFEs → epithelial neighborhood DEG/volcanoes | Secretory succession framing |
| 10 | `10_comprehensive_trajectory.R` | SFEs → trajectory DEG + figures | SecA→Intermediate→SecB trajectory |
| 11 | `11_metabolic_niche_analysis.R` | SFEs (+pathway_*) → metabolic niche tables | Fig 6 metabolic-niche narrative |

## Figures supported
Backend for **Fig 4** (composition/niche, 4G) and **Fig 6** (lymphocyte-excluded niche, 6B/6G); **Supp Data 6** (04).
Neighborhood assignments (02) and pathway scores (04) are upstream caches for `05_gradients_gams/`.
