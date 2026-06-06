# spatial/08_clinical_survival

Clinical feature construction, Cox survival, macrophage-niche survival, and
Xenium–protein correlation. All scripts source `../../config/config.R` +
`../00_setup/00_setup.R`; seed from `CFG$seed`. Epithelial label standardized to
**Intermediate**. (Original `setwd` lines removed.)

Two chains: clinical/Cox (01-02, ex-`10*`) and macrophage-niche survival
(03-07, ex-`29*`), plus protein correlation (08, ex-`41`).

## Ordered scripts

| # | Script | Inputs → Outputs | Panels |
|---|--------|------------------|--------|
| 01 | `01_clinical_v2.R` | SFEs + clinical_data_clean.csv → `10_clinical_v2/{per_patient_features_v2,cox_univariate_results}.csv` | **Fig 4J, 7A, 7B** |
| 02 | `02_clinical_comprehensive.R` | per_patient_features_v2 + clinical → comprehensive survival tables | Fig 7A/7B (supporting) |
| 03 | `03_compute_niche_metabolic_scores.R` | WT+TMA SFEs (pathway_hypoxia/glycolysis; dbscan frNN 50µm) → `29_macrophage_niche_survival/per_cell_niche_scores.rds` | Fig 6B/6E/6F (upstream) |
| 04 | `04_glmm_presence.R` | per_cell_niche_scores → GLMM presence results | Fig 6B |
| 05 | `05_paired_enrichment.R` | per_cell_niche_scores → paired enrichment tables | Fig 6E/6F |
| 06 | `06_functional_survival.R` | niche scores + clinical → survival results | Fig 6 survival support |
| 07 | `07_figures.R` | niche + survival outputs → `figures/29_macrophage_niche_survival/` | Fig 6B/6E/6F |
| 08 | `08_xenium_protein_correlation.R` | MFI xlsx + sfe_tma_filtered + per_core_proportions → `41_xenium_protein_correlation/per_core_xenium_protein.csv` + correlations | **Fig 7C, 7D** |

## Notes
- `05_paired_enrichment.R` keeps its original local `set.seed(29)` (a distinct
  enrichment-null seed, intentionally **not** converted to `CFG$seed` to preserve
  the exact published null distribution).

## Figures supported
**Fig 4J, 7A, 7B** (clinical/Cox), **Fig 6B/6E/6F** (niche survival), **Fig 7C, 7D** (protein).
