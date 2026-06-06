# spatial/06_spatial_stats

Spatial autocorrelation statistics for secretory polarization. Sources
`../../config/config.R` + `../00_setup/00_setup.R`; seed from `CFG$seed`.
Epithelial label standardized to **Intermediate**.

## Ordered scripts

| # | Script | Inputs → Outputs | Panels |
|---|--------|------------------|--------|
| 01 | `01_secretory_spatial_autocorrelation.R` | SFEs (cell_label, SecA_UCell/SecB_UCell/polarization_UCell, coords; 06f override) → `44_spatial_autocorrelation/{interpretation_summary,tma_patient_level_lee,bilisa_vs_label_crosstab}.csv` | **Fig 4E, 4F** |

Univariate Moran's I (LISA), bivariate Lee's L, and BiLISA HH/HL/LH/LL regime
classification of SecA vs SecB UCell scores (segregation vs co-occurrence).
