# Figure 7 — Survival; protein & TCGA validation

One reproduction script per panel. R scripts source `config/config.R` (paths/seed)
and `spatial/00_setup/00_setup.R` (`ref_palette`, `theme_lab`). Python scripts import
`config.config`. Inputs span the Xenium pipeline (`DATA_ROOT/2026_final_xenium_analysis`)
and the atlas TCGA deconvolution (`DATA_ROOT/2026_final_atlas`). Figures are written
under `OUTPUT_ROOT/figures/figure7/`.

| Order | Script | Panel(s) | Inputs → Outputs |
|---|---|---|---|
| 01 | `01_xenium_forest_cox.R` | 7A, 7B | xenium 10_clinical_v2 `cox_univariate_results.csv` → `xenium_forest_cox_combined.{pdf,png,svg}` |
| 02 | `02_xenium_protein_correlation_row.py` | 7C, 7D | xenium 41 `per_core_xenium_protein.csv` → `xenium_protein_correlation_row.{svg,png}` |
| 03 | `03_tcga_km_forest.R` | 7E, 7F, 7G | atlas 22 `22d_signature_scores.csv` → `tcga_km_{os,pfs}.{svg,pdf}`, `tcga_forest_stepwise.{svg,pdf}` (+ `_data.csv`) |

Notes:
- 03 computes the stepwise Cox in-script (unadjusted → + epi fraction → + stage + age;
  platinum intentionally omitted). The BayesPrism multivariate covariate set is a
  documented discrepancy.
