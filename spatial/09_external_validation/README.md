# spatial/09_external_validation

External validation of the SecA/SecB polarization signature in TCGA-OV
(cross-checks the atlas deconvolution-survival chain, atlas step 22d).
All scripts source `../../config/config.R` + `../00_setup/00_setup.R`; seed from
`CFG$seed`. **SecA/SecB gene lists are loaded from `../../shared/signatures.yml`**
(the noBCAM 7-gene set) rather than inlined.

## Ordered scripts

| # | Script | Inputs → Outputs | Role |
|---|--------|------------------|------|
| 01 | `01_tcga_external_validation.R` | TCGA-OV TPM + external clinical + CIBERSORTx; signatures.yml → `40_tcga_validation/` Cox tables + run log | Pre-registered validation (Cox × adjustment models) |
| 02 | `02_robustness.R` | `40_tcga_validation/` scores; signatures.yml → robustness tables (LOGO, GSVA/ssGSEA sensitivity, bootstrap) | Robustness |
| 03 | `03_figures_and_report.R` | `40_tcga_validation/` tables → `40_tcga_validation/figures/` + report | Figures + report |

## Notes
- `02_robustness.R` keeps its original local `set.seed(20260508)` for the bootstrap
  (distinct local seed; not converted to `CFG$seed`).

## Figures supported
Cross-validates **Fig 7E/7F/7G** (TCGA KM + forest, generated in the atlas tree from 22d).
