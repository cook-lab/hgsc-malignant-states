# spatial/00_setup — shared spatial setup

Single sourced setup for the entire Xenium spatial backend.

## Script
- `00_setup.R` — sources the central config (`config/config.R`), loads libraries,
  defines the canonical palettes (`ref_palette` with the standardized
  **"Intermediate epithelium"** label), the `theme_lab()` ggplot theme, the
  SecA/SecB signatures (`SECA_GENES`/`SECB_GENES`, loaded from
  `shared/signatures.yml` — the noBCAM 7-gene set), and the `load_sfe()` /
  `save_sfe()` helpers. Derives paths from config and resolves the cohort PIN.

## Inputs -> outputs
- IN: `config/config.{R,yml}`, `shared/signatures.yml`,
  `<data_root>/2026_final_xenium_analysis/metadata/samples.csv`
- OUT: none (defines objects/helpers in the caller; creates `output_root`
  subdirectories)

## Cohort PIN
Exposes `whole_tissue_samples` (the published 8 from `CFG$cohort$whole_tissue`),
`sfe_names_wt`, and `sfe_names_all` (= `sfe_tma` + the 8 whole tissues). The FTE
whole-tissue samples (`CFG$cohort$fte_exclude_wt` = FT1-1, EAOC-1-FTE) are NOT
in these vectors and are excluded from every downstream whole-tissue step.

## Usage
Every spatial stage script begins with `source("spatial/00_setup/00_setup.R")`
(run from the repo root with Rscript). The setup file locates `config/config.R`
two directories up via its own file path.

## Figures supported
None directly — supports the entire Xenium backend (Fig 4–7, SF10–SF14).
