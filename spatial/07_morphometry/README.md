# spatial/07_morphometry

Cell/nucleus and macrophage morphometry across the SecA→Intermediate→SecB axis.
All scripts source `../../config/config.R` + `../00_setup/00_setup.R`; seed from
`CFG$seed`. Epithelial label standardized to **Intermediate**. (The original
`setwd("/Volumes/CookLab/Sarah/...")` lines were removed.)

Two chains share this stage: cell/nucleus morphometrics (01-05, ex-`33*`) and
macrophage morphometrics (06-10, ex-`34*`).

## Ordered scripts

| # | Script | Inputs → Outputs | Panels |
|---|--------|------------------|--------|
| 01 | `01_compute_cell_nucleus_morphometrics.R` | SFEs (geometries; 06f) → `33_morphometrics/` per-cell + eligibility | Fig 5D/5E (upstream) |
| 02 | `02_pairwise_wt_morphometrics.R` | per-cell WT → `33_morphometrics/per_sample_summary_wt.csv` + stats | Fig 5D/5E (WT) |
| 03 | `03_pairwise_tma_morphometrics.R` | per-cell TMA → `33_morphometrics/per_patient_summary_tma.csv` + stats | Fig 5D/5E (TMA) |
| 04 | `04_cross_cohort_summary.R` | WT+TMA pairwise → `33_morphometrics/cross_cohort_summary.csv` | Fig 5D/5E (summary) |
| 05 | `05_morphometric_figures.R` | summaries → `figures/33_morphometrics/` | **Fig 5D, 5E** |
| 06 | `06_compute_macrophage_morphometrics_and_niche.R` | SFEs (RANN 50µm niche class) → `34_macrophage_morphometrics/per_cell_macrophage_morphometrics.rds` | Fig 6I/6J (upstream) |
| 07 | `07_pairwise_macrophage_wt.R` | per-cell rds → WT pairwise stats | Fig 6I/6J (WT) |
| 08 | `08_pairwise_macrophage_tma.R` | per-cell rds → TMA pairwise stats | Fig 6I/6J (TMA) |
| 09 | `09_cross_cohort_summary.R` | WT+TMA → `34_macrophage_morphometrics/cross_cohort_summary.csv` | Fig 6I/6J (summary) |
| 10 | `10_figures.R` | summaries → `figures/34_macrophage_morphometrics/` | **Fig 6I, 6J** |

## Figures supported
**Fig 5D, 5E** (nuclear area, N:C ratio) and **Fig 6I, 6J** (macrophage cell area, circularity).
