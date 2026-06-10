# Figure 6 — Glycolytic macrophages + SecB = lymphocyte-excluded niche

One reproduction script per panel. R scripts source `config/config.R` (paths/seed)
and `spatial/00_setup/00_setup.R` (`load_sfe`, `ref_palette`, `theme_lab`,
`log_session`). Python scripts import `config.config`. Pipeline inputs live under
`DATA_ROOT/2026_final_xenium_analysis`; figures are written under
`OUTPUT_ROOT/figures/figure6/`.

Whole-tissue cohort is taken from `CFG$cohort$whole_tissue` (the published 8 WT
samples; FTE whole-tissue samples excluded).
Epithelial labels standardized to SecA / Intermediate / SecB / Ciliated.

| Order | Script | Panel(s) | Inputs → Outputs |
|---|---|---|---|
| 01 | `01_tma_hypoxia_gradient_cores.R` | 6A, 6C | `sfe_tma_filtered` + 06f → `tma_hypoxia_gradient_cores/core{ID}_{pt}_{hypoxia,celltype}.png`, `gradient_metrics.csv` |
| 02 | `02_niche_immune_composition_gam.R` | 6B | 29 `per_cell_niche_scores.rds` → `niche_immune_composition_gam.{pdf,png,svg}` |
| 03 | `03_tma_gradient_cores_immune_highlight_svgs.R` | 6D | `sfe_tma_filtered` + 06f + `gradient_metrics.csv` (from 01) → `tma_hypoxia_gradient_cores/core{ID}_{pt}_{macrophage,tcell,plasmacell,nkcell,bcell,...}.{svg,png}` |
| 04 | `04_tma_gradient_cores_secb_highlight.R` | 6D | `sfe_tma_filtered` + 06f + `gradient_metrics.csv` (from 01) → `tma_hypoxia_gradient_cores/core{ID}_{pt}_secb.{svg,png}` |
| 05 | `05_fig_exhaustion_ucell.R` | 6E, 6F | 29 `per_cell_niche_scores.rds` + WT/TMA SFEs → `fig_{tcell,nkcell}_exhaustion_wt.{svg,png,pdf}` (+ TMA exclusion hist) |
| 06 | `06_gam_neighborhood_immunomod_genes_polarization.R` | 6G | 16b_v2 `neighborhood_features.rds` + WT SFEs (50µm frNN) → `gam_neighborhood_immunomod_genes_polarization.{pdf,png,svg}` |
| 07 | `07_gam_macrophage_focal_polarization.R` | 6H | 19e `tme_expression_polarization.rds` → `gam_macrophage_focal_polarization.{pdf,png,svg}` |
| 08 | `08_fig_macrophage_apoptosis_prolif_morphology.R` | 6I, 6J | 29 functional caches + 34 `per_cell_macrophage_morphometrics.rds` → `fig_macrophage_{apoptosis,proliferation,cell_area,circularity}_wt.{pdf,png,svg}` |
| 09 | `09_xenium_roi_whole_tissue_investigation.R` | 6K, 6L | `sfe_<sample>` (06f baked in), logcounts → `roi_whole_tissue_investigation/{sample}_{label}_{gene}.{png,svg}` |

Dependencies / quirks:
- 03 and 04 read `gradient_metrics.csv` produced by 01 — run 01 first.
- **09 requires a `<sample>` CLI argument** (and optional `xmin,xmax,ymin,ymax`,
  `<label>`, gene-subset). It renders nothing without a sample. For Fig 6K/6L use
  the SP24_24824 ROI calls (cell-type context + CTSL/MMP7/ICAM1); see script header.
- 6I/6J p-values differ between header text and the published PDF — code migrated
  faithfully; verify reproduced p (Fig 6J reproduces 6/8 vs published 8/8).
