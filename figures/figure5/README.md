# Figure 5 — Hypoxia/HIF/NF-κB program; adhesion remodelling

One reproduction script per panel. R scripts source `config/config.R` (paths/seed)
and `spatial/00_setup/00_setup.R` (`load_sfe`, `ref_palette`, `theme_lab`). Python
scripts import `config.config`. Pipeline inputs live under `DATA_ROOT/2026_final_*`;
figures are written under `OUTPUT_ROOT/figures/figure5/` (`cfg_path("figures_dir", ...)`).

Epithelial polarization labels are standardized to **SecA / Intermediate / SecB /
Ciliated** (was "Transitioning"). SecA/SecB 7-gene signatures live in
`shared/signatures.yml`; none of these generators inline that signature.

| Order | Script | Panel(s) | Inputs → Outputs |
|---|---|---|---|
| 01 | `01_gam_epithelial_pathways_polarization.R` | 5A | 19d `epithelial_expression_polarization.rds` + 9b pathway sets → `gam_epithelial_pathways_polarization.{pdf,png,svg}` |
| 02 | `02_roi_OTB_2384_roi06_svgs.R` | 5B, 5C | `sfe_OTB_2384` + 06f override + 9b sets → `roi_exploration_OTB_2384_secb/OTB_2384_roi06_*.{svg,png}` |
| 03 | `03_morphometrics_paired.py` | 5D, 5E | 33_morphometrics `per_sample_summary_wt.csv` / `per_patient_summary_tma.csv` (+ counts) → `xenium_{nuc_area,nuc_perimeter,nc_ratio}_paired.{png,svg}` |
| 04 | `04_atlas_seca_secb_autocrine_shift.py` | 5F | atlas 17 `17b_liana_global.csv` → `atlas_seca_secb_autocrine_shift.{png,svg}` |
| 05 | `05_gam_epithelial_genes_polarization.R` | 5G | 19d `epithelial_expression_polarization.rds` → `gam_epithelial_genes_polarization.{pdf,png,svg}` |
| 06 | `06_roi_SP24_24824_zoom_genes.R` | 5H, 5I | `sfe_SP24_24824` + 06f override → `ROI_figure_5/SP24_24824_roi_celltype_full.{svg,png}` (5H) + `SP24_24824_roi_<gene>.{svg,png}` (5I, full ROI: CTNNB1/ITGB5/MMP7/ICAM1) |

Notes:
- 5F: the generator reproduces the panel exactly, but thresholds on expression
  magnitude (`lrscore > 0.5`), not significance.
- 06 renders both panels over the published full ROI (x=[7800,9000], y=[-7100,-5900]):
  the cell-type map (5H, `SP24_24824_roi_celltype_full.*`) and the 4-gene expression
  maps (5I, CTNNB1/ITGB5/MMP7/ICAM1). The published 5H/5I overlay two ROI rectangles
  in Illustrator (an assembly-time annotation, not rendered here).
