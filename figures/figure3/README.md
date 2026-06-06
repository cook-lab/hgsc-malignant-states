# Figure 3 — Polarization is environmentally programmed (PDOs) + TCGA subtype

One reproduction script per panel. Config-resolved paths; `SEED`/`CFG$seed` from
config. Atlas NMF label standardized to **Intermediate** (was "Transitioning").

EXTERNAL DEPENDENCY: organoid (PDO) substrate lives outside the deposited
monorepo objects (`<organoids_root>`, default `2026_organoids`; override
`ORGANOIDS_ROOT`). Kept per the migration brief; see `docs/REPRODUCIBILITY.md`.

| # | Script | Panel | Inputs -> Output |
|---|---|---|---|
| 01 | `01_organoids_models_umap.R` | 3A | `<organoids_root>/output/01_.../seurat_untreated_baseline.rds` -> `organoids_models_umap.{png,svg}` |
| 02 | `02_organoids_secB_ucell_by_model.py` | 3B | `<organoids_root>/output/02_.../organoid_secB_classified_v5.csv` + `18_ucell_atlas/*` -> `organoids_secB_ucell_by_model.{png,svg}` |
| 03 | `03_organoids_secB_timecourse.py` | 3D | `<organoids_root>/output/08_OPTO98_.../opto98_ucell_scores_aligned.csv` (Growth) + atlas refs -> `organoids_secB_timecourse.{png,svg}` |
| 04 | `04_organoids_secB_perturbations.py` | 3E | same OPTO98 csv (Treatment) + atlas refs -> `organoids_secB_perturbations.{png,svg}` |
| 05 | `05_organoids_g2m_barplot.py` | 3F | `<organoids_root>/output/09_.../09b_extended_per_cell.parquet` + `metadata.csv` -> `organoids_g2m_barplot_secAB.{svg,png}` |
| 06 | `06_organoids_radar_characterization.py` | 3G | `<organoids_root>/output/09_.../09b_{progeny,hallmark,dorothea}_zscored.csv` -> `organoids_radar_characterization_secAB.{svg,png}` |
| 07 | `07_atlas_tcga_subtype_epitype.py` | 3H | `07_deconvolution_survival/20_consensusov/tables/20c_*` -> `atlas_tcga_subtype_epitype.{svg,png}` |

Not migrated: 3C (KRT19 flow-cytometry ridgeline) — wet-lab, plotted externally
(no committed generator anywhere). Outputs land in `figures_dir`.
