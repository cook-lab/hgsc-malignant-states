# Figure 3 — Polarization is environmentally programmed (PDOs) + TCGA subtype

One reproduction script per panel. Config-resolved paths; `SEED`/`CFG$seed` from
config. Atlas NMF label standardized to **Intermediate** (was "Transitioning").

EXTERNAL DEPENDENCY: organoid (PDO) substrate lives outside the deposited
monorepo objects (`<organoids_root>`, default `2026_organoids`; override
`ORGANOIDS_ROOT`). Kept per the migration brief.

| # | Script | Panel | Inputs -> Output |
|---|---|---|---|
| 01 | `01_organoids_models_umap.R` | 3A | `<organoids_root>/output/01_.../seurat_untreated_baseline.rds` -> `organoids_models_umap.{png,svg}` |
| 02 | `02_organoids_secB_ucell_by_model.py` | 3B | `<organoids_root>/output/02_.../organoid_secB_classified_v5.csv` + atlas cache¹ -> `organoids_secB_ucell_by_model.{png,svg}` |
| 03 | `03_organoids_secB_timecourse.py` | 3D | `<organoids_root>/output/08_OPTO98_.../opto98_ucell_scores_aligned.csv` (Growth) + atlas cache¹ -> `organoids_secB_timecourse.{png,svg}` |
| 04 | `04_organoids_secB_perturbations.py` | 3E | same OPTO98 csv (Treatment) + atlas cache¹ -> `organoids_secB_perturbations.{png,svg}` |
| 05 | `05_organoids_g2m_barplot.py` | 3F | `<organoids_root>/output/09_.../09b_extended_per_cell.parquet` + `metadata.csv` -> `organoids_g2m_barplot_secAB.{svg,png}` |
| 06 | `06_organoids_radar_characterization.py` | 3G | `<organoids_root>/output/09_.../09b_{progeny,hallmark,dorothea}_zscored.csv` -> `organoids_radar_characterization_secAB.{svg,png}` |
| 07 | `07_atlas_tcga_subtype_epitype.py` | 3H | `07_deconvolution_survival/20_consensusov/tables/20c_*` -> `atlas_tcga_subtype_epitype.{svg,png}` |

¹ Atlas reference cache (panels 3B/3D/3E dashed reference lines): the deposited
UCell cache at `data_root/2026_final_atlas/output/18_ucell_atlas/`
(`atlas_ucell_scores.csv` + `atlas_secretory_metadata.csv`). Per the cache
contract, **figures read the deposited cache directly from `data_root`** — they
do NOT re-run UCell or read `output_root`. The cache's `celltype_nmf` column
still carries the legacy **"Transitioning epithelium"** label; each script
renames it to **"Intermediate epithelium"** on read so the reference-line filter
(`SecA / Intermediate / SecB epithelium`) matches and the dashed `Atlas Int.`
line renders.

Not migrated: 3C (KRT19 flow-cytometry ridgeline) — wet-lab, plotted externally
(no committed generator anywhere). Outputs land in `figures_dir`.
