# Figure 1 — Integrated atlas; secretory polarization axis

One reproduction script per panel. All read config-resolved caches (no hardcoded
paths); seeds come from `config` (`SEED`). Epithelial label standardized to
**Intermediate** (was "Transitioning").

| # | Script | Panel | Inputs -> Output |
|---|---|---|---|
| 01 | `01_atlas_level1_umap.py` | 1A | `fig_data_fig1/meta.parquet` -> `atlas_level1_umap.{svg,png}` |
| 02 | `02_atlas_study_contribution.py` | 1B,1C,1D,1E,1I | `fig_data_dir/meta.parquet` + `fig_data_fig1/{meta.parquet, panel_b_cells_per_study.csv, panel_b_patients_per_study.csv, panel_g_composition_by_study.csv}` -> `atlas_study_contribution_with_epi_composition.{svg,png}` |
| 03 | `03_atlas_epithelial_secretory_ciliated_umap.py` | 1F | `fig_data_dir/meta.parquet` -> `atlas_epithelial_secretory_ciliated_umap.{svg,png}` |
| 04 | `04_atlas_SecA_nmf_factor_umap.py` | 1G | `fig_data_dir/meta.parquet` + `03_epithelial_nmf/11d_nmf_usage.csv` (Factor_3) -> `atlas_SecA_nmf_factor_umap.{svg,png}` |
| 05 | `05_atlas_SecB_nmf_factor_umap.py` | 1H | same + Factor_2 -> `atlas_SecB_nmf_factor_umap.{svg,png}` |
| 06 | `06_atlas_cnv_alluvial.py` | 1J | `05_cnv/tables/within_clone_coexistence.csv` -> `atlas_cnv_alluvial.{svg,png}` |

Upstream: `figures/_prep/fig_secretory_polarization_00_prepare_data.py` builds the
shared `meta.parquet` (schema_nmf). Outputs land in `figures_dir`.
