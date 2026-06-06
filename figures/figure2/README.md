# Figure 2 — Transcriptional programs, cycle, site/treatment distribution

One reproduction script per panel. Config-resolved paths; `SEED` from config.
Epithelial label standardized to **Intermediate** (was "Transitioning").
SecA/SecB signature markers (volcano) loaded from `shared/signatures.yml`.

| # | Script | Panel | Inputs -> Output |
|---|---|---|---|
| 01 | `01_atlas_radar_characterization.py` | 2A | `04_functional/21_{progeny,hallmark,dorothea,flux}_zscored.csv` -> `atlas_radar_characterization{,_1x4,_<set>}.{svg,png}` |
| 02 | `02_atlas_g2m_barplot.py` | 2B | `04_functional/nmf_characterization/cell_cycle_phase_proportions.csv` -> `atlas_g2m_barplot_secAB.{svg,png}` |
| 03 | `03_atlas_volcano_secA_secB.py` | 2C | `fig_data_dir/panel_i_deg_results.csv` + `shared/signatures.yml` -> `atlas_volcano_secA_secB.{svg,png}` |
| 04 | `04_atlas_site_secb_shift.py` | 2D,2E | `fig_data_dir/{panel_b_site_proportions,panel_c_paired_site}.csv` -> `atlas_site_secb_shift.{svg,png}` |
| 05 | `05_atlas_treatment_secb_shift.py` | 2F,2G | `fig_data_fig1/data_fig1i_treatment_proportions.csv` + `fig_data_dir/meta.parquet` -> `atlas_treatment_secb_shift.{svg,png}` |

Dependency note: 2B reads `cell_cycle_phase_proportions.csv` produced by the
Fig-2B NMF-characterization prep helper (`atlas_nmf_characterization_data` in the
source tree); migrate that helper into `figures/_prep/` if regenerating the cache.
