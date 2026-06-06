# figures/_prep — figure-data extraction helpers

Run-once data-prep helpers that produce the lightweight caches the per-figure
scripts read (so figures rebuild without re-loading multi-GB h5ads). Paths come
from `config/` (no hardcoded `/Volumes/...`). Set the global seed from config in
any stochastic step.

| Script | Inputs | Outputs | Feeds |
|---|---|---|---|
| `00_extract_atlas_obs.py` | `obj("atlas_scanvi")`; raw concat h5ad | `output_root/_prep_caches/atlas_obs_{prefilter,postfilter}.parquet` | SF1A-C |
| `00b_extract_integration_umaps.py` | integration h5ads (harmony/scvi/scanvi) + `obj("atlas_scanvi")` | `output_root/_prep_caches/integration_*_umap.parquet`, `atlas_final_umap.parquet` | SF2B/C |
| `fig_secretory_polarization_00_prepare_data.py` | `obj("atlas_epithelial")`; 11d usage; 04_functional score parquets | `fig_data_dir/meta.parquet` (schema_nmf 4-class) + `panel_{b,c,d,e,f,g,h,i}_*` | Fig 1F/G/H/I, Fig 2C/D/E/G, SF5/6/9 |

Notes
- `meta.parquet` is the shared `schema_nmf` cache. The epithelial polarization
  label is standardized to **Intermediate** (was "Transitioning"); downstream
  figure scripts select on "Intermediate".
- The whole-atlas `data_fig1/` extraction set (`fig_data_fig1/meta.parquet`,
  `panel_b_cells_per_study.csv`, `panel_b_patients_per_study.csv`,
  `panel_g_composition_by_study.csv`, `data_fig1i_treatment_proportions.csv`) is
  an upstream cache consumed by Fig 1A/B-E/I and Fig 2F; route it under
  `output_root/fig_data_fig1/` (see config `fig_data_dir`).
