# figures/_prep — figure-data extraction helpers

Run-once data-prep helpers that document **how** the lightweight figure caches
were derived (so the multi-GB h5ads need not be re-loaded to understand or
regenerate them). Paths come from `config/` (no hardcoded `/Volumes/...`). Set
the global seed from config in any stochastic step.

**Cache contract (read this before assuming a prep edit will change a figure).**
These helpers write their outputs to `output_root` (see the Outputs column).
The per-figure render scripts, however, read the **deposited** copies of these
caches from `data_root` — that is the canonical, published input set. So:

- Re-running a `_prep` helper regenerates an *inspection* copy under
  `output_root`; it lets you audit how a cache was built and compare it against
  the deposited version.
- It does **not** change what the figures render. The figures consume the
  `data_root` cache, which the `_prep` helpers do not (and must not) overwrite.
- Editing a `_prep` helper therefore does **not** propagate into the figures.
  To deliberately re-derive a published cache you would regenerate it here and
  then have the figure point at the new `output_root` copy — but the default,
  faithful path keeps the figures on the deposited `data_root` caches.

The exception is `figures/_prep/01_export_tma_barcode_patient_map.py`, whose
`output_root/metadata/tma_barcode_patient_map.csv` IS the live input read by
Fig 4A/B (a derived map deposited under `output_root`, not `data_root`).

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
