# atlas/01_preprocess_qc — preprocessing, QC, doublets, post-integration refilter

Raw-data ingestion through the QC'd, doublet-filtered, post-integration atlas
object that is the entry point for all downstream analysis.

## Ordered scripts

| # | Script | Inputs | Outputs |
|---|--------|--------|---------|
| 01 | `01_dataprep.py` | per-study raw matrices (`DATA_ROOT/2026_final_atlas/raw/<study>/`), `atlas_metadata.xlsx` | harmonised per-study h5ad → `processed/<study>.h5ad` |
| 02 | `02_concat_qc_doublets.py` | 13 per-study h5ad | `processed/atlas_concat_counts_only_X.h5ad`, then `processed/atlas_concatenated_filtered.h5ad` (QC: total_counts≥500, n_genes≥300, sample≥500 cells; Scrublet score < 0.3) + QC violins |
| — | **scVI → scANVI integration** | `atlas_concatenated_filtered.h5ad` | `integrated_scanvi.h5ad` → **TRUST BOUNDARY (not reproduced here)** |
| 03b | `03b_refilter_umap.py` | `integrated_scanvi.h5ad` | `obj("atlas_scanvi")` = `hgsc_atlas_scanvi.h5ad` (post-integration Scrublet < 0.25, neighbours, UMAP) + QC UMAPs |
| 04 | `04_atlas_qc_review.py` | `obj("atlas_scanvi")` | per-study QC summary CSVs + QC histograms/UMAPs |
| 05 | `05_umap_suite.py` | `obj("atlas_scanvi")` | 9 publication metadata UMAP SVGs |

All outputs are written under `output_root/01_preprocess_qc/`.

## scVI / scANVI integration — TRUST BOUNDARY

The scVI/scANVI integration was run on a compute cluster; the training code is not
part of this repository and the step is **not re-executed here**. It consumes the
concatenated, QC'd, doublet-filtered object produced by step 02
(`atlas_concatenated_filtered.h5ad`) and emits `integrated_scanvi.h5ad`.

Reproducible analysis resumes at **03b**, which applies the stricter post-integration
doublet filter and recomputes the embedding, writing
`hgsc_atlas_scanvi.h5ad`. This object is registered in `config/config.yml` as the
entry-point key **`atlas_scanvi`** (`obj("atlas_scanvi")`). Everything downstream of
03b runs locally and is reproducible from the deposited objects.

Steps 01 and 02 are retained as documented provenance for how the pre-integration
counts object was built; they require the raw per-study GEO/EGA inputs (not deposited)
and `mygene` network access. If you start from the deposited `hgsc_atlas_scanvi.h5ad`,
you do not need to run 01/02 or the integration.

## Figures supported

QC/overview substrate only. SF1 (per-study QC violins), SF2 (integration UMAPs),
Fig 1A (atlas UMAP) and SF3 (metadata UMAPs) are rendered by per-panel scripts in
`figures/` from `obj("atlas_scanvi")` / `obj("atlas_final")`.

## Conventions applied

Central config for all paths/objects; header docstrings; `np.random.seed(SEED)` and
seeded Scrublet (`random_state=SEED`) / UMAP (`random_state=SEED`); dead duplicate
notebook diagnostic cells removed; QC plots written to file instead of `plt.show()`.
No `Transitioning` labels and no inline SecA/SecB signatures occur in this stage.
