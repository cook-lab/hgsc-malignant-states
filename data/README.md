# Data

The large input objects are **not** stored in this repository (~1.4 TB total). They are
deposited at `<GEO/Zenodo accession — TBD>` and described in `docs/DATA_AVAILABILITY.md`.

To run the pipeline, download / point to them and set `DATA_ROOT` (or edit `config/config.yml`)
so these resolve:

| Key (`config.objects`) | Object | Role |
|---|---|---|
| `atlas_scanvi` | `hgsc_atlas_scanvi.h5ad` | scANVI integration output (deposited for convenience; reproducible via `atlas/01` step 03) |
| `atlas_final` | `hgsc_atlas_final.h5ad` | canonical atlas downstream entry-point |
| `atlas_epithelial` | `hgsc_atlas_epithelial.h5ad` | epithelial subset |
| `atlas_celltype_dir` | `celltype_h5ad/` | per-celltype subsets |
| `xenium_ref` | `xenium_celltype_downsampled.h5ad` | SingleR spatial reference |
| `sfe_tma_filtered` | `sfe/sfe_tma_filtered` | canonical TMA SpatialFeatureExperiment |
| `sfe_dir` | `sfe/` | per-whole-tissue SFEs |

The scVI/scANVI **integration is included in the repo** (the original cluster scripts `atlas/01_preprocess_qc/02_aggregate.py` … `07_process.py`):
the pipeline covers every step from raw data through final figures, including integration
(CellAssign → scVI → scANVI). The integrated object is **deposited for convenience** so you can
start downstream analysis from it directly, but it is **fully reproducible from raw data** via the
included integration code. Re-running integration is computationally expensive (originally a GPU
cluster job) and is not required to reproduce downstream results.

## The deposit also includes the analysis output caches

Because the repo is **deposit-driven** (figure/table scripts read their inputs from `DATA_ROOT`,
not from a freshly-run backend), the deposited bundle includes — alongside the entry objects —
the per-stage **analysis output caches** the figures consume, preserving the original directory
layout the scripts expect:

| Tree | Cache location (under DATA_ROOT) | Examples |
|---|---|---|
| atlas | `2026_final_atlas/output/<stage>/` | `11d_epithelial_nmf/`, `21_epitype_functional_characterization/`, `19_cnv/`, `18_ucell_atlas/`, `20_consensusov/`, `fig_secretory_polarization/data/` |
| atlas (fig data) | `202605_epitype_manuscript/final_publication_figures/data_fig1/` and `20260411_figures/`, `20260429_figures/data/` | per-panel CSV/parquet caches |
| spatial | `2026_final_xenium_analysis/output/<stage>/` | `06f_reclassification_polarization/`, `09_neighborhood/`, `10_clinical_v2/`, `16b_*`, `19d_*`, `22_vascular_proximity/`, `33_morphometrics/`, `44_spatial_autocorrelation/`, `06g_clean_split/` |

To **regenerate** these caches from the entry objects instead, run the `atlas/` and `spatial/`
backend stages (they write to `OUTPUT_ROOT`); see the repo README §4.

