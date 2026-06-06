# Data

The large input objects are **not** stored in this repository (~1.4 TB total). They are
deposited at `<GEO/Zenodo accession — TBD>` and described in `docs/DATA_AVAILABILITY.md`.

To run the pipeline, download / point to them and set `DATA_ROOT` (or edit `config/config.yml`)
so these resolve:

| Key (`config.objects`) | Object | Role |
|---|---|---|
| `atlas_scanvi` | `hgsc_atlas_scanvi.h5ad` | scANVI integration output (**trust boundary**) |
| `atlas_final` | `hgsc_atlas_final.h5ad` | canonical atlas downstream entry-point |
| `atlas_epithelial` | `hgsc_atlas_epithelial.h5ad` | epithelial subset |
| `atlas_celltype_dir` | `celltype_h5ad/` | per-celltype subsets |
| `xenium_ref` | `xenium_celltype_downsampled.h5ad` | SingleR spatial reference |
| `sfe_tma_filtered` | `sfe/sfe_tma_filtered` | canonical TMA SpatialFeatureExperiment |
| `sfe_dir` | `sfe/` | per-whole-tissue SFEs |

The scVI/scANVI **integration is a trust boundary** — it was run on a compute cluster and is
**not** reproduced here; its output object is the starting point for everything downstream.
