# atlas/02_annotation — cell-type annotation (level1 → level2), QC of labels, proportions

Resolves the integrated atlas into level-1 then level-2 cell types, validates the
epithelial labels, applies a metadata fix, and computes cell-type proportions /
counts. Entry object: `obj("atlas_scanvi")` (from `01_preprocess_qc`); the level-1
atlas (`hgsc_atlas_celltype_level1.h5ad`) and the final atlas (`obj("atlas_final")`)
are upstream/downstream entry-points.

## Ordered scripts

| # | Script | Inputs | Outputs |
|---|--------|--------|---------|
| 07 | `07_resolution_explorer.py` | `hgsc_atlas_celltype_level1.h5ad` | Leiden resolution sweep + scoring CSVs / HTML (resolution selection) |
| 07 | `07_epithelial_r04.py` | `hgsc_atlas_celltype_level1.h5ad`, `shared/signatures.yml` | epithelial Leiden r0.4 cluster stats / HTML |
| 08 | `08_finalize_celltype_level2.py` | `hgsc_atlas_celltype_level1.h5ad` | barcode→level2 maps + per-label marker CSVs; writes `obj("atlas_celltype_l2")` |
| 08 | `08_merge_celltype_level2.py` | level1 h5ad + barcode maps | `obj("atlas_celltype_l2")` (+ `is_low_complexity` flag) |
| 08c | `08c_rename_level2_labels.py` | `obj("atlas_celltype_l2")` + maps/markers | 19 reviewed level-2 label renames (in place, idempotent) |
| 09f | `09f_dotplot_canonical_markers.py` | `obj("atlas_final")` | `dotplot_stats.csv` + dotplot SVG/PDF → **SF4B** |
| 10c | `10c_epithelial_validation.py` | `obj("atlas_epithelial")` | epithelial validation CSVs + HTML report |
| 12a | `12a_fix_zhang2022_treatment_status.py` | level2 / epithelial / level1 h5ad | treatment_status fix (in place) — **must precede metadata export** |
| 12b | `12b_celltype_proportions.py` | `obj("atlas_celltype_l2")` | proportions figs + tables (supersedes step 12) |
| 12b | `12b_proportion_statistics.py` | proportion count tables | enrichment CSVs + statistical text report |
| 12c | `12c_celltype_counts.py` | `obj("atlas_celltype_l2")` | absolute counts figs + tables |
| 12d | `12d_epithelial_proportion_violins.py` | `obj("atlas_celltype_l2")` | epithelial composition violins by metadata |

All outputs are written under `output_root/02_annotation/`.

## Figures supported
SF4B (09f canonical-marker dotplot); composition substrate for Fig 1B-E and the
epithelial proportions feeding Fig 2D / SF9; metadata correctness for Fig 2F/2G and
Supp Data 1.

## Conventions applied
Central config for all paths/objects (`obj`, `path`); header docstrings; seeds set
from config (`np.random.seed(SEED)`, seeded UMAP/Leiden/silhouette subsample);
canonical noBCAM SecA/SecB sets loaded from `shared/signatures.yml` in the 07
scoring scripts (inline divergent lists removed); the epithelial polarization
**display** label "Transitioning" → "Intermediate" (12d short-label map).

### Label-naming note
The level-2 cell-type string **"Transitioning epithelial cell"** is a validated
annotation label (it keys the barcode maps and marker filenames) and is preserved
verbatim in 08c/12b/12c/10c. The manuscript "Transitioning → Intermediate"
standardisation applies to the epithelial **polarization** vocabulary
(SecA / Intermediate / SecB), which is where the rename was applied (12d display
labels; config `polarization.labels`).
