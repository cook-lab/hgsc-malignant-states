# atlas/08_xenium_reference: the SingleR reference for Xenium

This stage builds the atlas reference that SingleR uses to annotate the Xenium cells. It assigns 16 reference cell
types (`xenium_celltype`). The 13 level-1 lineages are kept, with three of them split:
- epithelial → secretory and ciliated
- DC → conventional and plasmacytoid
- T/NK → T and NK

The reference is downsampled to 1,000 cells per type. On Xenium, the secretory label is later split into SecA,
Intermediate and SecB by polarization, which gives the 18 Xenium cell types.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`. Both steps are the `as_run/` scripts that wrote the
deposited reference.

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `as_run/add_xenium_celltype.py` | `hgsc_atlas_celltype_level2.h5ad` (the level-2 object of `../02_annotation`, before its final filter) | `hgsc_atlas_xenium.h5ad`, which adds `xenium_celltype` | step 2 |
| 2 | `as_run/xenium_celltype_markers_and_report.py` | `hgsc_atlas_xenium.h5ad` | **`output/xenium_celltype/xenium_celltype_downsampled.h5ad`** (deposited; config `xenium_ref`): 1,000 cells per type, drawn with seed 42 | `spatial/03_annotation_polarization/01_annotation_singler.R` (Xenium SingleR) and `05_clean_split_rctd.R`; `figures/supplementary/19_figS14_annotation_validation_plot.py` (Figure S14); `tables/07_tableS7_xenium_gene_panel.py` (Table S7, gene names) |

The deposited reference is the entry point for all Xenium work. Steps 1–2 cannot be re-run on the deposit as written,
because their inputs are not deposited.

## `as_run/`

These are copies of the scripts that ran in the atlas project folder `202601_final_atlas`. They keep that project's
paths. `as_run/xenium_celltype_markers_and_report.py` is trimmed to the 1,000-cells-per-type downsample that writes
the reference: its one-vs-rest marker CSVs, panel-coverage analysis and HTML report were removed, and the copy opens
with a note on what was removed. `as_run/add_xenium_celltype.py` had nothing to remove and is unchanged. The untrimmed
script stays in `202601_final_atlas`.

## Missing code

None for this stage.
