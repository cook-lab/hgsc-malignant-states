# atlas/08_xenium_reference — atlas → Xenium SingleR reference bridge

Builds the cross-platform cell-type reference used by the Xenium spatial pipeline:
collapses the atlas level-2 annotation into a 16-type `xenium_celltype` schema,
downsamples to a balanced reference, computes per-type markers, and reports Xenium
panel coverage. The downsampled object is the canonical SingleR reference
(`obj("xenium_ref")`).

## Ordered scripts

| # | Script | Inputs | Outputs |
|---|--------|--------|---------|
| 1 | `add_xenium_celltype.py` | `obj("atlas_celltype_l2")` | `DATA_ROOT/2026_final_atlas/hgsc_atlas_xenium.h5ad` (adds `xenium_celltype`; 16 types via Epithelial→Ciliated/Secretory, DC→pDC/cDC, T/NK→T/NK) |
| 2 | `xenium_celltype_markers_and_report.py` | `hgsc_atlas_xenium.h5ad` + `output/xenium_panel_genes.txt` | `obj("xenium_ref")` = `output/xenium_celltype/xenium_celltype_downsampled.h5ad` (1000 cells/type), per-type marker CSVs, panel-coverage HTML |

Marker CSVs and the coverage report are written under `output_root/08_xenium_reference/`;
the downsampled reference object is written to its config-registered path
(`obj("xenium_ref")`), the canonical copy consumed by the spatial `06_annotation` step.

## Figures supported
Cross-platform bridge only. The SingleR reference produced here underpins every
Xenium cell-type panel (Fig 4–6, SF10/SF11) but renders no panel directly.

## Conventions applied
Central config for all paths/objects; header docstrings; seeds set from config
(`np.random.seed(SEED)`, `np.random.default_rng(SEED)` for the per-type downsample);
no `Transitioning` labels and no inline SecA/SecB signatures occur in this stage.
