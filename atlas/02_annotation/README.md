# atlas/02_annotation: cell-type annotation and the final atlas object

This stage annotates the integrated atlas at two levels: 13 lineages (level 1) and 51 cell types and states
(level 2). It then removes the clusters marked "Excluded" and corrects the Zhang 2022 treatment status and the
Vázquez-García metastatic site. The result is `hgsc_atlas_final.h5ad` (1,980,703 cells), the object every atlas
figure starts from. The stage also writes the per-lineage subsets and the per-lineage UMAPs, including the epithelial
UMAP that the figures read.

## Run order

Every step ran as the `as_run/` copy listed here; the steps are in dependency order. Paths without a root are under
`$DATA_ROOT/2026_final_atlas/`.

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `as_run/06_celltype_level1_leiden.py` | `hgsc_atlas_scanvi.h5ad` (`../01_preprocess_qc`) | `hgsc_atlas_celltype_level1.h5ad`: Leiden r = 0.3 on `X_scanvi`, level-1 labels from the CellAssign composition of each cluster, cluster 8 subclustered into smooth muscle and pericyte, two artifact clusters (38 cells) removed; `output/06_celltype_level1/leiden_r0.3_crosstab.csv` (Leiden × CellAssign, the table the map was read from) | step 2 |
| 2 | `as_run/06b_remap_level1.py` | level-1 object | the same object, with the cluster → lineage map corrected (13 lineages) | step 3 |
| 2b | `as_run/xenium_coverage_report_level1.py` | level-1 object | `output/06_celltype_level1/markers/level1_<type>.csv`: one-vs-rest Wilcoxon per level-1 type on at most 5,000 cells per type (seed 42), adjusted p < 0.05 and log fold change > 0 | Table S2 |
| 3 | `as_run/08_finalize_celltype_level2.py <compartment>` (or `--all`) | level-1 object | per compartment: subclustering on `X_scanvi`, level-2 labels, Wilcoxon markers (`output/08_celltype_level2/markers/<compartment>_<type>.csv`) and a barcode → level-2 map (`output/08_celltype_level2/barcode_maps/`) | step 4; Table S3 (`tables/02_tableS2-4_annotation_markers.py` reads the deposited `markers/`) |
| 4 | `as_run/08_finalize_celltype_level2.py --merge` | level-1 object; barcode maps | `hgsc_atlas_celltype_level2.h5ad` (adds `celltype_level2`) | step 5; `../08_xenium_reference` |
| 4b | `as_run/08b_rename_level2_labels.py` | level-2 object; barcode maps; marker files | renames after expert review of the level-2 report (excluded clusters rescued under retained names, the remaining excluded clusters renumbered), in place, in the object, the maps and the marker file names | step 5; Table S3 marker file names |
| 5 | `as_run/08c_rename_level2_labels_v2.py` | level-2 object; barcode maps; marker files | 19 level-2 renames, in place, in the object, the maps and the marker file names | step 6 |
| 6 | `as_run/12a_fix_zhang2022_treatment_status.py` | level-1 and level-2 objects; `hgsc_atlas_epithelial.h5ad` (the lab's epithelial object, not deposited) | Zhang 2022 `treatment_status` changed from "post-treatment" to "post-chemotherapy", in place | treatment panels (Fig 1d, Fig 2c,e, Figures S3a and S9a); Table S1 |
| 7 | `as_run/09d_publication_figures.py` | level-2 object | `hgsc_atlas_final.h5ad`, with the level-2 "Excluded" clusters removed (2,294,855 → 1,980,703 cells) | all atlas figures |
| 8 | `as_run/10c_rename_inflammatory_macrophage.py` | final and level-2 objects; the macrophage barcode map and marker file | "Inflammatory macrophage" renamed "Monocyte-derived macrophage", in place | level-2 labels (Table S3; Fig 6f,h) |
| 9 | `as_run/9a_atlas_celltypelevel1_v2.py` | `hgsc_atlas_final.h5ad` | `celltype_h5ad/hgsc_atlas_final_<lineage>.h5ad`, one per lineage (13) | stages 03, 04, 05 and 07; Fig 1g prep, Fig 2a (scFEA), 2b, 2g, Fig 6f,h, Figure S6, Table S6 |
| 10 | `as_run/9e_celltypelevel2_umaps.py` | lineage subsets | a per-lineage UMAP on `X_scanvi` (n = 15, min_dist 0.3), saved as `obsm["X_umap_local"]` in each subset | `../03_epithelial_nmf/as_run/` (`11d`, `11f`, `11g` and `fig_secretory_polarization_00_prepare_data.py` → `meta.parquet`, whose UMAP1/UMAP2 are the epithelial UMAP of every epithelial UMAP panel) |
| 11 | `as_run/fix_metastatic_site_vazquez.py --apply` | `hgsc_atlas_final.h5ad` | `metastatic_site` of the 7 Vázquez-García ascites samples (21,783 cells) changed from metastasis to ascites, in place, so that it matches `anatomic_site`; a backup of the original codes is written first | site panels (Fig 1c, Fig 2c,d, Figure S3c) |
| 12 | `as_run/fix_metastatic_site_docs.py` | `pre-processing and integration/atlas_metadata.xlsx`; `output/supplemental/T1_atlas_metadata.csv` | the same correction in both files | Table S1 (`tables/01_tableS1_atlas_studies.py`) |

**Not re-runnable on the deposit.** The level-1 and level-2 objects are not deposited, so steps 2–8 cannot be re-run
as written. The deposited `hgsc_atlas_final.h5ad` carries the labels they produced. The deposited `celltype_h5ad/`
subsets do not carry `X_umap_local`; the epithelial UMAP survives as UMAP1/UMAP2 in
`output/fig_secretory_polarization/data/meta.parquet`.

**How the figures use this stage.** The figures do not read `hgsc_atlas_final.h5ad` metadata directly. They read it
through `figures/_prep/atlas_final_tables.py`, which also applies the Nath 2021 sample relabelling
(`figures/_prep/nath2021_fix.py`; 158 patients, 362 samples).

## `as_run/`

These are trimmed copies of the scripts that ran in the atlas project folder `202601_final_atlas`, including the two
`fix_metastatic_site_*` in-place fixes. The code that wrote the objects and files above is kept as it ran. Figures,
summary prints, CSVs that nothing reads and a closing pointer to follow-up scripts (`12a`) were removed; the trimmed
copies (`06`, `06b`, `09d`, `12a`, `9a`, `9e`, `xenium_coverage_report_level1`) open with a note on what was removed. `08`, `08b`, `08c`, `10c`
and the two `fix_metastatic_site_*` fixes had nothing to remove and are unchanged. The untrimmed scripts stay in
`202601_final_atlas`. The copies keep that project's paths.

## Missing code

- The level-1 marker CSVs (`output/06_celltype_level1/markers/level1_*.csv`, the source of Table S2) are written by
  `as_run/xenium_coverage_report_level1.py`. `output/supplemental/T2_celltypelevel1_markers.csv.csv` is the same 13
  files stacked (51,592 gene–type pairs); the step that stacked them is not recorded.
  `tables/02_tableS2-4_annotation_markers.py` reads the deposited files.
- **[code missing]** The export of `output/supplemental/T3_epithelial_markeres.csv`, the ciliated-vs-secretory Wilcoxon
  table behind Table S4. Its 3,850 Ciliated-up rows equal the ciliated-cluster markers that the atlas project's
  `scripts/11v2_deg.py` and `11v4_deg.py` wrote (`output/11v2_DEG/`, `11v4_DEG/deg_ciliated_epithelial_cell.csv`); the
  export that added the 6,179 Secretory-up rows and the percentages is not recorded. The table builder reads the
  deposited file.
- **[code missing]** The first export of the per-cell metadata table `output/supplemental/T1_atlas_metadata.csv`
  (Table S1). Only its later correction (`as_run/fix_metastatic_site_docs.py`) is available.
