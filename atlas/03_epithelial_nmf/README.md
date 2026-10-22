# atlas/03_epithelial_nmf: epithelial NMF programs, gene classes and atlas UCell scores

This stage learns the epithelial gene programs with NMF (k = 10) on the epithelial compartment (575,366 cells).
Factor 2 is the SecB program and Factor 3 the SecA program. Factor 2 usage splits secretory cells into SecA,
Intermediate and SecB. The stage also:
- classifies genes as SecA-specific, SecB-specific or shared (177 / 177 / 124)
- checks the programs with per-patient consensus NMF
- exports the atlas secretory cells for UCell, which gives the atlas reference scores used by Xenium and the organoids
- writes the epithelial UMAP table (`meta.parquet`) that the figures read

## Run order

These are the `as_run/` scripts that wrote the deposited outputs, in dependency order. Paths without a root are under
`$DATA_ROOT/2026_final_atlas/`.

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `as_run/11d_epithelial_nmf.py` | `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` (counts, `celltype_level2`, `X_umap_local`) | `output/11d_epithelial_nmf/{11d_nmf_usage.csv, 11d_nmf_usage_raw.csv, 11d_nmf_loadings.csv, 11d_factor_scores.csv}` | `figures/_prep/atlas_final_tables.py` (usage, hence the states of every epithelial panel); Tables S5–S6; Figure S10; the KPCA.B alignment (`models/kpca_mouse/`); steps 2–4 and 7; stages 04 and 05 |
| 2 | `as_run/11e_nmf_characterization.py` | `11d_nmf_loadings.csv`, `11d_nmf_usage.csv` | `output/11e_nmf_characterization/11e_gene_classification.csv` (top-decile Factor 3 / Factor 2 loadings → SecA-specific, SecB-specific, shared) | Tables S5–S6; TCGA scores (`../07_deconvolution_survival`); `models/kpca_mouse/scripts/02_orthologs.py`; step 3b |
| 3a | `as_run/11f_consensus_nmf.py` | epithelial subset (including `X_umap_local`); `11d` usage and loadings | `output/11f_consensus_nmf/{11f_mp_summary.csv, 11f_mp_top_genes.csv, 11f_patient_program_counts.csv}` (per-patient NMF, k 5–10, 30 runs per k, clustered into metaprograms) | Table S5 |
| 3b | `as_run/11g_secAB_consensus_nmf.py` | epithelial subset; `11e` classes; `11d` outputs | `output/11g_secAB_consensus_nmf/{11g_mp_summary.csv, 11g_mp_top_genes.csv, 11g_mp_gene_weights.csv}` (the same on the 354 SecA- and SecB-specific genes, k 2–6) | Table S5 (summary and top genes) |
| 4 | `as_run/17_prepare_nmf_labels.py` | `hgsc_atlas_final.h5ad` obs; `11d_nmf_usage.csv` | `output/17_nmf_labels/celltype_nmf_mapping.csv` (per atlas cell: the level-2 label, or SecA / Transitioning (= Intermediate) / SecB / ciliated for epithelial cells) | step 5; `../07_deconvolution_survival` |
| 5 | `as_run/18a_ucell_atlas_export.py` | epithelial subset; step 4 mapping | `output/18_ucell_atlas/atlas_secretory_{counts.mtx.gz, barcodes.tsv, genes.tsv, metadata.csv}`, `atlas_gene_list.txt` (the SecA, Intermediate and SecB cells) | step 6; the Xenium UCell step `spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R`, which scores these cells in the Xenium gene space and sets the polarization thresholds; `atlas_secretory_metadata.csv` also → Fig 3b–e (`figures/figure3/`) and `models/organoids/scripts/` |
| 6 | `as_run/18b_ucell_atlas_scoring.R` | step 5 export; the PDO baseline object (`models/organoids/scripts/03_create_seurat_untreated_baseline.R`) for the shared gene space; the 7 + 7 genes | `output/18_ucell_atlas/{atlas_ucell_scores.csv, shared_gene_list.txt}` (14,271 genes shared with the organoids) | Fig 3b–e atlas reference (`figures/figure3/`); organoid UCell scores and comparisons (`models/organoids/scripts/`); Table S6 (gene space) |
| 7 | `as_run/fig_secretory_polarization_00_prepare_data.py` | epithelial subset with `X_umap_local` (`../02_annotation/as_run/9e_celltypelevel2_umaps.py`); `11d_nmf_usage.csv` | `output/fig_secretory_polarization/data/meta.parquet`: `UMAP1`/`UMAP2`, `Factor_2`, `schema_nmf` (SecA / Transitioning / SecB / Ciliated) and cell metadata | `figures/_prep/atlas_final_tables.py` (`UMAP1`/`UMAP2` only: the epithelial UMAP); `tables/06_tableS6_secA_secB_signatures.py` and `figures/figure6/10_fig6h_atlas_epithelial_pct.py` (`schema_nmf`: Table S6 and Fig 6h) |

Nothing in the repository reads `11d_nmf_usage_raw.csv` or `11g_mp_gene_weights.csv`.

**The NMF fit.** k = 10 on log-normalised expression of 3,007 genes: 3,000 seurat_v3 HVGs plus the forced SecA (9) and
SecB (8) NMF marker lists, which are the 7 + 7 genes of `shared/signatures.yml` plus BCAM and WT1 (SecA) and MMP7
(SecB). The fit uses a training draw of 49,997 epithelial cells (target 50,000), stratified by level-2 type with at
least 100 cells per type, and is then applied to all 575,366 epithelial cells. The loader turns Factor 2 usage into
states among non-ciliated cells: below p50 is SecA, p50–p75 Intermediate, and at or above p75 SecB.

**Re-running the NMF.** Steps 1, 3a, 3b and 7 read `X_umap_local`, which the deposited `celltype_h5ad/` subsets do
not carry. `01_epithelial_nmf.py` re-implements step 1 without it, on the deposited epithelial object, and writes to
`$OUTPUT_ROOT/03_epithelial_nmf/`:

| Script | Reads | Writes |
|---|---|---|
| `01_epithelial_nmf.py` | `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` | `11d_nmf_usage.csv`, `11d_nmf_loadings.csv`, `11d_factor_scores.csv`, `nmf_factor_mapping.json` |
| `01b_verify_against_deposited.py` | as above, plus the deposited `11d_nmf_loadings.csv` and `11d_nmf_usage.csv` | `verify_deposited/{verify_summary.json, loadings_by_factor.csv, feature_set_differences.csv}` |

`01b_verify_against_deposited.py` runs the feature selection, training draw and NMF fit of `01_epithelial_nmf.py` and
compares them with the deposited run: the 3,007 features are the same, and the loadings and the usage of the 49,997
training cells match to floating-point precision (largest difference below 10⁻¹²). It does not project all 575,366
cells. Figures read the deposited `11d` outputs.

## Environment

- **`ATLAS_UCELL_DIR`.** The Xenium step `spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R` reads the
  deposited `output/18_ucell_atlas/` export, which the final run used. Set `ATLAS_UCELL_DIR` only to point it at
  another copy.
- **Software.** Steps 1–5 and 7 ran in E2 (Python). Step 6 ran in E1 (R, UCell).

## `as_run/`

These are trimmed copies of the scripts that ran in the atlas project folder `202601_final_atlas`. The code that wrote
the deposited outputs that the paper, the tables and later steps read is kept as it ran. Figures, HTML reports, the
exploratory analyses of `11e` (the copy keeps only the gene classification) and the CSVs that nothing reads were
removed. The copy of `17` leaves out its label-summary CSV and an in-place write of a `celltype_nmf` column into
`hgsc_atlas_final.h5ad`; the deposited object has no such column, and nothing reads it. Each trimmed copy opens with
a note on what was removed. `18a_ucell_atlas_export.py` and `18b_ucell_atlas_scoring.R` had nothing to remove and are
unchanged. The untrimmed scripts stay in `202601_final_atlas`. The copies keep that project's paths; `11e` finds the
project folder from its own location (the parent of its `scripts/` folder).

## Missing code

None for this stage. The epithelial UMAP itself was computed in `../02_annotation/as_run/9e_celltypelevel2_umaps.py`.
