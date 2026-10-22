# `spatial/03_annotation_polarization/`: SingleR cell types, polarization and the label tables

Annotates every Xenium cell with SingleR against the scRNA-seq atlas, scores secretory cells for polarization
(SecB – SecA UCell score) in the gene space shared with the atlas, fixes the SecA / Intermediate / SecB thresholds
from the atlas, and writes one **SingleR label table** per sample. Those tables are what every figure and table script
reads; they are not the SFE `cell_label` column.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `01_annotation_singler.R` | normalised SFEs in `sfe_dir` (`sfe_names_all`); SingleR reference `$DATA_ROOT/2026_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad` (config `xenium_ref`); `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` | SFE `singler_label`, `singler_pruned`, `singler_score` (16 types; 429 genes: panel ∩ reference minus the 28 excluded genes); two console lines on the gene counts | 2, 3, 6; `spatial/02_qc/05b_core_qc.R`; `tables/07_tableS7_xenium_gene_panel.py` (parses the console log) |
| 2 | `02_adaptive_secretory_noBCAM.R` | SFE `singler_label`, `logcounts` | SFE `cell_label`, `secretory_subtype`, `score_SecA`, `score_SecB`, `ratio_SecB_SecA` | structural only: step 4 takes the secretory cells from `cell_label`; `spatial/02_qc/05b_core_qc.R` reads the column |
| 3 | `03_ucell_scoring_noBCAM.R` | SFEs; atlas UCell export in `ATLAS_UCELL_DIR` (default `$DATA_ROOT/2026_final_atlas/output/18_ucell_atlas/`: `atlas_gene_list.txt`, `atlas_secretory_{counts.mtx.gz,barcodes.tsv,genes.tsv,metadata.csv}`) | SFE `SecA_UCell`, `SecB_UCell`, `polarization_UCell` for SingleR "Secretory epithelium" cells (UCell on raw counts, 457 genes shared with the atlas, maxRank = 457); `$OUTPUT_ROOT/06d_annotation_noBCAM/{xenium,atlas}_ucell_scores.csv`, `checkpoint_{xenium,atlas}_scores.rds`, `shared_gene_list_xenium_atlas.txt` | 4 (atlas and SFE scores), 5, 6; `spatial/02_qc/05b_core_qc.R`; `figures/figure4/04_fig4d_bivariate_moran.R`, `figures/figure4/06_fig4e_bilisa.R` (SFE scores); `tables/06_tableS6_secA_secB_signatures.py` and `spatial/05_count_export/as_run/split_recompute_stage1_polarization.R` (gene list) |
| 4 | `04_reclassification_polarization.R` | `$OUTPUT_ROOT/06d_annotation_noBCAM/atlas_ucell_scores.csv`; SFE `cell_label`, `polarization_UCell` | `$OUTPUT_ROOT/06f_reclassification_polarization/threshold_summary.csv` (t_low = 0.109429, t_high = 0.210501, with the atlas AUC, accuracy, macro-F1 and class counts) and `.sfe_update_complete`; SFE `cell_label` rewritten for secretory cells (previous label kept as `cell_label_06e`) | 5 (reads the file); the values are typed into step 6, `spatial/02_qc/05b_core_qc.R` and the count export |
| 5 | `05_clean_split_rctd.R` | SFEs; SingleR reference; `genes_exclude_singler.txt`; `threshold_summary.csv` | `$OUTPUT_ROOT/06g_clean_split/<sfe>_{rctd_meta.csv,purified_counts.rds}`, `06g_summary<suffix>.csv` (per-SFE RCTD counts), `.sfe_update_complete`; SFE `cell_label` overwritten (previous label kept as `cell_label_06f`), `rctd_*` columns and a `purified` assay added | structural only: `spatial/05_count_export/` |
| 6 | `06_extract_singler.R <sfe_name>`, once per SFE | SFE `singler_label`, `singler_pruned`, `polarization_UCell`, coordinates, `patient_id`, `core_id`, `sample_type`, `cell_area` | `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` | every Xenium figure and table script (directly or through `spatial/00_setup/singler_labels.R`); `spatial/02_qc/06b_final_tma_cores.py` |

Steps 1–5 loop over every SFE in `sfe_names_all`. Steps 1, 3, 4 and 5 are heavy (step 5 peaks at about 21 GB per whole
tissue and 30–40 GB for a TMA).

```bash
mkdir -p $OUTPUT_ROOT/_logs
Rscript spatial/03_annotation_polarization/01_annotation_singler.R 2>&1 | tee $OUTPUT_ROOT/_logs/singler_parallel.log
Rscript spatial/03_annotation_polarization/02_adaptive_secretory_noBCAM.R
Rscript spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R
Rscript spatial/03_annotation_polarization/04_reclassification_polarization.R
Rscript spatial/03_annotation_polarization/05_clean_split_rctd.R --all      # RCTD on every SFE, then writeback
for s in sfe_WT01 sfe_WT02 sfe_WT03 sfe_WT04 sfe_WT05 sfe_WT06 sfe_WT07 sfe_WT08 sfe_tma sfe_COEUR_TMA_G sfe_COEUR_TMA_H; do
  Rscript spatial/03_annotation_polarization/06_extract_singler.R $s
done
```

`tables/07_tableS7_xenium_gene_panel.py` reads step 1's console log from `$OUTPUT_ROOT/_logs/singler_parallel.log`
and checks its lines "Shared genes (panel ∩ reference): 457" and "After removing 28 excluded genes: 429 genes for
SingleR" against Table S7; hence the `tee` above.

## The SingleR label tables

Columns: `barcode`, `x`, `y`, `singler_label`, `singler_pruned`, `polarization_UCell`, `cell_label_singler`,
`patient_id`, `core_id`, `sample_type`, `cell_area`.

`cell_label_singler` is `singler_label` with "Secretory epithelium" split on `polarization_UCell`: SecB if
≥ 0.210501, SecA if ≤ 0.109429, Intermediate otherwise (NA when the score is missing). The 16 SingleR types thus become
the 18 cell types of the manuscript. Step 4 itself classifies with SecA < t_low; step 6 and `02_qc/05b` use ≤ t_low.
The SecA and SecB genes (7 + 7, without BCAM) come from `shared/signatures.yml`.

The thresholds are the 75th percentile of polarization in atlas SecA cells (t_low) and the 25th percentile in atlas
SecB cells (t_high), with the atlas states taken from the NMF labels and atlas cells scored on the same 457 genes.

## Steps kept for structure only

- **Step 2** writes an initial ratio-based secretory subtype into `cell_label`. No panel uses these labels, but step 4
  takes its secretory cells from that column and `spatial/02_qc/05b_core_qc.R` reads it, so both fail without it.
- **Step 5** runs RCTD/SPLIT and is the final writer of the SFE `cell_label`. No panel uses these labels or the
  `purified` assay. The step is needed because the Xenium raw-count export (`spatial/05_count_export/`) reads its
  `06g_clean_split/` tables, and the export's cell set is the cells with a non-empty `cell_label` after this step.

## Environment

| Variable | Final run | Step |
|---|---|---|
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` | all |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | all |
| `SINGLER_WORKERS` | `12` (default) | 1 |
| `ATLAS_UCELL_DIR` | unset: the deposited export `$DATA_ROOT/2026_final_atlas/output/18_ucell_atlas/` | 3 |
| `RERUN_WRITEBACK_ONLY` | unset (a comma-separated SFE list limits the writeback to those SFEs) | 5 |

- Step 1 reads each SFE from `sfe_dir` and saves it under `$OUTPUT_ROOT/sfe`; an SFE that already has
  `singler_label` is skipped.
- Step 3 reuses `checkpoint_*_scores.rds` (or both score CSVs) in `06d_annotation_noBCAM/` if present, instead of
  re-running UCell.
- Steps 4 and 5 skip their SFE writeback when their `.sfe_update_complete` file exists.
- Step 5 also accepts a comma-separated SFE list with `--suffix <chunk>`, and `--writeback` alone; without arguments
  it runs the first whole tissue only.
