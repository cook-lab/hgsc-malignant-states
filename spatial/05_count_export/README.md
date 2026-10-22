# `spatial/05_count_export/`: Xenium raw-count export

Writes, per sample, the raw Xenium counts (all 477 panel genes), the cell centroids and a label column as plain files.
The Figure 6 Xenium analyses and Tables S9 and S10 read these files instead of the SFEs. Both scripts are `as_run/`
copies.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 0 | `spatial/03_annotation_polarization/05_clean_split_rctd.R` (RCTD/SPLIT; run in stage 03) | SFEs | `$OUTPUT_ROOT/06g_clean_split/<sfe>_{rctd_meta.csv,purified_counts.rds}`; SFE `cell_label` | steps 1 and 2 |
| 1 | `as_run/split_recompute_stage1_polarization.R` | every `06g_clean_split/<sfe>_purified_counts.rds` with its `_rctd_meta.csv`; `06d_annotation_noBCAM/shared_gene_list_xenium_atlas.txt` (from `spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R`) | `split_recompute/purified_polarization_scores.csv` (UCell SecA / SecB / polarization on the purified counts of RCTD secretory cells, split at the frozen thresholds) | step 2 |
| 2 | `as_run/02_deg_precompute.R <sfe_name>`, once per sample | SFE raw `counts`, `cellSeg` centroids and `cell_label`; `$OUTPUT_ROOT/split_recompute/purified_polarization_scores.csv` | `$OUTPUT_ROOT/fig6_viewer/degdata/<sample>_counts.mtx.gz`, `<sample>_cells.csv.gz` (`barcode`, `x`, `y`, `cell_label`), `<sample>_genes.txt` | `figures/figure6/04_fig6fhi_xenium_niche_macrophages.R`, `figures/figure6/05_fig6fh_xenium_reciprocal_effects.R`, `figures/figure6/06_fig6g_macrophage_secB_neighbours.py`, `figures/figure6/07_fig6h_xenium_pct_expressing.py`, `figures/figure6/fig6_roi_picker/roi_server.py` (6i); `tables/09_tableS9_ucell_gene_sets.py` (`WT01_genes.txt`); Table S10 through `figures/figure6/04_fig6fhi_xenium_niche_macrophages.R` |

Samples exported in the final run: `sfe_WT01` … `sfe_WT08`, `sfe_tma` (TMA A, both slides; output stem `tma`),
`sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H`.

```bash
Rscript spatial/05_count_export/as_run/split_recompute_stage1_polarization.R      # edit its paths first (see below)
for s in sfe_WT01 sfe_WT02 sfe_WT03 sfe_WT04 sfe_WT05 sfe_WT06 sfe_WT07 sfe_WT08 sfe_tma sfe_COEUR_TMA_G sfe_COEUR_TMA_H; do
  Rscript spatial/05_count_export/as_run/02_deg_precompute.R $s
done
```

`02_deg_precompute.R` takes two optional arguments, a `sample_id` filter and an output stem; the final export did not
use them.

## What the export contains

- **Counts:** raw integer counts of the full panel (477 genes, including MMP7 and LGALS3), genes × cells.
- **Cell set:** cells with a non-empty SFE `cell_label`, the column last written by the RCTD/SPLIT step 0. For the
  whole tissues and TMA A this is every cell of the SingleR label tables. For TMA B it is not: 139,233 on-core tumour
  cells of block G and 196,586 of block H that are in the SingleR label tables are absent from the export.
- **Labels:** `cell_label` is the SFE label with the step 1 SecA / Intermediate / SecB overlay. The Figure 6 scripts do
  not use this column. They take cell identity from the SingleR label tables (`$OUTPUT_ROOT/singler_recompute/`) and
  keep only cells present in the export.

## `as_run/`

| File | Copied from |
|---|---|
| `as_run/split_recompute_stage1_polarization.R` | the recovered scratch area of the workspace project folder `202609_final_hgsc_spatial` |
| `as_run/02_deg_precompute.R` | the interactive explorer project folder `hgsc-spatial-explorer` (it served an ROI differential-expression viewer) |

Both are unchanged copies of the scripts that ran: all of their code builds or reports the export. They keep their
original paths:
- `split_recompute_stage1_polarization.R` names the final run's output folder for its inputs (`G`: `06g_clean_split`;
  the `06d_annotation_noBCAM` gene list) and its output (`OUT`: `split_recompute`). Point them at `$OUTPUT_ROOT/...`
  before running. Its SecA and SecB genes and the frozen thresholds are written inline.
- `02_deg_precompute.R` resolves its paths through `config/` and `$OUTPUT_ROOT`, and sources
  `spatial/00_setup/00_setup.R` by relative path. Run it from this repository's root, with
  `RERUN_SFE_DIR=$OUTPUT_ROOT/sfe` and `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8`.
