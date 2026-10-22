# `spatial/02_qc/`: cell QC, probe QC, core QC and the final TMA core set

Removes technical artefacts at the cell level, records which panel genes are kept out of annotation (probe QC), and
builds the one TMA core set that every TMA analysis uses (378 patients, 503 cores, 1,542,448 cells). The stage runs in
two parts: cell QC before annotation, core QC after `spatial/03_annotation_polarization/` (it needs the SingleR labels
and polarization scores).

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `04_qc_filter_normalize.R` | raw SFEs in `sfe_dir` (`RERUN_SFE_DIR`; all of `sfe_names_all`) | SFEs rewritten in place: cells kept with ≥ 10 counts, ≥ 5 genes, ≤ 5% control-probe counts and area 10–500 µm²; 64 control-probe features removed (477 genes left); `logcounts` and `logcounts_area` added; `$OUTPUT_ROOT/03_04_qc/filtering_summary.csv` | `spatial/03_annotation_polarization/` and everything after |
| 2 | `as_run/explore_probe_qc.Rmd` (sources `as_run/00_setup.R`) | the SFEs of the project folder `2026_final_xenium_analysis` (TMA A and the 8 whole tissues, before the final re-run), 200,000-cell stratified subsample; that project's copy of the SingleR reference (16,000 cells, 16 types); `metadata/panel_markers.csv` | in `2026_final_xenium_analysis/output/05_probe_qc/`: `probe_qc_full.csv`, `flag_qc.csv`, `gene_exclusion_decisions.csv`, `genes_exclude_singler.txt` (28 genes), `genes_monitor_singler.txt` (read by no repo script) | step 3; `spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R` (`flag_qc.csv`); `tables/07_tableS7_xenium_gene_panel.py` (the three CSVs); `tables/as_run/generate_xenium_gene_panel_table.R` (`probe_qc_full.csv`, `gene_exclusion_decisions.csv`) |
| 3 | manual copy | `$DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/genes_exclude_singler.txt` | `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (identical copy) | `spatial/03_annotation_polarization/01_annotation_singler.R`, `05_clean_split_rctd.R`; `figures/supplementary/19_figS14_annotation_validation_plot.py`; `tables/07_tableS7_xenium_gene_panel.py` |
| — | *run stage 03 here* | | | |
| 4 | `05b_core_qc.R`, once per array | SFE `singler_label`, `singler_score`, `singler_pruned`, `polarization_UCell`, `core_id`, `patient_id`, `sample_type` (and the `cell_label` column) of `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` | `$OUTPUT_ROOT/07_core_qc/core_qc_summary_singler_{tma,COEUR_TMA_G,COEUR_TMA_H}.csv` | step 5 |
| 5 | `06b_final_tma_cores.py` | step 4 summaries; `$OUTPUT_ROOT/singler_recompute/{tma,COEUR_TMA_G,COEUR_TMA_H}_singler_cells.csv.gz`; clinical workbooks `chum_tma_clinical_data_updated.xlsx` (TMA A) and `coeur_clinical_tmagh.xlsx` (TMA B) in `data/clinical data/` of the repository's parent folder | `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv`, `final_tma_patients.csv`, `final_tma_counts.csv` | `spatial/00_setup/final_tma_cores.R`, `figures/_prep/final_tma_cores.py` → TMA arms of Figs 4a, 4d, 4g, 4i, 5b, 5j, 6g, 7a–f, Figures S12–S15; Tables S8, S11 (Table S8 also checks `final_tma_counts.csv`) |

```bash
Rscript spatial/02_qc/04_qc_filter_normalize.R
mkdir -p $OUTPUT_ROOT/05_probe_qc
cp $DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/genes_exclude_singler.txt $OUTPUT_ROOT/05_probe_qc/
# ... stage 03 ...
for a in tma COEUR_TMA_G COEUR_TMA_H; do
  CORE_QC_LABELS=singler CORE_QC_TMAS=sfe_$a CORE_QC_OUT_TAG=_singler_$a \
    Rscript -e 'source("spatial/02_qc/05b_core_qc.R")'
done
python spatial/02_qc/06b_final_tma_cores.py
```

## Probe QC

The probe-QC notebook (step 2) compares each panel gene with the scRNA-seq reference: correlation of gene–gene
co-expression profiles, top-20-partner Jaccard overlap, mean and variance ratios, detection-rate and non-zero-mean
differences, and overdispersion. Genes with ≥ 3 flags fail (31 genes). A 4-tier rule, which weighs each failing gene
against how many reliable markers its cell type keeps, excludes 28 of them from SingleR annotation and keeps 3 markers
of cell types with fewer than 5 reliable markers. The final run did not re-run the notebook: it read the notebook's
outputs from
`$DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/`, and its copy of the exclusion list (step 3) is
md5-identical to the one there.

## Core QC and the final core set

`05b_core_qc.R` flags each tumour core of one array:
- fewer than 200 cells;
- less than 10% epithelial cells;
- low SingleR confidence (median score below the array's 5th percentile, or more than 10% pruned labels);
- composition outlier (Mahalanobis distance on the centred log-ratio of the 18-type composition, χ² p < 0.01);
- discordant pair (within-patient Jensen–Shannon divergence above the 90th percentile, or Spearman r < 0.5).

Cores with 2 or more flags are excluded. With `CORE_QC_LABELS=singler` the composition uses the SingleR labels split
at the frozen thresholds (as `06_extract_singler.R`).

`06b_final_tma_cores.py` then keeps tumour cores that pass core QC, from patients with primary HGSC on central review
(unreviewed histology kept) who were treatment-naive at sampling, and keeps a patient's cores only if they hold
≥ 500 cells. TMA B pools blocks G and H (a patient lies in one block; asserted). Counts by step
(`final_tma_counts.csv`):

| Step | TMA A patients / cores | TMA B patients / cores | Total patients / cores / cells |
|---|---|---|---|
| profiled tumour cores | 101 / 200 | 317 / 376 | 418 / 576 / 1,703,244 |
| after core QC | 98 / 183 | 307 / 355 | 405 / 538 / 1,626,163 |
| + primary HGSC, treatment-naive | 95 / 177 | 287 / 330 | 382 / 507 / 1,543,538 |
| **final** (≥ 500 cells per patient) | **95 / 177** | **283 / 326** | **378 / 503 / 1,542,448** |
| final patients with ≥ 2 final cores | 82 / 164 | 43 / 86 | 125 / 250 / 791,653 |

Panel-specific minimums (for example ≥ 500 secretory cells per core in Fig 4d) are applied on top of this set by each
panel script.

## Environment

| Variable | Final run | Step |
|---|---|---|
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` (step 1 rewrites the SFEs in this folder) | 1, 4 |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | 1, 4 |
| `CORE_QC_LABELS` | `singler` (default; any other value keeps the SFE `cell_label`) | 4 |
| `CORE_QC_TMAS` | one array per run: `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` (default `sfe_tma`; a comma-separated list pools arrays) | 4 |
| `CORE_QC_OUT_TAG` | `_singler_tma`, `_singler_COEUR_TMA_G`, `_singler_COEUR_TMA_H` | 4 |
| `OUTPUT_ROOT` | output root | 1–5 |

`05b_core_qc.R` finds `config/` relative to its own folder, which R only knows when the file is `source()`d. Run it as
above from the repository root, or with `Rscript 05b_core_qc.R` from inside `spatial/02_qc/` (then `OUTPUT_ROOT` must
be an absolute path). `Rscript spatial/02_qc/05b_core_qc.R` from the root stops at `config/config.R`.

`06b_final_tma_cores.py` reads the clinical workbooks from `data/clinical data/` in the folder that contains the
repository (`CLIN`); in the lab's layout that is `$DATA_ROOT/202609_final_hgsc_spatial/data/clinical data/`, the folder
that Figure 7 and Tables S8 and S11 read.

## `as_run/`

Trimmed copies of the scripts that ran in the lab project folder `2026_final_xenium_analysis`: the code that wrote the
outputs the paper uses is kept as it ran, and the notebook's plots, tables and summaries were removed. Each copy opens
with a note on what was removed; the untrimmed scripts stay in `2026_final_xenium_analysis`. The copies keep that
project's paths (the notebook sets `root.dir` to the project folder, and `00_setup.R` sets the project paths).

| File | Copy of (in `2026_final_xenium_analysis`) | Kept |
|---|---|---|
| `as_run/explore_probe_qc.Rmd` | `scripts/sandbox/explore_probe_qc.Rmd` | the chunks that compute the probe metrics and flags, the level-1 marker coverage, the 4-tier decisions, and write `probe_qc_full.csv`, `flag_qc.csv`, `gene_exclusion_decisions.csv`, `genes_exclude_singler.txt` and `genes_monitor_singler.txt` |
| `as_run/00_setup.R` | `scripts/00_setup.R` | libraries, project paths, output folders, `load_sfe()`, `log_session()`; sourced by this notebook and by `spatial/06_clinical/as_run/explore_clinical_associations.Rmd` |
