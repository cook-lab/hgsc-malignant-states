# `spatial/`: Xenium backend

These stages take the raw Xenium output to the per-cell tables that every Xenium panel reads: Figures 4–7, Figures
S11–S18 and Tables S7–S11. Run them in the order below. Each stage README lists its scripts with inputs and outputs.

## Cohort

| Arm | Samples | Cells |
|---|---|---|
| Whole tissues | 8 sections, WT01–WT08 (4 treatment-naive, 4 post-chemotherapy) | 1,898,772 |
| TMA A (CHUM cohort) | 101 patients / 200 tumour cores profiled; two slides merged into one object (`sfe_tma`) | |
| TMA B (COEUR cohort) | 317 patients / 376 tumour cores profiled; blocks G and H (`sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H`) | |
| **Final TMA core set** | **378 patients, 503 cores** (TMA A 95 / 177; TMA B 283 / 326) | **1,542,448** |
| Fallopian tube epithelium | 15 FTE cores on the TMA A slides | 34,281 |

The panel has 477 genes. Every TMA analysis uses the final core set (`02_qc/06b_final_tma_cores.py`). The backend
steps loop over 11 SFEs: `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` and `sfe_WT01`–`sfe_WT08` (`sfe_names_all`
in `00_setup/00_setup.R`, from `config.yml` `cohort`).

## Run order

```
raw Xenium output (one folder per run)
 01_build_sfe        SFE per sample: 8 WT, TMA A, TMA B blocks G and H
 02_qc   04          cell QC and normalisation (in place)
         probe QC    28 genes excluded from annotation (as_run notebook; its outputs are read from $DATA_ROOT,
                     and the exclusion list is copied by hand to $OUTPUT_ROOT/05_probe_qc/)
 03_annotation_polarization
         01  SingleR (429 genes, 16 types)
         02  initial secretory subtypes             (structural only)
         03  UCell SecA / SecB / polarization (457 genes shared with the atlas)
         04  frozen thresholds: SecA <= 0.109429, SecB >= 0.210501
         05  RCTD/SPLIT                             (structural only)
         06  SingleR label tables   -> $OUTPUT_ROOT/singler_recompute/
 02_qc   05b, 06b    per-array core QC -> final TMA core set (needs stage 03)
 04_pathway_scoring  37 UCell pathway modules -> SFE pathway_* (all cells)
 05_count_export     raw-count export for Figure 6f–i and Tables S9, S10 (as_run)
 06_clinical         TMA A survival table (as_run)
```

| Stage | What it does |
|---|---|
| [`00_setup/`](00_setup/README.md) | Shared R setup (`load_sfe()`, palettes, signatures) and the two helpers the figure scripts source: the final core set and the SingleR labels |
| [`01_build_sfe/`](01_build_sfe/README.md) | Builds one SpatialFeatureExperiment (SFE) per sample from the raw Xenium folders; attaches TMA core and patient labels |
| [`02_qc/`](02_qc/README.md) | Cell QC and normalisation; probe QC; per-array core QC; the final TMA core set |
| [`03_annotation_polarization/`](03_annotation_polarization/README.md) | SingleR cell types, UCell polarization, frozen thresholds, the SingleR label tables |
| [`04_pathway_scoring/`](04_pathway_scoring/README.md) | The 37 UCell pathway modules (Table S9), scored on every cell |
| [`05_count_export/`](05_count_export/README.md) | Per-sample raw-count export read by the Figure 6 Xenium analyses and Tables S9, S10 |
| [`06_clinical/`](06_clinical/README.md) | TMA A clinical and survival table read by Figure 7a–d and Tables S8, S11 |

Stage 02 runs in two parts: cell QC (`04`) before annotation, core QC (`05b`, `06b`) after it, because core QC needs
the SingleR labels and polarization scores.

## Cell labels: use the SingleR label tables

Figure and table scripts take cell identity from the **SingleR label tables**,
`$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` (`03_annotation_polarization/06_extract_singler.R`), not
from `colData(sfe)$cell_label`.

- `cell_label_singler` is the SingleR label, with "Secretory epithelium" split on `polarization_UCell`
  (SecB – SecA UCell score): SecA at ≤ 0.109429, SecB at ≥ 0.210501, Intermediate in between. This gives the 18 cell
  types of the manuscript.
- The thresholds are the 75th percentile of atlas SecA cells and the 25th percentile of atlas SecB cells on the same
  score (`03_annotation_polarization/04_reclassification_polarization.R`, `threshold_summary.csv`). They are frozen and
  typed into `06_extract_singler.R`, `02_qc/05b_core_qc.R` and the count export.
- The SFE column `cell_label` is last written by the RCTD/SPLIT step (`03_annotation_polarization/05_clean_split_rctd.R`).
  No panel takes cell identity from it. It stays in the SFEs because `03_annotation_polarization/04`, `02_qc/05b` and
  the count export read it; the count export takes its cell set from it.
- R scripts read the tables through `00_setup/singler_labels.R` (`singler_labels()`); Python scripts read the CSVs
  directly. There are 11 tables in use: WT01–WT08, `tma`, `COEUR_TMA_G`, `COEUR_TMA_H`.

## Main outputs

| Output | Path | Written by | Read by |
|---|---|---|---|
| SFEs (HDF5) | `$OUTPUT_ROOT/sfe/sfe_<sample>` | `01_build_sfe`; filtered and normalised by `02_qc/04`; columns added by `03_annotation_polarization/01`–`05` and `04_pathway_scoring/04c` | every Xenium step |
| SingleR label tables | `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` | `03_annotation_polarization/06_extract_singler.R` | every Xenium figure and table script |
| Polarization thresholds | `$OUTPUT_ROOT/06f_reclassification_polarization/threshold_summary.csv` | `03_annotation_polarization/04` | `03_annotation_polarization/05`; the values are typed into later scripts |
| Final TMA core set | `$OUTPUT_ROOT/07_core_qc/final_tma_{cores,patients,counts}.csv` | `02_qc/06b` | `00_setup/final_tma_cores.R`, `figures/_prep/final_tma_cores.py`, Table S8 |
| Pathway gene sets and scores | `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv`; SFE `pathway_*` | `04_pathway_scoring` | Figs 4f, 4j and 5, Figures S16–S18 (SFE scores); Fig 6f,h (a hand copy of the gene-set file); Tables S7, S9 |
| Raw-count export | `$OUTPUT_ROOT/fig6_viewer/degdata/` | `05_count_export/as_run` | Figure 6f–i; Tables S9, S10 |
| TMA A clinical table | `per_patient_features.csv`, read from `$DATA_ROOT/2026_final_xenium_analysis/output/10_clinical/` | `06_clinical/as_run` | Figure 7a–d; Tables S8, S11 |

On the laboratory's data, the final core tables (`final_tma_*.csv`) and the label tables carry raw core and patient
identifiers. Published tables use study codes (WT01–WT08, TMA-P###).

## Environment

Run every script from the repository root (R scripts `source()` helpers by repo-relative path). Two scripts locate
`config/` from their own folder instead and must be `source()`d: `02_qc/05b_core_qc.R` and
`04_pathway_scoring/04_pathway_scoring.R` (`Rscript -e 'source("spatial/02_qc/05b_core_qc.R")'`).

| Variable | Value in the final run (default) | Needed by |
|---|---|---|
| `DATA_ROOT`, `OUTPUT_ROOT` | deposit root; output root (`./output` if unset) | everything |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | every R script (otherwise `config.yml` is read incorrectly) |
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` | `01_build_sfe` (where SFEs are written), `02_qc/04`, `03_annotation_polarization/01` and `load_sfe()` (`sfe_dir`) |
| `RAW_XENIUM_DIR` | `$DATA_ROOT/202609_final_hgsc_spatial/data/HGSC xenium ` (the folder name ends with a space) | `01_build_sfe` |
| `DEARRAY_DIR` | `$DATA_ROOT/202609_final_hgsc_spatial/data/dearray` | `01_build_sfe/03_build_sfe_new_coeur.R` |
| `SINGLER_WORKERS` | `12` | `03_annotation_polarization/01` |
| `ATLAS_UCELL_DIR` | unset: the deposited atlas export `$DATA_ROOT/2026_final_atlas/output/18_ucell_atlas/` | `03_annotation_polarization/03` |
| `RERUN_WRITEBACK_ONLY` | unset | `03_annotation_polarization/05` |
| `CORE_QC_LABELS`, `CORE_QC_TMAS`, `CORE_QC_OUT_TAG` | `singler`; one array per run; `_singler_<array>` | `02_qc/05b` |
| `FINAL_TMA_CORES` | `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` | both final-core loaders |

**Where SFEs are read and written.** `load_sfe()` reads `$OUTPUT_ROOT/sfe/<name>` first, then `sfe_dir`. `sfe_dir`
is `RERUN_SFE_DIR`, or, if that is unset, `$DATA_ROOT/2026_final_xenium_analysis/output/sfe` (config
`objects.sfe_dir`), which on the lab's data holds SFEs from before the final re-run. `save_sfe()` writes only under
`$OUTPUT_ROOT/sfe`, but `02_qc/04_qc_filter_normalize.R` rewrites each SFE in place in `sfe_dir`. Set
`RERUN_SFE_DIR=$OUTPUT_ROOT/sfe` for the whole backend.

The `as_run/` copies keep the paths of the folders they ran in; each stage README says what to change.

## Missing code

- TMA A core → patient labels. `01_build_sfe/05_rebuild_existing_cohort_rerun.R` copies core, patient and
  sample-type labels from an earlier processed object,
  `$DATA_ROOT/2026_final_xenium_analysis/data/previously processed/xenium_tma_final_sfe_hdf5/se.rds`. The 2025 code
  that made them is in `01_build_sfe/as_run/tmaA_*`; its input core map `TMA_xenium_map.xlsx` is **[not on disk]**.
- TMA B de-array (blocks G and H): the core circles were set by hand in `01_build_sfe/as_run/dearray_editor.html`;
  the cell → core → patient tables are written by `01_build_sfe/as_run/dearray_cell_core_assignment.py`.
- The TMA B clinical extract `coeur_clinical_perpatient_extracted.csv` is written by
  `06_clinical/as_run/coeur_clinical_extract.py`.
- Manual step: the copy of the probe-QC exclusion list `genes_exclude_singler.txt` into `$OUTPUT_ROOT/05_probe_qc/`
  (see `02_qc/`).
- `06_clinical/as_run/explore_clinical_associations.Rmd` reads a neighbourhood table that no longer exists, so it
  cannot be re-run as is; the readers use its output table.
