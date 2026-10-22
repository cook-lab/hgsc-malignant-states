# `spatial/01_build_sfe/`: SpatialFeatureExperiment build

Reads the raw Xenium output folders and writes one HDF5-backed SpatialFeatureExperiment (SFE) per sample, with raw
counts, cell and nucleus segmentation and, for the TMAs, core, patient and sample-type labels. All cells are kept;
filtering and normalisation happen in `spatial/02_qc/04_qc_filter_normalize.R`.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `05_rebuild_existing_cohort_rerun.R` | `$RAW_XENIUM_DIR/<run>/` (6 single-tissue runs, the dual-tissue run `WT05_WT01`, TMA A slides `CHUM_TMA_1` and `CHUM_TMA_2`); TMA A labels from `$DATA_ROOT/2026_final_xenium_analysis/data/previously processed/xenium_tma_final_sfe_hdf5/se.rds` | `$RERUN_SFE_DIR/sfe_WT01` … `sfe_WT08`, `sfe_tma` | `spatial/02_qc/04_qc_filter_normalize.R` |
| – | `as_run/dearray_cell_core_assignment.py` (after the core circles were set in `as_run/dearray_editor.html`) | `$RAW_XENIUM_DIR/COEUR_TMA_G/cells.parquet`, `COEUR_TMA_H/cells.parquet`; the editor's core circles; `data/clinical data/coeur_clinical_tmagh.xlsx` (rank punch 1 / 2 → COEUR patient) | `$DEARRAY_DIR/{G,H}_cell_core_assignment.parquet`, `{G,H}_core_summary.csv` | step 2 |
| 2 | `03_build_sfe_new_coeur.R` | `$RAW_XENIUM_DIR/COEUR_TMA_G/`, `COEUR_TMA_H/`; `$DEARRAY_DIR/{G,H}_cell_core_assignment.parquet` | `$RERUN_SFE_DIR/sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` | `spatial/02_qc/04_qc_filter_normalize.R` |

Both scripts are heavy (each reads full Xenium runs). A sample whose output folder already exists is skipped.

```bash
Rscript spatial/01_build_sfe/05_rebuild_existing_cohort_rerun.R
Rscript spatial/01_build_sfe/03_build_sfe_new_coeur.R
```

Details:
- **Whole tissues.** WT05 and WT01 share one slide. Step 1 splits it with DBSCAN on cell coordinates (eps 200 µm,
  minPts 50) and assigns the upper cluster to WT05 and the lower to WT01.
- **TMA A.** The two slides are merged into `sfe_tma`, the second shifted in x so the slides do not overlap. Core,
  patient and sample-type labels (tumour or fallopian tube) are copied from the earlier processed object by cell id;
  cells without a match take the labels of their nearest matched cell.
- **TMA B.** Blocks G and H are separate objects. Core and patient labels are joined by cell id from the de-array
  tables; cells on a core get `core_flag` "On core" and `sample_type` tumour, all other cells "Off core".
- Morphology images are not read (the reader cannot open the newer multi-channel image folder); the scripts stage each
  run without them.

## Environment

| Variable | Default | Notes |
|---|---|---|
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` (`./output/sfe` if `OUTPUT_ROOT` is unset) | the final run used `$OUTPUT_ROOT/sfe` |
| `RAW_XENIUM_DIR` | `$DATA_ROOT/202609_final_hgsc_spatial/data/HGSC xenium ` (the folder name ends with a space) | one subfolder per run, named as above |
| `DEARRAY_DIR` | `$DATA_ROOT/202609_final_hgsc_spatial/data/dearray` | holds the TMA B de-array tables |
| `DATA_ROOT` | deposit root | also gives the TMA A labels object |

## `as_run/`

`as_run/dearray_editor.html` is a copy of the interactive page used to review and correct the TMA B core grid of
blocks G and H. It was recovered from the scratch area of the workspace project folder `202609_final_hgsc_spatial`.
The page loads an automatic core-detection seed per block and saves edited cores to the database of its hosting page;
opened on its own it shows the layout but saves nothing. The only change from the page that ran is the example text of
one input field, which named a real patient number.

`as_run/dearray_cell_core_assignment.py` assigns each cell to the core circle it falls in (nearest of the 3 closest
centres) and each core to its COEUR patient through the core rank. It ran inline in the analysis session (2026-09-13)
and is recovered from the session log, without a third TMA that is not in the paper. It reads the final circles and the
block offset from the editor's files in a temporary session folder, which is gone. The same circles (centre, radius,
rank) are in `{G,H}_core_summary.csv`: rebuilt from those, with a zero offset, the script reproduces both
`*_cell_core_assignment.parquet` tables exactly (checked 2026-10-22, every cell).

`as_run/tmaA_01_core_id_into_sfe.R`, `tmaA_02_core_flag_sample_type.R` and `tmaA_03_shift_and_merge.R` are the
2025 scripts (run interactively in the home folder) that made the TMA A labels of the earlier processed object `se.rds`:
core and patient from the core boxes of the Xenium core map, the QC flag and sample type (cores 203–217 fallopian
tube), and the merge of the two slides with slide-2 cell ids prefixed `tma2_`. The flag and sample-type rules and the
prefix hold for every cell of `se.rds`, with one patient per core (checked 2026-10-22).

The TMA B core circles themselves were placed by hand in the editor. They started from an automatic detection and rank
labelling run as a series of inline helper steps in the same session; those helpers are not in this repository.

## Missing code

- TMA A labels: step 1 inherits them from `previously processed/xenium_tma_final_sfe_hdf5/se.rds` (read from
  `$DATA_ROOT/2026_final_xenium_analysis/data/`). The code that made them is in `as_run/` (`tmaA_01`–`tmaA_03`, see
  above). **[not on disk]** its input, the Xenium core map `TMA_xenium_map.xlsx` (per core: box coordinates, core and
  patient), so the labels are read from `se.rds` rather than recomputed.
- TMA B de-array: nothing is missing. The cell → core → patient step is `as_run/dearray_cell_core_assignment.py`; the
  core circles are a manual step (see `as_run/` above).
