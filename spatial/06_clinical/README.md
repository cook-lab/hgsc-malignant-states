# `spatial/06_clinical/`: TMA clinical survival tables

Produces the per-patient clinical tables (overall and progression-free survival censored at 5 years, FIGO stage,
residual disease, age) that the TMA survival analyses join to the Xenium features: Figure 7a–d and Tables S8 and S11.
The TMA A tables come from two `as_run/` copies from the lab project folder `2026_final_xenium_analysis`; the TMA B
extract comes from `as_run/coeur_clinical_extract.py`.

The treatment and histology rules of the final TMA core set are applied separately, from the clinical workbooks, by
`spatial/02_qc/06b_final_tma_cores.py`.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `as_run/clean_clinical_data.R` | `data/clinical_data_updated.xlsx` (TMA A, in `2026_final_xenium_analysis`) | `data/clinical_data_clean.csv` (same project): empty rows dropped; treatment status, CA-125, residual disease and chemotherapy response recoded | step 2 |
| 2 | `as_run/explore_clinical_associations.Rmd` | that project's filtered TMA A SFE (`sfe_tma_filtered`, from before the final re-run); `output/09_neighborhood/neighborhood_assignments_k10.csv`; `data/clinical_data_clean.csv` | `output/10_clinical/per_patient_features.csv`: per-patient cell-type densities and neighbourhood proportions, the cleaned clinical columns, and OS and PFS censored at 60 months (`os_time_5y`, `os_event_5y`, `pfs_time_5y`, `pfs_event_5y`) | `figures/figure7/01_fig7abcd_tma_cox_km_tables.R`; `tables/08_tableS8_xenium_cohort.py`; `tables/11_tableS11_xenium_tma_cox.py` |
| 3 | `as_run/coeur_clinical_extract.py <out_dir>`, then copied to `exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv` | `coeur_clinical_tmagh.xlsx` (sheet "clinical data") | `coeur_clinical_perpatient_extracted.csv` (`patient`, `os_time`, `os_event`, `pfs_time`, `pfs_event`, `figo_code`, `debulking`, `age`) | same readers as step 2 |

The readers load:
- `per_patient_features.csv` from `$DATA_ROOT/2026_final_xenium_analysis/output/10_clinical/`, and use its columns
  `patient_id`, `os_time_5y`, `os_event_5y`, `pfs_time_5y`, `pfs_event_5y`, `stage`, `residual_binary` and `age`;
- the TMA B extract from `$DATA_ROOT/202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/`;
- the clinical workbooks `chum_tma_clinical_data_updated.xlsx` (TMA A) and `coeur_clinical_tmagh.xlsx` (TMA B) from
  `$DATA_ROOT/202609_final_hgsc_spatial/data/clinical data/`, for the treatment and histology exclusions.

These tables carry the laboratory's patient identifiers; for access to the clinical tables see `data/README.md`
("Clinical").

Checks in the repository:
- `tables/08_tableS8_xenium_cohort.py` re-reads the TMA A workbook and asserts that the OS months, age and residual
  disease of `per_patient_features.csv` match it.
- The same script re-derives the TMA B fields from `coeur_clinical_tmagh.xlsx` and asserts that they equal the TMA B
  extract.

## `as_run/`

| File | Copy of (in `2026_final_xenium_analysis`) | Kept |
|---|---|---|
| `as_run/clean_clinical_data.R` | `scripts/sandbox/clean_clinical_data.R` | all the code that writes `data/clinical_data_clean.csv` |
| `as_run/explore_clinical_associations.Rmd` | `scripts/sandbox/explore_clinical_associations.Rmd` | sections 1–2, which write `output/10_clinical/per_patient_features.csv`, and the session info |
| `as_run/coeur_clinical_extract.py` | no file: it ran inline in the analysis session on 2026-09-15 and is recovered from the session log | all of it (unchanged) |

Both are trimmed copies of the scripts that ran: the code that wrote the outputs the paper uses is kept as it ran, and
the notebook's exploratory sections (clinical summary, Kaplan–Meier, Cox, ratio, forest plots, clustering) and the
cleaning script's closing console summary were removed. Each copy opens with a note on what was removed; the untrimmed
scripts stay in `2026_final_xenium_analysis`.

Both keep their original paths: `clean_clinical_data.R` reads and writes relative to the project folder, and the
notebook sets `root.dir` to the project folder and sources that project's `scripts/00_setup.R` (copied as
`spatial/02_qc/as_run/00_setup.R`).

## Missing code

- None for the TMA B extract: `as_run/coeur_clinical_extract.py` (recovered from the analysis session log; it ran
  inline on 2026-09-15) writes a file byte-identical to `coeur_clinical_perpatient_extracted.csv` (checked 2026-10-22).
  `tables/08_tableS8_xenium_cohort.py` also reproduces its fields from the workbook.
- `explore_clinical_associations.Rmd` cannot be re-run as is: it reads
  `output/09_neighborhood/neighborhood_assignments_k10.csv` of `2026_final_xenium_analysis`, which no longer exists.
  Its patients are those present in all three inputs: the tumour cells of `sfe_tma_filtered`, that neighbourhood table
  and the cleaned clinical table.
