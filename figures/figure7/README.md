# Figure 7: polarization, survival and bevacizumab benefit

Secretory polarization and outcome in three cohorts. Panels a–d: Kaplan–Meier curves and multivariable Cox models in
the Xenium TMA patients represented by two tumour cores. Panels e,f: Xenium transcripts and polarization against
published immunofluorescence (IF) protein on TMA A cores. Panels g–j: TCGA-OV. Panels k–n: bevacizumab benefit by
polarization group in the ICON7 trial.

Each drawing script writes one combined row: a–d, e–f, g–j or k–n. Three output stems carry other letters than the
published panels: `xenium_tma_survival_row_ABC` is panels a–d, `tcga_survival_row_EFG` is panels g–j and
`icon7_survival_row_HIJ` is panels k–n.

## Requires

| Upstream step | Writes | Used by |
|---|---|---|
| Xenium backend: `spatial/01_build_sfe/` → `spatial/02_qc/04_qc_filter_normalize.R` → `spatial/03_annotation_polarization/01_annotation_singler.R` … `04_reclassification_polarization.R` → `06_extract_singler.R sfe_<array>` (TMA core and patient assignment: see `spatial/01_build_sfe/`). SingleR reads the probe-QC exclusion list `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (from `spatial/02_qc/as_run/explore_probe_qc.Rmd`, copied by hand; see [`spatial/02_qc/README.md`](../../spatial/02_qc/README.md)) | `$OUTPUT_ROOT/singler_recompute/{tma,COEUR_TMA_G,COEUR_TMA_H}_singler_cells.csv.gz`; `$OUTPUT_ROOT/sfe/sfe_tma` (logcounts) | 01, 03 |
| Final TMA core set: `spatial/02_qc/05b_core_qc.R` (once per array) → `spatial/02_qc/06b_final_tma_cores.py`; loader `spatial/00_setup/final_tma_cores.R` | `$OUTPUT_ROOT/07_core_qc/final_tma_{cores,patients,counts}.csv` (378 patients, 503 cores) | 01, 03 |
| TMA A clinical table: `spatial/06_clinical/as_run/clean_clinical_data.R` → `spatial/06_clinical/as_run/explore_clinical_associations.Rmd` | `$DATA_ROOT/2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv` (5-year OS and PFS, stage, residual disease, age) | 01 |
| TMA B clinical table: `spatial/06_clinical/as_run/coeur_clinical_extract.py` | `$DATA_ROOT/202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv` | 01 |
| Raw TMA clinical workbooks (clinical data; see `data/README.md`) | `$DATA_ROOT/202609_final_hgsc_spatial/data/clinical data/{chum_tma_clinical_data_updated.xlsx, coeur_clinical_tmagh.xlsx}` | 01 (neoadjuvant and histology exclusions) |
| TMA A immunofluorescence (published) | `$DATA_ROOT/2026_final_xenium_analysis/data/Copy of MFI all tissues TMA valid all markers_.xlsx`, sheet `Pour SPSS` | 03 |
| TCGA-OV scores: `atlas/07_deconvolution_survival/as_run/22a_cibersort_reference_v2.py` → `atlas/07_deconvolution_survival/as_run/22b_bayesprism_deconv.R` (epithelial fraction); `atlas/03_epithelial_nmf/as_run/11e_nmf_characterization.py` (177 SecA- and 177 SecB-specific genes) → `atlas/07_deconvolution_survival/as_run/22d_signature_survival.py` | deposited `$DATA_ROOT/2026_final_atlas/output/22_tcga_deconvolution/22d_signature_scores.csv` (score, tertile, epithelial fraction, clinical per tumour) | 05 |
| ICON7 cohort: `icon7/scripts/00_setup.R` … `04_descriptive.R` (see [`icon7/README.md`](../../icon7/README.md)) | `cohort_filtered.tsv` under `$DATA_ROOT` (config `objects.icon7_cohort`) | 06 |

The TCGA-OV steps are trimmed copies of the atlas project scripts (`202601_final_atlas`) that wrote the deposited
`22d_signature_scores.csv`: the code behind that table is kept as it ran, and the plots, analyses and outputs the paper
does not use were removed. The ICON7 copies are described in [`icon7/README.md`](../../icon7/README.md). Steps 05 and 06 read the
deposited tables.

## Run order

Outputs are under `$OUTPUT_ROOT/figures/figure7/` unless written otherwise.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 01 | `01_fig7abcd_tma_cox_km_tables.R` | SingleR cell tables of the three TMA arrays (final cores only); TMA A clinical table; TMA B clinical table; raw clinical workbooks | `tma_2core_patients.csv` (per-patient measures and clinical data; lab-internal identifiers), `cox_2core_results.csv` (6 models) | a–d |
| 02 | `02_fig7abcd_tma_survival_row.R` (+ `_km_forest_helpers.R`) | `tma_2core_patients.csv`, `cox_2core_results.csv` | `xenium_tma_survival_row_ABC.{svg,png}`, `xenium_tma_km_logrank.csv`, `xenium_tma_forest_fdr.csv` | a–d |
| 03 | `03_fig7ef_protein_correlation.R` | `sfe_tma` (logcounts), `$OUTPUT_ROOT/singler_recompute/tma_singler_cells.csv.gz`, final TMA cores, IF workbook | `$OUTPUT_ROOT/41_xenium_protein_correlation/per_core_xenium_protein_final.csv` (lab-internal core numbers), `correlations_final.csv` | e, f |
| 04 | `04_fig7ef_protein_correlation_plot.py` | `per_core_xenium_protein_final.csv` (`PROTEIN_CSV` overrides) | `xenium_protein_correlation_row.{svg,png}` | e, f |
| 05 | `05_fig7ghij_tcga_km_forest.R` (+ `_km_forest_helpers.R`) | `22d_signature_scores.csv` | `tcga_survival_row_EFG.{svg,png}`, `tcga_forest_stepwise_data.csv` | g–j |
| 06 | `06_fig7klmn_icon7_survival_row.R` (+ `_km_forest_helpers.R`) | `cohort_filtered.tsv` | `icon7_survival_row_HIJ.{svg,png}`, `icon7_km_pfs_within_group_hr.csv`, `icon7_timevarying_forest_15mo_data.csv` | k–n |

`_km_forest_helpers.R` holds the shared Kaplan–Meier, forest and row-layout code of 02, 05 and 06. Tables S8 and S11
check against the outputs of 01–02 (through the copies below), and Table S12 against those of 06, so run this folder
before `tables/build_all_tables.sh`.

**Copies for the tables.** Tables S8 and S11 read three outputs of 01–02 from `$OUTPUT_ROOT/final_panels/figure7/`
under other names. Copy them there after step 02 (the files are identical copies):

```bash
mkdir -p $OUTPUT_ROOT/final_panels/figure7
F=$OUTPUT_ROOT/figures/figure7; P=$OUTPUT_ROOT/final_panels/figure7
cp $F/cox_2core_results.csv     $P/Fig7_ABCD_cox_2core_results.csv   # Tables S8, S11
cp $F/xenium_tma_forest_fdr.csv $P/Fig7_ABCD_forest_fdr.csv          # Table S11
cp $F/xenium_tma_km_logrank.csv $P/Fig7_ABCD_km_logrank.csv          # Table S11
```

## Environment

| Setting | Step | Value in the final run |
|---|---|---|
| Working directory | all R steps | Repository root (01 and 03 `source()` helpers by repo-relative path) |
| `LANG`, `LC_ALL` | all R steps | `en_US.UTF-8` |
| `OUTPUT_ROOT` | all | Output root (default `./output`). 01 writes to `$OUTPUT_ROOT/figures/figure7/`, the folder 02, 05 and 06 use (config `figures_dir`); 01 does not create it, so run `mkdir -p $OUTPUT_ROOT/figures/figure7` before 01 |
| `DATA_ROOT` | 01, 03, 05, 06 | Deposit root: the clinical tables and workbooks (01), the IF workbook (03), `22d_signature_scores.csv` (05), `cohort_filtered.tsv` (06) |
| `FINAL_TMA_CORES` | 01, 03 | default, `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` |
| `RERUN_SFE_DIR` | 03 | `$OUTPUT_ROOT/sfe`. `load_sfe("sfe_tma")` looks in `$OUTPUT_ROOT/sfe` first; with `RERUN_SFE_DIR` unset, a missing object is loaded without a message from the earlier SFEs under `$DATA_ROOT` (config `objects.sfe_dir`) |
| `PROTEIN_CSV` | 04 | default, `$OUTPUT_ROOT/41_xenium_protein_correlation/per_core_xenium_protein_final.csv` |
| Interpreters | all | R 4.4.3 (E1), Python 3.12.9 (E2); see `env/README.md` |

## Panel → script

| Panel | Drawn by | Numbers computed by |
|---|---|---|
| 7a | `02` | `01`: cohort = patients with primary, treatment-naive HGSC and ≥ 2 final tumour cores; patients given neoadjuvant chemotherapy or with non-HGSC central histology are removed; complete stage, residual disease and age; follow-up censored at 5 years. `02`: lowest vs highest tertile of per-patient mean polarization, 5-year OS, log-rank (n = 104; 35 vs 35; `xenium_tma_km_logrank.csv`) |
| 7b | `02` | As 7a, 5-year PFS (n = 103; 35 vs 34) |
| 7c | `02` | `01`: Cox `Surv ~ z(measure) + stage + residual disease + age` for SecA proportion, SecB proportion and mean polarization (UCell); HR per 1 SD; Benjamini–Hochberg q across the 6 tests of 7c,d (`cox_2core_results.csv`; `xenium_tma_forest_fdr.csv`). Same models as Table S11 |
| 7d | `02` | As 7c, PFS |
| 7e | `04` | `03`: per-core mean log counts of KRT7 and KRT19 over all cells of the core vs IF MFI, final TMA A cores with IF (n = 157); Pearson r in `correlations_final.csv` (`04` recomputes r and ρ for the plot) |
| 7f | `04` | `03`: per-core mean polarization of the SingleR secretory cells vs IF MFI of KRT7, KRT19, KRT18, VIM and E-cadherin; Pearson r and Spearman ρ |
| 7g | `05` | Tertiles of the TCGA polarization score (`polar_tertile` in `22d_signature_scores.csv`); highest vs lowest, 5-year OS, log-rank (SecA n = 104, SecB n = 106). The log-rank p is drawn on the panel only; no CSV holds it |
| 7h | `05` | As 7g, 5-year PFS (SecA n = 104, SecB n = 105) |
| 7i | `05` | Cox on the polarization score in the full cohort: unadjusted, + epithelial fraction (BayesPrism), + stage and age (OS n = 314; 311 with stage and age); `tcga_forest_stepwise_data.csv` |
| 7j | `05` | As 7i, PFS (n = 312; 309) |
| 7k | `06` | FIGO III/IV patients with complete PFS, polarization, stage, debulking status and age (n = 191), split at the median polarization. SecA-polarized (n = 96): PFS by arm, log-rank; bevacizumab HR adjusted for FIGO stage, debulking status and age (`icon7_km_pfs_within_group_hr.csv`). Also Table S12, model B |
| 7l | `06` | As 7k, SecB-polarized (n = 95) |
| 7m | `06` | SecB-polarized: bevacizumab HR for PFS before and after 15 months (split at 456 days), adjusted for FIGO stage, debulking status and age; Benjamini–Hochberg q across the 4 estimates of 7m,n; period interaction (`icon7_timevarying_forest_15mo_data.csv`). Also Table S12, model C |
| 7n | `06` | As 7m, SecA-polarized |

### Numbers in the text

- TMA survival (OS n = 104, 40 events; PFS n = 103, 71 events; log-rank p = 0.011 and 0.11; the six HRs with 95% CI
  and q): `01` → `02` (`xenium_tma_km_logrank.csv`, `xenium_tma_forest_fdr.csv`); also Table S11.
- Methods "Survival", 5 neoadjuvant and 4 non-HGSC patients excluded: these patients are removed upstream by
  `spatial/02_qc/06b_final_tma_cores.py`, so `01` logs no exclusions. The counts are reproduced in Table S8, sheet
  `TMA_patients`.
- KRT7 r = 0.64 and KRT19 r = 0.63; Methods "177 final TMA A cores from 95 patients, 157–163 with IF for each
  marker": `03` (`correlations_final.csv`).
- TCGA-OV (314 tumours; log-rank p = 0.010 for OS and PFS; OS HR 2.00 and PFS HR 1.79 unadjusted; OS HR 1.84 and PFS HR
  1.65 after adjustment for epithelial fraction, stage and age): `05` (`tcga_forest_stepwise_data.csv`; the log-rank p
  is on the panel only). `22d_signature_scores.csv` scores 323 tumours; 314 is the OS-evaluable cohort.
- ICON7 (n = 191; SecA 96, SecB 95; bevacizumab PFS HR 0.36 in SecB-polarized and 0.80 in SecA-polarized tumours,
  with KM medians; time-varying HRs 0.23 and 0.87, period interaction p = 0.022 and 0.19): `06`; Table S12 models B
  and C, sheet `KM_medians`.
- ICON7 treatment × polarization interaction (HR 0.44, likelihood-ratio p = 0.025): Table S12 model D
  (`tables/12_tableS12_icon7_cox.py`).
- ICON7 standard-arm OS HR 1.38 per 1 SD: Table S12 model A, the per-arm continuous-polarization Cox
  (`icon7/scripts/15_chemo_arm_analysis.R` section B.2, rebuilt in `tables/12_tableS12_icon7_cox.py`).

## Not in this folder

- The TMA B (COEUR) clinical table `coeur_clinical_perpatient_extracted.csv` (per-patient OS and PFS, events, FIGO
  stage, debulking and age, extracted from `coeur_clinical_tmagh.xlsx`) is written by
  `spatial/06_clinical/as_run/coeur_clinical_extract.py` and read by 01 from
  `$DATA_ROOT/202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/`. Table S8 re-derives the same
  fields from the workbook and asserts that they are equal.
- **[code missing]** Download and preparation of the TCGA-OV inputs, `tcga_ecotyper.txt` (expression) and
  `tcga_hla_clinical.csv` (clinical), in the atlas project folder `202601_final_atlas` (`data/cibersort_data_prev/`).
- Upstream of 7a–f: the TMA A and TMA B de-array code is in `spatial/01_build_sfe/as_run/`; the TMA A core map it
  read is **[not on disk]** (see [`spatial/README.md`](../../spatial/README.md)).
- **TMA A clinical table.** `spatial/06_clinical/as_run/explore_clinical_associations.Rmd` is a trimmed copy of the
  notebook that ran. It reads an intermediate neighbourhood table that is no longer on disk, so it does not re-run
  unchanged. Table S8 checks the table's OS months, age and residual disease against the raw workbook.
- **Manual:** the copies of three 01–02 outputs into `$OUTPUT_ROOT/final_panels/figure7/` (see Run order) and the
  probe-QC exclusion list (see Requires).
- **Assembly.** Panel letters, and the OS / PFS labels over the Kaplan–Meier panels, were added in Affinity Designer
  when the figure was assembled from the four rows.
