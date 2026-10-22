# icon7: ICON7 trial transcriptomes (GSE140082)

This folder turns the public ICON7 microarray data (GEO GSE140082, 380 arrays) into one table of per-patient SecA,
SecB and polarization scores with clinical data, `cohort_filtered.tsv`. That table is the entry object of Figure 7k–n
and Table S12 (config `objects.icon7_cohort`, under `$DATA_ROOT`).

## As-run copies

`scripts/` holds copies of the scripts that ran in the ICON7 folder of the Xenium project folder
`2026_final_xenium_analysis` (`davids side quests/ICON7/scripts/`), on a second machine with R 4.5.2 (environment E5
in `env/README.md`; lock file `env/renv_R4.5.2_bioc3.21_previous.lock`; the UCell version used was not recorded).
`02`, `04` and `15` are trimmed copies: the code that wrote the outputs the paper uses is kept as it ran, and the
tables and plots the paper does not use were removed. Each trimmed copy opens with a note on what was removed. `00`,
`01` and `03` are unchanged. The untrimmed scripts stay in the lab project folder. In the copies:

- Each script except `00_setup.R` sets `ROOT` at the top to the ICON7 folder on that machine. Set `ROOT` to your own
  ICON7 folder; the scripts read from and write to `data/raw/`, `data/processed/` and `results/` below it.
- The header of `00_setup.R` lists scripts `05`–`07`. Only the scripts listed below fed the manuscript; the other
  scripts of the module are not included.
- The deposited ICON7 folder is `$DATA_ROOT/2026_final_xenium_analysis/davids side quests/ICON7/`: the config keys
  `objects.icon7_cohort`, `objects.icon7_results_dir` and `objects.icon7_expr` point into it.

## Run order

Paths are relative to `ROOT`.

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 00 | `scripts/00_setup.R` | – | installs and checks packages (CRAN: data.table, survival, survminer, ggplot2, dplyr, tidyr, patchwork, readr, stringr; Bioconductor: GEOquery, illuminaHumanv4.db, UCell, AnnotationDbi, Biobase) | all |
| 01 | `scripts/01_parse_metadata.R` | `data/raw/GSE140082_series_matrix.txt.gz` | `data/processed/sample_metadata.tsv` (380 arrays: treatment arm, clinical characteristics, survival) | 02, 03 |
| 02 | `scripts/02_build_gene_matrix.R` | `data/raw/GSE140082_geo.normdata.gz`, `sample_metadata.tsv` | `data/processed/expr_gene_symbol.rds` (gene × array, log2 normalized; detection p-value columns dropped; the probe with the highest mean kept per gene), `data/processed/probe_annotation.tsv` | 03; Table S6 (gene presence; it reads the deposited file at the path of config `objects.icon7_expr`) |
| 03 | `scripts/03_score_signatures.R` | `expr_gene_symbol.rds`, `sample_metadata.tsv` | `data/processed/sample_metadata_scored.tsv`: UCell scores of the 7 SecA and 7 SecB genes, each z-scored across the 380 arrays; `polarization_ucell` = z(SecB) − z(SecA). Also a z-scored mean-expression version, used only as a check | 04; Table S6 parses its gene lists from the copy in the deposited ICON7 folder (`scripts/03_score_signatures.R`) |
| 04 | `scripts/04_descriptive.R` | `sample_metadata_scored.tsv` | `data/processed/cohort_filtered.tsv` (serous histology, high grade and `manuscript_analysis359 == 1`) | `figures/figure7/06_fig7klmn_icon7_survival_row.R`, `tables/12_tableS12_icon7_cox.py`, step 15 |
| 15 | `scripts/15_chemo_arm_analysis.R` | `cohort_filtered.tsv` | `results/tables/17_chemo_arm_continuous_polarization.tsv`, `results/tables/17b_continuous_polarization_per_arm.tsv` (section B.2, the only section kept) | `tables/12_tableS12_icon7_cox.py` (check only) |

Section B.2 of `15_chemo_arm_analysis.R` is the per-arm continuous-polarization Cox model, Table S12 model A.
`tables/12_tableS12_icon7_cox.py` rebuilds model A itself and checks it against the two B.2 output tables, which it
reads from config `objects.icon7_results_dir`.

## Downstream filters

Step 04 keeps all patients with high-grade serous histology (`cohort_filtered.tsv`). Step 15 and the figure and
table scripts then keep FIGO stage III/IV patients with complete survival, polarization, FIGO stage, debulking status
and age (n = 191). The figure and table scripts split them at the median polarization (SecA-polarized n = 96,
SecB-polarized n = 95). Every model adjusts for FIGO stage, debulking status and age.

## Inputs not produced by code

- `data/raw/GSE140082_series_matrix.txt.gz` and `data/raw/GSE140082_geo.normdata.gz` were downloaded from GEO by hand.
  No download script is included.
