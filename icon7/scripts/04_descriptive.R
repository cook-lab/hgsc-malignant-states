# Trimmed copy (2026-10-21) of 2026_final_xenium_analysis/davids side quests/ICON7/scripts/04_descriptive.R, the
# script that ran. Kept: the HGSC filter, which writes cohort_filtered.tsv (Figure 7k-n, Table S12, script 15).
# Removed: the cohort summary table, the subtype score table and tests, and the subtype boxplot.
# The untrimmed script is in the lab project folder.
# ============================================================================
# 04_descriptive.R
# Cohort filtering + descriptive panel:
#   - Apply HGSC filter: histology.serous == "serous" & newgrade == "high.grade"
#     & manuscript_analysis359 == 1
#
# Outputs:
#   data/processed/cohort_filtered.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

ROOT <- "/Users/dpcook/My Drive/Lab/Manuscripts/hgsc_epitypes/supplemental_analysis/ICON7"
IN_FILE <- file.path(ROOT, "data/processed/sample_metadata_scored.tsv")
OUT_COH <- file.path(ROOT, "data/processed/cohort_filtered.tsv")
OUT_T1  <- file.path(ROOT, "results/tables/02_cohort_table1.tsv")
OUT_ST  <- file.path(ROOT, "results/tables/03_subtype_scores.tsv")
OUT_FIG <- file.path(ROOT, "results/figs/descriptive_subtype_scores.pdf")

scored <- fread(IN_FILE)

message(sprintf("Loaded %d scored samples", nrow(scored)))
message(sprintf("histology.serous values: %s",
                paste(sort(unique(scored$histology.serous)), collapse = " | ")))
message(sprintf("newgrade values: %s",
                paste(sort(unique(scored$newgrade)), collapse = " | ")))
message(sprintf("manuscript_analysis359 values: %s",
                paste(sort(unique(scored$manuscript_analysis359)), collapse = " | ")))
message(sprintf("treatment values: %s",
                paste(sort(unique(scored$treatment)), collapse = " | ")))

# ---- HGSC filter ---------------------------------------------------------
cohort <- scored[
  histology.serous == "serous" &
    newgrade == "high.grade" &
    manuscript_analysis359 == 1
]
message(sprintf("After HGSC filter: n=%d (expected ~280-310)", nrow(cohort)))

# Treatment arms
message("Treatment arm counts:")
print(table(cohort$treatment, useNA = "ifany"))
message("TCGA subtype counts:")
print(table(cohort$t1_cluster_name, useNA = "ifany"))

fwrite(cohort, OUT_COH, sep = "\t")
