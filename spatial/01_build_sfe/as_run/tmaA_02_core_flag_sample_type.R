# Trimmed copy (2026-10-22) of ~/core_qc_viz.R, the script that ran (2025-10-14, interactively, in the home folder).
# Kept: the per-core QC summary and its lowest-5% lists, from which the flagged cores were chosen by hand (flag_cores),
# the core_flag column, and sample_type (cores 203-217 = fallopian tube, all others tumour), saved into tma_1 / tma_2.
# Removed: the QC histograms and box plots before it and two closing colnames() prints. The untrimmed script is in the
# home folder. Both rules hold for every cell of the earlier object se.rds (checked 2026-10-22).

library(dplyr)

# Function to get core summary from a SFE object
get_core_summary <- function(sfe, name) {
  cell_meta <- as.data.frame(colData(sfe))
  cell_meta %>%
    group_by(core_id) %>%
    summarise(
      median_nCounts = median(nCounts, na.rm = TRUE),
      median_nGenes = median(nGenes, na.rm = TRUE),
      median_nCounts_normed = median(nCounts_normed, na.rm = TRUE),
      median_nGenes_normed = median(nGenes_normed, na.rm = TRUE),
      dataset = name,
      .groups = "drop"
    )
}

# Core summaries for tma_1 and tma_2
core_summary_1 <- get_core_summary(tma_1, "tma_1")
core_summary_2 <- get_core_summary(tma_2, "tma_2")

# Combine both datasets
core_summary_all <- bind_rows(core_summary_1, core_summary_2)

# Function to get lowest 5% cores
lowest_5pct <- function(x, core_ids) {
  cutoff <- quantile(x, 0.05, na.rm = TRUE)
  cores <- core_ids[x <= cutoff]
  list(cores = cores, cutoff = cutoff)
}

# Calculate lowest 5% for each metric
low_nCounts <- lowest_5pct(core_summary_all$median_nCounts, core_summary_all$core_id)
low_nGenes <- lowest_5pct(core_summary_all$median_nGenes, core_summary_all$core_id)
low_nCounts_normed <- lowest_5pct(core_summary_all$median_nCounts_normed, core_summary_all$core_id)
low_nGenes_normed <- lowest_5pct(core_summary_all$median_nGenes_normed, core_summary_all$core_id)

# Show results
low_nCounts
low_nGenes
low_nCounts_normed
low_nGenes_normed


library(dplyr)
# Vector of core IDs to flag
flag_cores <- c(6, 13, 46, 114, 143, 161, 172, 26, 75, 107, 111, 118, 203, 206, 208, 209)

# Add core_flag column based on core_id
colData(tma_1)$core_flag <- ifelse(as.integer(colData(tma_1)$core_id) %in% flag_cores, "FLAG", "")
colData(tma_2)$core_flag <- ifelse(as.integer(colData(tma_2)$core_id) %in% flag_cores, "FLAG", "")

fallopian_cores <- 203:217

# Assign sample_type based on core_id
colData(tma_1)$sample_type <- ifelse(as.integer(colData(tma_1)$core_id) %in% fallopian_cores,
                                     "fallopian", "tumour")
colData(tma_2)$sample_type <- ifelse(as.integer(colData(tma_2)$core_id) %in% fallopian_cores,
                                     "fallopian", "tumour")

saveHDF5SummarizedExperiment(tma_2, dir = "tma_2_sfe_hdf5", replace = TRUE)
saveHDF5SummarizedExperiment(tma_1, dir = "tma_1_sfe_hdf5", replace = TRUE)
