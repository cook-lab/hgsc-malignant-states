# ============================================================================
# 04_vascular_proximity.R
# ----------------------------------------------------------------------------
# PURPOSE: Distance-to-vasculature: for each SecA / Intermediate / SecB cell compute distance to nearest vascular cell (RANN nn2 k=1); per-sample medians + Wilcoxon.
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label, spatial coords
#
# OUTPUTS:
#   - output/22_vascular_proximity/vascular_distance_summary.csv
#   - vascular_distance_all_cells.csv
#
# MANUSCRIPT PANEL(S): Fig 4H, SF12.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

library(RANN)

out_dir <- file.path("output", "22_vascular_proximity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
wt_names <- c("sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432", "sfe_OTB_2454",
              "sfe_OTB_2457", "sfe_OTB_2461", "sfe_SP24_24824", "sfe_SP24_25573")

vascular_types <- c("Pericyte", "Endothelial")
query_types    <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

# ---------------------------------------------------------------------------
# 1. Compute distances per WT sample
# ---------------------------------------------------------------------------
results <- list()

for (nm in wt_names) {
  label <- sub("sfe_", "", nm)
  message(sprintf("[%s] Loading %s ...", Sys.time(), label))

  sfe <- load_sfe(nm)
  cd <- data.table(
    cell_label = sfe$cell_label,
    x = spatialCoords(sfe)[, 1],
    y = spatialCoords(sfe)[, 2]
  )

  vasc <- cd[cell_label %in% vascular_types]
  message(sprintf("  %s vascular cells (%s Pericyte, %s Endothelial)",
                  format(nrow(vasc), big.mark = ","),
                  sum(vasc$cell_label == "Pericyte"),
                  sum(vasc$cell_label == "Endothelial")))

  if (nrow(vasc) < 10) {
    message("  Skipping — too few vascular cells")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  vasc_coords <- as.matrix(vasc[, .(x, y)])

  for (ct in query_types) {
    query <- cd[cell_label == ct]
    if (nrow(query) < 10) next

    nn <- nn2(vasc_coords, as.matrix(query[, .(x, y)]), k = 1)

    results[[length(results) + 1]] <- data.table(
      sample_id        = label,
      cell_label       = ct,
      dist_to_vascular = as.numeric(nn$nn.dists)
    )
  }

  message(sprintf("  Queried %s epithelial cells",
                  format(sum(cd$cell_label %in% query_types), big.mark = ",")))
  rm(sfe, cd); gc(verbose = FALSE)
}

dist_dt <- rbindlist(results)
message(sprintf("\nTotal measurements: %s", format(nrow(dist_dt), big.mark = ",")))

# ---------------------------------------------------------------------------
# 2. Also run on TMA
# ---------------------------------------------------------------------------
message(sprintf("[%s] Loading TMA ...", Sys.time()))
sfe_tma <- load_sfe("sfe_tma_filtered")
cd_tma <- data.table(
  cell_label  = sfe_tma$cell_label,
  patient_id  = sfe_tma$patient_id,
  core_id     = sfe_tma$core_id,
  sample_type = sfe_tma$sample_type,
  x = spatialCoords(sfe_tma)[, 1],
  y = spatialCoords(sfe_tma)[, 2]
)

# Per-core within TMA (tumour cores only).
# Each TMA core is a physically separate tissue piece with its own coordinate
# system — pooling across a patient's cores would create artificial nearest
# neighbours across unrelated coordinates. Distances are computed per core and
# aggregated to patient level downstream.
tma_results <- list()
tumour_core_ids <- unique(cd_tma[sample_type == "tumour", .(core_id, patient_id)])
for (i in seq_len(nrow(tumour_core_ids))) {
  cid <- tumour_core_ids$core_id[i]
  pid <- tumour_core_ids$patient_id[i]
  sub <- cd_tma[core_id == cid & sample_type == "tumour"]
  vasc <- sub[cell_label %in% vascular_types]
  if (nrow(vasc) < 5) next

  vasc_coords <- as.matrix(vasc[, .(x, y)])
  for (ct in query_types) {
    query <- sub[cell_label == ct]
    if (nrow(query) < 5) next
    nn <- nn2(vasc_coords, as.matrix(query[, .(x, y)]), k = 1)
    tma_results[[length(tma_results) + 1]] <- data.table(
      sample_id        = paste0("TMA_", pid),
      patient_id       = as.character(pid),
      core_id          = as.character(cid),
      cell_label       = ct,
      dist_to_vascular = as.numeric(nn$nn.dists)
    )
  }
}
rm(sfe_tma, cd_tma); gc(verbose = FALSE)

tma_dist <- rbindlist(tma_results, fill = TRUE)
message(sprintf("TMA measurements: %s across %s cores / %s patients",
                format(nrow(tma_dist), big.mark = ","),
                length(unique(tma_dist$core_id)),
                length(unique(tma_dist$patient_id))))

# Combine (TMA has extra patient_id / core_id columns; WT fills with NA)
dist_all <- rbind(
  dist_dt[, tissue := "whole_tissue"],
  tma_dist[, tissue := "TMA"],
  use.names = TRUE, fill = TRUE
)

# ---------------------------------------------------------------------------
# 3. Summary statistics
# ---------------------------------------------------------------------------
# WT: straight per-sample summary (each WT = one tissue, single coord frame)
wt_summary <- dist_all[tissue == "whole_tissue",
  .(median_dist = median(dist_to_vascular),
    mean_dist   = mean(dist_to_vascular),
    sd_dist     = sd(dist_to_vascular),
    q25         = quantile(dist_to_vascular, 0.25),
    q75         = quantile(dist_to_vascular, 0.75),
    n           = .N),
  by = .(cell_label, sample_id, tissue)]

# TMA: first median per core (each core is a physically independent tissue),
# then average per-core medians to patient level (one value per patient).
tma_core_dt <- dist_all[tissue == "TMA",
  .(median_core = median(dist_to_vascular),
    mean_core   = mean(dist_to_vascular),
    n_core      = .N),
  by = .(cell_label, sample_id, patient_id, core_id)]

tma_patient_summary <- tma_core_dt[,
  .(median_dist = mean(median_core),   # patient = mean of core medians
    mean_dist   = mean(mean_core),
    sd_dist     = sd(median_core),
    q25         = quantile(median_core, 0.25),
    q75         = quantile(median_core, 0.75),
    n           = sum(n_core),
    n_cores     = .N),
  by = .(cell_label, sample_id, patient_id)]
tma_patient_summary[, tissue := "TMA"]

# Combine (drop patient_id / n_cores columns for WT; they don't apply there)
summary_dt <- rbind(
  wt_summary[, `:=`(patient_id = NA_character_, n_cores = NA_integer_)][],
  tma_patient_summary[, .(cell_label, sample_id, tissue, median_dist, mean_dist,
                           sd_dist, q25, q75, n, patient_id, n_cores)],
  use.names = TRUE, fill = TRUE
)

# Save the per-core TMA layer for transparency / re-use
fwrite(tma_core_dt, file.path(out_dir, "vascular_distance_tma_per_core.csv"))

# Per-sample paired comparison (WT only)
wt_summary <- summary_dt[tissue == "whole_tissue"]
paired <- merge(
  wt_summary[cell_label == "SecA epithelium", .(sample_id, secA_med = median_dist)],
  wt_summary[cell_label == "SecB epithelium", .(sample_id, secB_med = median_dist)],
  by = "sample_id"
)

# Statistical tests
sink(file.path(out_dir, "vascular_distance_stats.txt"))
cat("=== Vascular Proximity Analysis ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("Per-sample median distances (WT):\n")
print(paired)
cat(sprintf("\nOverall median — SecA: %.1f µm, SecB: %.1f µm\n",
            median(paired$secA_med), median(paired$secB_med)))

if (nrow(paired) >= 4) {
  wt <- wilcox.test(paired$secA_med, paired$secB_med, paired = TRUE)
  cat(sprintf("\nPaired Wilcoxon signed-rank test (n=%d samples):\n", nrow(paired)))
  cat(sprintf("  V = %.0f, p = %.4g\n", wt$statistic, wt$p.value))
}

# Pooled cell-level test (WT)
wt_pooled <- wilcox.test(
  dist_dt[cell_label == "SecA epithelium"]$dist_to_vascular,
  dist_dt[cell_label == "SecB epithelium"]$dist_to_vascular
)
cat(sprintf("\nPooled Wilcoxon rank-sum (all WT cells):\n"))
cat(sprintf("  W = %.0f, p = %.4g\n", wt_pooled$statistic, wt_pooled$p.value))

# Effect size
cat(sprintf("\nEffect size (median difference): %.1f µm\n",
            median(paired$secA_med) - median(paired$secB_med)))
sink()

# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------
fwrite(dist_all, file.path(out_dir, "vascular_distance_all_cells.csv"))
fwrite(summary_dt, file.path(out_dir, "vascular_distance_summary.csv"))

message(sprintf("\n[%s] Done — saved to %s", Sys.time(), out_dir))
message("  vascular_distance_all_cells.csv")
message("  vascular_distance_summary.csv")
message("  vascular_distance_stats.txt")
