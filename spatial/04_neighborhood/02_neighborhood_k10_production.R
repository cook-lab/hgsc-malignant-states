# ============================================================================
# 02_neighborhood_k10_production.R
# ----------------------------------------------------------------------------
# PURPOSE: Final k=10 k-means neighborhood assignment; labels each cell's 50um niche and writes neighborhood + neighborhood_name to every SFE.
#
# INPUTS:
#   - output/09_neighborhood/neighborhood_feature_matrix.rds
#   - SFEs (load_sfe): sfe_tma_filtered + 8 whole-tissue
#
# OUTPUTS:
#   - output/09_neighborhood/neighborhood_assignments.csv (canonical k=10 assignments; read by stages 04-08 + clinical)
#   - kmeans_centers_k10.csv, neighborhood_composition_k10.csv, neighborhood_name_mapping_k10.csv
#   - SFEs updated with neighborhood, neighborhood_name
#
# MANUSCRIPT PANEL(S): Upstream cache for Fig 4G/6G niche GAMs (neighborhood assignments).
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

library(FNN)

message("\n=== 09c: k=10 Neighborhood Production Assignment ===")

# --- Load feature matrix ---
nb_dir <- file.path(out_dir, "09_neighborhood")
dat <- readRDS(file.path(nb_dir, "neighborhood_feature_matrix.rds"))
nb_matrix_filt <- dat$nb_matrix_filt
meta_filt <- dat$meta_filt
all_features_mat <- dat$all_features_mat
empty_mask <- dat$empty_mask

message("Loaded feature matrix: ", format(nrow(nb_matrix_filt), big.mark = ","), " cells")

# --- Stratified subsample ---
set.seed(CFG$seed)
sub_idx <- unlist(tapply(seq_len(nrow(nb_matrix_filt)), meta_filt$sample_id, function(idx) {
  n_take <- max(1, round(length(idx) * 100000 / nrow(nb_matrix_filt)))
  sample(idx, min(n_take, length(idx)))
}))
nb_sub <- nb_matrix_filt[sub_idx, ]
message("Subsample: ", format(nrow(nb_sub), big.mark = ","), " cells")

# --- k=10 clustering ---
message("Running k=10 k-means...")
set.seed(CFG$seed)
km <- kmeans(nb_sub, centers = 10, nstart = 50, iter.max = 200)

# Assign all cells
message("Assigning all cells to nearest center...")
assigns <- FNN::get.knnx(km$centers, nb_matrix_filt, k = 1)$nn.index[, 1]

# --- Neighborhood naming ---
# Use the SINGLE canonical `nb_names` map sourced from 00_setup.R (keyed
# nb_1..nb_10; identities verified against neighborhood_composition.csv). A
# divergent local map here previously shadowed it and disagreed with 00_setup /
# 03 / 05-08 on which cluster is SecA vs SecB (audit B2). Do NOT redefine locally.
# NOTE: this integer->biology map assumes the canonical k-means numbering; if a
# fresh re-run renumbers clusters, verify identities against the
# neighborhood_composition_k10.csv written below before trusting the labels.

# Create labeled assignments
nb_id    <- paste0("nb_", assigns)
nb_named <- nb_names[paste0("nb_", assigns)]

# Build full assignment table (including unassigned)
all_features_mat$neighborhood    <- NA_character_
all_features_mat$neighborhood_name <- NA_character_
all_features_mat$neighborhood[!empty_mask]      <- nb_id
all_features_mat$neighborhood_name[!empty_mask]  <- nb_named
all_features_mat$neighborhood[empty_mask]        <- "nb_unassigned"
all_features_mat$neighborhood_name[empty_mask]   <- "Unassigned"

# --- Summary ---
message("\n=== k=10 Neighborhood Summary ===")
sizes <- table(factor(assigns, levels = 1:10))
for (i in 1:10) {
  message(sprintf("  nb_%d  %-30s  %s cells (%.1f%%)",
                  i, nb_names[paste0("nb_", i)],
                  format(sizes[i], big.mark = ","),
                  100 * sizes[i] / sum(sizes)))
}
message("  Unassigned: ", sum(empty_mask))

# --- Save assignments CSV ---
assign_df <- all_features_mat[, c("cell_id", "sample_id", "neighborhood", "neighborhood_name")]
write.csv(assign_df, file.path(nb_dir, "neighborhood_assignments.csv"), row.names = FALSE)
message("\nSaved: neighborhood_assignments.csv")

# Save centers
write.csv(km$centers, file.path(nb_dir, "kmeans_centers_k10.csv"))

# Save composition
comp_mat <- t(sapply(1:10, function(i) colMeans(nb_matrix_filt[assigns == i, , drop = FALSE])))
rownames(comp_mat) <- paste0("nb_", 1:10)
write.csv(comp_mat, file.path(nb_dir, "neighborhood_composition_k10.csv"))

# Save name mapping
name_map <- data.frame(
  neighborhood = paste0("nb_", 1:10),
  neighborhood_name = nb_names[paste0("nb_", 1:10)],
  n_cells = as.numeric(sizes),
  pct = round(100 * as.numeric(sizes) / sum(sizes), 1),
  stringsAsFactors = FALSE
)
write.csv(name_map, file.path(nb_dir, "neighborhood_name_mapping_k10.csv"), row.names = FALSE)

# --- Write to SFEs ---
message("\n=== Writing neighborhoods to SFEs ===")

sfe_names <- c("sfe_tma_filtered", sfe_names_wt)

for (sname in sfe_names) {
  message("\n--- Processing ", sname, " ---")
  sfe <- load_sfe(sname)
  message("  Loaded: ", format(ncol(sfe), big.mark = ","), " cells")

  # Archive old neighborhood if it exists
  if ("neighborhood" %in% colnames(colData(sfe))) {
    sfe$neighborhood_old <- sfe$neighborhood
    message("  Archived old neighborhood → neighborhood_old")
  }
  if ("neighborhood_name" %in% colnames(colData(sfe))) {
    sfe$neighborhood_name_old <- sfe$neighborhood_name
    message("  Archived old neighborhood_name → neighborhood_name_old")
  }

  # Match cell_ids
  match_idx <- match(colnames(sfe), assign_df$cell_id)
  n_matched <- sum(!is.na(match_idx))
  n_na <- sum(is.na(match_idx))

  sfe$neighborhood      <- assign_df$neighborhood[match_idx]
  sfe$neighborhood_name <- assign_df$neighborhood_name[match_idx]

  message("  Matched: ", format(n_matched, big.mark = ","),
          " | NA: ", format(n_na, big.mark = ","))

  # Show distribution
  tab <- table(sfe$neighborhood_name, useNA = "ifany")
  for (nm in sort(names(tab))) {
    message("    ", nm, ": ", format(tab[nm], big.mark = ","))
  }

  # Realize and save
  message("  Realizing assays...")
  for (a in assayNames(sfe)) {
    assay(sfe, a) <- as(assay(sfe, a), "dgCMatrix")
  }

  message("  Saving...")
  save_sfe(sfe, sname)
  rm(sfe); gc(verbose = FALSE)
  message("  Done.")
}

message("\n=== 09c Complete ===")
message("All 9 SFEs updated with: neighborhood, neighborhood_name")
message("Old assignments preserved in: neighborhood_old, neighborhood_name_old")
message("Neighborhood identities (nb_1..nb_10) come from 00_setup nb_names; verify against neighborhood_composition_k10.csv.")

log_session()
