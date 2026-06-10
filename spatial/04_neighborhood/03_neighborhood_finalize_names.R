# ============================================================================
# 03_neighborhood_finalize_names.R
# ----------------------------------------------------------------------------
# PURPOSE: Companion to k10 production: refresh neighborhood_name on every SFE by mapping nb_1..nb_10 through nb_names. Name-only; does not recompute clusters.
#
# INPUTS:
#   - SFEs (load_sfe) with existing neighborhood colData (nb_1..nb_10)
#
# OUTPUTS:
#   - same SFEs with neighborhood_name refreshed
#
# MANUSCRIPT PANEL(S): Supports Fig 4G/6G niche naming.
# RUNTIME TIER: moderate
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

sfe_names <- c("sfe_tma_filtered", sfe_names_wt)

message("=== 09c finalize neighborhood names ===")
message("nb_names mapping in use:")
for (i in seq_along(nb_names)) {
  message(sprintf("  %-5s -> %s", names(nb_names)[i], nb_names[i]))
}

for (nm in sfe_names) {
  message(sprintf("\n[%s] Loading %s ...", format(Sys.time()), nm))
  sfe <- load_sfe(nm)
  if (!"neighborhood" %in% colnames(colData(sfe))) {
    message("  ! missing 'neighborhood' colData — skipping")
    next
  }

  old_counts <- if ("neighborhood_name" %in% colnames(colData(sfe)))
    sort(table(sfe$neighborhood_name), decreasing = TRUE) else NULL

  # Map nb_i -> readable name; unassigned/NA stays "Unassigned"
  new_name <- ifelse(sfe$neighborhood %in% names(nb_names),
                     unname(nb_names[sfe$neighborhood]),
                     "Unassigned")
  sfe$neighborhood_name <- new_name

  new_counts <- sort(table(sfe$neighborhood_name), decreasing = TRUE)

  message("  updated neighborhood_name counts:")
  print(new_counts)

  save_sfe(sfe, nm)
  rm(sfe); gc(verbose = FALSE)
}

message("\n=== Done — neighborhood_name refreshed on all SFEs ===")
