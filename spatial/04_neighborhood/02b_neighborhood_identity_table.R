# ============================================================================
# 02b_neighborhood_identity_table.R
# ----------------------------------------------------------------------------
# PURPOSE: Document the neighborhood integer -> biology mapping BY DATA. Reads
#   the k=10 neighborhood composition (mean cell-type proportions per nb_) and
#   emits a compact identity table that pairs each nb_ integer with (a) its
#   canonical readable name from 00_setup `nb_names` and (b) the empirically
#   dominant cell types that ground that name. This makes the nb_1..nb_10 ->
#   biology mapping auditable from the deposited composition rather than from a
#   hand-curated comment.
#
# INPUTS (first that exists wins):
#   - <output_root>/09_neighborhood/neighborhood_composition_k10.csv  (regenerated)
#   - <data_root>/2026_final_xenium_analysis/output/09_neighborhood/neighborhood_composition.csv
#       (deposited fallback; uses the legacy "Transitioning epithelium" column
#        label, which is normalized to "Intermediate epithelium" on read)
#   Both are nb_ (rows) x cell-type (cols) matrices of mean proportions.
#
# OUTPUTS:
#   - <output_root>/09_neighborhood/neighborhood_identity_table.csv
#       columns: neighborhood (nb_1..nb_10), name (00_setup nb_names),
#                top1/top2/top3 (dominant cell types) + their mean proportions.
#
# MANUSCRIPT PANEL(S): documentation/QC for Fig 4G / 6G niche naming
#   (not a panel itself). Grounds the nb_names nomenclature used throughout
#   spatial/04_neighborhood.
# RUNTIME TIER: fast
#
# Paths routed through central config; epithelial label "Transitioning" ->
# "Intermediate" normalized on read; nb_names sourced from 00_setup (single
# source of truth for the integer->biology nomenclature).
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # nb_names, palettes, paths
set.seed(CFG$seed)

# --- Resolve composition input (output_root first, deposited data_root fallback) ---

nb_dir   <- file.path(out_dir, "09_neighborhood")
comp_out <- file.path(nb_dir, "neighborhood_composition_k10.csv")
comp_dep <- file.path(xen_root, "output", "09_neighborhood", "neighborhood_composition.csv")

comp_path <- if (file.exists(comp_out)) comp_out else comp_dep
if (!file.exists(comp_path)) {
  stop("Neighborhood composition CSV not found in either location:\n  ",
       comp_out, "\n  ", comp_dep)
}
message("=== 02b neighborhood identity table ===")
message("Reading composition: ", comp_path)

# row.names = 1 -> nb_ rownames; columns are cell-type mean proportions.
comp <- read.csv(comp_path, row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)

# Normalize the legacy epithelial label so the deposited fallback matches the
# standardized nomenclature used everywhere else ("Transitioning" -> "Intermediate").
colnames(comp) <- sub("^Transitioning epithelium$", "Intermediate epithelium",
                      colnames(comp))

# --- Build the identity table -------------------------------------------------
# For each nb_ in canonical order, take the top ~3 cell types by mean proportion
# and pair them with the canonical readable name from 00_setup `nb_names`.

n_top   <- 3
nb_ids  <- names(nb_names)                  # nb_1 .. nb_10 (canonical order)
present <- intersect(nb_ids, rownames(comp))
missing <- setdiff(nb_ids, rownames(comp))
if (length(missing) > 0) {
  warning("Composition CSV missing rows for: ", paste(missing, collapse = ", "))
}

rows <- lapply(present, function(nb) {
  props  <- sort(unlist(comp[nb, , drop = TRUE]), decreasing = TRUE)
  top    <- head(props, n_top)
  # Pad to n_top if a neighborhood somehow has fewer cell-type columns.
  ct     <- c(names(top), rep(NA_character_, n_top - length(top)))
  pr     <- c(round(unname(top), 4), rep(NA_real_, n_top - length(top)))
  data.frame(
    neighborhood   = nb,
    name           = unname(nb_names[nb]),
    top1_celltype  = ct[1], top1_prop = pr[1],
    top2_celltype  = ct[2], top2_prop = pr[2],
    top3_celltype  = ct[3], top3_prop = pr[3],
    stringsAsFactors = FALSE
  )
})

identity_tbl <- do.call(rbind, rows)

# Preserve canonical nb_ ordering (nb_1, nb_2, ..., nb_10) rather than lexical.
identity_tbl <- identity_tbl[match(present, identity_tbl$neighborhood), ]
rownames(identity_tbl) <- NULL

# --- Write -------------------------------------------------------------------

if (!dir.exists(nb_dir)) dir.create(nb_dir, recursive = TRUE)
out_csv <- file.path(nb_dir, "neighborhood_identity_table.csv")
write.csv(identity_tbl, out_csv, row.names = FALSE)

message("\nNeighborhood identity table (data-grounded nb_ -> biology):")
for (i in seq_len(nrow(identity_tbl))) {
  r <- identity_tbl[i, ]
  message(sprintf("  %-5s %-32s %s (%.3f), %s (%.3f), %s (%.3f)",
                  r$neighborhood, r$name,
                  r$top1_celltype, r$top1_prop,
                  r$top2_celltype, r$top2_prop,
                  r$top3_celltype, r$top3_prop))
}
message("\nSaved: ", out_csv)
