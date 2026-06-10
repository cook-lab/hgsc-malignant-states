# ============================================================================
# 01b_neighborhood_feature_matrix.R
# ----------------------------------------------------------------------------
# PURPOSE: Build the 50um neighborhood cell-type-proportion feature matrix that
#   feeds the k=10 niche clustering. For every cell, count its neighbours within
#   a 50um radius (dbscan::frNN, self excluded) and express the neighbourhood as
#   a proportion vector over the 18 canonical cell types. Computed per TMA core
#   and per whole-tissue sample so neighbours never cross core/sample edges.
#
# INPUTS:
#   - SFEs (load_sfe): sfe_tma_filtered + 8 whole-tissue, with cell_label (06f/06g),
#     core_id + sample_type on the TMA.
#
# OUTPUTS:
#   - output/09_neighborhood/neighborhood_feature_matrix.rds
#       list(nb_matrix_filt, meta_filt, all_features_mat, empty_mask)
#       (same element names + shape that 02_neighborhood_k10_production.R reads)
#
# PIPELINE ROLE: producer for 02_neighborhood_k10_production.R (upstream cache
#   for the Fig 4G / Fig 6G niche GAMs). Number 01b so it runs before 02.
#
# RUNTIME TIER: heavy (per-cell frNN over all TMA cores + 8 whole tissues).
#
# Ported from 2026_final_xenium_analysis/scripts/09b_k_comparison.Rmd
#   (compute_neighborhood_features / per-core + per-sample feature build) and
#   09_neighborhood.R. Analytical logic preserved; paths routed through central
#   config, seed from CFG$seed, epithelial label "Transitioning" -> "Intermediate".
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- local({                       # robust: works via `Rscript <path>`, source(), or from-dir
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1]))
  else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
})
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, theme_lab, palettes
set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(dbscan)
  library(data.table)
})

message("\n=== 09b: Neighborhood feature matrix (50um frNN proportions) ===")

# --- Constants ---------------------------------------------------------------
RADIUS <- 50   # um

# 18 canonical cell types (lineage order). Epithelial polarization label
# standardized to "Intermediate epithelium" (was "Transitioning epithelium").
celltype_order <- c(
  "Ciliated epithelium", "SecA epithelium", "Intermediate epithelium",
  "SecB epithelium", "Mesothelial",
  "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
  "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

WT_SAMPLES <- sfe_names_wt

# --- Helpers -----------------------------------------------------------------

# Per-cell 50um neighbourhood cell-type proportions (self excluded).
compute_neighborhood_features <- function(coords, labels, radius = RADIUS,
                                          all_types = celltype_order) {
  n <- nrow(coords)
  fr <- dbscan::frNN(coords, eps = radius)
  from_idx <- rep(seq_along(fr$id), lengths(fr$id))
  to_idx   <- unlist(fr$id)
  keep <- from_idx != to_idx
  from_idx <- from_idx[keep]
  to_idx   <- to_idx[keep]
  dt <- data.table::data.table(from = from_idx, nb_type = labels[to_idx])
  counts <- dt[, .N, by = .(from, nb_type)]
  nb_mat <- matrix(0, nrow = n, ncol = length(all_types),
                   dimnames = list(NULL, all_types))
  for (tp in all_types) {
    sub <- counts[nb_type == tp]
    if (nrow(sub) > 0) nb_mat[sub$from, tp] <- sub$N
  }
  row_totals <- rowSums(nb_mat)
  row_totals[row_totals == 0] <- 1
  nb_mat <- nb_mat / row_totals
  nb_mat
}

extract_coords_labels <- function(sfe, label_col = "cell_label") {
  xy <- spatialCoords(sfe)
  labs <- as.character(colData(sfe)[[label_col]])
  # Repo label standardization: the deposited SFEs carry the original
  # "Transitioning epithelium" label; rename to "Intermediate epithelium" so the
  # neighbour labels match celltype_order (else those counts would silently drop).
  labs[labs == "Transitioning epithelium"] <- "Intermediate epithelium"
  data.frame(x = xy[, 1], y = xy[, 2], label = labs,
             cell_id = colnames(sfe), stringsAsFactors = FALSE)
}

# --- Extract coords + labels -------------------------------------------------
message("\n=== Extracting data ===")

# TMA (per-core neighbourhoods; drop off-core cells)
sfe_tma <- load_sfe("sfe_tma_filtered")
tma_dat <- extract_coords_labels(sfe_tma)
tma_dat$core_id <- as.character(sfe_tma$core_id)
tma_dat <- tma_dat[!is.na(tma_dat$core_id) & tma_dat$core_id != "Off core", ]
rm(sfe_tma); gc(verbose = FALSE)
message("TMA: ", format(nrow(tma_dat), big.mark = ","), " on-core cells")

# Whole tissue
wt_dat <- list()
for (sname in WT_SAMPLES) {
  sfe <- load_sfe(sname)
  df  <- extract_coords_labels(sfe)
  df$sample_id <- sname
  wt_dat[[sname]] <- df
  rm(sfe); gc(verbose = FALSE)
}
message("WT: ", format(sum(sapply(wt_dat, nrow)), big.mark = ","), " cells")

# --- Compute features --------------------------------------------------------
message("\n=== Computing 50um neighborhood proportions ===")

# TMA per-core (sample_id tag = "sfe_tma_filtered" so the consumer matches the
# merged TMA SFE by cell_id).
core_ids <- sort(unique(tma_dat$core_id))
tma_features_list <- list()
for (cid in core_ids) {
  sub <- tma_dat[tma_dat$core_id == cid, ]
  if (nrow(sub) < 20) next
  nb_mat <- compute_neighborhood_features(as.matrix(sub[, c("x", "y")]),
                                          sub$label, radius = RADIUS)
  tma_features_list[[cid]] <- data.frame(
    cell_id = sub$cell_id, sample_id = "sfe_tma_filtered",
    nb_mat, check.names = FALSE, stringsAsFactors = FALSE
  )
}
tma_features <- do.call(rbind, tma_features_list)
message("  TMA cores processed: ", length(tma_features_list))

# WT per-sample
wt_features_list <- list()
for (sname in WT_SAMPLES) {
  df <- wt_dat[[sname]]
  nb_mat <- compute_neighborhood_features(as.matrix(df[, c("x", "y")]),
                                          df$label, radius = RADIUS)
  wt_features_list[[sname]] <- data.frame(
    cell_id = df$cell_id, sample_id = sname,
    nb_mat, check.names = FALSE, stringsAsFactors = FALSE
  )
}
wt_features <- do.call(rbind, wt_features_list)

# --- Assemble + filter zero-neighbour cells ----------------------------------
all_features_mat <- rbind(tma_features, wt_features)
nb_matrix <- as.matrix(all_features_mat[, celltype_order])
empty_mask <- rowSums(nb_matrix) == 0
nb_matrix_filt <- nb_matrix[!empty_mask, ]
meta_filt <- all_features_mat[!empty_mask, c("cell_id", "sample_id")]

message("\nFeature matrix: ", format(nrow(nb_matrix_filt), big.mark = ","),
        " x ", ncol(nb_matrix_filt),
        "  (dropped ", format(sum(empty_mask), big.mark = ","),
        " zero-neighbour cells)")

# --- Save --------------------------------------------------------------------
out_rds <- cfg_path("output_root", "09_neighborhood", "neighborhood_feature_matrix.rds")
saveRDS(list(
  nb_matrix_filt   = nb_matrix_filt,
  meta_filt        = meta_filt,
  all_features_mat = all_features_mat,
  empty_mask       = empty_mask
), out_rds)
message("Saved: ", out_rds)

log_session()
