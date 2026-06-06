# ============================================================================
# 01_niche_succession_gams.R
# ----------------------------------------------------------------------------
# PURPOSE: Niche succession GAMs: model 50um neighborhood composition + pathway activity along the SecA->SecB polarization gradient (WT-primary, TMA-validation). Base version (supersedes v1/v2).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label, pathway_*, polarization_UCell
#
# OUTPUTS:
#   - output/16b_niche_succession_gams/neighborhood_features.rds (key cache)
#   - GAM tables + inflection points + figures
#
# MANUSCRIPT PANEL(S): Fig 4G, Fig 6G, SF12.
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

library(dbscan)
library(data.table)
library(mgcv)
library(circlize)

set.seed(CFG$seed)

message("\n", strrep("=", 70))
message("Script 16b: Niche Succession GAMs Along Polarization Gradient")
message("  Design: WT-primary / TMA-validation split")
message(strrep("=", 70))

# --- Output directories ------------------------------------------------------

out_path <- file.path(out_dir, "16b_niche_succession_gams")
fig_path <- file.path(fig_dir, "16b_niche_succession_gams")
for (d in c(out_path, fig_path)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- SFE names ----------------------------------------------------------------

sfe_names <- c(
  "sfe_tma_filtered",
  "sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
  "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
  "sfe_SP24_24824", "sfe_SP24_25573"
)

# --- Cell type order (18 types) -----------------------------------------------

celltype_order <- c(
  "Ciliated epithelium", "SecA epithelium", "Intermediate epithelium",
  "SecB epithelium",
  "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
  "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

secretory_types <- c("SecA epithelium", "Intermediate epithelium",
                     "SecB epithelium")

immune_types <- c("T cell", "NK cell", "B cell", "Plasma cell",
                  "Macrophage", "Conventional dendritic cell",
                  "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell")

# --- Pathway column names (22 from 9b) ----------------------------------------

pathway_names <- c(
  "pathway_proliferation", "pathway_apoptosis", "pathway_hypoxia",
  "pathway_rtk_ras", "pathway_tgfb", "pathway_pi3k_akt_mtor",
  "pathway_wnt", "pathway_myc", "pathway_p53", "pathway_hippo",
  "pathway_notch", "pathway_jak_stat", "pathway_immune_checkpoint",
  "pathway_type_i_ifn", "pathway_type_ii_ifn", "pathway_cytotoxicity",
  "pathway_antigen_presentation", "pathway_emt", "pathway_angiogenesis",
  "pathway_complement", "pathway_chemokine", "pathway_nfkb"
)

RADIUS <- 50  # micrometers

# --- Delete ALL old cached RDS files (schema changed: per-cell-type features) --

cell_data_path <- file.path(out_path, "all_cell_data_16b.rds")
nb_features_path <- file.path(out_path, "neighborhood_features.rds")
gam_results_path <- file.path(out_path, "gam_results.rds")

# Cache validation: only delete if schema is wrong (missing nb_ct_ columns)
# Once correct schema is cached, reuse it.
if (file.exists(nb_features_path)) {
  test_nb <- tryCatch(readRDS(nb_features_path), error = function(e) NULL)
  if (!is.null(test_nb) && any(grepl("^nb_ct_", names(test_nb)))) {
    message("Cached neighborhood_features.rds has per-cell-type columns — reusing.")
  } else {
    message("Cached neighborhood_features.rds missing per-cell-type columns — deleting.")
    file.remove(nb_features_path)
    if (file.exists(gam_results_path)) file.remove(gam_results_path)
  }
  rm(test_nb); gc(verbose = FALSE)
}

# ============================================================================
# PART 1: Data Extraction
# ============================================================================

message("\n", strrep("=", 70))
message("PART 1: Data Extraction")
message(strrep("=", 70))

cell_data_list <- list()

for (sname in sfe_names) {
  message("\n  Loading ", sname, " ...")
  t0 <- Sys.time()

  sfe <- load_sfe(sname)
  cd <- as.data.frame(colData(sfe))
  xy <- spatialCoords(sfe)

  # Sample ID and core_id

  is_tma <- grepl("tma", sname, ignore.case = TRUE)
  if (is_tma) {
    sample_id <- "TMA"
    keep <- !is.na(cd$core_id) & !grepl("Off", cd$core_id, ignore.case = TRUE)
  } else {
    sample_id <- sname
    keep <- rep(TRUE, nrow(cd))
  }

  cd <- cd[keep, ]
  xy <- xy[keep, , drop = FALSE]

  # Extract pathway scores
  pathway_cols <- intersect(pathway_names, colnames(cd))
  missing_pw <- setdiff(pathway_names, colnames(cd))

  # Build data.table
  dt <- data.table(
    cell_id    = colnames(sfe)[keep],
    sample_id  = sample_id,
    cell_label = cd$cell_label,
    x          = xy[, 1],
    y          = xy[, 2]
  )

  # Add core_id for TMA cells (NA for WT samples)
  if (is_tma) {
    dt[, core_id := cd$core_id]
  } else {
    dt[, core_id := NA_character_]
  }

  # Pathway scores
  for (pw in pathway_cols) {
    set(dt, j = pw, value = cd[[pw]])
  }
  for (pw in missing_pw) {
    set(dt, j = pw, value = NA_real_)
  }

  # Polarization score (only valid for secretory cells)
  if ("polarization_UCell" %in% colnames(cd)) {
    dt[, polarization_UCell := cd$polarization_UCell]
  } else {
    dt[, polarization_UCell := NA_real_]
  }

  cell_data_list[[sname]] <- dt

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("    %s cells (%d secretory) in %.1fs",
                  format(nrow(dt), big.mark = ","),
                  sum(dt$cell_label %in% secretory_types, na.rm = TRUE),
                  elapsed))

  rm(sfe, cd, xy, dt)
  gc(verbose = FALSE)
}

all_cells <- rbindlist(cell_data_list, fill = TRUE)
rm(cell_data_list)
gc(verbose = FALSE)

message(sprintf("\nTotal cells: %s", format(nrow(all_cells), big.mark = ",")))
message(sprintf("Total secretory: %s",
                format(sum(all_cells$cell_label %in% secretory_types, na.rm = TRUE),
                       big.mark = ",")))

# Report TMA core counts
tma_cores <- all_cells[sample_id == "TMA", uniqueN(core_id)]
message(sprintf("TMA cores: %d", tma_cores))

saveRDS(all_cells, cell_data_path)
message("  Saved all_cell_data_16b.rds")

# Summary
message("\nCell label distribution:")
print(all_cells[, .N, by = cell_label][order(-N)])

message("\nSample distribution:")
print(all_cells[, .N, by = sample_id][order(-N)])

if (tma_cores > 0) {
  message("\nTMA core distribution:")
  print(all_cells[sample_id == "TMA", .N, by = core_id][order(-N)])
}

# ============================================================================
# PART 2: Neighborhood Feature Computation (50um)
# ============================================================================

message("\n", strrep("=", 70))
message("PART 2: Neighborhood Feature Computation (50um radius)")
message("  WT: per-sample frNN on whole tissue")
message("  TMA: per-core frNN to respect core boundaries")
message("  Features: proportions, counts, all-neighbor pathways, per-cell-type pathways")
message(strrep("=", 70))

# Helper function: compute neighborhood features for a set of cells
compute_nb_features <- function(dt_cells, group_label) {
  labels <- dt_cells$cell_label
  sec_idx <- which(labels %in% secretory_types)
  n_sec <- length(sec_idx)

  if (n_sec < 10) {
    message(sprintf("    Skipping %s: only %d secretory cells (< 10)", group_label, n_sec))
    return(NULL)
  }

  coords <- as.matrix(dt_cells[, .(x, y)])
  message(sprintf("    %s: %s total cells, %d secretory query cells",
                  group_label, format(nrow(dt_cells), big.mark = ","), n_sec))

  # Build frNN tree
  message(sprintf("    Building frNN tree for %s ...", group_label))
  fr <- dbscan::frNN(coords, eps = RADIUS)
  message("    frNN done.")

  # Pathway matrix for all cells
  pw_mat <- as.matrix(dt_cells[, ..pathway_names])

  # Pre-allocate result matrices
  n_query <- length(sec_idx)

  # Cell type proportions
  ct_prop_mat <- matrix(0, nrow = n_query, ncol = length(celltype_order),
                        dimnames = list(NULL, celltype_order))
  # Cell type counts (raw)
  ct_count_mat <- matrix(0L, nrow = n_query, ncol = length(celltype_order),
                         dimnames = list(NULL, celltype_order))
  # All-neighbor pathway means
  pw_all_mat  <- matrix(NA_real_, nrow = n_query, ncol = length(pathway_names),
                        dimnames = list(NULL, pathway_names))
  # Per-cell-type pathway means: 18 cell types x 22 pathways = 396 features
  # Stored as a list of matrices, one per cell type
  ct_pw_mats <- lapply(celltype_order, function(ct) {
    matrix(NA_real_, nrow = n_query, ncol = length(pathway_names),
           dimnames = list(NULL, pathway_names))
  })
  names(ct_pw_mats) <- celltype_order

  n_nb_total  <- integer(n_query)

  # Process each secretory cell
  for (ii in seq_len(n_query)) {
    qi <- sec_idx[ii]
    nb_ids   <- fr$id[[qi]]
    nb_dists <- fr$dist[[qi]]

    # Remove self
    not_self <- nb_ids != qi
    nb_ids   <- nb_ids[not_self]
    nb_dists <- nb_dists[not_self]

    if (length(nb_ids) == 0) next

    nb_labels <- labels[nb_ids]

    # Cell type counts and proportions
    tab <- table(factor(nb_labels, levels = celltype_order))
    total <- sum(tab)
    ct_count_mat[ii, ] <- as.integer(tab)
    ct_prop_mat[ii, ] <- as.numeric(tab) / total
    n_nb_total[ii] <- total

    # Mean pathway scores of ALL neighbors
    if (total > 0) {
      pw_nb <- pw_mat[nb_ids, , drop = FALSE]
      pw_all_mat[ii, ] <- colMeans(pw_nb, na.rm = TRUE)
    }

    # Per-cell-type pathway means
    for (ct in celltype_order) {
      ct_mask <- nb_labels == ct
      n_ct <- sum(ct_mask)
      if (n_ct > 0) {
        pw_ct <- pw_mat[nb_ids[ct_mask], , drop = FALSE]
        ct_pw_mats[[ct]][ii, ] <- colMeans(pw_ct, na.rm = TRUE)
      }
      # If n_ct == 0, stays NA_real_ (already initialized)
    }
  }

  # Assemble data.table
  dt_nb <- data.table(
    cell_id            = dt_cells$cell_id[sec_idx],
    sample_id          = dt_cells$sample_id[sec_idx],
    core_id            = dt_cells$core_id[sec_idx],
    cell_label         = labels[sec_idx],
    polarization_UCell = dt_cells$polarization_UCell[sec_idx],
    x                  = coords[sec_idx, 1],
    y                  = coords[sec_idx, 2],
    n_neighbors        = n_nb_total
  )

  # Add cell type proportions (prefixed)
  for (ct in celltype_order) {
    col_name <- paste0("prop_", gsub(" ", "_", ct))
    set(dt_nb, j = col_name, value = ct_prop_mat[, ct])
  }

  # Add cell type counts (prefixed)
  for (ct in celltype_order) {
    col_name <- paste0("count_", gsub(" ", "_", ct))
    set(dt_nb, j = col_name, value = ct_count_mat[, ct])
  }

  # Add all-neighbor pathway means
  for (pw in pathway_names) {
    set(dt_nb, j = paste0("nb_mean_", pw), value = pw_all_mat[, pw])
  }

  # Add per-cell-type pathway means
  for (ct in celltype_order) {
    ct_safe <- gsub(" ", "_", ct)
    for (pw in pathway_names) {
      col_name <- paste0("nb_ct_", ct_safe, "_", pw)
      set(dt_nb, j = col_name, value = ct_pw_mats[[ct]][, pw])
    }
  }

  return(dt_nb)
}

nb_list <- list()

# --- WT samples: process per-sample (whole tissue frNN) ---
wt_sample_ids <- setdiff(unique(all_cells$sample_id), "TMA")

for (sid in wt_sample_ids) {
  message(sprintf("\n  Processing WT sample %s ...", sid))
  t0 <- Sys.time()

  dt_sample <- all_cells[sample_id == sid]
  dt_nb <- compute_nb_features(dt_sample, group_label = sid)

  if (!is.null(dt_nb)) {
    nb_list[[sid]] <- dt_nb
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("    Done in %.1fs", elapsed))

  rm(dt_sample, dt_nb)
  gc(verbose = FALSE)
}

# --- TMA: process PER-CORE (frNN within each core boundary) ---
dt_tma <- all_cells[sample_id == "TMA"]
tma_core_ids <- unique(dt_tma$core_id)
tma_core_ids <- tma_core_ids[!is.na(tma_core_ids)]

message(sprintf("\n  Processing TMA: %d cores individually ...", length(tma_core_ids)))

for (cid in tma_core_ids) {
  message(sprintf("\n  Processing TMA core %s ...", cid))
  t0 <- Sys.time()

  dt_core <- dt_tma[core_id == cid]
  dt_nb <- compute_nb_features(dt_core, group_label = paste0("TMA_core_", cid))

  if (!is.null(dt_nb)) {
    nb_list[[paste0("TMA_", cid)]] <- dt_nb
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("    Done in %.1fs", elapsed))

  rm(dt_core, dt_nb)
  gc(verbose = FALSE)
}

rm(dt_tma)
gc(verbose = FALSE)

nb_features <- rbindlist(nb_list, fill = TRUE)
rm(nb_list)
gc(verbose = FALSE)

message(sprintf("\nTotal secretory cells with neighborhood features: %s",
                format(nrow(nb_features), big.mark = ",")))
message(sprintf("  WT: %s", format(nrow(nb_features[sample_id != "TMA"]), big.mark = ",")))
message(sprintf("  TMA: %s", format(nrow(nb_features[sample_id == "TMA"]), big.mark = ",")))

saveRDS(nb_features, nb_features_path)
message("  Saved neighborhood_features.rds")

# Remove cells with NA polarization
nb_features <- nb_features[!is.na(polarization_UCell)]
message(sprintf("Secretory cells with valid polarization: %s",
                format(nrow(nb_features), big.mark = ",")))

# Split into WT and TMA
nb_wt  <- nb_features[sample_id != "TMA"]
nb_tma <- nb_features[sample_id == "TMA"]

message(sprintf("  WT with valid polarization: %s", format(nrow(nb_wt), big.mark = ",")))
message(sprintf("  TMA with valid polarization: %s", format(nrow(nb_tma), big.mark = ",")))

# Save neighborhood features as CSV (summary)
fwrite(nb_features[, .(cell_id, sample_id, core_id, cell_label, polarization_UCell,
                        n_neighbors)],
       file.path(out_path, "neighborhood_summary.csv"))

# ============================================================================
# PART 3: GAM Fitting (WT Primary Analysis Only)
# ============================================================================

message("\n", strrep("=", 70))
message("PART 3: GAM Fitting — Primary Analysis (WT only, 8 independent samples)")
message("  Families: betar for proportions, nb for counts, gaussian for pathways")
message("  Per-sample GAMs as primary evidence of consistency")
message(strrep("=", 70))

# Subsample for GAM fitting (stratified by WT sample ONLY)
TARGET_N <- 50000
n_wt <- nrow(nb_wt)

if (n_wt > TARGET_N) {
  message(sprintf("Subsampling WT %d -> %d (stratified by sample)", n_wt, TARGET_N))
  nb_wt[, .N, by = sample_id] |> print()

  # Proportional stratified sampling
  sample_props <- nb_wt[, .N, by = sample_id]
  sample_props[, n_draw := pmax(50, round(TARGET_N * N / n_wt))]

  idx_list <- list()
  for (i in seq_len(nrow(sample_props))) {
    sid <- sample_props$sample_id[i]
    n_draw <- sample_props$n_draw[i]
    pool <- which(nb_wt$sample_id == sid)
    idx_list[[sid]] <- pool[sample(length(pool), min(n_draw, length(pool)))]
  }
  sub_idx <- unlist(idx_list)
  gam_data <- nb_wt[sub_idx]
  message(sprintf("  Subsampled to %d WT cells", nrow(gam_data)))
} else {
  gam_data <- copy(nb_wt)
  message(sprintf("  Using all %d WT cells (below target)", nrow(gam_data)))
}

# Define feature columns by type
ct_prop_cols  <- paste0("prop_", gsub(" ", "_", celltype_order))
ct_count_cols <- paste0("count_", gsub(" ", "_", celltype_order))
pw_all_cols   <- paste0("nb_mean_", pathway_names)

# Per-cell-type pathway columns
ct_pw_cols <- character(0)
for (ct in celltype_order) {
  ct_safe <- gsub(" ", "_", ct)
  for (pw in pathway_names) {
    ct_pw_cols <- c(ct_pw_cols, paste0("nb_ct_", ct_safe, "_", pw))
  }
}

all_feature_cols <- c(ct_prop_cols, ct_count_cols, pw_all_cols, ct_pw_cols)

# Feature type lookup
feature_type_map <- c(
  setNames(rep("proportion", length(ct_prop_cols)), ct_prop_cols),
  setNames(rep("count", length(ct_count_cols)), ct_count_cols),
  setNames(rep("pathway", length(pw_all_cols)), pw_all_cols),
  setNames(rep("ct_pathway", length(ct_pw_cols)), ct_pw_cols)
)

# Labels for features
feature_labels <- c(
  setNames(celltype_order, ct_prop_cols),
  setNames(paste0(celltype_order, " (count)"), ct_count_cols),
  setNames(gsub("^pathway_", "", pathway_names) |>
             gsub("_", " ", x = _) |>
             tools::toTitleCase(),
           pw_all_cols)
)
# Per-cell-type pathway labels
for (ct in celltype_order) {
  ct_safe <- gsub(" ", "_", ct)
  for (pw in pathway_names) {
    col_name <- paste0("nb_ct_", ct_safe, "_", pw)
    pw_label <- gsub("^pathway_", "", pw) |> gsub("_", " ", x = _) |> tools::toTitleCase()
    feature_labels[col_name] <- paste0(ct, " - ", pw_label)
  }
}

# Determine GAM family for each feature type
get_gam_family <- function(feat_type) {
  switch(feat_type,
    proportion = betar(link = "logit"),
    count      = nb(),
    pathway    = gaussian(),
    ct_pathway = gaussian(),
    gaussian()
  )
}

# Fit GAMs
message(sprintf("\nFitting GAMs for %d features (WT primary) ...", length(all_feature_cols)))
message(sprintf("  %d proportion, %d count, %d pathway, %d per-cell-type pathway",
                length(ct_prop_cols), length(ct_count_cols),
                length(pw_all_cols), length(ct_pw_cols)))

gam_results <- list()

for (feat in all_feature_cols) {
  y <- gam_data[[feat]]
  x <- gam_data$polarization_UCell
  feat_type <- feature_type_map[feat]

  # Skip if all NA or zero variance
  valid <- !is.na(y) & !is.na(x)
  n_valid <- sum(valid)
  frac_na <- 1 - n_valid / length(y)

  # Skip features with >50% NA

  if (frac_na > 0.5) {
    message(sprintf("  Skipping %s (%.0f%% NA, >50%%)", feat, frac_na * 100))
    next
  }

  if (n_valid < 200 || var(y[valid], na.rm = TRUE) < 1e-12) {
    message(sprintf("  Skipping %s (insufficient data or zero variance)", feat))
    next
  }

  fit_y <- y[valid]
  fit_x <- x[valid]

  # For betar family: clamp proportions to (0,1) exclusive
  if (feat_type == "proportion") {
    fit_y <- pmin(pmax(fit_y, 0.001), 0.999)
  }

  fit_df <- data.frame(y = fit_y, x = fit_x)
  fam <- get_gam_family(feat_type)

  tryCatch({
    fit <- gam(y ~ s(x, bs = "tp", k = 20), data = fit_df,
               family = fam, method = "REML")

    # Extract smooth term significance
    sm <- summary(fit)
    p_val <- sm$s.table[1, "p-value"]
    edf   <- sm$s.table[1, "edf"]
    r_sq  <- sm$r.sq
    dev_explained <- sm$dev.expl

    # Predict along gradient
    x_seq <- seq(min(fit_df$x), max(fit_df$x), length.out = 200)
    pred <- predict(fit, newdata = data.frame(x = x_seq), se.fit = TRUE,
                    type = "response")

    gam_results[[feat]] <- list(
      feature     = feat,
      label       = feature_labels[feat],
      feature_type = feat_type,
      family_used = as.character(fam$family),
      fit         = fit,
      p_value     = p_val,
      edf         = edf,
      r_sq        = r_sq,
      dev_expl    = dev_explained,
      n_cells     = nrow(fit_df),
      pred_df     = data.frame(
        polarization = x_seq,
        fitted       = as.numeric(pred$fit),
        se           = as.numeric(pred$se.fit),
        lower        = as.numeric(pred$fit - 1.96 * pred$se.fit),
        upper        = as.numeric(pred$fit + 1.96 * pred$se.fit)
      )
    )

    message(sprintf("  %s [%s/%s]: edf=%.1f, p=%.2e, dev=%.1f%%",
                    feat, feat_type, as.character(fam$family),
                    edf, p_val, dev_explained * 100))
  }, error = function(e) {
    message(sprintf("  Error fitting %s: %s", feat, e$message))
  })
}

message(sprintf("\nSuccessfully fit %d / %d GAMs (WT primary)", length(gam_results),
                length(all_feature_cols)))

# --- Per-sample GAMs as primary test (ALL successfully fit features) ----------
message("\nFitting per-sample GAMs for all fitted features (WT only) ...")
message("  Per-sample consistency is the primary evidence.")

wt_sample_ids_gam <- unique(gam_data$sample_id)
per_sample_results <- list()

for (feat in names(gam_results)) {
  feat_type <- feature_type_map[feat]
  fam <- get_gam_family(feat_type)

  # Get direction of pooled fit (sign of slope at midpoint vs endpoints)
  pooled_pred <- gam_results[[feat]]$pred_df
  pooled_direction <- sign(pooled_pred$fitted[nrow(pooled_pred)] - pooled_pred$fitted[1])

  for (sid in wt_sample_ids_gam) {
    sub <- gam_data[sample_id == sid]
    y <- sub[[feat]]
    x <- sub$polarization_UCell
    valid <- !is.na(y) & !is.na(x)
    if (sum(valid) < 100) next

    fit_y <- y[valid]
    fit_x <- x[valid]

    # Clamp proportions for betar
    if (feat_type == "proportion") {
      fit_y <- pmin(pmax(fit_y, 0.001), 0.999)
    }

    fit_df <- data.frame(y = fit_y, x = fit_x)

    tryCatch({
      fit_s <- gam(y ~ s(x, bs = "tp", k = 10), data = fit_df,
                   family = fam, method = "REML")
      sm_s <- summary(fit_s)

      # Determine direction from per-sample fit
      x_seq_s <- seq(min(fit_df$x), max(fit_df$x), length.out = 50)
      pred_s <- predict(fit_s, newdata = data.frame(x = x_seq_s), type = "response")
      sample_direction <- sign(pred_s[length(pred_s)] - pred_s[1])

      per_sample_results[[paste0(feat, "__", sid)]] <- data.table(
        feature          = feat,
        feature_type     = feat_type,
        sample_id        = sid,
        p_value          = sm_s$s.table[1, "p-value"],
        edf              = sm_s$s.table[1, "edf"],
        r_sq             = sm_s$r.sq,
        dev_expl         = sm_s$dev.expl,
        direction        = sample_direction,
        pooled_direction = pooled_direction,
        same_direction   = (sample_direction == pooled_direction)
      )
    }, error = function(e) NULL)
  }
}

per_sample_dt <- data.table()
if (length(per_sample_results) > 0) {
  per_sample_dt <- rbindlist(per_sample_results)
  fwrite(per_sample_dt, file.path(out_path, "gam_per_sample_wt.csv"))
  message(sprintf("  Saved per-sample WT GAMs: %d fits across %d features",
                  nrow(per_sample_dt), uniqueN(per_sample_dt$feature)))
}

# --- Aggregate per-sample consistency stats ---
per_sample_summary <- data.table()
if (nrow(per_sample_dt) > 0) {
  per_sample_summary <- per_sample_dt[, .(
    n_samples_tested      = .N,
    n_samples_significant = sum(p_value < 0.05),
    n_samples_same_direction = sum(same_direction, na.rm = TRUE)
  ), by = feature]
}

# Save GAM results summary with new columns
gam_summary <- data.table(
  feature      = sapply(gam_results, `[[`, "feature"),
  feature_type = sapply(gam_results, `[[`, "feature_type"),
  family_used  = sapply(gam_results, `[[`, "family_used"),
  label        = sapply(gam_results, `[[`, "label"),
  p_value      = sapply(gam_results, `[[`, "p_value"),
  edf          = sapply(gam_results, `[[`, "edf"),
  r_sq         = sapply(gam_results, `[[`, "r_sq"),
  dev_expl     = sapply(gam_results, `[[`, "dev_expl"),
  n_cells      = sapply(gam_results, `[[`, "n_cells"),
  significant  = sapply(gam_results, `[[`, "p_value") < 0.05,
  p_adj        = p.adjust(sapply(gam_results, `[[`, "p_value"), method = "BH")
)

# Merge per-sample consistency
if (nrow(per_sample_summary) > 0) {
  gam_summary <- merge(gam_summary, per_sample_summary, by = "feature", all.x = TRUE)
} else {
  gam_summary[, n_samples_tested := NA_integer_]
  gam_summary[, n_samples_significant := NA_integer_]
  gam_summary[, n_samples_same_direction := NA_integer_]
}

gam_summary <- gam_summary[order(p_value)]

fwrite(gam_summary, file.path(out_path, "gam_summary.csv"))
message("  Saved gam_summary.csv")

# Save full GAM results
saveRDS(gam_results, gam_results_path)
message("  Saved gam_results.rds")

# Print per-sample consistency summary
message("\nPer-sample consistency (top features by deviance explained):")
print(gam_summary[order(-dev_expl)][1:min(20, nrow(gam_summary)),
  .(feature_type, label, dev_expl, p_adj, n_samples_significant, n_samples_same_direction)])


# ============================================================================
# PART 3b: TMA Per-Core Validation (Spearman Correlations)
# ============================================================================

message("\n", strrep("=", 70))
message("PART 3b: TMA Per-Core Validation (Spearman Correlations)")
message(strrep("=", 70))

MIN_CORE_CELLS <- 100

# Features to validate: significant features from WT GAMs
sig_features <- gam_summary[p_adj < 0.05]$feature
message(sprintf("Validating %d significant WT features across TMA cores", length(sig_features)))

tma_core_ids_valid <- nb_tma[, .N, by = core_id][N >= MIN_CORE_CELLS]$core_id
message(sprintf("TMA cores with >= %d secretory cells: %d / %d",
                MIN_CORE_CELLS, length(tma_core_ids_valid),
                uniqueN(nb_tma$core_id)))

validation_results <- list()

for (cid in tma_core_ids_valid) {
  dt_core <- nb_tma[core_id == cid]

  for (feat in sig_features) {
    y <- dt_core[[feat]]
    x <- dt_core$polarization_UCell

    valid <- !is.na(y) & !is.na(x)
    n_valid <- sum(valid)
    if (n_valid < 10) next

    # Spearman correlation
    cor_test <- tryCatch(
      cor.test(x[valid], y[valid], method = "spearman", exact = FALSE),
      error = function(e) NULL
    )

    if (!is.null(cor_test)) {
      validation_results[[length(validation_results) + 1]] <- data.table(
        core_id   = cid,
        feature   = feat,
        label     = feature_labels[feat],
        rho       = cor_test$estimate,
        p_value   = cor_test$p.value,
        n_cells   = n_valid,
        direction = fifelse(cor_test$estimate > 0, "positive", "negative")
      )
    }
  }
}

if (length(validation_results) > 0) {
  validation_dt <- rbindlist(validation_results)

  # Adjust p-values within each feature
  validation_dt[, p_adj := p.adjust(p_value, method = "BH"), by = feature]
  validation_dt[, significant := p_adj < 0.05]

  # Save per-core results
  fwrite(validation_dt, file.path(out_path, "tma_validation_per_core.csv"))
  message(sprintf("  Saved tma_validation_per_core.csv (%d feature-core tests)",
                  nrow(validation_dt)))

  # Summary: replication rate per feature
  # Determine expected direction from WT GAM (sign of overall Spearman on WT data)
  wt_direction <- data.table(feature = character(), wt_direction = character(), wt_rho = numeric())
  for (feat in sig_features) {
    y <- nb_wt[[feat]]
    x <- nb_wt$polarization_UCell
    valid <- !is.na(y) & !is.na(x)
    if (sum(valid) > 100) {
      rho_wt <- cor(x[valid], y[valid], method = "spearman")
      wt_direction <- rbind(wt_direction, data.table(
        feature = feat,
        wt_direction = fifelse(rho_wt > 0, "positive", "negative"),
        wt_rho = rho_wt
      ))
    }
  }

  validation_summary <- validation_dt[, .(
    n_cores_tested   = .N,
    n_significant    = sum(significant),
    n_concordant     = NA_integer_,
    mean_rho         = mean(rho, na.rm = TRUE),
    median_rho       = median(rho, na.rm = TRUE)
  ), by = .(feature, label)]

  # Add concordance (sig + same direction as WT)
  validation_summary <- merge(validation_summary, wt_direction, by = "feature", all.x = TRUE)

  for (i in seq_len(nrow(validation_summary))) {
    feat <- validation_summary$feature[i]
    expected_dir <- validation_summary$wt_direction[i]
    if (!is.na(expected_dir)) {
      validation_summary$n_concordant[i] <- validation_dt[
        feature == feat & significant == TRUE & direction == expected_dir, .N
      ]
    }
  }

  validation_summary[, replication_rate := n_concordant / n_cores_tested]
  validation_summary <- validation_summary[order(-replication_rate)]

  fwrite(validation_summary, file.path(out_path, "tma_validation_summary.csv"))
  message("  Saved tma_validation_summary.csv")

  # Print summary
  message("\nTMA Validation Summary:")
  message(sprintf("  Features tested: %d", nrow(validation_summary)))
  message(sprintf("  Cores used: %d", length(tma_core_ids_valid)))

  well_replicated <- validation_summary[replication_rate >= 0.5]
  message(sprintf("  Features replicated in >= 50%% of cores: %d / %d",
                  nrow(well_replicated), nrow(validation_summary)))

  message("\nTop replicated features:")
  print(validation_summary[1:min(10, nrow(validation_summary)),
                            .(label, n_cores_tested, n_concordant,
                              replication_rate, mean_rho, wt_rho)])
} else {
  validation_dt <- data.table()
  validation_summary <- data.table()
  message("  No TMA cores with sufficient secretory cells for validation.")
}


# ============================================================================
# PART 4: Inflection Point Detection (WT Primary)
# ============================================================================

message("\n", strrep("=", 70))
message("PART 4: Inflection Point Detection (WT Primary)")
message(strrep("=", 70))

inflection_list <- list()

for (feat in names(gam_results)) {
  res <- gam_results[[feat]]
  pred_df <- res$pred_df

  x_vals <- pred_df$polarization
  y_vals <- pred_df$fitted

  # Numerical second derivative (central differences)
  dx <- diff(x_vals)
  dy <- diff(y_vals)
  first_deriv <- dy / dx

  # Second derivative
  dx2 <- (dx[-length(dx)] + dx[-1]) / 2
  d2y <- diff(first_deriv) / dx2

  x_mid <- x_vals[-c(1, length(x_vals))]

  # Find zero-crossings of second derivative
  sign_changes <- which(diff(sign(d2y)) != 0)

  if (length(sign_changes) > 0) {
    # Interpolate exact zero-crossing positions
    for (sc in sign_changes) {
      # Linear interpolation
      x1 <- x_mid[sc]
      x2 <- x_mid[sc + 1]
      y1 <- d2y[sc]
      y2 <- d2y[sc + 1]
      x_inflection <- x1 - y1 * (x2 - x1) / (y2 - y1)

      # Get fitted value at inflection
      idx_closest <- which.min(abs(x_vals - x_inflection))
      y_inflection <- y_vals[idx_closest]

      inflection_list[[length(inflection_list) + 1]] <- data.table(
        feature         = feat,
        label           = res$label,
        inflection_x    = x_inflection,
        inflection_y    = y_inflection,
        p_value         = res$p_value,
        feature_type    = feature_type_map[feat]
      )
    }
  }
}

if (length(inflection_list) > 0) {
  inflections <- rbindlist(inflection_list)
  # Keep only inflection points for significant features
  sig_features_infl <- gam_summary[p_adj < 0.05]$feature
  inflections <- inflections[feature %in% sig_features_infl]

  # For features with multiple inflection points, keep the one with largest
  # second-derivative magnitude change (most prominent)
  inflections[, rank := frank(-abs(inflection_x - median(inflection_x))),
              by = feature]
  inflections_primary <- inflections[rank == 1]

  fwrite(inflections, file.path(out_path, "inflection_points_all.csv"))
  fwrite(inflections_primary, file.path(out_path, "inflection_points_primary.csv"))
  message(sprintf("  Found %d inflection points (%d primary) across %d features",
                  nrow(inflections), nrow(inflections_primary),
                  uniqueN(inflections$feature)))
} else {
  inflections <- data.table()
  inflections_primary <- data.table()
  message("  No inflection points detected.")
}


# ============================================================================
# PART 5: Figures
# ============================================================================

message("\n", strrep("=", 70))
message("PART 5: Figures")
message(strrep("=", 70))

# --------------------------------------------------------------------------
# Figure 5a: Succession Ribbon (stacked area of cell type proportions)
#            Primary (WT) data only
# --------------------------------------------------------------------------

message("\nFigure 5a: Succession ribbon (Primary — WT) ...")

# Bin polarization into 100 bins (WT only)
nb_wt[, pol_bin := cut(polarization_UCell,
                        breaks = 100,
                        labels = FALSE,
                        include.lowest = TRUE)]
pol_breaks <- seq(min(nb_wt$polarization_UCell, na.rm = TRUE),
                  max(nb_wt$polarization_UCell, na.rm = TRUE),
                  length.out = 101)
pol_mids <- (pol_breaks[-1] + pol_breaks[-101]) / 2

# Compute mean proportions per bin
ribbon_data <- list()
for (ct in celltype_order) {
  col_name <- paste0("prop_", gsub(" ", "_", ct))
  bin_means <- nb_wt[, .(mean_prop = mean(get(col_name), na.rm = TRUE)),
                     by = pol_bin][order(pol_bin)]
  ribbon_data[[ct]] <- data.table(
    celltype    = ct,
    pol_bin     = bin_means$pol_bin,
    polarization = pol_mids[bin_means$pol_bin],
    proportion  = bin_means$mean_prop
  )
}
ribbon_dt <- rbindlist(ribbon_data)
ribbon_dt[, celltype := factor(celltype, levels = rev(celltype_order))]

p_5a <- ggplot(ribbon_dt, aes(x = polarization, y = proportion, fill = celltype)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = ref_palette, name = "Cell type") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Polarization score (SecA \u2192 SecB)",
       y = "Neighbor proportion",
       title = "Niche succession along secretory polarization",
       subtitle = "Primary (WT, n=8 samples)") +
  theme_lab() +
  theme(legend.key.size = unit(0.6, "lines"),
        legend.text = element_text(size = 6))

ggsave(file.path(fig_path, "fig5a_succession_ribbon.pdf"),
       p_5a, width = 10, height = 5)
message("  Saved fig5a_succession_ribbon.pdf")


# --------------------------------------------------------------------------
# Figure 5b: GAM Curve Panels (one per cell type) — Primary (WT)
# --------------------------------------------------------------------------

message("\nFigure 5b: GAM curve panels (Primary — WT) ...")

ct_gam_plots <- list()
for (ct in celltype_order) {
  feat <- paste0("prop_", gsub(" ", "_", ct))
  if (!feat %in% names(gam_results)) next

  res <- gam_results[[feat]]
  pred_df <- res$pred_df
  sig_label <- ifelse(res$p_value < 0.001, "***",
                      ifelse(res$p_value < 0.01, "**",
                             ifelse(res$p_value < 0.05, "*", "ns")))

  ct_color <- ref_palette[ct]

  ct_gam_plots[[ct]] <- ggplot(pred_df, aes(x = polarization)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ct_color, alpha = 0.25) +
    geom_line(aes(y = fitted), color = ct_color, linewidth = 0.8) +
    labs(title = ct,
         subtitle = sprintf("p=%s %s", formatC(res$p_value, format = "e", digits = 1),
                            sig_label),
         x = NULL, y = "Proportion") +
    theme_lab(base_size = 7) +
    theme(plot.title = element_text(size = 7),
          plot.subtitle = element_text(size = 5.5, color = "gray40"))
}

if (length(ct_gam_plots) > 0) {
  p_5b <- wrap_plots(ct_gam_plots, ncol = 6) +
    plot_annotation(
      title = "GAM smooth fits: cell type proportions along polarization",
      subtitle = "Primary (WT, n=8 samples)",
      theme = theme_lab()
    )

  ggsave(file.path(fig_path, "fig5b_gam_celltype_panels.pdf"),
         p_5b, width = 14, height = 9)
  message("  Saved fig5b_gam_celltype_panels.pdf")
}


# --------------------------------------------------------------------------
# Figure 5c: Pathway GAM Heatmap — Primary (WT)
# --------------------------------------------------------------------------

message("\nFigure 5c: Pathway GAM heatmap (Primary — WT) ...")

# Build matrix: pathways x 100 polarization bins
pw_gam_feats <- intersect(pw_all_cols, names(gam_results))

if (length(pw_gam_feats) > 0) {
  # Use 100 evenly spaced points along gradient
  x_range <- range(gam_data$polarization_UCell, na.rm = TRUE)
  x_grid <- seq(x_range[1], x_range[2], length.out = 100)

  pw_heatmap_mat <- matrix(NA_real_, nrow = length(pw_gam_feats), ncol = 100,
                           dimnames = list(pw_gam_feats, NULL))

  for (feat in pw_gam_feats) {
    pred <- predict(gam_results[[feat]]$fit,
                    newdata = data.frame(x = x_grid),
                    type = "response")
    pw_heatmap_mat[feat, ] <- as.numeric(pred)
  }

  # Z-score per row
  pw_z <- t(scale(t(pw_heatmap_mat)))

  # Clean row names
  row_labels <- feature_labels[rownames(pw_z)]

  # Column annotation: SecA/Transition/SecB zones
  pol_tertiles <- quantile(x_grid, probs = c(1/3, 2/3))
  zone <- ifelse(x_grid <= pol_tertiles[1], "SecA-like",
                 ifelse(x_grid >= pol_tertiles[2], "SecB-like", "Transition"))

  col_ha <- HeatmapAnnotation(
    Zone = zone,
    col = list(Zone = c("SecA-like" = "#E6A141",
                        "Transition" = "#C08E48",
                        "SecB-like" = "#9A7D55")),
    annotation_name_side = "left",
    annotation_legend_param = list(title = "Gradient zone"),
    show_legend = TRUE
  )

  col_fun <- colorRamp2(
    seq(-2, 2, length.out = 11),
    rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )

  ht <- Heatmap(
    pw_z,
    name = "Z-score",
    col = col_fun,
    cluster_rows = TRUE,
    clustering_method_rows = "ward.D2",
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_labels = row_labels,
    row_names_gp = gpar(fontsize = 7),
    top_annotation = col_ha,
    column_title = "Pathway activity along SecA \u2192 SecB gradient (Primary — WT)",
    column_title_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(title = "Z-score"),
    use_raster = FALSE,
    raster_quality = 3,
    width = unit(12, "cm"),
    height = unit(10, "cm"),
    rect_gp = gpar(col = NA)
  )

  pdf(file.path(fig_path, "fig5c_pathway_gam_heatmap.pdf"),
      width = 8, height = 7)
  draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  message("  Saved fig5c_pathway_gam_heatmap.pdf")
}


# --------------------------------------------------------------------------
# Figure 5d: Inflection Point Summary — Primary (WT)
# --------------------------------------------------------------------------

message("\nFigure 5d: Inflection point summary (Primary — WT) ...")

if (nrow(inflections_primary) > 0) {
  infl_plot <- copy(inflections_primary)
  infl_plot[, label := factor(label, levels = label[order(inflection_x)])]

  # Color by feature type
  infl_plot[, color := fifelse(feature_type == "proportion",
                               ref_palette[as.character(label)],
                               "#4A4A4A")]
  # For pathways, use a neutral color
  infl_colors <- infl_plot$color
  names(infl_colors) <- as.character(infl_plot$label)

  p_5d <- ggplot(infl_plot, aes(x = inflection_x, y = label, color = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(size = 2.5) +
    scale_color_manual(values = infl_colors, guide = "none") +
    labs(x = "Polarization score at inflection point",
         y = NULL,
         title = "Inflection points along polarization gradient",
         subtitle = "Primary (WT) — where each neighborhood feature changes most rapidly") +
    theme_lab() +
    theme(axis.text.y = element_text(size = 6))

  ggsave(file.path(fig_path, "fig5d_inflection_summary.pdf"),
         p_5d, width = 7, height = max(4, nrow(infl_plot) * 0.25 + 1))
  message("  Saved fig5d_inflection_summary.pdf")
}


# --------------------------------------------------------------------------
# Figure 5e: Pioneer vs Climax comparison — Primary (WT)
# --------------------------------------------------------------------------

message("\nFigure 5e: Pioneer vs Climax comparison (Primary — WT) ...")

# Define SecA extreme (bottom 10%) and SecB extreme (top 10%) — WT only
pol_q10 <- quantile(nb_wt$polarization_UCell, probs = 0.10, na.rm = TRUE)
pol_q90 <- quantile(nb_wt$polarization_UCell, probs = 0.90, na.rm = TRUE)

pioneer <- nb_wt[polarization_UCell <= pol_q10]
climax  <- nb_wt[polarization_UCell >= pol_q90]

pioneer[, zone := "SecA extreme\n(bottom 10%)"]
climax[, zone := "SecB extreme\n(top 10%)"]

pc_data <- rbind(pioneer, climax)

# Reshape to long format for cell type proportions
ct_long <- list()
for (ct in celltype_order) {
  col_name <- paste0("prop_", gsub(" ", "_", ct))
  ct_long[[ct]] <- data.table(
    zone       = pc_data$zone,
    sample_id  = pc_data$sample_id,
    celltype   = ct,
    proportion = pc_data[[col_name]]
  )
}
ct_long_dt <- rbindlist(ct_long)
ct_long_dt[, celltype := factor(celltype, levels = celltype_order)]

# Summarize per celltype x zone
ct_summary <- ct_long_dt[, .(
  mean_prop = mean(proportion, na.rm = TRUE),
  sd_prop   = sd(proportion, na.rm = TRUE),
  n         = .N
), by = .(celltype, zone)]

p_5e <- ggplot(ct_long_dt, aes(x = celltype, y = proportion, fill = zone)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.8,
               position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("SecA extreme\n(bottom 10%)" = "#E6A141",
                                "SecB extreme\n(top 10%)"     = "#9A7D55"),
                    name = "Gradient zone") +
  labs(x = NULL,
       y = "Neighbor proportion",
       title = "Pioneer (SecA) vs Climax (SecB) microenvironment",
       subtitle = "Primary (WT, n=8 samples)") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(fig_path, "fig5e_pioneer_vs_climax.pdf"),
       p_5e, width = 10, height = 5)
message("  Saved fig5e_pioneer_vs_climax.pdf")

# Save comparison statistics
pc_stats <- list()
for (ct in celltype_order) {
  col_name <- paste0("prop_", gsub(" ", "_", ct))
  x1 <- pioneer[[col_name]]
  x2 <- climax[[col_name]]
  if (length(x1) > 5 && length(x2) > 5) {
    wt <- wilcox.test(x1, x2, conf.int = FALSE)
    pc_stats[[ct]] <- data.table(
      celltype      = ct,
      pioneer_mean  = mean(x1, na.rm = TRUE),
      climax_mean   = mean(x2, na.rm = TRUE),
      diff          = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
      p_value       = wt$p.value
    )
  }
}
if (length(pc_stats) > 0) {
  pc_stats_dt <- rbindlist(pc_stats)
  pc_stats_dt[, p_adj := p.adjust(p_value, method = "BH")]
  fwrite(pc_stats_dt, file.path(out_path, "pioneer_vs_climax_stats.csv"))
  message("  Saved pioneer_vs_climax_stats.csv")
}


# --------------------------------------------------------------------------
# Figure 5f: Spatial polarization map (representative WT sample)
# --------------------------------------------------------------------------

message("\nFigure 5f: Spatial polarization map (Primary — WT) ...")

# Pick the whole-tissue sample with the most secretory cells
sec_counts <- nb_wt[, .N, by = sample_id][order(-N)]
if (nrow(sec_counts) > 0) {
  rep_sample <- sec_counts$sample_id[1]
  message(sprintf("  Representative sample: %s (%d secretory cells)",
                  rep_sample, sec_counts$N[1]))

  # Get all cells from this sample for background
  bg_cells <- all_cells[sample_id == rep_sample]

  # Secretory cells from this sample with polarization
  sec_cells <- nb_wt[sample_id == rep_sample]

  # Background: non-secretory cells in light grey
  p_5f <- ggplot() +
    spatial_point_layer(
      data = bg_cells[!cell_label %in% secretory_types],
      mapping = aes(x = x, y = y),
      color = "grey90", size = 0.05, alpha = 0.3
    ) +
    spatial_point_layer(
      data = sec_cells,
      mapping = aes(x = x, y = y, color = polarization_UCell),
      size = 0.15, alpha = 0.8
    ) +
    scale_color_gradientn(
      colors = expr_spatial,
      name = "Polarization\n(SecA \u2192 SecB)"
    ) +
    coord_fixed() +
    labs(title = sprintf("Secretory polarization gradient (%s)",
                         gsub("sfe_", "", rep_sample)),
         subtitle = "Primary (WT)",
         x = "X (\u00b5m)", y = "Y (\u00b5m)") +
    theme_lab() +
    theme(axis.text = element_text(size = 5))

  ggsave(file.path(fig_path, "fig5f_spatial_polarization.pdf"),
         p_5f, width = 8, height = 7)
  message("  Saved fig5f_spatial_polarization.pdf")

  # Also save version with contours showing gradient direction
  p_5f_contour <- p_5f +
    geom_density_2d(
      data = sec_cells[polarization_UCell > median(sec_cells$polarization_UCell,
                                                    na.rm = TRUE)],
      mapping = aes(x = x, y = y),
      color = "#9A7D55", linewidth = 0.3, alpha = 0.6, bins = 6
    ) +
    geom_density_2d(
      data = sec_cells[polarization_UCell < median(sec_cells$polarization_UCell,
                                                    na.rm = TRUE)],
      mapping = aes(x = x, y = y),
      color = "#E6A141", linewidth = 0.3, alpha = 0.6, bins = 6,
      linetype = "dashed"
    ) +
    labs(subtitle = "Primary (WT) — Contours: solid = SecB-enriched, dashed = SecA-enriched")

  ggsave(file.path(fig_path, "fig5f_spatial_polarization_contours.pdf"),
         p_5f_contour, width = 8, height = 7)
  message("  Saved fig5f_spatial_polarization_contours.pdf")
}


# --------------------------------------------------------------------------
# Figure 5g: TMA Validation Summary
# --------------------------------------------------------------------------

message("\nFigure 5g: TMA validation summary ...")

if (nrow(validation_summary) > 0) {

  # 5g-i: Barplot of replication rate across features
  val_plot <- copy(validation_summary)
  val_plot[, label := factor(label, levels = label[order(replication_rate)])]

  # Color by feature type
  val_plot[, feature_type := feature_type_map[feature]]

  p_5g_bar <- ggplot(val_plot, aes(x = replication_rate, y = label, fill = feature_type)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("proportion" = "#E6A141",
                                  "count" = "#CC79A7",
                                  "pathway" = "#56B4E9",
                                  "ct_pathway" = "#009E73"),
                      name = "Feature type") +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(),
                       expand = c(0, 0.01)) +
    labs(x = "Replication rate (concordant cores / total cores)",
         y = NULL,
         title = "TMA per-core validation of WT niche succession trends",
         subtitle = sprintf("Validation (TMA, %d independent cores, >= %d secretory cells each)",
                            length(tma_core_ids_valid), MIN_CORE_CELLS)) +
    theme_lab() +
    theme(axis.text.y = element_text(size = 6))

  ggsave(file.path(fig_path, "fig5g_tma_validation_replication.pdf"),
         p_5g_bar, width = 8, height = max(4, nrow(val_plot) * 0.25 + 1))
  message("  Saved fig5g_tma_validation_replication.pdf")

  # 5g-ii: Distribution of per-core Spearman rho for top features
  top_val_features <- validation_summary[replication_rate >= 0.5]$feature
  if (length(top_val_features) > 0) {
    top_val_dt <- validation_dt[feature %in% top_val_features[1:min(12, length(top_val_features))]]
    top_val_dt[, label := feature_labels[feature]]

    # Add WT reference direction
    top_val_dt <- merge(top_val_dt, wt_direction[, .(feature, wt_rho)],
                        by = "feature", all.x = TRUE)

    p_5g_rho <- ggplot(top_val_dt, aes(x = rho, y = label)) +
      geom_vline(xintercept = 0, linetype = "solid", color = "gray70") +
      geom_boxplot(outlier.size = 0.5, linewidth = 0.3, alpha = 0.6,
                   fill = "grey90") +
      geom_point(aes(color = significant), size = 1, alpha = 0.7,
                 position = position_jitter(height = 0.15, width = 0)) +
      # Mark WT direction
      geom_point(aes(x = wt_rho), shape = 18, size = 3, color = "#D55E00") +
      scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "gray60"),
                         name = "Core sig.\n(FDR < 0.05)") +
      labs(x = "Spearman rho (per TMA core)",
           y = NULL,
           title = "Per-core Spearman correlations (top replicated features)",
           subtitle = "Diamond = WT pooled rho; Validation (TMA)") +
      theme_lab() +
      theme(axis.text.y = element_text(size = 6.5))

    ggsave(file.path(fig_path, "fig5g_tma_validation_rho_dist.pdf"),
           p_5g_rho, width = 8, height = max(4, length(top_val_features) * 0.4 + 1))
    message("  Saved fig5g_tma_validation_rho_dist.pdf")
  }
} else {
  message("  Skipping validation figures (no validation results)")
}


# ============================================================================
# Final saves and summary
# ============================================================================

message("\n", strrep("=", 70))
message("Final summary")
message(strrep("=", 70))

# Save key objects
saveRDS(inflections_primary, file.path(out_path, "inflection_points.rds"))
saveRDS(gam_summary, file.path(out_path, "gam_summary.rds"))

# Print summary
message("\nPrimary analysis (WT):")
message(sprintf("  Total features tested: %d", nrow(gam_summary)))
message(sprintf("  Significant (FDR < 0.05): %d", sum(gam_summary$p_adj < 0.05)))
message(sprintf("  Inflection points found: %d", nrow(inflections_primary)))

message("\nTop 10 features by deviance explained (WT primary):")
print(gam_summary[order(-dev_expl)][1:min(10, nrow(gam_summary)),
                                     .(feature_type, label, dev_expl, p_adj, edf,
                                       n_samples_significant, n_samples_same_direction)])

if (nrow(validation_summary) > 0) {
  message("\nValidation analysis (TMA):")
  message(sprintf("  Cores used: %d", length(tma_core_ids_valid)))
  message(sprintf("  Features validated: %d", nrow(validation_summary)))
  message(sprintf("  Features replicated (>= 50%% cores): %d",
                  sum(validation_summary$replication_rate >= 0.5, na.rm = TRUE)))
}

message("\nOutputs saved to:")
message(sprintf("  Tables: %s", out_path))
message(sprintf("  Figures: %s", fig_path))

message("\n=== 16b_niche_succession_gams.R complete ===\n")
