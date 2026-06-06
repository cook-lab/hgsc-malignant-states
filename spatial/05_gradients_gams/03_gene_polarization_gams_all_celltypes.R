# ============================================================================
# 03_gene_polarization_gams_all_celltypes.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-gene polarization GAMs for all (non-epithelial) cell types: each TME cell inherits the polarization score of its nearest epithelial neighbor (RANN nn2), then GAMs are fit per cell type.
#
# INPUTS:
#   - SFEs (cfg sfe dir): 8 WT, cell_label + polarization_UCell + logcounts
#
# OUTPUTS:
#   - output/19e_gene_gams_all_celltypes/tme_expression_polarization.rds
#   - gene_gams_<celltype>.rds, gene_gam_summary_<celltype>.csv
#
# MANUSCRIPT PANEL(S): Fig 6H, SF14.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config (replaces bare relative output/ + sfe paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
set.seed(CFG$seed)

library(data.table)
library(SpatialFeatureExperiment)
library(SummarizedExperiment)
library(HDF5Array)
library(mgcv)
library(RANN)

out_dir <- cfg_path("output_root", "19e_gene_gams_all_celltypes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sfe_dir <- dirname(cfg_obj("sfe_tma_filtered"))
wt_names <- c("sfe_OTB_2384","sfe_OTB_2417","sfe_OTB_2432","sfe_OTB_2454",
              "sfe_OTB_2457","sfe_OTB_2461","sfe_SP24_24824","sfe_SP24_25573")

epi_types <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")
tme_types <- c("Macrophage","T cell","NK cell","Fibroblast","Endothelial",
               "Pericyte","B cell","Plasma cell","Smooth muscle")

# ---- Step 1: For each TME cell, assign nearest epithelial polarization ----
message("=== Step 1: Assigning polarization scores to TME cells ===")

all_tme <- list()

for (nm in wt_names) {
  message("  Loading ", nm, "...")
  sfe <- loadHDF5SummarizedExperiment(file.path(sfe_dir, nm))
  
  cl <- colData(sfe)$cell_label
  coords <- spatialCoords(sfe)
  pol <- colData(sfe)$polarization_UCell
  
  epi_mask <- cl %in% epi_types & !is.na(pol)
  tme_mask <- cl %in% tme_types
  
  if (sum(epi_mask) < 100 || sum(tme_mask) < 100) {
    rm(sfe); gc(verbose = FALSE)
    next
  }
  
  # Find nearest epithelial neighbor for each TME cell
  epi_coords <- coords[epi_mask, ]
  epi_pol <- pol[epi_mask]
  tme_coords <- coords[tme_mask, ]
  
  nn <- nn2(epi_coords, tme_coords, k = 1)
  
  # Only keep TME cells within 50um of an epithelial cell
  close_mask <- nn$nn.dists[, 1] <= 50
  
  # Get expression for TME cells
  tme_idx <- which(tme_mask)[close_mask]
  
  lc <- as.matrix(logcounts(sfe[, tme_idx]))
  
  tme_dt <- data.table(
    cell_id = colnames(sfe)[tme_idx],
    cell_label = cl[tme_idx],
    sample_id = sub("sfe_", "", nm),
    nearest_epi_polarization = epi_pol[nn$nn.idx[close_mask, 1]],
    nn_dist = nn$nn.dists[close_mask, 1]
  )
  
  expr_dt <- as.data.table(t(lc))
  expr_dt[, cell_id := colnames(sfe)[tme_idx]]
  
  tme_dt <- merge(tme_dt, expr_dt, by = "cell_id")
  all_tme[[nm]] <- tme_dt
  
  message("    ", sum(close_mask), " TME cells within 50um of epithelium")
  rm(sfe, lc, expr_dt, nn); gc(verbose = FALSE)
}

tme_all <- rbindlist(all_tme, fill = TRUE)
message("Total TME cells: ", nrow(tme_all))
message("Cell types: ", paste(unique(tme_all$cell_label), collapse=", "))

# Save checkpoint
message("Saving checkpoint...")
saveRDS(tme_all, file.path(out_dir, "tme_expression_polarization.rds"))

# ---- Step 2: Fit GAMs per cell type ----
all_genes <- setdiff(names(tme_all), c("cell_id","cell_label","sample_id",
                                         "nearest_epi_polarization","nn_dist"))

pol_range <- range(tme_all$nearest_epi_polarization, na.rm = TRUE)
pred_grid <- data.frame(nearest_epi_polarization = seq(pol_range[1], pol_range[2], length.out = 200))

# Only fit for cell types with enough cells
ct_counts <- tme_all[, .N, by = cell_label][order(-N)]
message("\nCell type counts:")
print(ct_counts)

tme_mtime <- file.info(file.path(out_dir, "tme_expression_polarization.rds"))$mtime
for (ct in ct_counts[N >= 1000]$cell_label) {
  safe_ct <- tolower(gsub(" ", "_", ct))
  out_rds <- file.path(out_dir, paste0("gene_gams_", safe_ct, ".rds"))
  if (file.exists(out_rds) && !is.na(tme_mtime) &&
      file.info(out_rds)$mtime >= tme_mtime) {
    message("\n=== Skipping ", ct, " (fresh gene_gams file already exists) ===")
    next
  }
  message("\n=== Fitting GAMs for ", ct, " (", ct_counts[cell_label == ct]$N, " cells) ===")

  ct_data <- tme_all[cell_label == ct]
  
  # Subsample equally along polarization axis
  set.seed(CFG$seed)
  n_target <- min(50000, nrow(ct_data))
  n_bins <- 50
  ct_data[, pol_bin := cut(nearest_epi_polarization, breaks = n_bins, labels = FALSE)]
  per_bin <- ceiling(n_target / n_bins)
  
  ct_sub <- ct_data[, {
    if (.N <= per_bin) .SD
    else .SD[sample(.N, per_bin)]
  }, by = pol_bin]
  ct_sub[, pol_bin := NULL]
  ct_data[, pol_bin := NULL]
  
  message("  Subsampled to ", nrow(ct_sub), " cells")
  
  results <- list()
  
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    if (i %% 100 == 0) message("  Gene ", i, "/", length(all_genes))
    
    y <- ct_sub[[gene]]
    if (is.null(y) || all(is.na(y)) || mean(y > 0, na.rm=TRUE) < 0.01) next
    
    tryCatch({
      fit <- gam(y ~ s(nearest_epi_polarization, k = 8), 
                 data = ct_sub, method = "REML")
      
      pred <- predict(fit, newdata = pred_grid, se.fit = TRUE)
      pred_df <- data.frame(
        polarization = pred_grid$nearest_epi_polarization,
        fitted = as.numeric(pred$fit),
        se = as.numeric(pred$se.fit),
        lower = as.numeric(pred$fit - 1.96 * pred$se.fit),
        upper = as.numeric(pred$fit + 1.96 * pred$se.fit)
      )
      
      delta <- pred_df$fitted[200] - pred_df$fitted[1]
      
      results[[gene]] <- list(
        gene = gene,
        cell_type = ct,
        pred_df = pred_df,
        direction = ifelse(delta > 0, "positive", "negative"),
        dev_expl = summary(fit)$dev.expl,
        delta = delta,
        mean_expr = mean(y, na.rm = TRUE),
        pct_expressing = mean(y > 0, na.rm = TRUE) * 100
      )
    }, error = function(e) NULL)
  }
  
  # Save per cell type
  safe_ct <- gsub(" ", "_", tolower(ct))
  saveRDS(results, file.path(out_dir, paste0("gene_gams_", safe_ct, ".rds")))
  
  # Summary CSV
  summary_rows <- lapply(results, function(r) {
    data.table(gene = r$gene, cell_type = r$cell_type, direction = r$direction,
               dev_expl = r$dev_expl, delta = r$delta, 
               mean_expr = r$mean_expr, pct_expressing = r$pct_expressing)
  })
  summary_dt <- rbindlist(summary_rows)[order(-dev_expl)]
  fwrite(summary_dt, file.path(out_dir, paste0("gene_gam_summary_", safe_ct, ".csv")))
  
  message("  Saved ", length(results), " gene GAMs for ", ct)
}

message("\n=== DONE ===")
