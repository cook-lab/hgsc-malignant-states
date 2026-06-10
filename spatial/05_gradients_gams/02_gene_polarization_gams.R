# ============================================================================
# 02_gene_polarization_gams.R
# ----------------------------------------------------------------------------
# PURPOSE: Precompute GAMs for all panel genes across the UCell polarization axis on epithelial cells (all WT samples). Per-gene predictions for instant querying.
#
# INPUTS:
#   - output/16b_niche_succession_gams/neighborhood_features.rds (polarization_UCell)
#   - SFEs (cfg sfe dir) logcounts
#
# OUTPUTS:
#   - output/19d_gene_polarization_gams/epithelial_expression_polarization.rds
#   - gene_gam_results.rds, gene_gam_summary.csv, query_gene_gam.R
#
# MANUSCRIPT PANEL(S): Fig 5A, Fig 5G, SF13.
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

out_dir <- cfg_path("output_root", "19d_gene_polarization_gams")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sfe_dir <- dirname(cfg_obj("sfe_tma_filtered"))
wt_names <- paste0("sfe_", CFG$cohort$whole_tissue)

epi_types <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")

# ---- Step 1: Collect all epithelial expression + polarization ----
message("=== Step 1: Loading expression from all WT samples ===")

# Load polarization scores
# NOTE: deliberate path rename. The original source read from
# "16b_niche_succession_gams_v2"; in this repo the in-repo producer
# (01_niche_succession_gams.R) writes the non-`_v2` name, so reading the
# non-`_v2` directory here is internally consistent (not a bug).
nf <- readRDS(cfg_path("output_root", "16b_niche_succession_gams", "neighborhood_features.rds"))
nf <- as.data.table(nf)
# Idempotent legacy-label rename: the deposited cache may still carry the
# legacy epithelial value "Transitioning epithelium"; downstream filtering
# keys on the canonical "Intermediate epithelium". Rename at the read point
# (harmless if already canonical).
nf[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
pol_dt <- nf[, .(cell_id, polarization_UCell, sample_id, cell_label)]
pol_dt <- pol_dt[cell_label %in% epi_types]
pol_dt <- pol_dt[!grepl("TMA", sample_id)]
message("  Epithelial cells with polarization: ", nrow(pol_dt))

rm(nf); gc()

# Get gene list from first sample
sfe_tmp <- loadHDF5SummarizedExperiment(file.path(sfe_dir, wt_names[1]))
all_genes <- rownames(sfe_tmp)
n_genes <- length(all_genes)
message("  Panel: ", n_genes, " genes")
rm(sfe_tmp); gc()

# Load expression sample by sample and combine
all_expr <- list()
for (nm in wt_names) {
  message("  Loading ", nm, "...")
  sfe <- loadHDF5SummarizedExperiment(file.path(sfe_dir, nm))

  cl <- colData(sfe)$cell_label
  # Idempotent legacy-label rename: deposited SFEs still carry the legacy
  # "Transitioning epithelium"; epi_mask keys on "Intermediate epithelium".
  cl[cl == "Transitioning epithelium"] <- "Intermediate epithelium"
  epi_mask <- cl %in% epi_types
  sfe_epi <- sfe[, epi_mask]
  
  # Get logcounts for all genes
  lc <- as.matrix(logcounts(sfe_epi))
  
  expr_dt <- as.data.table(t(lc))
  expr_dt[, cell_id := colnames(sfe_epi)]
  
  all_expr[[nm]] <- expr_dt
  rm(sfe, sfe_epi, lc); gc(verbose = FALSE)
}

expr_all <- rbindlist(all_expr, fill = TRUE)
message("  Combined: ", nrow(expr_all), " cells x ", ncol(expr_all) - 1, " genes")

# Merge with polarization
expr_all <- merge(expr_all, pol_dt, by = "cell_id")
message("  After merge: ", nrow(expr_all), " cells")

rm(all_expr, pol_dt); gc()

# Save checkpoint
message("  Saving expression checkpoint...")
saveRDS(expr_all, file.path(out_dir, "epithelial_expression_polarization.rds"))

# ---- Step 2: Fit GAMs for all genes ----
message("\n=== Step 2: Fitting GAMs for ", n_genes, " genes ===")

# Prediction grid
pol_range <- range(expr_all$polarization_UCell, na.rm = TRUE)
pred_grid <- data.frame(polarization_UCell = seq(pol_range[1], pol_range[2], length.out = 200))

# Subsample EQUALLY along polarization axis (200K total, uniform bins)
set.seed(CFG$seed)
n_target <- 200000
n_bins <- 100
expr_all[, pol_bin := cut(polarization_UCell, breaks = n_bins, labels = FALSE)]
per_bin <- ceiling(n_target / n_bins)

expr_sub <- expr_all[, {
  if (.N <= per_bin) .SD
  else .SD[sample(.N, per_bin)]
}, by = pol_bin]
expr_sub[, pol_bin := NULL]
expr_all[, pol_bin := NULL]

message("  Subsampled to ", nrow(expr_sub), " cells (uniform across ", n_bins, " polarization bins)")

results <- list()
summary_rows <- list()

for (i in seq_along(all_genes)) {
  gene <- all_genes[i]
  
  if (i %% 50 == 0) message("  Gene ", i, "/", n_genes, ": ", gene)
  
  y <- expr_sub[[gene]]
  if (is.null(y) || all(is.na(y))) next
  
  fit_dt <- data.table(y = y, polarization_UCell = expr_sub$polarization_UCell,
                        sample_id = expr_sub$sample_id)
  
  # Global GAM
  tryCatch({
    fit <- gam(y ~ s(polarization_UCell, k = 10), data = fit_dt, method = "REML")
    
    pred <- predict(fit, newdata = pred_grid, se.fit = TRUE)
    pred_df <- data.frame(
      polarization = pred_grid$polarization_UCell,
      fitted = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      lower = as.numeric(pred$fit - 1.96 * pred$se.fit),
      upper = as.numeric(pred$fit + 1.96 * pred$se.fit)
    )
    
    # Direction and deviance
    delta <- pred_df$fitted[200] - pred_df$fitted[1]
    direction <- ifelse(delta > 0, "positive", "negative")
    dev_expl <- summary(fit)$dev.expl
    r_sq <- summary(fit)$r.sq
    p_val <- summary(fit)$s.table[1, "p-value"]
    
    # Per-sample consistency
    n_same <- 0
    n_tested <- 0
    for (sid in unique(fit_dt$sample_id)) {
      sub <- fit_dt[sample_id == sid]
      if (nrow(sub) < 100) next
      n_tested <- n_tested + 1
      tryCatch({
        fit_s <- gam(y ~ s(polarization_UCell, k = 5), data = sub, method = "REML")
        pred_s <- predict(fit_s, newdata = data.frame(
          polarization_UCell = c(pol_range[1], pol_range[2])))
        delta_s <- pred_s[2] - pred_s[1]
        if (sign(delta_s) == sign(delta)) n_same <- n_same + 1
      }, error = function(e) NULL)
    }
    
    # Mean expression and percent expressing
    mean_expr <- mean(y, na.rm = TRUE)
    pct_expr <- mean(y > 0, na.rm = TRUE) * 100
    
    results[[gene]] <- list(
      gene = gene,
      pred_df = pred_df,
      direction = direction,
      dev_expl = dev_expl,
      r_sq = r_sq,
      p_value = p_val,
      mean_expr = mean_expr,
      pct_expressing = pct_expr,
      n_samples_same_direction = n_same,
      n_samples_tested = n_tested
    )
    
    summary_rows[[gene]] <- data.table(
      gene = gene, direction = direction, dev_expl = dev_expl, r_sq = r_sq,
      p_value = p_val, mean_expr = mean_expr, pct_expressing = pct_expr,
      delta = delta, n_same_dir = n_same, n_tested = n_tested
    )
    
  }, error = function(e) {
    message("    Error on ", gene, ": ", conditionMessage(e))
  })
}

# ---- Step 3: Save results ----
message("\n=== Step 3: Saving results ===")

saveRDS(results, file.path(out_dir, "gene_gam_results.rds"))

summary_dt <- rbindlist(summary_rows)
summary_dt[, p_adj := p.adjust(p_value, method = "BH")]
summary_dt <- summary_dt[order(-dev_expl)]
fwrite(summary_dt, file.path(out_dir, "gene_gam_summary.csv"))

message("  Saved ", length(results), " gene GAMs")
message("  Saved gene_gam_summary.csv")

# ---- Step 4: Quick query function ----
# Save a helper script for plotting any gene
query_script <- '
# Quick gene polarization plot
# Usage: Rscript query_gene_gam.R GENE_NAME [output.png]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript query_gene_gam.R GENE_NAME [output.png]")

gene <- toupper(args[1])
out_file <- if (length(args) >= 2) args[2] else paste0(gene, "_polarization.png")

library(data.table)
library(ggplot2)

res <- readRDS("__GAM_RDS__")

if (!gene %in% names(res)) stop(gene, " not found. Available: ", paste(head(names(res), 10), collapse=", "), "...")

r <- res[[gene]]
pred <- as.data.table(r$pred_df)

p <- ggplot(pred, aes(x = polarization, y = fitted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#7B2D3B", alpha = 0.15) +
  geom_line(color = "#7B2D3B", linewidth = 1) +
  theme_classic(base_size = 12) +
  theme(text = element_text(color = "black"), axis.text = element_text(color = "black")) +
  labs(title = paste0(gene, " expression along SecA to SecB polarization"),
       subtitle = paste0("Dev.expl=", round(r$dev_expl*100, 2), "% | ", r$direction,
                         " | ", r$n_samples_same_direction, "/", r$n_samples_tested, " samples"),
       x = "Polarization score (SecA to SecB)",
       y = paste0(gene, " (logcounts)"))

ggsave(out_file, p, width = 7, height = 4, dpi = 200, bg = "white")
message("Saved: ", out_file)
'
# Interpolate the RESOLVED absolute path at generation time so the emitted
# helper is self-contained (it never sources config, so cfg_path would be
# undefined when a user runs `Rscript query_gene_gam.R GENE`). Use a literal
# placeholder substitution (not sprintf) because the script body may contain
# `%` characters.
query_script <- gsub("__GAM_RDS__", file.path(out_dir, "gene_gam_results.rds"),
                     query_script, fixed = TRUE)
writeLines(query_script, file.path(out_dir, "query_gene_gam.R"))

message("\n=== DONE ===")
message("Results: ", out_dir)
message("Query any gene: Rscript ", file.path(out_dir, "query_gene_gam.R"), " GENE_NAME")
