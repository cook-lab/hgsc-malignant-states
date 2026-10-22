# Fig 4f GAM fit from the per-cell SingleR inputs (08_fig4f_gam_features.R).
# ALL secretory cells, no subsampling. Each tissue contributes equally via cell weights
# w = (N_total / n_tissues) / n_tissue (mean weight 1). Per-feature weighted z-score; mgcv::bam(y ~ s(x, k = 10), weights,
# fREML, discrete), 200-point grid over the 1st-99th percentile of polarization (fit on all cells), +-1.96 SE, label min_gap 0.20.
# Outputs: figures_dir/figure4/gam_curves.csv, gam_labels.csv (drawn by 10_fig4f_gam_plot.py; gam_curves.csv is also the
# pooled-vs-main check of Figure S16). Manuscript panel: Fig 4f.
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
set.seed(CFG$seed)
SP <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute")
fig_dir <- cfg_path("figures_dir","figure4")
N_GRID <- 200; K_POOL <- 10
features <- c("dist_to_vascular","nb_matrix","cell_density","nb_hypoxia")
labels   <- c("Distance to vasculature","Matrix remodeling","Cell density","Hypoxia")
colors   <- c("#B87A7A","#E07850","#E6A141","#D14E6C")

nf <- fread(file.path(SP,"singler_gam_percell.csv"))
cat("per-cell rows:", nrow(nf), " samples:", uniqueN(nf$sample), "\n")
nf_sub <- nf   # all cells
nf_sub[, w := (nrow(nf_sub) / uniqueN(nf_sub$sample)) / .N, by = sample]   # equal total weight per tissue
cat("cells used:", nrow(nf_sub), " per-tissue weights:\n"); print(nf_sub[, .(n = .N, w = w[1]), by = sample])
for (f in features) { v <- nf_sub[[f]]; ok <- is.finite(v); m <- weighted.mean(v[ok], nf_sub$w[ok]); s <- sqrt(weighted.mean((v[ok]-m)^2, nf_sub$w[ok]))
  if (is.finite(s) && s > 0) nf_sub[, (f) := (v-m)/s] }

pooled_all <- list()
for (i in seq_along(features)) {
  d <- data.frame(y = nf_sub[[features[i]]], x = nf_sub$polarization, w = nf_sub$w)
  d <- d[!is.na(d$y) & !is.na(d$x) & is.finite(d$y), ]
  g <- seq(quantile(d$x, 0.01), quantile(d$x, 0.99), length.out = N_GRID)   # display 1st-99th percentile of polarization (fit uses all cells)
  fit <- bam(y ~ s(x, k = K_POOL), data = d, weights = w, method = "fREML", discrete = TRUE)
  pr <- predict(fit, newdata = data.frame(x = g), se.fit = TRUE)
  pooled_all[[i]] <- data.table(polarization = g, fitted = as.numeric(pr$fit),
    lower = as.numeric(pr$fit - 1.96*pr$se.fit),
    upper = as.numeric(pr$fit + 1.96*pr$se.fit), feature = labels[i])
  cat(sprintf("  %-26s OK  (n=%d)\n", labels[i], nrow(d)))
}
pooled <- rbindlist(pooled_all, fill = TRUE)
pooled[, feature := factor(feature, levels = labels)]
cmap <- setNames(colors, labels)
label_df <- pooled[, .SD[which.max(polarization)], by = feature]
label_df[, color := cmap[as.character(feature)]]
label_df <- label_df[order(fitted)]
min_gap <- 0.20; label_df[, fitted_orig := fitted]
for (i in 2:nrow(label_df))
  if (label_df$fitted[i] - label_df$fitted[i-1] < min_gap)
    label_df$fitted[i] <- label_df$fitted[i-1] + min_gap
fwrite(pooled,   file.path(fig_dir,"gam_curves.csv"))
fwrite(label_df, file.path(fig_dir,"gam_labels.csv"))
cat("\nendpoint (right edge) fitted values:\n")
print(label_df[, .(feature, fitted_orig = round(fitted_orig,3))])
cat("\nwrote gam_curves.csv / gam_labels.csv\nDONE\n")
