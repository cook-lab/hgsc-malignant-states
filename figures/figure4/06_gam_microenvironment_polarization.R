#!/usr/bin/env Rscript
# ============================================================================
# Figure 4G — Pooled GAM: 50um neighborhood z-score vs focal polarization
# ----------------------------------------------------------------------------
# PURPOSE: Pooled GAM (mgcv REML) of four neighborhood features (distance to
#   vasculature, matrix remodeling, cell density, hypoxia), each z-scored, vs
#   focal-cell polarization score, across whole-tissue samples.
#
# INPUTS:
#   - output_root/16b_niche_succession_gams_v2/neighborhood_features.rds
#   - SFE objects : <data_root>/<sfe_dir>/sfe_<sample>  (NN distance to vasculature)
#   - load_sfe / theme_lab from spatial/00_setup/00_setup.R
#
# OUTPUTS:
#   - figures_dir/gam_microenvironment_polarization.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 4G.
#
# RUNTIME TIER: heavy (loads each whole-tissue SFE for NN distance; fits GAMs).
# ============================================================================

# --- Shared spatial setup (provides config, load_sfe, theme_lab, seed) -------
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure4"
.setup_path <- normalizePath(file.path(.fig_dir, "..", "..", "spatial", "00_setup", "00_setup.R"),
                             mustWork = FALSE)
if (!file.exists(.setup_path)) .setup_path <- "spatial/00_setup/00_setup.R"
source(.setup_path)

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2); library(scales); library(RANN)
})

set.seed(CFG$seed)

# --- Paths -------------------------------------------------------------------
fig_dir <- path.expand(CFG$paths$figures_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
OUT_STEM <- file.path(fig_dir, "gam_microenvironment_polarization")

# --- Parameters --------------------------------------------------------------
N_CAP_PER_SAMPLE <- 5000
N_GRID           <- 200
K_POOL           <- 10
vascular_types   <- c("Pericyte", "Endothelial")

# --- Load neighborhood features ----------------------------------------------
message("Loading neighborhood_features.rds ...")
nf <- readRDS(cfg_path("output_root", "16b_niche_succession_gams_v2",
                       "neighborhood_features.rds"))
setDT(nf)
nf <- nf[sample_id != "TMA" & !is.na(polarization_UCell)]
message("  WT cells with polarization: ", nrow(nf))

nf_sub <- nf[, { if (.N <= N_CAP_PER_SAMPLE) .SD else .SD[sample(.N, N_CAP_PER_SAMPLE)] },
             by = sample_id]
message("  equal-cap total: ", nrow(nf_sub))
rm(nf); gc(verbose = FALSE)

# --- Distance to vasculature from SFE objects --------------------------------
message("\nComputing distance to vasculature ...")
dist_list <- list()
for (sfe_name in unique(nf_sub$sample_id)) {
  message(sprintf("  Loading %s ...", sfe_name))
  sfe <- load_sfe(sfe_name)
  coords <- spatialCoords(sfe)
  cd <- data.table(cell_id = colnames(sfe), cell_label = sfe$cell_label,
                   x = coords[, 1], y = coords[, 2])
  vasc <- cd[cell_label %in% vascular_types]
  if (nrow(vasc) < 10) { rm(sfe); gc(verbose = FALSE); next }
  vasc_coords <- as.matrix(vasc[, .(x, y)])
  epi <- cd[cell_id %in% nf_sub[sample_id == sfe_name]$cell_id]
  if (nrow(epi) > 0) {
    nn <- nn2(vasc_coords, as.matrix(epi[, .(x, y)]), k = 1)
    dist_list[[length(dist_list) + 1]] <- data.table(
      cell_id = epi$cell_id, sample_id = sfe_name,
      dist_to_vascular = as.numeric(nn$nn.dists))
  }
  rm(sfe, cd); gc(verbose = FALSE)
}
dist_dt <- rbindlist(dist_list)
nf_sub <- merge(nf_sub, dist_dt[, .(cell_id, sample_id, dist_to_vascular)],
                by = c("cell_id", "sample_id"), all.x = TRUE)

# --- Derived features --------------------------------------------------------
area_mm2 <- pi * 50^2 / 1e6
nf_sub[, cell_density := n_neighbors / area_mm2]

features <- c("dist_to_vascular", "nb_mean_pathway_matrix_remodeling",
              "cell_density", "nb_mean_pathway_hypoxia")
labels   <- c("Distance to vasculature", "Matrix remodeling",
              "Cell density", "Hypoxia")
colors   <- c("#B87A7A", "#E07850", "#E6A141", "#D14E6C")

for (feat in features) {
  v <- nf_sub[[feat]]; mu <- mean(v, na.rm = TRUE); sd_ <- sd(v, na.rm = TRUE)
  if (is.finite(sd_) && sd_ > 0) nf_sub[, (feat) := (v - mu) / sd_]
}

# --- GAM fitting -------------------------------------------------------------
fit_feature <- function(dt, feat, label) {
  d <- data.frame(y = dt[[feat]], x = dt$polarization_UCell, sid = factor(dt$sample_id))
  d <- d[!is.na(d$y) & !is.na(d$x) & is.finite(d$y), ]
  if (nrow(d) < 200 || sd(d$y) == 0) return(NULL)
  grid <- seq(min(d$x), max(d$x), length.out = N_GRID)
  fit <- tryCatch(gam(y ~ s(x, k = K_POOL), data = d, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pr <- predict(fit, newdata = data.frame(x = grid), se.fit = TRUE)
  list(pooled = data.table(polarization = grid, fitted = as.numeric(pr$fit),
                           lower = as.numeric(pr$fit - 1.96 * pr$se.fit),
                           upper = as.numeric(pr$fit + 1.96 * pr$se.fit),
                           feature = label))
}

message("\nFitting GAMs ...")
pooled_all <- list()
for (i in seq_along(features)) {
  res <- fit_feature(nf_sub, features[i], labels[i])
  if (!is.null(res)) { pooled_all[[i]] <- res$pooled; message(sprintf("  %s: OK", labels[i])) }
}
pooled <- rbindlist(pooled_all, fill = TRUE)
pooled[, feature := factor(feature, levels = labels)]
color_map <- setNames(colors, labels)

label_df <- pooled[, .SD[which.max(polarization)], by = feature]
label_df[, color := color_map[as.character(feature)]]
label_df <- label_df[order(fitted)]
min_gap <- 0.20
label_df[, fitted_orig := fitted]
for (i in 2:nrow(label_df)) {
  if (label_df$fitted[i] - label_df$fitted[i - 1] < min_gap)
    label_df$fitted[i] <- label_df$fitted[i - 1] + min_gap
}

# --- Plot --------------------------------------------------------------------
p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = pooled, aes(x = polarization, ymin = lower, ymax = upper, fill = feature),
              alpha = 0.12, show.legend = FALSE) +
  geom_line(data = pooled, aes(x = polarization, y = fitted, color = feature),
            linewidth = 0.8, show.legend = FALSE) +
  geom_segment(data = label_df[abs(fitted - fitted_orig) > 0.02],
               aes(x = polarization, xend = polarization + 0.004,
                   y = fitted_orig, yend = fitted, color = feature),
               linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = label_df, aes(x = polarization + 0.006, y = fitted,
                                 label = feature, color = feature),
            hjust = 0, size = 2.0, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  labs(x = "Focal cell polarization score (SecA to SecB)",
       y = "50µm neighborhood\nz-score") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.22))) +
  theme_lab() +
  theme(plot.margin = margin(4, 60, 4, 4),
        axis.title.y = element_text(size = 6, angle = 90, margin = margin(r = 4)),
        axis.title.x = element_text(size = 6, margin = margin(t = 4)),
        axis.text = element_text(size = 5.5))

w <- 4.25; h <- 2.5
ggsave(paste0(OUT_STEM, ".pdf"), p, width = w, height = h, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), p, width = w, height = h, dpi = 450, bg = "white")
ggsave(paste0(OUT_STEM, ".svg"), p, width = w, height = h, bg = "white")
message("\nSaved: ", OUT_STEM, ".{pdf,png,svg}\nDONE")
