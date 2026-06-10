#!/usr/bin/env Rscript
# ============================================================================
# Figure 6G — GAM: 50µm neighborhood immunomodulatory genes vs polarization
# ----------------------------------------------------------------------------
# PURPOSE
#   Pooled GAM (mgcv REML) of 50µm neighborhood mean expression (z-scored) of
#   selected immunomodulatory genes vs the focal epithelial cell polarization
#   score. WT samples only. For each epithelial cell, the neighborhood mean is
#   the average logcounts of each gene across cells within 50µm (dbscan::frNN).
#   Genes: C7, CD55, TGFB1, INHBA, IL32.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/16b_niche_succession_gams_v2/
#     neighborhood_features.rds  (per-cell polarization_UCell)
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_<sample> (load_sfe)
#   Cohort: CFG$cohort$whole_tissue (published 8 WT samples).
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe, theme_lab).
#
# OUTPUTS
#   figures_dir/figure6/gam_neighborhood_immunomod_genes_polarization.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 6G
# RUNTIME TIER: heavy (frNN neighborhood means across 8 WT SFEs)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2)
  library(dbscan); library(ragg)
})

set.seed(CFG$seed)

OUT_STEM <- cfg_path("figures_dir", "figure6",
                     "gam_neighborhood_immunomod_genes_polarization")

N_GRID  <- 200
K_POOL  <- 10
RADIUS  <- 50   # µm
N_CAP_PER_SAMPLE <- 5000

GENES  <- c("C7", "CD55", "TGFB1", "INHBA", "IL32")
COLORS <- c("C7"    = "#87CEFA",
             "CD55"  = "#D14E6C",
             "TGFB1" = "#B87A7A",
             "INHBA" = "#8A5DAF",
             "IL32"  = "#5665B6")

EPI_TYPES <- c("SecA epithelium", "Intermediate epithelium",
               "SecB epithelium", "Ciliated epithelium")

WT_SAMPLES <- paste0("sfe_", CFG$cohort$whole_tissue)

# ── Load neighborhood features (for polarization scores) ─────────────────
message("[1] Loading neighborhood features for polarization ...")
nf <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                       "16b_niche_succession_gams_v2",
                       "neighborhood_features.rds"))
setDT(nf)
nf <- nf[sample_id != "TMA" & !is.na(polarization_UCell)]
nf_key <- nf[, .(cell_id, sample_id, polarization_UCell)]
rm(nf); gc(verbose = FALSE)

# ── Compute neighborhood gene means per SFE ──────────────────────────────
message("\n[2] Computing 50µm neighborhood gene means per sample ...")

all_results <- list()

for (sname in WT_SAMPLES) {
  message(sprintf("\n  === %s ===", sname))
  sfe <- load_sfe(sname)

  avail <- intersect(GENES, rownames(sfe))
  message(sprintf("    Genes available: %s / %s", length(avail), length(GENES)))
  if (length(avail) == 0) { rm(sfe); gc(verbose = FALSE); next }

  coords <- spatialCoords(sfe)
  labels <- as.character(sfe$cell_label)
  # Rename-mismatch fix (idempotent): deposited SFE still carries the legacy
  # "Transitioning epithelium"; standardize on the SAME 'labels' vector that
  # the EPI_TYPES mask below keys on, so Intermediate cells are not dropped.
  labels[labels == "Transitioning epithelium"] <- "Intermediate epithelium"
  barcodes <- colnames(sfe)

  epi_mask <- labels %in% EPI_TYPES
  epi_ids  <- barcodes[epi_mask]

  nf_sample <- nf_key[sample_id == sname]
  epi_with_pol <- intersect(epi_ids, nf_sample$cell_id)
  message(sprintf("    Epithelial cells with polarization: %s",
                  format(length(epi_with_pol), big.mark = ",")))

  if (length(epi_with_pol) < 200) { rm(sfe); gc(verbose = FALSE); next }

  if (length(epi_with_pol) > N_CAP_PER_SAMPLE) {
    epi_with_pol <- sample(epi_with_pol, N_CAP_PER_SAMPLE)
  }

  lc <- assay(sfe, "logcounts")
  gene_mat <- as.matrix(lc[avail, , drop = FALSE])

  epi_idx <- match(epi_with_pol, barcodes)
  all_coords <- as.matrix(coords)

  message("    Running frNN ...")
  nn <- frNN(all_coords, eps = RADIUS, query = all_coords[epi_idx, , drop = FALSE],
             sort = FALSE)

  message("    Computing neighborhood means ...")
  nb_means <- matrix(NA_real_, nrow = length(epi_with_pol), ncol = length(avail),
                     dimnames = list(epi_with_pol, avail))

  for (i in seq_along(epi_with_pol)) {
    nb_idx <- nn$id[[i]]
    if (length(nb_idx) < 3) next
    nb_means[i, ] <- rowMeans(gene_mat[avail, nb_idx, drop = FALSE], na.rm = TRUE)
  }

  dt <- data.table(cell_id = epi_with_pol, sample_id = sname)
  for (g in avail) {
    dt[, (g) := nb_means[, g]]
  }

  dt <- merge(dt, nf_sample[, .(cell_id, polarization_UCell)], by = "cell_id")

  all_results[[sname]] <- dt
  message(sprintf("    Result: %d cells", nrow(dt)))

  rm(sfe, lc, gene_mat, nn, nb_means); gc(verbose = FALSE)
}

combined <- rbindlist(all_results, fill = TRUE)
message(sprintf("\n[3] Combined: %s cells across %d samples",
                format(nrow(combined), big.mark = ","), length(all_results)))

# ── Melt + z-score ───────────────────────────────────────────────────────
message("[4] Melting and z-scoring ...")
gene_long <- melt(combined[, c("cell_id", "polarization_UCell", "sample_id", GENES),
                            with = FALSE],
                   id.vars = c("cell_id", "polarization_UCell", "sample_id"),
                   variable.name = "gene", value.name = "nb_mean")
gene_long[, gene := as.character(gene)]
gene_long <- gene_long[is.finite(nb_mean)]

gene_long[, nb_z := (nb_mean - mean(nb_mean, na.rm = TRUE)) /
                      sd(nb_mean, na.rm = TRUE),
           by = gene]

# ── GAM fitting ──────────────────────────────────────────────────────────
pol_range <- range(gene_long$polarization_UCell, na.rm = TRUE)
grid_pol  <- seq(pol_range[1], pol_range[2], length.out = N_GRID)

message("\n[5] Fitting pooled GAMs ...")
pooled_all <- list()
for (g in GENES) {
  d <- gene_long[gene == g & !is.na(nb_z)]
  if (nrow(d) < 200) { message("  SKIP: ", g); next }
  fit <- tryCatch(
    gam(nb_z ~ s(polarization_UCell, k = K_POOL),
        data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) { message("  FAILED: ", g); next }
  pp <- predict(fit, newdata = data.frame(polarization_UCell = grid_pol),
                se.fit = TRUE)
  pooled_all[[g]] <- data.table(
    polarization = grid_pol,
    fitted = as.numeric(pp$fit),
    lower  = as.numeric(pp$fit - 1.96 * pp$se.fit),
    upper  = as.numeric(pp$fit + 1.96 * pp$se.fit),
    feature = g
  )
  message(sprintf("  %s: OK (n=%d)", g, nrow(d)))
}
pooled <- rbindlist(pooled_all)
pooled[, feature := factor(feature, levels = GENES)]

# ── End-of-line label positions ──────────────────────────────────────────
label_df <- pooled[, .SD[which.max(polarization)], by = feature]
label_df[, color := COLORS[as.character(feature)]]
label_df <- label_df[order(fitted)]
min_gap <- 0.20
label_df[, fitted_orig := fitted]
if (nrow(label_df) > 1) {
  for (i in 2:nrow(label_df)) {
    if (label_df$fitted[i] - label_df$fitted[i - 1] < min_gap) {
      label_df$fitted[i] <- label_df$fitted[i - 1] + min_gap
    }
  }
}

# ── Plot ─────────────────────────────────────────────────────────────────
message("\n[6] Plotting ...")

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = pooled,
              aes(x = polarization, ymin = lower, ymax = upper,
                  fill = feature),
              alpha = 0.12, show.legend = FALSE) +
  geom_line(data = pooled,
            aes(x = polarization, y = fitted, color = feature),
            linewidth = 0.8, show.legend = FALSE) +
  geom_segment(data = label_df[abs(fitted - fitted_orig) > 0.02],
               aes(x = polarization, xend = polarization + 0.004,
                   y = fitted_orig, yend = fitted, color = feature),
               linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = label_df,
            aes(x = polarization + 0.006, y = fitted,
                label = feature, color = feature),
            hjust = 0, size = 2.0, fontface = "bold.italic",
            show.legend = FALSE) +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Focal cell polarization score (SecA to SecB)",
       y = expression(atop("Neighborhood expression", "z-score (50 µm)"))) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.22))) +
  theme_lab() +
  theme(
    plot.margin = margin(4, 60, 4, 4),
    axis.title.y = element_text(size = 6, angle = 90, margin = margin(r = 4)),
    axis.title.x = element_text(size = 6, margin = margin(t = 4)),
    axis.text    = element_text(size = 5.5)
  )

# ── Save ─────────────────────────────────────────────────────────────────
w <- 5.1
h <- 2.5

ggsave(paste0(OUT_STEM, ".pdf"), p, width = w, height = h, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), p, width = w, height = h, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(OUT_STEM, ".svg"), p, width = w, height = h, bg = "white")

message("\nDONE")
