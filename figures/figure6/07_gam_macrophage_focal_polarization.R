#!/usr/bin/env Rscript
# ============================================================================
# Figure 6H — Macrophage GAM: expression vs nearest-epithelial polarization
# ----------------------------------------------------------------------------
# PURPOSE
#   Pooled GAM (mgcv REML) of macrophage-specific features (z-scored) against
#   the focal (nearest-epithelial) polarization score. Features: Glycolysis
#   UCell, NF-kB UCell, CXCL10, C1QC, A2M, ICAM1, CTSL. If pathway scores are
#   absent from the 19e cache, glycolysis / NF-kB proxies are computed from
#   gene means.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/19e_gene_gams_all_celltypes/
#     tme_expression_polarization.rds  (per-cell expr + nearest_epi_polarization)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (theme_lab).
#
# OUTPUTS
#   figures_dir/figure6/gam_macrophage_focal_polarization.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 6H
# RUNTIME TIER: moderate (subsample + GAM fits)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2)
  library(ragg); library(svglite)
})

set.seed(CFG$seed)

OUT_STEM <- cfg_path("figures_dir", "figure6", "gam_macrophage_focal_polarization")

N_GRID  <- 200
K_POOL  <- 10

# ── Load tme_expression_polarization (19e cache) ────────────────────────
message("[1] Loading tme_expression_polarization.rds ...")
tme <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                        "19e_gene_gams_all_celltypes",
                        "tme_expression_polarization.rds"))
setDT(tme)

mac <- tme[cell_label == "Macrophage" & !is.na(nearest_epi_polarization)]
message(sprintf("    Macrophages with polarization: %s",
                format(nrow(mac), big.mark = ",")))
rm(tme); gc(verbose = FALSE)

# Subsample uniformly across polarization bins
n_bins    <- 50
N_PER_BIN <- 3000
mac[, pol_bin := cut(nearest_epi_polarization, breaks = n_bins, labels = FALSE)]
mac <- mac[, if (.N <= N_PER_BIN) .SD else .SD[sample(.N, N_PER_BIN)],
            by = pol_bin]
message(sprintf("[subsample] kept %s macrophages", format(nrow(mac), big.mark = ",")))

# ── Pathway scores in cache? ─────────────────────────────────────────────
has_glyc <- "pathway_glycolysis" %in% colnames(mac)
has_nfkb <- "pathway_nfkb" %in% colnames(mac)
message(sprintf("    pathway_glycolysis in data: %s", has_glyc))
message(sprintf("    pathway_nfkb in data: %s", has_nfkb))

GLYC_GENES <- c("ENO1", "LDHA", "PDK1", "PGK1", "SLC2A1", "SLC16A3")
NFKB_GENES <- c("BIRC3", "CD40", "FAS", "ICAM1", "NFKB1", "NFKB2",
                 "NFKBIA", "RELA", "TNF", "TRAF6", "VCAM1")

if (!has_glyc) {
  avail_glyc <- intersect(GLYC_GENES, colnames(mac))
  message(sprintf("    Computing glycolysis proxy from %d/%d genes: %s",
                  length(avail_glyc), length(GLYC_GENES),
                  paste(avail_glyc, collapse = ", ")))
  if (length(avail_glyc) >= 3) {
    mac[, pathway_glycolysis := rowMeans(.SD, na.rm = TRUE),
         .SDcols = avail_glyc]
  }
}

if (!has_nfkb) {
  avail_nfkb <- intersect(NFKB_GENES, colnames(mac))
  message(sprintf("    Computing NF-kB proxy from %d/%d genes: %s",
                  length(avail_nfkb), length(NFKB_GENES),
                  paste(avail_nfkb, collapse = ", ")))
  if (length(avail_nfkb) >= 3) {
    mac[, pathway_nfkb := rowMeans(.SD, na.rm = TRUE),
         .SDcols = avail_nfkb]
  }
}

# ── Feature specification ────────────────────────────────────────────────
SINGLE_GENES <- c("CXCL10", "C1QC", "A2M", "ICAM1", "CTSL")

avail <- intersect(SINGLE_GENES, colnames(mac))
message(sprintf("\n[2] Genes available: %s", paste(avail, collapse = ", ")))

setnames(mac, "pathway_glycolysis", "Glycolysis", skip_absent = TRUE)
setnames(mac, "pathway_nfkb", "NF-kB", skip_absent = TRUE)

FEATURES <- c("Glycolysis", "NF-kB", avail)
COLORS   <- c(
  "Glycolysis" = "#D14E6C",
  "NF-kB"      = "#8A5DAF",
  "CXCL10"     = "#87CEFA",
  "C1QC"       = "#6B8E23",
  "A2M"        = "#8B9B6B",
  "ICAM1"      = "#5665B6",
  "CTSL"       = "#B87A7A"
)

# ── Melt + z-score ───────────────────────────────────────────────────────
message("[3] Melting and z-scoring ...")

mac[, polarization := nearest_epi_polarization]

gene_long <- melt(mac[, c("cell_id", "polarization", "sample_id", FEATURES),
                       with = FALSE],
                   id.vars = c("cell_id", "polarization", "sample_id"),
                   variable.name = "feature", value.name = "value")
gene_long[, feature := as.character(feature)]
gene_long <- gene_long[is.finite(value)]

gene_long[, value_z := (value - mean(value, na.rm = TRUE)) /
                         sd(value, na.rm = TRUE),
           by = feature]

# ── GAM fitting ──────────────────────────────────────────────────────────
pol_range <- range(gene_long$polarization, na.rm = TRUE)
grid_pol  <- seq(pol_range[1], pol_range[2], length.out = N_GRID)

message("\n[4] Fitting pooled GAMs ...")
pooled_all <- list()
for (f in FEATURES) {
  d <- gene_long[feature == f & !is.na(value_z)]
  if (nrow(d) < 200) { message("  SKIP: ", f); next }
  fit <- tryCatch(
    gam(value_z ~ s(polarization, k = K_POOL),
        data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) { message("  FAILED: ", f); next }
  pp <- predict(fit, newdata = data.frame(polarization = grid_pol),
                se.fit = TRUE)
  pooled_all[[f]] <- data.table(
    polarization = grid_pol,
    fitted = as.numeric(pp$fit),
    lower  = as.numeric(pp$fit - 1.96 * pp$se.fit),
    upper  = as.numeric(pp$fit + 1.96 * pp$se.fit),
    feature = f
  )
  message(sprintf("  %s: OK (n=%d, dev.expl=%.3f)",
                  f, nrow(d), summary(fit)$dev.expl))
}
pooled <- rbindlist(pooled_all)
pooled[, feature := factor(feature, levels = FEATURES)]

# ── End-of-line label positions ──────────────────────────────────────────
label_df <- pooled[, .SD[which.max(polarization)], by = feature]
label_df[, color := COLORS[as.character(feature)]]
label_df <- label_df[order(fitted)]
min_gap <- diff(range(pooled$fitted)) * 0.055
label_df[, fitted_orig := fitted]
if (nrow(label_df) > 1) {
  for (i in 2:nrow(label_df)) {
    if (label_df$fitted[i] - label_df$fitted[i - 1] < min_gap) {
      label_df$fitted[i] <- label_df$fitted[i - 1] + min_gap
    }
  }
}

x_range <- range(pooled$polarization)
x_nudge <- diff(x_range) * 0.012

# ── Plot ─────────────────────────────────────────────────────────────────
message("\n[5] Plotting ...")

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = pooled,
              aes(x = polarization, ymin = lower, ymax = upper,
                  fill = feature),
              alpha = 0.10, show.legend = FALSE) +
  geom_line(data = pooled,
            aes(x = polarization, y = fitted, color = feature),
            linewidth = 0.8, show.legend = FALSE) +
  geom_segment(data = label_df[abs(fitted - fitted_orig) > 0.02],
               aes(x = polarization, xend = polarization + x_nudge * 0.6,
                   y = fitted_orig, yend = fitted, color = feature),
               linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = label_df,
            aes(x = polarization + x_nudge, y = fitted,
                label = feature, color = feature),
            hjust = 0, size = 2.2, fontface = "bold.italic",
            show.legend = FALSE) +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Focal cell polarization score (SecA to SecB)",
       y = "Macrophage expression\n(z-scored)") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.25),
                     expand = expansion(mult = c(0.02, 0.22))) +
  theme_lab() +
  theme(
    plot.margin  = margin(4, 60, 4, 4),
    axis.title.y = element_text(size = 6, angle = 90, margin = margin(r = 4)),
    axis.title.x = element_text(size = 6, margin = margin(t = 4)),
    axis.text    = element_text(size = 5.5)
  )

# ── Save ─────────────────────────────────────────────────────────────────
w <- 5.1
h <- 2.8

ggsave(paste0(OUT_STEM, ".pdf"), p, width = w, height = h, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), p, width = w, height = h, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(OUT_STEM, ".svg"), p, width = w, height = h, bg = "white")

message("\nDONE")
