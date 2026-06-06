#!/usr/bin/env Rscript
# ============================================================================
# Figure 5G — Pooled GAM: epithelial gene expression vs polarization
# ----------------------------------------------------------------------------
# PURPOSE
#   Whole-tissue, epithelial-only pooled GAM (mgcv REML) of z-scored single-
#   gene expression against the focal-cell polarization score (SecA -> SecB).
#   Genes: CDH2, CTNNB1, ITGB5, MMP7, TGM2, ICAM1.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/19d_gene_polarization_gams/
#     epithelial_expression_polarization.rds   (19d cache)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (theme_lab).
#
# OUTPUTS
#   figures_dir/figure5/gam_epithelial_genes_polarization.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 5G
# RUNTIME TIER: moderate (subsample + GAM fits)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2)
})

set.seed(CFG$seed)

OUT_STEM <- cfg_path("figures_dir", "figure5", "gam_epithelial_genes_polarization")

N_GRID  <- 200
K_POOL  <- 10

# ── Load epithelial expression + polarization ────────────────────────────
message("[load] epithelial expression + polarization (cached, 19d)")
expr <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                         "19d_gene_polarization_gams",
                         "epithelial_expression_polarization.rds"))
setDT(expr)
expr <- expr[!is.na(polarization_UCell)]

# Subsample uniformly across polarization bins
n_bins    <- 50
N_PER_BIN <- 4000
expr[, pol_bin := cut(polarization_UCell, breaks = n_bins, labels = FALSE)]
expr <- expr[, if (.N <= N_PER_BIN) .SD else .SD[sample(.N, N_PER_BIN)],
              by = pol_bin]
message("[subsample] kept ", nrow(expr), " epithelial cells")

# ── Gene feature specification ───────────────────────────────────────────
GENES  <- c("CDH2", "CTNNB1", "ITGB5", "MMP7", "TGM2", "ICAM1")
LABELS <- c("CDH2"   = "CDH2",
             "CTNNB1" = "CTNNB1",
             "ITGB5"  = "ITGB5",
             "MMP7"   = "MMP7",
             "TGM2"   = "TGM2",
             "ICAM1"  = "ICAM1")
COLORS <- c("CDH2"   = "#E6A141",
             "CTNNB1" = "#8FBC8F",
             "ITGB5"  = "#5665B6",
             "MMP7"   = "#8A5DAF",
             "TGM2"   = "#B87A7A",
             "ICAM1"  = "#D14E6C")

stopifnot(all(GENES %in% colnames(expr)))

# ── Melt to long format + z-score ────────────────────────────────────────
message("\n[melt] reshaping to long format ...")
gene_long <- melt(expr[, c("cell_id", "polarization_UCell", "sample_id", GENES),
                        with = FALSE],
                   id.vars = c("cell_id", "polarization_UCell", "sample_id"),
                   variable.name = "gene", value.name = "expression")
gene_long[, gene := as.character(gene)]

message("[zscore] z-scoring per gene ...")
gene_long[, expr_z := (expression - mean(expression, na.rm = TRUE)) /
                        sd(expression, na.rm = TRUE),
           by = gene]

# ── GAM fitting ──────────────────────────────────────────────────────────
pol_range <- range(gene_long$polarization_UCell, na.rm = TRUE)
grid_pol  <- seq(pol_range[1], pol_range[2], length.out = N_GRID)

message("\n[gam] fitting pooled GAMs ...")
pooled_all <- list()
for (g in GENES) {
  d <- gene_long[gene == g & !is.na(expr_z)]
  label <- LABELS[g]
  fit <- tryCatch(
    gam(expr_z ~ s(polarization_UCell, k = K_POOL),
        data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) { message("  FAILED: ", label); next }
  pp <- predict(fit, newdata = data.frame(polarization_UCell = grid_pol),
                se.fit = TRUE)
  pooled_all[[g]] <- data.table(
    polarization = grid_pol,
    fitted = as.numeric(pp$fit),
    lower  = as.numeric(pp$fit - 1.96 * pp$se.fit),
    upper  = as.numeric(pp$fit + 1.96 * pp$se.fit),
    feature = label
  )
  message(sprintf("  %s: OK (n=%d)", label, nrow(d)))
}
pooled <- rbindlist(pooled_all)
pooled[, feature := factor(feature, levels = unname(LABELS))]

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
message("\n[plot] generating figure ...")

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
       y = "Epithelial expression\nz-score") +
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
