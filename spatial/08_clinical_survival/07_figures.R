# ============================================================================
# 07_figures.R
# ----------------------------------------------------------------------------
# PURPOSE: Render macrophage-niche survival figures (immune composition, exhaustion).
#
# INPUTS:
#   - output/29_macrophage_niche_survival/ niche + survival outputs
#
# OUTPUTS:
#   - output/figures/29_macrophage_niche_survival/ figures
#
# MANUSCRIPT PANEL(S): Fig 6B/6E/6F
# RUNTIME TIER: moderate
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

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(mgcv)
})

# READ the published (frozen, deposited) niche-survival caches for paper-fidelity:
# these hold the as-published values (e.g. spatial-distance n=88). Re-running the
# backend (03/05/06/06b) writes regenerated copies under output_root that drift
# slightly from the deposited snapshot (unseeded original labeling), so prefer the frozen deposit when present.
.frozen_29 <- cfg_path("data_root", "2026_final_xenium_analysis", "output", "29_macrophage_niche_survival")
OUT_DIR <- if (dir.exists(.frozen_29)) .frozen_29 else file.path(out_dir, "29_macrophage_niche_survival")
FIG_DIR <- file.path(out_dir, "29_macrophage_niche_survival", "figures")  # figures always written under output_root
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LYMPH_LBLS <- c("T cell", "NK cell", "B cell", "Plasma cell")

COL_MAC   <- "#D65146"
COL_LYMPH <- "#4A8BC9"
COL_T     <- "#2E7D4F"
COL_NK    <- "#E6A141"
COL_B     <- "#9467BD"
COL_PLASMA <- "#E377C2"
COL_UP    <- "#D65146"   # red = macrophage-enriched
COL_DOWN  <- "#4A8BC9"   # blue = lymphocyte-enriched

# ---------------------------------------------------------------------------
# 29A — density-vs-stress GAM
# ---------------------------------------------------------------------------
niche <- readRDS(file.path(OUT_DIR, "per_cell_niche_scores.rds"))
wt <- niche$wt; tma <- niche$tma
wt[, is_mac    := cell_label == "Macrophage"]
wt[, is_lymph  := cell_label %in% LYMPH_LBLS]
tma[, is_mac   := cell_label == "Macrophage"]
tma[, is_lymph := cell_label %in% LYMPH_LBLS]

fit_gam_density <- function(dt, is_col, xvar = "niche_metabolic_stress_z",
                              group_col = NULL) {
  d <- dt[, .(y = get(is_col), x = get(xvar),
              grp = if (is.null(group_col)) "all" else get(group_col))]
  d <- d[is.finite(x)]
  # Pooled
  fit_pool <- tryCatch(
    mgcv::gam(y ~ s(x, k = 10), data = d, family = binomial,
               method = "REML"),
    error = function(e) NULL)
  if (is.null(fit_pool)) return(NULL)
  x_grid <- seq(quantile(d$x, 0.02), quantile(d$x, 0.98), length.out = 200)
  pr <- predict(fit_pool, data.frame(x = x_grid), type = "link", se.fit = TRUE)
  pooled_dt <- data.table(x = x_grid,
                          fit = plogis(pr$fit) * 100,
                          lo  = plogis(pr$fit - 1.96 * pr$se.fit) * 100,
                          hi  = plogis(pr$fit + 1.96 * pr$se.fit) * 100)
  # Per-group
  per_dt <- NULL
  if (!is.null(group_col)) {
    grp_list <- lapply(unique(d$grp), function(g) {
      dg <- d[grp == g]
      if (nrow(dg) < 500 || sum(dg$y) < 20) return(NULL)
      qs <- quantile(dg$x, c(0.05, 0.95))
      fg <- tryCatch(
        mgcv::gam(y ~ s(x, k = 5), data = dg, family = binomial,
                   method = "REML"),
        error = function(e) NULL)
      if (is.null(fg)) return(NULL)
      gx <- x_grid[x_grid >= qs[1] & x_grid <= qs[2]]
      if (length(gx) == 0) return(NULL)
      p <- predict(fg, data.frame(x = gx), type = "response")
      data.table(grp = g, x = gx, fit = p * 100)
    })
    per_dt <- rbindlist(grp_list[!vapply(grp_list, is.null, logical(1))])
  }
  list(pooled = pooled_dt, per_group = per_dt)
}

message("[29A] fitting density GAMs...")
mac_gam   <- fit_gam_density(wt, "is_mac",   group_col = "sample_key")
lymph_gam <- fit_gam_density(wt, "is_lymph", group_col = "sample_key")

panel_A_wt <- ggplot() +
  geom_line(data = mac_gam$per_group,
            aes(x = x, y = fit, group = grp),
            colour = COL_MAC, linetype = "dotted", alpha = 0.4,
            linewidth = 0.4) +
  geom_line(data = lymph_gam$per_group,
            aes(x = x, y = fit, group = grp),
            colour = COL_LYMPH, linetype = "dotted", alpha = 0.4,
            linewidth = 0.4) +
  geom_ribbon(data = mac_gam$pooled,
              aes(x = x, ymin = lo, ymax = hi),
              fill = COL_MAC, alpha = 0.15) +
  geom_ribbon(data = lymph_gam$pooled,
              aes(x = x, ymin = lo, ymax = hi),
              fill = COL_LYMPH, alpha = 0.15) +
  geom_line(data = mac_gam$pooled,
            aes(x = x, y = fit), colour = COL_MAC, linewidth = 1.1) +
  geom_line(data = lymph_gam$pooled,
            aes(x = x, y = fit), colour = COL_LYMPH, linewidth = 1.1) +
  labs(x = "Niche metabolic stress (z-score, K=50\u00b5m)",
       y = "Cell-type density (% of TME cells)",
       title = "Macrophage rises, lymphocyte falls across niche hypoxia (WT)",
       subtitle = "Pooled GAM + 95% CI. Per-sample GAMs dotted. n=8 WT samples.") +
  annotate("text", x = 2.0, y = max(mac_gam$pooled$fit) * 1.02,
            label = "Macrophage", colour = COL_MAC,
            hjust = 1, fontface = "bold", size = 3.5) +
  annotate("text", x = 2.0, y = min(lymph_gam$pooled$fit) * 1.02,
            label = "Lymphocyte\n(T+NK+B+plasma)", colour = COL_LYMPH,
            hjust = 1, fontface = "bold", size = 3.2) +
  theme_lab(base_size = 10)
# TMA panel
message("[29A] fitting TMA density GAM...")
tma_ok <- tma[!is.na(patient_id) & patient_id != ""]
mac_gam_tma   <- fit_gam_density(tma_ok, "is_mac",   group_col = "patient_id")
lymph_gam_tma <- fit_gam_density(tma_ok, "is_lymph", group_col = "patient_id")
panel_A_tma <- ggplot() +
  geom_ribbon(data = mac_gam_tma$pooled,
              aes(x = x, ymin = lo, ymax = hi),
              fill = COL_MAC, alpha = 0.2) +
  geom_ribbon(data = lymph_gam_tma$pooled,
              aes(x = x, ymin = lo, ymax = hi),
              fill = COL_LYMPH, alpha = 0.2) +
  geom_line(data = mac_gam_tma$pooled,
            aes(x = x, y = fit), colour = COL_MAC, linewidth = 1.1) +
  geom_line(data = lymph_gam_tma$pooled,
            aes(x = x, y = fit), colour = COL_LYMPH, linewidth = 1.1) +
  labs(x = "Niche metabolic stress (z-score, K=50\u00b5m)",
       y = "Cell-type density (% of TME cells)",
       title = "TMA validation",
       subtitle = "Pooled GAM + 95% CI") +
  theme_lab(base_size = 10)

ggsave(file.path(FIG_DIR, "fig29A_density_vs_stress.png"),
       panel_A_wt | panel_A_tma,
       width = 12, height = 4.5, dpi = 400, bg = "white",
       device = ragg::agg_png)
message("  saved: fig29A_density_vs_stress.png")

# ---------------------------------------------------------------------------
# 29B — paired enrichment forest
# ---------------------------------------------------------------------------
message("[29B] paired enrichment forest...")
pe_wt  <- fread(file.path(OUT_DIR, "paired_enrichment_per_sample_wt.csv"))
pe_tma <- fread(file.path(OUT_DIR, "paired_enrichment_per_patient_tma.csv"))
pe_all <- rbind(pe_wt, pe_tma, fill = TRUE)
# Primary test: stress_decile + T+NK+B+plasma
primary <- pe_all[bin == "stress_decile" & lymph_set == "T+NK+B+plasma" &
                    eligible == TRUE]
# Use group_col-agnostic id
primary[, id := ifelse(!is.na(sample_key) & sample_key != "",
                         sample_key, patient_id)]
if (!"id" %in% colnames(primary)) primary[, id := "x"]
# Handle per_patient columns depending on how fread typed them
# sample_key present in WT rows only, patient_id in TMA rows
primary[, cohort := fifelse(is.na(patient_id) | patient_id == "", "WT", "TMA")]
# sort within cohort by LR
primary[, group_plot := factor(id, levels = primary[order(LR)]$id)]

panel_B <- ggplot(primary,
                  aes(x = LR, y = group_plot,
                      colour = ifelse(LR > 0, "mac-enriched", "lymph-enriched"))) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = LR, yend = group_plot),
               linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("mac-enriched" = COL_MAC,
                                   "lymph-enriched" = COL_LYMPH),
                       name = NULL) +
  facet_wrap(~ cohort, scales = "free_y", ncol = 2) +
  labs(x = "log2( Enrichment_mac / Enrichment_lymph ) across top-vs-bottom niche stress decile",
       y = NULL,
       title = "Paired per-sample/patient log-enrichment: macrophages vs lymphocytes in top hypoxic decile",
       subtitle = sprintf("WT: median LR=%+.2f (n=%d samples, %.0f%% positive, p=%.3g); TMA: median LR=%+.2f (n=%d patients, %.0f%% positive, p=%.3g)",
                           median(primary[cohort == "WT"]$LR),
                           sum(primary$cohort == "WT"),
                           100 * mean(primary[cohort == "WT"]$LR > 0),
                           wilcox.test(primary[cohort == "WT"]$LR, mu = 0)$p.value,
                           median(primary[cohort == "TMA"]$LR),
                           sum(primary$cohort == "TMA"),
                           100 * mean(primary[cohort == "TMA"]$LR > 0),
                           wilcox.test(primary[cohort == "TMA"]$LR, mu = 0)$p.value)) +
  theme_lab(base_size = 10) +
  theme(axis.text.y = element_text(size = rel(0.55)),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig29B_paired_enrichment_forest.png"),
       panel_B, width = 11, height = 9, dpi = 400, bg = "white",
       device = ragg::agg_png)
message("  saved: fig29B_paired_enrichment_forest.png")

# ---------------------------------------------------------------------------
# 29C — functional survival grid
# ---------------------------------------------------------------------------
message("[29C] functional survival grid...")
fs <- fread(file.path(OUT_DIR, "functional_survival_summary.csv"))
fs[, sig := p_wilcox < 0.05]
fs[, score_label := gsub("func_|_score", "", score)]
fs[, score_label := gsub("_", " ", score_label)]
# Put major cell types in a row
fs[, cell_type := factor(cell_type,
                          levels = c("Macrophage", "T cell", "NK cell",
                                      "B cell", "Plasma cell"))]

panel_C <- ggplot(fs,
                   aes(x = score_label, y = cohort, fill = median_delta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f\n(%d/%d%%%s)",
                                  median_delta,
                                  round(pct_positive),
                                  100,
                                  ifelse(sig, "*", ""))),
             size = 2.3, colour = "grey10") +
  scale_fill_gradient2(low = COL_LYMPH, mid = "white", high = COL_MAC,
                        midpoint = 0, name = "\u0394 (top-bot)") +
  facet_wrap(~ cell_type, nrow = 1, scales = "free_x") +
  labs(x = "Functional score", y = NULL,
       title = "Functional survival in top-decile niche hypoxia (Δ top-vs-bottom)",
       subtitle = "Red = higher in top (macrophage-favoring) / Blue = lower in top (lymphocyte-suppressing). * p<0.05 paired Wilcoxon.") +
  theme_lab(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = rel(0.8)),
        legend.position = "right",
        strip.text = element_text(face = "bold", size = rel(0.9)))

ggsave(file.path(FIG_DIR, "fig29C_functional_survival_grid.png"),
       panel_C, width = 14, height = 4, dpi = 400, bg = "white",
       device = ragg::agg_png)
message("  saved: fig29C_functional_survival_grid.png")

# ---------------------------------------------------------------------------
# 29D — spatial distance test
# ---------------------------------------------------------------------------
message("[29D] spatial distance figure...")
sd_wt  <- fread(file.path(OUT_DIR, "spatial_distance_per_sample_wt.csv"))
sd_tma <- fread(file.path(OUT_DIR, "spatial_distance_per_patient_tma.csv"))

sd_wt_long <- melt(sd_wt, id.vars = c("sample_id", "cohort"),
                    measure.vars = c("Macrophage", "Lymphocyte"),
                    variable.name = "cell_type", value.name = "med_dist")
sd_tma_long <- melt(sd_tma[!is.na(delta)], id.vars = c("patient_id", "cohort"),
                     measure.vars = c("Macrophage", "Lymphocyte"),
                     variable.name = "cell_type", value.name = "med_dist")
sd_wt_long[, id := sample_id]
sd_tma_long[, id := patient_id]

panel_D_wt <- ggplot(sd_wt_long,
                      aes(x = cell_type, y = med_dist,
                          group = id, colour = cell_type)) +
  geom_line(aes(group = id), colour = "grey70", alpha = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(Macrophage = COL_MAC,
                                   Lymphocyte = COL_LYMPH),
                       guide = "none") +
  labs(x = NULL, y = "Median distance to top-hypoxia anchor (\u00b5m)",
       title = sprintf("WT spatial distance (n=%d samples)", nrow(sd_wt)),
       subtitle = sprintf("\u0394 = %.1f \u00b5m; %.0f%% mac-closer; p=%.3g (paired, 1-sided)",
                           median(sd_wt$delta), 100 * mean(sd_wt$delta < 0),
                           wilcox.test(sd_wt$Macrophage, sd_wt$Lymphocyte,
                                       paired = TRUE, alternative = "less")$p.value)) +
  theme_lab(base_size = 10)

panel_D_tma <- ggplot(sd_tma_long,
                       aes(x = cell_type, y = med_dist,
                           group = id, colour = cell_type)) +
  geom_line(aes(group = id), colour = "grey80", alpha = 0.3, linewidth = 0.3) +
  geom_point(size = 1, alpha = 0.8) +
  scale_colour_manual(values = c(Macrophage = COL_MAC,
                                   Lymphocyte = COL_LYMPH),
                       guide = "none") +
  labs(x = NULL, y = "Median distance to top-hypoxia anchor (\u00b5m)",
       title = sprintf("TMA spatial distance (n=%d patients)", nrow(sd_tma[!is.na(delta)])),
       subtitle = sprintf("\u0394 = %.1f \u00b5m; %.0f%% mac-closer; p=%.3g (paired, 1-sided)",
                           median(sd_tma$delta, na.rm = TRUE),
                           100 * mean(sd_tma$delta < 0, na.rm = TRUE),
                           wilcox.test(sd_tma[!is.na(delta)]$Macrophage,
                                       sd_tma[!is.na(delta)]$Lymphocyte,
                                       paired = TRUE, alternative = "less")$p.value)) +
  theme_lab(base_size = 10)

ggsave(file.path(FIG_DIR, "fig29D_spatial_distance.png"),
       panel_D_wt | panel_D_tma,
       width = 10, height = 5, dpi = 400, bg = "white",
       device = ragg::agg_png)
message("  saved: fig29D_spatial_distance.png")

# ---------------------------------------------------------------------------
# 29E — summary manuscript figure (stack A top-wide, then B and C)
# ---------------------------------------------------------------------------
combo <- (panel_A_wt | panel_A_tma) /
          (panel_D_wt | panel_D_tma) /
          panel_C /
          panel_B +
          plot_layout(heights = c(1, 1, 1, 2.8))
ggsave(file.path(FIG_DIR, "fig29E_manuscript_summary.png"),
       combo, width = 14, height = 18, dpi = 400, bg = "white",
       device = ragg::agg_png)
message("  saved: fig29E_manuscript_summary.png")

message("\nAll figures saved to: ", FIG_DIR)
