#!/usr/bin/env Rscript
# ============================================================================
# Figure 6E,6F — T cell / NK cell exhaustion (UCell) low vs high stress
# ----------------------------------------------------------------------------
# PURPOSE
#   Re-scored T cell and NK cell exhaustion with UCell (consistent with all
#   other signature scoring) across stress deciles. Per-sample paired medians
#   of bottom (decile 1) vs top (decile 10) stress, with cell-level violins.
#   Wilcoxon signed-rank on per-sample deltas.
#   Exhaustion gene sets (trimmed to panel):
#     T cell: PDCD1, HAVCR2, LAG3, TIGIT, CTLA4
#     NK cell: TIGIT, HAVCR2, LAG3
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/29_macrophage_niche_survival/
#     per_cell_niche_scores.rds  (stress_decile)
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_<sample> (load_sfe)
#   Cohort: CFG$cohort$whole_tissue (published 8 WT samples) + sfe_tma_filtered.
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe, ref_palette, theme_lab).
#
# OUTPUTS
#   figures_dir/figure6/fig_tcell_exclusion_tma.{svg,png,pdf}
#   figures_dir/figure6/fig_tcell_exhaustion_wt.{svg,png,pdf}    (Fig 6E)
#   figures_dir/figure6/fig_nkcell_exhaustion_wt.{svg,png,pdf}   (Fig 6F)
#
# MANUSCRIPT PANEL(S): Fig 6E (T cell exhaustion), Fig 6F (NK cell exhaustion)
# RUNTIME TIER: heavy (loads 8 WT SFEs + TMA; per-cell UCell scoring)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(UCell)
  library(Matrix); library(ragg); library(svglite)
})

set.seed(CFG$seed)

FIG_DIR <- cfg_path("figures_dir", "figure6")
OUT_29  <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                    "29_macrophage_niche_survival")

MIN_CELLS <- 15

# Published whole-tissue cohort (config single source of truth), sfe_ prefixed.
WT_SAMPLES <- paste0("sfe_", CFG$cohort$whole_tissue)

# Trimmed exhaustion sets (panel-restricted; receptors only).
TCELL_EXHAUST <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4")
NK_EXHAUST    <- c("TIGIT", "HAVCR2", "LAG3")

# ── Colours ───────────────────────────────────────────────────────────────
COL_T     <- unname(ref_palette["T cell"])
COL_NK    <- unname(ref_palette["NK cell"])
COL_T_LO  <- "#BDE0FD"
COL_T_HI  <- "#4A9BD9"
COL_NK_LO <- "#9DD4DE"
COL_NK_HI <- "#2E8A9E"

# ── Load niche cache (for stress_decile) ──────────────────────────────────
message("[1/4] Loading niche cache ...")
niche <- readRDS(file.path(OUT_29, "per_cell_niche_scores.rds"))

# ── UCell scoring per SFE ─────────────────────────────────────────────────
message("[2/4] UCell exhaustion scoring ...")

pw_sets <- list(tcell_exhaustion = TCELL_EXHAUST,
                nk_exhaustion    = NK_EXHAUST)

score_sfe <- function(sfe, cell_ids, label) {
  keep <- intersect(cell_ids, colnames(sfe))
  if (length(keep) < 50) return(NULL)
  sfe_sub <- sfe[, keep]

  mat <- as(assay(sfe_sub, "counts"), "dgCMatrix")

  for (nm in names(pw_sets)) {
    avail <- intersect(pw_sets[[nm]], rownames(mat))
    message(sprintf("    %s: %d/%d genes available", nm,
                    length(avail), length(pw_sets[[nm]])))
  }

  uc <- ScoreSignatures_UCell(matrix = mat, features = pw_sets,
                               maxRank = nrow(mat), name = "",
                               chunk.size = 5000)
  dt <- as.data.table(uc, keep.rownames = "cell_id")
  dt
}

# ── WT samples ────────────────────────────────────────────────────────────
wt_list <- list()
for (s in WT_SAMPLES) {
  samp <- sub("^sfe_", "", s)
  message(sprintf("\n  === %s ===", samp))
  sfe <- load_sfe(s)

  sample_niche <- niche$wt[sample_key == samp]
  immune_ids <- sample_niche[cell_label %in% c("T cell", "NK cell")]$cell_id

  scores <- score_sfe(sfe, immune_ids, samp)
  if (!is.null(scores)) {
    merged <- merge(sample_niche[cell_id %in% scores$cell_id],
                    scores, by = "cell_id")
    wt_list[[samp]] <- merged
    message(sprintf("    merged: %d cells", nrow(merged)))
  }
  rm(sfe); gc(verbose = FALSE)
}
wt <- rbindlist(wt_list, fill = TRUE)
message(sprintf("\n  WT total: %s immune cells scored",
                format(nrow(wt), big.mark = ",")))

# ── TMA ───────────────────────────────────────────────────────────────────
message("\n  === TMA ===")
sfe_t <- load_sfe("sfe_tma_filtered")
tma_niche <- niche$tma[!is.na(patient_id) & patient_id != ""]
tma_immune_ids <- tma_niche[cell_label %in% c("T cell", "NK cell")]$cell_id
tma_scores <- score_sfe(sfe_t, tma_immune_ids, "TMA")
tma <- merge(tma_niche[cell_id %in% tma_scores$cell_id],
             tma_scores, by = "cell_id")
message(sprintf("    TMA total: %s immune cells scored",
                format(nrow(tma), big.mark = ",")))
rm(sfe_t); gc(verbose = FALSE)

# ── Compute per-group deltas ──────────────────────────────────────────────
message("\n[3/4] Computing deltas ...")

compute_deltas_ct <- function(dt, group_col, ct, score_col) {
  ct_dt <- dt[cell_label == ct]
  d <- ct_dt[, .(
    n_top   = sum(stress_decile == 10, na.rm = TRUE),
    n_bot   = sum(stress_decile == 1,  na.rm = TRUE),
    med_top = median(get(score_col)[stress_decile == 10], na.rm = TRUE),
    med_bot = median(get(score_col)[stress_decile == 1],  na.rm = TRUE)
  ), by = group_col]
  d[, delta := med_top - med_bot]
  d
}

# w_t / w_nk are only well-defined when >= 2 samples are eligible (Wilcoxon
# signed-rank needs a paired sample). Initialize NA-valued sentinels so the
# downstream figure annotations (which reference w_t$p.value / w_nk$p.value)
# degrade gracefully instead of erroring on an undefined symbol when the
# eligible-sample count drops below 2.
w_t  <- list(p.value = NA_real_)
w_nk <- list(p.value = NA_real_)

wt_t_per <- compute_deltas_ct(wt, "sample_key", "T cell", "tcell_exhaustion")
wt_t_elig <- wt_t_per[n_top >= MIN_CELLS & n_bot >= MIN_CELLS & is.finite(delta)]
if (nrow(wt_t_elig) >= 2) {
  w_t <- wilcox.test(wt_t_elig$delta, mu = 0)
  message(sprintf("  T cell: %d/%d samples eligible, median delta=%.4f, p=%.4f",
                  nrow(wt_t_elig), nrow(wt_t_per), median(wt_t_elig$delta), w_t$p.value))
} else {
  message(sprintf("  T cell: %d/%d samples eligible (<2) — skipping Wilcoxon test",
                  nrow(wt_t_elig), nrow(wt_t_per)))
}

wt_nk_per <- compute_deltas_ct(wt, "sample_key", "NK cell", "nk_exhaustion")
wt_nk_elig <- wt_nk_per[n_top >= MIN_CELLS & n_bot >= MIN_CELLS & is.finite(delta)]
if (nrow(wt_nk_elig) >= 2) {
  w_nk <- wilcox.test(wt_nk_elig$delta, mu = 0)
  message(sprintf("  NK cell: %d/%d samples eligible, median delta=%.4f, p=%.4f",
                  nrow(wt_nk_elig), nrow(wt_nk_per), median(wt_nk_elig$delta), w_nk$p.value))
} else {
  message(sprintf("  NK cell: %d/%d samples eligible (<2) — skipping Wilcoxon test",
                  nrow(wt_nk_elig), nrow(wt_nk_per)))
}

tma_t_per <- compute_deltas_ct(tma, "patient_id", "T cell", "tcell_exhaustion")
n_tma_elig <- sum(tma_t_per$n_top >= MIN_CELLS & tma_t_per$n_bot >= MIN_CELLS,
                   na.rm = TRUE)
message(sprintf("  TMA T cell: %d/%d patients eligible",
                n_tma_elig, nrow(tma_t_per)))

# =========================================================================
# [4/4] FIGURES — each as its own SVG
# =========================================================================
message("\n[4/4] Generating figures ...")

# ── Figure A: TMA T cell exclusion histogram ─────────────────────────────
message("  TMA T cell exclusion ...")

tma_long <- melt(tma_t_per[, .(patient_id, n_top, n_bot)],
                  id.vars = "patient_id",
                  variable.name = "bin", value.name = "n_cells")
tma_long[, bin := factor(fifelse(bin == "n_bot",
                                  "Bottom decile\n(low stress)",
                                  "Top decile\n(high stress)"),
                          levels = c("Bottom decile\n(low stress)",
                                     "Top decile\n(high stress)"))]

n_total <- nrow(tma_t_per)

p_excl <- ggplot(tma_long, aes(x = n_cells, fill = bin)) +
  geom_histogram(binwidth = 2, colour = "white", linewidth = 0.15,
                 position = "identity", alpha = 0.7) +
  geom_vline(xintercept = MIN_CELLS, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  annotate("text", x = MIN_CELLS + 1, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("min. threshold\n(n=%d)", MIN_CELLS),
           size = 1.8, colour = "grey30") +
  scale_fill_manual(values = c("Bottom decile\n(low stress)" = COL_T_LO,
                                "Top decile\n(high stress)" = COL_T_HI),
                     name = NULL) +
  labs(x = "T cells per patient per decile bin",
       y = "Number of TMA patients") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           label = sprintf("TMA (n=%d patients)\n%d/%d eligible for\npaired comparison",
                           n_total, n_tma_elig, n_total),
           size = 1.8, colour = "grey30") +
  theme_lab() +
  theme(
    legend.position   = c(0.72, 0.65),
    legend.text       = element_text(size = 5),
    legend.key.size   = unit(0.3, "cm"),
    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
    axis.title        = element_text(size = 6),
    axis.text         = element_text(size = 5.5),
    plot.margin       = margin(4, 8, 4, 4)
  )

stem_excl <- file.path(FIG_DIR, "fig_tcell_exclusion_tma")
ggsave(paste0(stem_excl, ".svg"), p_excl, width = 3.2, height = 2.5, bg = "white")
ggsave(paste0(stem_excl, ".png"), p_excl, width = 3.2, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_excl, ".pdf"), p_excl, width = 3.2, height = 2.5, bg = "white")
message("    saved: ", basename(stem_excl))

# ── Figure B: WT T cell exhaustion (UCell) — Fig 6E ───────────────────────
message("  WT T cell exhaustion ...")

tc_cells <- wt[cell_label == "T cell" & stress_decile %in% c(1, 10)]
tc_cells[, stress_bin := factor(
  fifelse(stress_decile == 1,
          "Low stress\n(bottom decile)",
          "High stress\n(top decile)"),
  levels = c("Low stress\n(bottom decile)",
             "High stress\n(top decile)")
)]

pct_pos_t <- round(100 * mean(wt_t_elig$delta > 0), 1)

tc_paired <- rbind(
  wt_t_elig[, .(sample = sample_key, x = 1, y = med_bot)],
  wt_t_elig[, .(sample = sample_key, x = 2, y = med_top)]
)

tc_cells[, x_num := fifelse(stress_decile == 1, 1, 2)]
tc_paired[, direction := fifelse(
  y[x == 2] > y[x == 1], "up", "down"), by = sample]

n_up_t <- sum(wt_t_elig$delta > 0)

# NA-safe annotation: include the p-value line only when the Wilcoxon test ran
# (>= 2 eligible samples); otherwise show counts without a "p = NA" line.
lab_tc <- if (is.finite(w_t$p.value)) {
  sprintf("n=%d samples\n%d/%d increase\np = %.3f",
          nrow(wt_t_elig), n_up_t, nrow(wt_t_elig), w_t$p.value)
} else {
  sprintf("n=%d samples\n%d/%d increase",
          nrow(wt_t_elig), n_up_t, nrow(wt_t_elig))
}

p_tc <- ggplot() +
  geom_violin(data = tc_cells,
              aes(x = x_num, y = tcell_exhaustion,
                  fill = factor(x_num), group = x_num),
              colour = NA, alpha = 0.15,
              scale = "width", width = 0.6) +
  geom_line(data = tc_paired,
            aes(x = x, y = y, group = sample, colour = direction),
            linewidth = 0.6, alpha = 0.8) +
  geom_point(data = tc_paired,
             aes(x = x, y = y, fill = factor(x)),
             shape = 21, size = 2.5, colour = "grey30", stroke = 0.3) +
  scale_fill_manual(values = c("1" = COL_T_LO, "2" = COL_T_HI),
                    guide = "none") +
  scale_colour_manual(values = c("up" = COL_T_HI, "down" = "grey70"),
                      guide = "none") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Low stress\n(bottom decile)",
                                "High stress\n(top decile)"),
                     expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL,
       y = "T cell exhaustion\n(UCell score)") +
  annotate("text", x = 1.5, y = Inf, vjust = 1.3,
           label = lab_tc,
           size = 1.8, colour = "grey30") +
  theme_lab() +
  theme(
    axis.title      = element_text(size = 6),
    axis.text.x     = element_text(size = 5.5),
    axis.text.y     = element_text(size = 5.5),
    plot.margin     = margin(4, 4, 4, 4)
  )

stem_tc <- file.path(FIG_DIR, "fig_tcell_exhaustion_wt")
ggsave(paste0(stem_tc, ".svg"), p_tc, width = 2.5, height = 2.5, bg = "white")
ggsave(paste0(stem_tc, ".png"), p_tc, width = 2.5, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_tc, ".pdf"), p_tc, width = 2.5, height = 2.5, bg = "white")
message("    saved: ", basename(stem_tc))

# ── Figure C: WT NK cell exhaustion (UCell) — Fig 6F ──────────────────────
message("  WT NK cell exhaustion ...")

nk_cells <- wt[cell_label == "NK cell" & stress_decile %in% c(1, 10)]
nk_cells[, stress_bin := factor(
  fifelse(stress_decile == 1,
          "Low stress\n(bottom decile)",
          "High stress\n(top decile)"),
  levels = c("Low stress\n(bottom decile)",
             "High stress\n(top decile)")
)]

pct_pos_nk <- round(100 * mean(wt_nk_elig$delta > 0), 1)

nk_paired <- rbind(
  wt_nk_elig[, .(sample = sample_key, x = 1, y = med_bot)],
  wt_nk_elig[, .(sample = sample_key, x = 2, y = med_top)]
)

nk_cells[, x_num := fifelse(stress_decile == 1, 1, 2)]
nk_paired[, direction := fifelse(
  y[x == 2] > y[x == 1], "up", "down"), by = sample]

n_up_nk <- sum(wt_nk_elig$delta > 0)

# NA-safe annotation: include the p-value line only when the Wilcoxon test ran
# (>= 2 eligible samples); otherwise show counts without a "p = NA" line.
lab_nk <- if (is.finite(w_nk$p.value)) {
  sprintf("n=%d samples\n%d/%d increase\np = %.3f",
          nrow(wt_nk_elig), n_up_nk, nrow(wt_nk_elig), w_nk$p.value)
} else {
  sprintf("n=%d samples\n%d/%d increase",
          nrow(wt_nk_elig), n_up_nk, nrow(wt_nk_elig))
}

p_nk <- ggplot() +
  geom_violin(data = nk_cells,
              aes(x = x_num, y = nk_exhaustion,
                  fill = factor(x_num), group = x_num),
              colour = NA, alpha = 0.15,
              scale = "width", width = 0.6) +
  geom_line(data = nk_paired,
            aes(x = x, y = y, group = sample, colour = direction),
            linewidth = 0.6, alpha = 0.8) +
  geom_point(data = nk_paired,
             aes(x = x, y = y, fill = factor(x)),
             shape = 21, size = 2.5, colour = "grey30", stroke = 0.3) +
  scale_fill_manual(values = c("1" = COL_NK_LO, "2" = COL_NK_HI),
                    guide = "none") +
  scale_colour_manual(values = c("up" = COL_NK_HI, "down" = "grey70"),
                      guide = "none") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Low stress\n(bottom decile)",
                                "High stress\n(top decile)"),
                     expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL,
       y = "NK cell exhaustion\n(UCell score)") +
  annotate("text", x = 1.5, y = Inf, vjust = 1.3,
           label = lab_nk,
           size = 1.8, colour = "grey30") +
  theme_lab() +
  theme(
    axis.title      = element_text(size = 6),
    axis.text.x     = element_text(size = 5.5),
    axis.text.y     = element_text(size = 5.5),
    plot.margin     = margin(4, 4, 4, 4)
  )

stem_nk <- file.path(FIG_DIR, "fig_nkcell_exhaustion_wt")
ggsave(paste0(stem_nk, ".svg"), p_nk, width = 2.5, height = 2.5, bg = "white")
ggsave(paste0(stem_nk, ".png"), p_nk, width = 2.5, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_nk, ".pdf"), p_nk, width = 2.5, height = 2.5, bg = "white")
message("    saved: ", basename(stem_nk))

message("\nAll figures saved.")
message("DONE")
