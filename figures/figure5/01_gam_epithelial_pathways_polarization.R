#!/usr/bin/env Rscript
# ============================================================================
# Figure 5A — Pooled GAM: epithelial UCell pathway scores vs polarization
# ----------------------------------------------------------------------------
# PURPOSE
#   Whole-tissue, epithelial-only pooled GAM (mgcv REML) of z-scored UCell
#   pathway scores against the focal-cell polarization score (SecA -> SecB).
#   Features: Proliferation, Oncogenic growth (mean of z-scored Myc/PI3K/Wnt/
#   Notch activating), NF-kB, Glycolysis, DNA damage repair, Epigenetic
#   silencing, and LDHA (z-scored expression).
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/19d_gene_polarization_gams/
#     epithelial_expression_polarization.rds   (19d cache; per-cell expr + polarization)
#   data_root/2026_final_xenium_analysis/output/9b_scoring/pathway_gene_sets_v2.csv
#   Custom sigs (inline): ddr (BRCA1,BRCA2,ATM,CHEK2),
#     epi_silencing (DNMT3A,HDAC1,KDM5A,NCOR1,NCOR2).
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (theme_lab).
#
# OUTPUTS
#   figures_dir/figure5/gam_epithelial_pathways_polarization.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 5A
# RUNTIME TIER: moderate (UCell scoring on subsampled epithelial cells)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# --- Config + shared setup (script is 2 levels deep: figures/figure5/) -------
.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2); library(UCell)
})

set.seed(CFG$seed)

OUT_STEM <- cfg_path("figures_dir", "figure5", "gam_epithelial_pathways_polarization")
XEN_OUT  <- function(...) cfg_path("data_root", "2026_final_xenium_analysis", "output", ...)

N_GRID  <- 200
K_POOL  <- 10

# ── Load epithelial expression + polarization ────────────────────────────
message("[load] epithelial expression + polarization (cached, 19d)")
expr <- readRDS(XEN_OUT("19d_gene_polarization_gams",
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

# ── Pathway gene sets ────────────────────────────────────────────────────
pw_dt <- fread(XEN_OUT("9b_scoring", "pathway_gene_sets_v2.csv"))

PW_SCORE <- c("glycolysis", "proliferation", "myc_activating",
              "pi3k_activating", "wnt_activating", "notch_activating", "nfkb")

pw_sets <- split(pw_dt[pathway %in% PW_SCORE]$gene,
                  pw_dt[pathway %in% PW_SCORE]$pathway)

CUSTOM_SETS <- list(
  ddr           = c("BRCA1", "BRCA2", "ATM", "CHEK2"),
  epi_silencing = c("DNMT3A", "HDAC1", "KDM5A", "NCOR1", "NCOR2")
)
pw_sets <- c(pw_sets, CUSTOM_SETS)

ALL_PW <- c(PW_SCORE, "ddr", "epi_silencing")

META_COLS <- c("cell_id", "polarization_UCell", "sample_id",
                "cell_label", "pol_bin")
gene_cols <- setdiff(colnames(expr), META_COLS)
pw_sets <- lapply(pw_sets, function(g) intersect(g, gene_cols))

message("[pathway gene coverage]")
for (nm in names(pw_sets)) {
  message(sprintf("  %s: %d genes", nm, length(pw_sets[[nm]])))
}

# ── UCell scoring ────────────────────────────────────────────────────────
mat <- Matrix::Matrix(t(as.matrix(expr[, ..gene_cols])), sparse = TRUE)
colnames(mat) <- expr$cell_id
message("[ucell] scoring ", ncol(mat), " cells across ",
        length(pw_sets), " pathways")

ucell_scores <- ScoreSignatures_UCell(matrix = mat,
                                       features = pw_sets,
                                       maxRank  = nrow(mat),
                                       name     = "",
                                       chunk.size = 500)
ucell_scores <- as.data.table(ucell_scores, keep.rownames = "cell_id")

pw_long <- melt(ucell_scores, id.vars = "cell_id",
                variable.name = "pathway", value.name = "ucell")
pw_long[, pathway := as.character(pathway)]
pw_long <- merge(pw_long,
                  expr[, .(cell_id, polarization_UCell, sample_id)],
                  by = "cell_id")

# ── Define final features for the plot ───────────────────────────────────
PLOT_PW_DIRECT <- c("proliferation", "nfkb", "glycolysis", "ddr", "epi_silencing")
PLOT_PW <- c("proliferation", "oncogenic_growth", "nfkb",
             "glycolysis", "ddr", "epi_silencing",
             "LDHA")

nfkb_label <- "NF-κB"
PLOT_LABELS <- c("proliferation"    = "Proliferation",
                  "oncogenic_growth" = "Oncogenic growth",
                  "nfkb"             = nfkb_label,
                  "glycolysis"       = "Glycolysis",
                  "ddr"              = "DNA damage repair",
                  "epi_silencing"    = "Epigenetic silencing",
                  "LDHA"             = "LDHA")
PLOT_COLORS <- setNames(
  c("#E6A141", "#5665B6", "#8A5DAF", "#D14E6C", "#8FBC8F", "#56AFC4", "#B87A7A"),
  c("Proliferation", "Oncogenic growth", nfkb_label,
    "Glycolysis", "DNA damage repair", "Epigenetic silencing", "LDHA")
)

# ── Z-score per pathway (including individual activating pathways) ───────
message("\n[zscore] z-scoring per pathway ...")
onco_pw <- c("myc_activating", "pi3k_activating",
             "wnt_activating", "notch_activating")
pw_plot <- pw_long[pathway %in% c(PLOT_PW_DIRECT, onco_pw)]
pw_plot[, ucell_z := (ucell - mean(ucell, na.rm = TRUE)) / sd(ucell, na.rm = TRUE),
         by = pathway]

# ── Oncogenic growth (mean of z-scored activating pathways) ──────────────
message("\n[oncogenic] averaging z-scored Myc + PI3K + Wnt + Notch ...")
onco_wide <- dcast(pw_plot[pathway %in% onco_pw],
                    cell_id + polarization_UCell + sample_id ~ pathway,
                    value.var = "ucell_z", fun.aggregate = mean)
onco_wide[, oncogenic_growth := rowMeans(.SD, na.rm = TRUE),
           .SDcols = onco_pw]
onco_long <- onco_wide[, .(cell_id, polarization_UCell, sample_id,
                            pathway = "oncogenic_growth",
                            ucell_z = oncogenic_growth)]

# ── LDHA z-scored expression ─────────────────────────────────────────────
message("\n[LDHA] adding z-scored gene expression ...")
ldha_vals <- expr[["LDHA"]]
ldha_z <- (ldha_vals - mean(ldha_vals, na.rm = TRUE)) / sd(ldha_vals, na.rm = TRUE)
ldha_long <- expr[, .(cell_id, polarization_UCell, sample_id)]
ldha_long[, pathway := "LDHA"]
ldha_long[, ucell_z := ldha_z]

pw_plot <- rbind(pw_plot[pathway %in% PLOT_PW_DIRECT,
                          .(cell_id, polarization_UCell, sample_id,
                            pathway, ucell_z)],
                  onco_long,
                  ldha_long)

# ── GAM fitting ──────────────────────────────────────────────────────────
pol_range <- range(pw_plot$polarization_UCell, na.rm = TRUE)
grid_pol  <- seq(pol_range[1], pol_range[2], length.out = N_GRID)

message("\n[gam] fitting pooled GAMs ...")
pooled_all <- list()
for (pw in PLOT_PW) {
  d <- pw_plot[pathway == pw & !is.na(ucell_z)]
  label <- PLOT_LABELS[pw]
  fit <- tryCatch(
    gam(ucell_z ~ s(polarization_UCell, k = K_POOL),
        data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) { message("  FAILED: ", label); next }
  pp <- predict(fit, newdata = data.frame(polarization_UCell = grid_pol),
                se.fit = TRUE)
  pooled_all[[pw]] <- data.table(
    polarization = grid_pol,
    fitted = as.numeric(pp$fit),
    lower  = as.numeric(pp$fit - 1.96 * pp$se.fit),
    upper  = as.numeric(pp$fit + 1.96 * pp$se.fit),
    feature = label
  )
  message(sprintf("  %s: OK (n=%d)", label, nrow(d)))
}
pooled <- rbindlist(pooled_all)
pooled[, feature := factor(feature, levels = unname(PLOT_LABELS))]

# ── End-of-line label positions ──────────────────────────────────────────
label_df <- pooled[, .SD[which.max(polarization)], by = feature]
label_df[, color := PLOT_COLORS[as.character(feature)]]
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
  geom_text(data = label_df[feature != "LDHA"],
            aes(x = polarization + 0.006, y = fitted,
                label = feature, color = feature),
            hjust = 0, size = 2.0, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = label_df[feature == "LDHA"],
            aes(x = polarization + 0.006, y = fitted,
                label = feature, color = feature),
            hjust = 0, size = 2.0, fontface = "bold.italic",
            show.legend = FALSE) +
  scale_color_manual(values = PLOT_COLORS) +
  scale_fill_manual(values = PLOT_COLORS) +
  labs(x = "Focal cell polarization score (SecA to SecB)",
       y = "Epithelial UCell\nz-score") +
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
w <- 4.25
h <- 2.5

ggsave(paste0(OUT_STEM, ".pdf"), p, width = w, height = h, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), p, width = w, height = h, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(OUT_STEM, ".svg"), p, width = w, height = h, bg = "white")

message("\nDONE")
