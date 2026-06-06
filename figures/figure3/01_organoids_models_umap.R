#!/usr/bin/env Rscript
# ============================================================================
# Figure 3A — Baseline organoid (PDO) UMAP coloured by model
# ----------------------------------------------------------------------------
# PURPOSE: Baseline PDO UMAP showing where each of the 8 untreated baseline
#   models sits in the integrated embedding, coloured by biological origin
#   (ascites vs primary). Manuscript-ready (no in-figure title / panel letter).
#
# INPUTS:
#   - <organoids_root>/output/01_Data_Processing_and_QC/
#       seurat_untreated_baseline.rds   (Seurat v5; pdo_model, umap reduction)
#     EXTERNAL DEPENDENCY: PDO data are not part of the deposited monorepo
#     objects. Override ORGANOIDS_ROOT to point at the 2026_organoids tree.
#
# OUTPUTS:
#   - figures_dir/organoids_models_umap.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 3A.
#
# RUNTIME TIER: moderate (loads a Seurat object; renders rasterized points).
# ============================================================================

# --- Central config (this file is 2 levels under the repo root) -------------
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure3"
.config_path <- normalizePath(file.path(.fig_dir, "..", "..", "config", "config.R"),
                              mustWork = FALSE)
if (!file.exists(.config_path)) .config_path <- "config/config.R"
source(.config_path)

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table)
})

set.seed(CFG$seed)

# --- Paths -------------------------------------------------------------------
INPUT   <- file.path(cfg_obj("organoids_root"),
                     "output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds")
FIGDIR  <- file.path(path.expand(CFG$paths$figures_dir))
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)
OUT_PNG <- file.path(FIGDIR, "organoids_models_umap.png")
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)

FA <- 6;  FK <- 5.5;  FN <- 5

MODEL_LEVELS   <- c("OCAD106", "OCAD93", "OCAD96", "OCAD97",
                    "OPTO112", "OPTO129", "OPTO98", "PDO66")
ASCITES_MODELS <- c("OCAD106", "OCAD93", "OCAD96", "OCAD97")
ASCITES_FILL   <- "#5665B6"   # B cell
PRIMARY_FILL   <- "#8FBC8F"   # Macrophage

# --- Load and extract --------------------------------------------------------
cat("Loading: ", INPUT, "\n", sep = "")
stopifnot(file.exists(INPUT))
s <- readRDS(INPUT)
cat("  cells: ", ncol(s), "  genes: ", nrow(s), "\n", sep = "")

stopifnot("umap"      %in% names(s@reductions))
stopifnot("pdo_model" %in% colnames(s@meta.data))

emb <- as.data.table(Embeddings(s, "umap"), keep.rownames = "cell")
setnames(emb, c("cell", "UMAP1", "UMAP2"))
emb[, pdo_model := factor(s@meta.data$pdo_model[match(cell, colnames(s))],
                          levels = MODEL_LEVELS)]
stopifnot(!any(is.na(emb$pdo_model)))

emb[, origin := factor(
  ifelse(as.character(pdo_model) %in% ASCITES_MODELS, "ascites", "primary"),
  levels = c("ascites", "primary"))]

emb <- emb[sample(.N)]

labels_dt <- emb[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = pdo_model]

# --- Plot --------------------------------------------------------------------
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = origin)) +
  geom_point(size = 0.18, alpha = 0.85, shape = 16, stroke = 0) +
  scale_colour_manual(values = c(ascites = ASCITES_FILL, primary = PRIMARY_FILL),
                      labels = c(ascites = "PDO — ascites", primary = "PDO — primary"),
                      name = NULL) +
  geom_text(data = labels_dt, aes(x = UMAP1, y = UMAP2, label = pdo_model),
            inherit.aes = FALSE, size = FN / .pt, colour = "black",
            family = "Helvetica", fontface = "plain") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1),
                               keyheight = unit(7, "pt"), keywidth = unit(7, "pt"))) +
  coord_fixed() +
  theme_classic(base_size = FA) +
  theme(
    text             = element_text(family = "Helvetica", colour = "black"),
    axis.line        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    axis.title       = element_text(size = FA, colour = "black"),
    axis.title.x     = element_text(margin = margin(t = 2)),
    axis.title.y     = element_text(margin = margin(r = 2)),
    axis.text        = element_text(size = FK, colour = "black"),
    legend.text      = element_text(size = FN, colour = "black"),
    legend.title     = element_blank(),
    legend.key       = element_blank(),
    legend.position  = "right",
    legend.margin    = margin(0, 0, 0, 2),
    legend.box.spacing = unit(2, "pt"),
    plot.margin      = margin(2, 2, 2, 2),
    plot.title       = element_blank())

ggsave(OUT_PNG, p, width = 88, height = 62, units = "mm", dpi = 600)
ggsave(OUT_SVG, p, width = 88, height = 62, units = "mm")
cat("\nSaved:\n  ", basename(OUT_PNG), "\n  ", basename(OUT_SVG), "\n", sep = "")
