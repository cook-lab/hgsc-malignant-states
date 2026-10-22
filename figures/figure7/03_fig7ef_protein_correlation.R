#!/usr/bin/env Rscript
# ============================================================================
# Figure 7e,f input -- per-core Xenium transcripts and polarization vs IF protein, final TMA A cores
# ----------------------------------------------------------------------------
# PURPOSE
#   Builds the per-core table behind Fig 7e,f on the SingleR annotation and the final TMA core set:
#     - KRT7_xenium, KRT19_xenium: mean logcounts over all cells of the core;
#     - polarization_mean: mean polarization_UCell over the core's SingleR secretory cells (SecA, Intermediate, SecB),
#       as polar_mean in the Fig 7a-d survival models;
#     - joined by core id to the published immunofluorescence MFI (workbook below, sheet "Pour SPSS").
#   Cores: TMA A cores of the final set (spatial/00_setup/final_tma_cores.R) that have IF data.
#
# INPUTS
#   sfe_tma (load_sfe; logcounts), output_root/singler_recompute/tma_singler_cells.csv.gz,
#   data_dir/Copy of MFI all tissues TMA valid all markers_.xlsx (sheet "Pour SPSS"), the final TMA core set
#
# OUTPUTS (output_root/41_xenium_protein_correlation/)
#   per_core_xenium_protein_final.csv   core_id, IF MFI (VIM, ECAD, KRT7, KRT18, KRT19), KRT7_xenium, KRT19_xenium,
#                                       polarization_mean, n_cells_core, n_secretory_core (raw ids; lab-internal)
#   correlations_final.csv               n, Pearson and Spearman for every pair drawn in Fig 7e,f
#
# MANUSCRIPT PANEL(S): Fig 7e, 7f (plot: 04_fig7ef_protein_correlation_plot.py with PROTEIN_CSV set to the table above)
# RUNTIME TIER: moderate (loads sfe_tma)
# USAGE (repo root): DATA_ROOT=<deposit> OUTPUT_ROOT=<outputs> Rscript figures/figure7/03_fig7ef_protein_correlation.R
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/final_tma_cores.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment); library(readxl); library(data.table)
})

out_p <- file.path(out_dir, "41_xenium_protein_correlation")
dir.create(out_p, showWarnings = FALSE, recursive = TRUE)
SEC <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
MARKERS <- c("VIM", "ECAD", "KRT7", "KRT18", "KRT19")

# ---- 1. IF protein (published TMA A immunofluorescence) ----
mfi <- as.data.table(suppressWarnings(read_excel(file.path(data_dir, "Copy of MFI all tissues TMA valid all markers_.xlsx"),
                                                 sheet = "Pour SPSS")))
setnames(mfi, c("ID", "E-CADH"), c("core_id", "ECAD"))
mfi[, core_id := as.character(core_id)]
mfi <- mfi[!is.na(core_id), c("core_id", MARKERS), with = FALSE]

# ---- 2. Final TMA A cores: cells, labels, polarization ----
keep <- final_core_ids("sfe_tma")
cells <- fread(file.path(out_dir, "singler_recompute", "tma_singler_cells.csv.gz"),
               select = c("barcode", "core_id", "cell_label_singler", "polarization_UCell"),
               colClasses = list(character = "core_id"))
cells <- cells[core_id %in% keep]

sfe <- suppressWarnings(load_sfe("sfe_tma"))
lc <- assay(sfe, "logcounts")[c("KRT7", "KRT19"), , drop = FALSE]
idx <- match(cells$barcode, colnames(sfe))
stopifnot(!anyNA(idx))
cells[, KRT7_lc := as.numeric(lc["KRT7", idx])]
cells[, KRT19_lc := as.numeric(lc["KRT19", idx])]
rm(sfe, lc); gc(verbose = FALSE)

core <- cells[, .(n_cells_core = .N, n_secretory_core = sum(cell_label_singler %in% SEC),
                  KRT7_xenium = mean(KRT7_lc), KRT19_xenium = mean(KRT19_lc),
                  polarization_mean = mean(polarization_UCell[cell_label_singler %in% SEC], na.rm = TRUE)),
              by = core_id]

# ---- 3. Merge + correlations ----
m <- merge(mfi, core, by = "core_id")
fwrite(m, file.path(out_p, "per_core_xenium_protein_final.csv"))
message("final TMA A cores: ", length(keep), " | with IF data: ", nrow(m))

pairs <- c(list(c("KRT7_xenium", "KRT7"), c("KRT19_xenium", "KRT19")),
           lapply(MARKERS, function(g) c("polarization_mean", g)))
cor_tbl <- rbindlist(lapply(pairs, function(p) {
  ok <- complete.cases(m[[p[1]]], m[[p[2]]]); x <- m[[p[1]]][ok]; y <- m[[p[2]]][ok]
  pe <- cor.test(x, y, method = "pearson"); sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  data.table(x = p[1], y = p[2], n = sum(ok), pearson = unname(pe$estimate), p_pearson = pe$p.value,
             spearman = unname(sp$estimate), p_spearman = sp$p.value)
}))
fwrite(cor_tbl, file.path(out_p, "correlations_final.csv"))
print(cor_tbl[, .(x, y, n, pearson = round(pearson, 2), p_pearson = signif(p_pearson, 2),
                  spearman = round(spearman, 2), p_spearman = signif(p_spearman, 2))])
log_session()
