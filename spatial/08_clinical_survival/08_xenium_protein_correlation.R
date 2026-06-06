# ============================================================================
# 08_xenium_protein_correlation.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-TMA-core correlation of Xenium transcript / SecB% / mean polarization against immunofluorescence protein MFI.
#
# INPUTS:
#   - data/ MFI all tissues TMA xlsx + sfe_tma_filtered (logcounts, polarization_UCell)
#   - output/38_FTE_baseline_TMA/per_core_proportions_wide.csv
#
# OUTPUTS:
#   - output/41_xenium_protein_correlation/per_core_xenium_protein.csv
#   - correlations.csv + figures
#
# MANUSCRIPT PANEL(S): Fig 7C, Fig 7D.
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
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

out_p <- file.path(out_dir, "41_xenium_protein_correlation")
dir.create(out_p, showWarnings = FALSE, recursive = TRUE)

# ---- 1. MFI protein -------------------------------------------------------
mfi_path <- file.path(data_dir, "Copy of MFI all tissues TMA valid all markers_.xlsx")
mfi <- suppressWarnings(read_excel(mfi_path, sheet = "Pour SPSS"))
mfi <- mfi |>
  rename(core_id = ID, ECAD = `E-CADH`) |>
  mutate(core_id = as.character(core_id)) |>
  filter(!is.na(core_id))

protein_markers <- c("VIM", "ECAD", "KRT7", "KRT18", "KRT19")
message("MFI rows: ", nrow(mfi),
        "  | non-missing per marker: ",
        paste0(protein_markers, "=", colSums(!is.na(mfi[, protein_markers])),
               collapse = ", "))

# ---- 2. Per-core SecB% (precomputed) --------------------------------------
prop_path <- file.path(out_dir, "38_FTE_baseline_TMA", "per_core_proportions_wide.csv")
prop <- read.csv(prop_path, check.names = FALSE) |>
  transmute(core_id      = as.character(core_id),
            patient_id   = as.character(patient_id),
            sample_type  = sample_type,
            n_cells_core = total,
            secb_pct     = `SecB epithelium`,
            seca_pct     = `SecA epithelium`,
            trans_pct    = `Transitioning epithelium`)

# ---- 3. Per-core mean Xenium transcripts + polarization -------------------
sfe <- load_sfe("sfe_tma_filtered")
xen_genes <- c("KRT7", "KRT19")
stopifnot(all(xen_genes %in% rownames(sfe)))

logc  <- as(assay(sfe, "logcounts")[xen_genes, , drop = FALSE], "CsparseMatrix")
core  <- as.character(sfe$core_id)
polar <- sfe$polarization_UCell

# Mean logcounts per core for each gene (sum / n_cells)
core_factor <- factor(core)
n_per_core  <- as.integer(table(core_factor))
sum_per_core <- as.matrix(logc %*% Matrix::sparse.model.matrix(~ 0 + core_factor))
colnames(sum_per_core) <- levels(core_factor)
mean_per_core <- sweep(sum_per_core, 2, n_per_core, "/")

xen_expr <- as.data.frame(t(mean_per_core)) |>
  tibble::rownames_to_column("core_id") |>
  rename(KRT7_xenium = KRT7, KRT19_xenium = KRT19)

# Mean polarization (na.rm — non-secretory cells lack scores)
polar_mean <- tapply(polar, core_factor, mean, na.rm = TRUE)
polar_df <- data.frame(core_id = names(polar_mean),
                       polarization_mean = as.numeric(polar_mean))

# ---- 4. Merge -------------------------------------------------------------
merged <- mfi |>
  inner_join(prop,     by = "core_id") |>
  inner_join(xen_expr, by = "core_id") |>
  inner_join(polar_df, by = "core_id")

message("Cores merged (MFI ∩ Xenium): ", nrow(merged))
write.csv(merged,
          file.path(out_p, "per_core_xenium_protein.csv"),
          row.names = FALSE)

# ---- 5. Correlations + scatter helper -------------------------------------
cor_pair <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(c(n = sum(ok), spearman = NA, p_spearman = NA,
                            pearson = NA, p_pearson = NA))
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  p <- cor.test(x[ok], y[ok], method = "pearson")
  c(n = sum(ok),
    spearman = unname(s$estimate),
    p_spearman = s$p.value,
    pearson = unname(p$estimate),
    p_pearson = p$p.value)
}

scatter_cor <- function(df, x_col, y_col, x_lab, y_lab, title) {
  d <- df[, c(x_col, y_col)]
  names(d) <- c("x", "y")
  d <- d[complete.cases(d), ]
  if (nrow(d) < 5) return(ggplot() + theme_void() +
                            ggtitle(paste0(title, "\n(n<5)")))
  s  <- suppressWarnings(cor.test(d$x, d$y, method = "spearman"))
  pe <- cor.test(d$x, d$y, method = "pearson")
  sub <- sprintf("n=%d  rho=%.2f (p=%.2g)  r=%.2f (p=%.2g)",
                 nrow(d), s$estimate, s$p.value, pe$estimate, pe$p.value)
  ggplot(d, aes(x, y)) +
    geom_point(shape = 16, size = 1.2, alpha = 0.75, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4,
                color = "#D14E6C", formula = y ~ x) +
    labs(x = x_lab, y = y_lab, title = title, subtitle = sub) +
    theme_lab()
}

# ---- 6. Build correlation table ------------------------------------------
cor_rows <- list()
add_cor <- function(label, x_col, y_col) {
  v <- cor_pair(merged[[x_col]], merged[[y_col]])
  cor_rows[[label]] <<- data.frame(
    comparison = label,
    x = x_col, y = y_col,
    n = v["n"], spearman = v["spearman"], p_spearman = v["p_spearman"],
    pearson = v["pearson"], p_pearson = v["p_pearson"],
    row.names = NULL
  )
}

# Transcript vs protein (KRT7, KRT19 only)
for (g in c("KRT7", "KRT19"))
  add_cor(paste0(g, "_xenium_vs_", g, "_protein"),
          paste0(g, "_xenium"), g)

# SecB% vs protein (5 markers)
for (g in protein_markers)
  add_cor(paste0("SecB_pct_vs_", g, "_protein"), "secb_pct", g)

# Polarization mean vs protein (5 markers)
for (g in protein_markers)
  add_cor(paste0("polarization_mean_vs_", g, "_protein"),
          "polarization_mean", g)

cor_tbl <- do.call(rbind, cor_rows)
write.csv(cor_tbl, file.path(out_p, "correlations.csv"), row.names = FALSE)
print(cor_tbl)

# ---- 7. Figures -----------------------------------------------------------
# Fig 1: transcript vs protein (KRT7, KRT19)
p_krt7  <- scatter_cor(merged, "KRT7_xenium",  "KRT7",
                       "KRT7 mean logcounts (Xenium)",
                       "KRT7 protein MFI",  "KRT7")
p_krt19 <- scatter_cor(merged, "KRT19_xenium", "KRT19",
                       "KRT19 mean logcounts (Xenium)",
                       "KRT19 protein MFI", "KRT19")
fig1 <- (p_krt7 | p_krt19) +
  plot_annotation(title = "Xenium transcript vs IF protein (per TMA core)")
ggsave(file.path(out_p, "fig_transcript_vs_protein.pdf"),
       fig1, width = 7, height = 3.4)
ggsave(file.path(out_p, "fig_transcript_vs_protein.png"),
       fig1, width = 7, height = 3.4, dpi = 300)

# Fig 2: SecB% vs each protein
mk_protein_panels <- function(x_col, x_lab, fname, fig_title) {
  pls <- lapply(protein_markers, function(g) {
    scatter_cor(merged, x_col, g,
                x_lab,
                paste0(g, " protein MFI"),
                g)
  })
  fig <- wrap_plots(pls, ncol = 5) +
    plot_annotation(title = fig_title)
  ggsave(file.path(out_p, paste0(fname, ".pdf")),
         fig, width = 18, height = 4)
  ggsave(file.path(out_p, paste0(fname, ".png")),
         fig, width = 18, height = 4, dpi = 300)
}

mk_protein_panels("secb_pct",
                  "SecB epithelium (% of core cells, Xenium)",
                  "fig_secb_vs_protein",
                  "SecB% (Xenium) vs IF protein MFI per core")

mk_protein_panels("polarization_mean",
                  "polarization_UCell mean (Xenium)",
                  "fig_polarization_vs_protein",
                  "Mean SecA-to-SecB polarization (Xenium) vs IF protein MFI per core")

message("Done. Outputs in: ", out_p)
