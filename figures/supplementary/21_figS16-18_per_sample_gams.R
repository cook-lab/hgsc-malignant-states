#!/usr/bin/env Rscript
# =============================================================================
# Figures S16 / S17 / S18 — per-tissue GAM consistency (final data, SingleR labels)
# -----------------------------------------------------------------------------
# PURPOSE (internal figure ids SF14 / SF15 / SF16 = Figures S16 / S17 / S18; the plot script
# 22_figS16-18_per_sample_gams_plot.py writes the final numbers)
#   SF14 = S16  environment vs focal polarization (as Fig 4f)          secretory cells
#   SF15 = S17  epithelium-intrinsic vs focal polarization (as Fig 5a) secretory cells
#   SF16 = S18  macrophage-intrinsic vs local epithelial polarization (as Fig 5i, x = nb_pol)
#
# NO SUBSAMPLING. The pooled curve is fit exactly as the final main-figure fits
# (figures/figure4/09_fig4f_gam_fit.R, figures/figure5/02_fig5a_epithelial_programs_fit.R,
# figures/figure5/09_fig5i_macrophage_programs_fit.R): all cells, tissue-equal weights
# w = (N / n_tissues) / n_tissue, weighted per-feature z-score, mgcv::bam(y ~ s(x, k = 10),
# weights, fREML, discrete = TRUE), 200-point grid over the pooled 1st-99th percentile of x, +-1.96 SE.
# Per-sample curves: the SAME z-scaled y (one scale), bam(y ~ s(x, k = 10), fREML, discrete)
# on all cells of that sample, 200-point grid over that sample's own 1st-99th percentile of x.
#
# Per-cell inputs (all cells):
#   singler_recompute/singler_gam_percell.csv            SF14 (Fig 4f features)
#   split_recompute/fig5a_epi_program_scores_presplit.csv SF15 programs (proliferation, ddr, nfkb, glycolysis, LDHA)
#   split_recompute/fig5e_macrophage_percell_presplit.csv SF16 programs (glycolytic, nfkb_activated) + nb_pol
#   SFE colData pathway_myc/pi3k/wnt_activating (UCell, 477-gene panel, all cells; 04c_pathway_scoring_allcells.R)
#   single genes not in the files: log1p(CP10k) from counts(sfe) (library size = all 477 panel genes)
#
# Per-panel consistency (stats CSV): k_rise = n tissues whose per-sample curve rises
# (fitted at right edge - left edge > 0); paired two-sided Wilcoxon signed-rank of the per-tissue
# mean y in the top vs bottom 20% of x (within tissue), across the 8 tissues (as the Fig 5i stats).
# SF15 also fits TGFB1, INHBA, C7, CDH2 and ITGB5, which the plot script drops before drawing
# Figure S17 and writing the S16-S18 statistics table.
#
# OUTPUTS (output_root/supplementary/sf14_16_per_sample_gams/):
#   extra_epi_percell.csv.gz, extra_mac_percell.csv.gz (cache), curves_pooled.csv, curves_per_sample.csv,
#   per_sample_consistency.csv, pooled_vs_main.csv, sample_meta.csv (de-identified codes only);
#   all read by 22_figS16-18_per_sample_gams_plot.py
# MANUSCRIPT PANEL(S): Figures S16a-d, S17a-n, S18a-h (curves and per-panel statistics).
# RUNTIME TIER: slow (bam fits on all cells; the first run reads counts from the 8 whole-tissue SFEs).
# =============================================================================
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({
  library(data.table); library(mgcv); library(SummarizedExperiment); library(Matrix)
})
set.seed(CFG$seed)
R   <- cfg_path("output_root")
OUT <- file.path(R, "supplementary", "sf14_16_per_sample_gams"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_GRID <- 200; K <- 10

# ---------------------------------------------------------------- sample key --
source(file.path("config", "study_codes.R"))   # study codes via the private key (never written out)
key <- data.table(source_id = CFG$cohort$whole_tissue)[, study_code := study_code_of("Whole tissue", source_id)]
st8 <- fread(file.path(Sys.getenv("OUTPUT_ROOT", "output"), "final_panels", "supplementary", "tables",
                       "SuppTable8_xenium_clinical_cohort_Whole_tissues.csv"))
key <- merge(key, st8[, .(study_code = `Tissue code`, treatment = `Treatment status at sampling`)], by = "study_code")
key[, sample := paste0("sfe_", source_id)]
stopifnot(nrow(key) == 8)
code_of <- setNames(key$study_code, key$sample)
fwrite(key[order(study_code), .(code = study_code, treatment)], file.path(OUT, "sample_meta.csv"))

# ------------------------------------------------------ per-cell input tables --
env <- fread(file.path(R, "singler_recompute", "singler_gam_percell.csv"))
epi <- fread(file.path(R, "split_recompute", "fig5a_epi_program_scores_presplit.csv"))
mac <- fread(file.path(R, "split_recompute", "fig5e_macrophage_percell_presplit.csv"))
cat("env", nrow(env), " epi", nrow(epi), " mac", nrow(mac), "\n")

EPI_PATH  <- c("pathway_myc_activating", "pathway_pi3k_activating", "pathway_wnt_activating")
EPI_GENES <- c("TGFB1","INHBA","IL32","CD55","C7","CDH2","CTNNB1","ITGB5","ICAM1","MMP7","TGM2")
MAC_GENES <- c("A2M","C1QC","CXCL10","CXCL11","CTSL","ICAM1")

f_epi <- file.path(OUT, "extra_epi_percell.csv.gz"); f_mac <- file.path(OUT, "extra_mac_percell.csv.gz")
if (!file.exists(f_epi) || !file.exists(f_mac)) {
  xe <- list(); xm <- list()
  for (nm in key$sample) {
    t0 <- Sys.time(); message("[", format(t0), "] ", nm)
    sfe <- load_sfe(nm)
    C   <- as(counts(sfe), "CsparseMatrix")
    lib <- Matrix::colSums(C)                                   # all 477 panel genes
    stopifnot(all(c(EPI_GENES, MAC_GENES) %in% rownames(C)), all(EPI_PATH %in% colnames(colData(sfe))))
    cp  <- function(genes, j) { M <- as.matrix(t(C[genes, j, drop = FALSE])); log1p(sweep(M, 1, lib[j], "/") * 1e4) }
    # epithelium (Fig 5a cell set)
    b  <- epi[sample == nm, barcode]; j <- match(b, colnames(sfe)); stopifnot(!anyNA(j))
    de <- data.table(sample = nm, barcode = b,
                     as.data.table(as.data.frame(colData(sfe)[j, EPI_PATH, drop = FALSE])),
                     as.data.table(cp(EPI_GENES, j)))
    stopifnot(!anyNA(de$pathway_myc_activating))
    xe[[nm]] <- de
    # macrophages (Fig 5e/h cell set)
    b  <- mac[sample == nm, barcode]; j <- match(b, colnames(sfe)); stopifnot(!anyNA(j))
    xm[[nm]] <- data.table(sample = nm, barcode = b, as.data.table(cp(MAC_GENES, j)))
    message(sprintf("   epi %d  mac %d  (%.1f min)", nrow(de), length(j),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    rm(sfe, C, lib); gc(verbose = FALSE)
  }
  fwrite(rbindlist(xe), f_epi); fwrite(rbindlist(xm), f_mac)
}
xe <- fread(f_epi); xm <- fread(f_mac)
setnames(xe, EPI_PATH, c("myc","pi3k","wnt"))
setnames(xm, MAC_GENES, paste0("mac_", MAC_GENES))
epi <- merge(epi, xe, by = c("sample","barcode"), sort = FALSE); stopifnot(nrow(epi) == nrow(xe))
mac <- merge(mac, xm, by = c("sample","barcode"), sort = FALSE); stopifnot(nrow(mac) == nrow(xm))

# ------------------------------------------------------------- panel specs --
# id, column, label, x column; order = panel letter order
SPEC <- list(
  SF14 = list(dt = "env", x = "polarization", feats = data.table(
    col   = c("cell_density","dist_to_vascular","nb_matrix","nb_hypoxia"),
    label = c("Cell density","Distance to nearest vascular cell (µm)","Matrix remodelling","Hypoxia"),
    main  = c("Cell density","Distance to vasculature","Matrix remodeling","Hypoxia"))),
  SF15 = list(dt = "epi", x = "polarization", feats = data.table(
    col   = c("proliferation","myc","pi3k","wnt","ddr","nfkb","glycolysis",
              "TGFB1","INHBA","IL32","CD55","C7","LDHA","CDH2","CTNNB1","ITGB5","ICAM1","MMP7","TGM2"),
    label = c("Proliferation","Myc","PI3K","Wnt","DNA damage repair","NF-κB","Glycolysis",
              "TGFB1","INHBA","IL32","CD55","C7","LDHA","CDH2","CTNNB1","ITGB5","ICAM1","MMP7","TGM2"),
    main  = c("Proliferation",NA,NA,NA,"DNA damage repair","NF-κB","Glycolysis",
              NA,NA,NA,NA,NA,"LDHA",NA,NA,NA,NA,NA,NA))),
  SF16 = list(dt = "mac", x = "nb_pol", feats = data.table(
    col   = c("glycolytic","nfkb_activated", paste0("mac_", MAC_GENES)),
    label = c("Glycolytic","NF-κB activated", MAC_GENES),
    main  = c("Glycolytic","NF-κB activated", rep(NA, 6))))
)
SOURCE <- c(cell_density = "all cells within 50 µm (per mm²)", dist_to_vascular = "µm to nearest pericyte/endothelial",
            nb_matrix = "neighbourhood UCell", nb_hypoxia = "neighbourhood UCell",
            proliferation = "UCell (SFE pathway_proliferation)", myc = "UCell (SFE pathway_myc_activating)",
            pi3k = "UCell (SFE pathway_pi3k_activating)", wnt = "UCell (SFE pathway_wnt_activating)",
            ddr = "UCell custom BRCA1/BRCA2/ATM/CHEK2 (01_fig5a_epithelial_programs_build.R)", nfkb = "UCell (SFE pathway_nfkb)",
            glycolysis = "UCell (SFE pathway_glycolysis)", LDHA = "logcounts(sfe) (as Fig 5a)",
            glycolytic = "UCell (Fig 5e macrophage program)", nfkb_activated = "UCell (Fig 5e macrophage program)")

wz <- function(v, w) { ok <- is.finite(v); m <- weighted.mean(v[ok], w[ok]); s <- sqrt(weighted.mean((v[ok]-m)^2, w[ok]))
  if (is.finite(s) && s > 0) (v - m) / s else v }
fitgrid <- function(y, x, w = NULL) {
  d <- data.frame(y = y, x = x); if (!is.null(w)) d$w <- w
  d <- d[is.finite(d$y) & is.finite(d$x), ]
  g <- seq(quantile(d$x, 0.01), quantile(d$x, 0.99), length.out = N_GRID)
  fit <- if (is.null(w)) bam(y ~ s(x, k = K), data = d, method = "fREML", discrete = TRUE)
         else bam(y ~ s(x, k = K), data = d, weights = w, method = "fREML", discrete = TRUE)
  pr <- predict(fit, newdata = data.frame(x = g), se.fit = TRUE)
  list(tab = data.table(x = g, fitted = as.numeric(pr$fit), lower = as.numeric(pr$fit - 1.96*pr$se.fit),
                        upper = as.numeric(pr$fit + 1.96*pr$se.fit)), n = nrow(d))
}

curves_pool <- list(); curves_samp <- list(); stats <- list()
for (fig in names(SPEC)) {
  sp <- SPEC[[fig]]
  D  <- copy(get(sp$dt))
  if (sp$dt == "mac") D <- D[is.finite(nb_pol)]                     # as 09_fig5i_macrophage_programs_fit.R (nb_pol axis)
  D[, w := (nrow(D) / uniqueN(D$sample)) / .N, by = sample]        # tissue-equal weights (on the fit set)
  for (i in seq_len(nrow(sp$feats))) {
    col <- sp$feats$col[i]; lab <- sp$feats$label[i]; letter <- LETTERS[i]
    raw <- D[[col]]
    z <- wz(raw, D$w)
    t0 <- Sys.time()
    P <- fitgrid(z, D[[sp$x]], D$w)
    curves_pool[[paste(fig, col)]] <- cbind(data.table(figure = fig, panel = letter, feature = lab), P$tab)
    ks <- c(); hi <- c(); lo <- c()
    for (s in key[order(study_code), sample]) {
      ii <- which(D$sample == s)
      S  <- fitgrid(z[ii], D[[sp$x]][ii])
      curves_samp[[paste(fig, col, s)]] <- cbind(data.table(figure = fig, panel = letter, feature = lab,
                                                            code = code_of[[s]]), S$tab)
      ks[code_of[[s]]] <- tail(S$tab$fitted, 1) - S$tab$fitted[1]
      xs <- D[[sp$x]][ii]; ys <- z[ii]; ok <- is.finite(xs) & is.finite(ys); xs <- xs[ok]; ys <- ys[ok]
      q <- quantile(xs, c(0.2, 0.8))
      lo[code_of[[s]]] <- mean(ys[xs <= q[1]]); hi[code_of[[s]]] <- mean(ys[xs >= q[2]])
    }
    wt <- wilcox.test(hi, lo, paired = TRUE, exact = TRUE)
    pe <- P$tab$fitted
    stats[[paste(fig, col)]] <- data.table(figure = fig, panel = letter, feature = lab,
      source = if (!is.na(SOURCE[col])) SOURCE[[col]] else paste0("log1p CP10k, counts(sfe)",
                                                                   if (fig == "SF16") ", macrophages" else ", secretory epithelium"),
      x_axis = if (fig == "SF16") "local epithelial polarization (nb_pol, z)" else "focal polarization",
      n_cells = P$n, pooled_left = pe[1], pooled_right = tail(pe, 1), pooled_delta = tail(pe, 1) - pe[1],
      k_rise = sum(ks > 0), k_fall = sum(ks < 0),
      k_top_gt_bottom = sum(hi > lo), median_top_minus_bottom = median(hi - lo),
      wilcoxon_p = wt$p.value,
      per_sample_delta = paste(sprintf("%s:%+.2f", names(ks), ks), collapse = ";"))
    cat(sprintf("%s %s %-40s n=%8d pooled %+.2f  rise %d/8  top>bot %d/8  p=%.4f  (%.1fs)\n", fig, letter, lab, P$n,
                tail(pe,1) - pe[1], sum(ks > 0), sum(hi > lo), wt$p.value, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
CP <- rbindlist(curves_pool); CS <- rbindlist(curves_samp); ST <- rbindlist(stats)
fwrite(CP, file.path(OUT, "curves_pooled.csv")); fwrite(CS, file.path(OUT, "curves_per_sample.csv"))
fwrite(ST, file.path(OUT, "per_sample_consistency.csv"))

# -------------------------------------------- pooled vs main-figure curves --
m4 <- fread(file.path(R, "figures", "figure4", "gam_curves.csv"))[, .(feature, x = polarization, fitted)]
m5 <- fread(file.path(R, "figures", "figure5", "gam_curves_5a.csv"))[, .(feature, x = polarization, fitted)]
m5h <- fread(file.path(R, "figures", "figure5", "fig5e_mac_gam_curves.csv"))[axis == "nb_pol", .(feature, x, fitted)]
MAIN <- list(SF14 = m4, SF15 = m5, SF16 = m5h)
cmp <- list()
for (fig in names(SPEC)) for (i in seq_len(nrow(SPEC[[fig]]$feats))) {
  mf <- SPEC[[fig]]$feats$main[i]; if (is.na(mf)) next
  a <- CP[figure == fig & feature == SPEC[[fig]]$feats$label[i]]; b <- MAIN[[fig]][feature == mf]
  stopifnot(nrow(b) == N_GRID)
  cmp[[length(cmp)+1]] <- data.table(figure = fig, feature = SPEC[[fig]]$feats$label[i], main_feature = mf,
    max_abs_dx = max(abs(a$x - b$x)), max_abs_dfit = max(abs(a$fitted - b$fitted)))
}
CMP <- rbindlist(cmp); fwrite(CMP, file.path(OUT, "pooled_vs_main.csv")); print(CMP)
cat("\nDONE\n")
