# Fig 5a GAM fit -- epithelium-intrinsic program scores vs focal polarization.
# Spec identical to the Fig 4f fit (figures/figure4/09_fig4f_gam_fit.R):
# ALL secretory cells (no subsampling), tissue-equal cell weights w = (N/n_tissues)/n_tissue, weighted per-feature z-score,
# mgcv::bam(y ~ s(x, k=10), weights, fREML, discrete), 200-point grid over the 1st-99th percentile of polarization,
# +-1.96 SE. Output: figures_dir/figure5/gam_curves_5a.csv, drawn by 11_fig5_compose_panels.py (Fig 5a) and used as the
# pooled-vs-main check of Figure S17 (figures/supplementary/21_figS16-18_per_sample_gams.R).
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
set.seed(CFG$seed)
IN  <- cfg_path("output_root","split_recompute","fig5a_epi_program_scores_presplit.csv")
fig_dir <- cfg_path("figures_dir","figure5"); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
N_GRID <- 200; K_POOL <- 10

features <- c("nfkb","glycolysis","LDHA","TNFA",
              "oncogenic_growth","ddr","epi_silencing","proliferation")
labels   <- c("NF-κB","Glycolysis","LDHA","TNFA",
              "Oncogenic growth","DNA damage repair","Epigenetic silencing","Proliferation")

nf <- fread(IN)
cat("secretory epithelial cells:", nrow(nf), " samples:", uniqueN(nf$sample), "\n")
nf_sub <- nf   # all cells
nf_sub[, w := (nrow(nf_sub) / uniqueN(nf_sub$sample)) / .N, by = sample]   # equal total weight per tissue
cat("cells used:", nrow(nf_sub), "\n")
for (f in features) { v<-nf_sub[[f]]; ok<-is.finite(v); m<-weighted.mean(v[ok],nf_sub$w[ok]); s<-sqrt(weighted.mean((v[ok]-m)^2,nf_sub$w[ok]))
  if (is.finite(s) && s>0) nf_sub[, (f) := (v-m)/s] }

pooled_all <- list()
for (i in seq_along(features)) {
  d <- data.frame(y = nf_sub[[features[i]]], x = nf_sub$polarization, w = nf_sub$w)
  d <- d[!is.na(d$y) & !is.na(d$x) & is.finite(d$y), ]
  g <- seq(quantile(d$x, 0.01), quantile(d$x, 0.99), length.out = N_GRID)
  fit <- bam(y ~ s(x, k = K_POOL), data = d, weights = w, method = "fREML", discrete = TRUE)
  pr  <- predict(fit, newdata = data.frame(x = g), se.fit = TRUE)
  pooled_all[[i]] <- data.table(polarization=g, fitted=as.numeric(pr$fit),
    lower=as.numeric(pr$fit-1.96*pr$se.fit), upper=as.numeric(pr$fit+1.96*pr$se.fit),
    feature=labels[i])
  cat(sprintf("  %-22s OK (n=%d)  edge=%+.2f\n", labels[i], nrow(d),
              tail(as.numeric(pr$fit),1)))
}
pooled <- rbindlist(pooled_all, fill=TRUE)
pooled[, feature := factor(feature, levels=labels)]
fwrite(pooled,   file.path(fig_dir,"gam_curves_5a.csv"))
cat("\nwrote gam_curves_5a.csv\nDONE\n")
