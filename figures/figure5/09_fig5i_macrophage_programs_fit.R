# Fig 5i — macrophage-intrinsic program GAMs vs the local epithelial niche.
# Same fit spec as 5a/4f: ALL macrophages with a defined niche value (no subsampling),
# tissue-equal cell weights w = (N/n_tissues)/n_tissue, weighted z-score, mgcv::bam(y ~ s(x, k=10), weights, fREML, discrete),
# 200-point grid over the 1st-99th percentile, +-1.96 SE. Fits both axes; Fig 5i and the Figure S18 check use axis nb_pol.
# Output: figures_dir/figure5/fig5e_mac_gam_curves.csv (file name keeps an earlier letter), drawn by 11_fig5_compose_panels.py.
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
set.seed(CFG$seed)
d <- fread(cfg_path("output_root","split_recompute","fig5e_macrophage_percell_presplit.csv"))
fig_dir <- cfg_path("figures_dir","figure5")
N_GRID <- 200; K <- 10
feats  <- c("glycolytic","nfkb_activated","chemokine","tissue_resident")
labels <- c("Glycolytic","NF-κB activated","Chemokine","Tissue-resident identity")
cat("macrophages:", nrow(d), " samples:", uniqueN(d$sample), "\n")
out <- list()
for (ax in c("nb_hyp","nb_pol")) {
  sub <- d[is.finite(get(ax))]
  sub[, w := (nrow(sub) / uniqueN(sub$sample)) / .N, by = sample]   # equal total weight per tissue
  for (f in feats) { v<-sub[[f]]; ok<-is.finite(v); m<-weighted.mean(v[ok],sub$w[ok]); s<-sqrt(weighted.mean((v[ok]-m)^2,sub$w[ok]))
    if (is.finite(s)&&s>0) sub[, (f) := (v-m)/s] }
  cat("\naxis", ax, " n (all cells) =", nrow(sub), "\n")
  for (i in seq_along(feats)) {
    dd <- data.frame(y=sub[[feats[i]]], x=sub[[ax]], w=sub$w)
    dd <- dd[is.finite(dd$y) & is.finite(dd$x), ]
    g <- seq(quantile(dd$x,0.01), quantile(dd$x,0.99), length.out=N_GRID)
    fit <- bam(y ~ s(x, k=K), data=dd, weights=w, method="fREML", discrete=TRUE)
    pr <- predict(fit, newdata=data.frame(x=g), se.fit=TRUE)
    out[[paste(ax,feats[i])]] <- data.table(axis=ax, feature=labels[i], x=g,
      fitted=as.numeric(pr$fit), lower=as.numeric(pr$fit-1.96*pr$se.fit),
      upper=as.numeric(pr$fit+1.96*pr$se.fit))
    cat(sprintf("  %-26s edge %+.2f -> %+.2f\n", labels[i],
                as.numeric(pr$fit)[1], tail(as.numeric(pr$fit),1)))
  }
}
fwrite(rbindlist(out), file.path(fig_dir,"fig5e_mac_gam_curves.csv"))
cat("\nDONE\n")
