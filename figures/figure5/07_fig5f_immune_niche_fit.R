# Fig 5f -- binomial GAM: P(TME cell is type T | local epithelial hypoxia), fitted per
# sample (each whole tissue and each TMA array) and immune type, k = 6; a fit needs >= 50
# cells of the type. x = z-scored mean epithelial hypoxia within 50um (06_fig5f_immune_niche_build.R);
# the grid spans the 1st-99th percentile of x over all rows (whole tissues and TMA).
# Output: figures_dir/figure5/fig5c_gam_persample.csv (file name keeps an earlier letter). 11_fig5_compose_panels.py
# draws the mean of the whole-tissue fits (dotted lines = the tissues). Manuscript panel: Fig 5f.
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(data.table); library(mgcv)})
set.seed(CFG$seed)
d <- fread(cfg_path("output_root","split_recompute","fig5c_tme_percell_presplit.csv"))
fig_dir <- cfg_path("figures_dir","figure5")
# FIG5C_AXIS = hyp (default; Fig 5f) | pol (writes _pol files that no panel uses). Output names keep hyp unsuffixed.
AXIS <- Sys.getenv("FIG5C_AXIS", "hyp"); XCOL <- paste0("x_", AXIS)
SUF <- if (AXIS == "hyp") "" else paste0("_", AXIS)
d[, x_hyp := get(XCOL)]; d <- d[is.finite(x_hyp)]
TYPES <- c("Macrophage","T cell","NK cell","B cell","Plasma cell")
K <- 6; NG <- 100
XLIM <- quantile(d$x_hyp, c(0.01, 0.99), na.rm = TRUE)
grid <- seq(XLIM[1], XLIM[2], length.out = NG)

fit_one <- function(sub, ct) {
  y <- as.integer(sub$cell_label == ct)
  if (sum(y) < 50) return(NULL)
  m <- tryCatch(gam(y ~ s(x_hyp, k = K), data = data.frame(y = y, x_hyp = sub$x_hyp),
                    family = binomial, method = "REML"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  p <- predict(m, newdata = data.frame(x_hyp = grid), se.fit = TRUE)
  data.table(x = grid, fit = plogis(p$fit),
             lo = plogis(p$fit - 1.96*p$se.fit), hi = plogis(p$fit + 1.96*p$se.fit))
}
persample <- list()
for (co in c("WT","TMA")) {
  sub <- d[cohort == co]
  for (ct in TYPES) {
    for (s in unique(sub$sample)) {
      rs <- fit_one(sub[sample == s], ct)
      if (!is.null(rs)) persample[[paste(co,ct,s)]] <- cbind(cohort=co, cell_type=ct, sample=s, rs)
    }
  }
}
ps <- rbindlist(persample)
fwrite(ps, file.path(fig_dir,paste0("fig5c_gam_persample",SUF,".csv")))
cat("DONE\n")
