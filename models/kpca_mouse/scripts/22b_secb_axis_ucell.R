# Recovered copy (2026-10-22). This code ran inline (an `Rscript -` heredoc, with KMP_DUPLICATE_LIB_OK=TRUE) from
# the mouse project folder 202608_kpca_secb_mouse/ on 2026-09-05 00:05 UTC and was never saved as a file; it is
# recovered from the analysis session log. Paths are relative to the mouse project folder, so run it from there.
# Writes output/22_secb_axis/22_ucell_scores.csv, the tumour scores that step 24 checks its own scores against. The
# code below is unchanged.
# The run printed: "UCell done on 7752 tumour cells".

suppressPackageStartupMessages({library(UCell); library(Matrix)})
X <- as(Matrix::readMM("output/22_secb_axis/tumour_X.mtx"),"CsparseMatrix")   # genes x cells
rownames(X) <- read.csv("output/22_secb_axis/genes.csv")$gene
colnames(X) <- read.csv("output/22_secb_axis/cells.csv")$cell
seca <- read.csv("output/22_secb_axis/seca_genes.csv")$gene
secb <- read.csv("output/22_secb_axis/secb_genes.csv")$gene
sets <- list(SecA=seca, SecB=secb)
# UCell expects genes x cells; X is lognorm (UCell ranks internally, robust to that)
sc <- ScoreSignatures_UCell(X, features=sets, ncores=1)
res <- data.frame(cell=rownames(sc), SecA_UCell=sc[,"SecA_UCell"], SecB_UCell=sc[,"SecB_UCell"])
res$polarization <- res$SecB_UCell - res$SecA_UCell
write.csv(res, "output/22_secb_axis/22_ucell_scores.csv", row.names=FALSE)
cat("UCell done on", nrow(res), "tumour cells\n")
cat(sprintf("SecA_UCell mean=%.3f  SecB_UCell mean=%.3f  polarization mean=%.3f\n",
            mean(res$SecA_UCell), mean(res$SecB_UCell), mean(res$polarization)))
cat("polarization quantiles:\n"); print(round(quantile(res$polarization, c(.05,.25,.5,.75,.95)),3))
