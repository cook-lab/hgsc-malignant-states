#!/usr/bin/env Rscript
# Fig 2f prep — multi-seed consensusOV calls on whole-sample pseudobulk
# =====================================================================
# PURPOSE   consensusOV trains a random forest on each call, so single-seed calls vary
#           (~5% of samples between seeds). This repeats the identical pipeline of
#           fig2f_consensusov.R (TMM -> log2(CPM + 1) -> Entrez -> get.subtypes) for
#           N seeds, averages the four random-forest probabilities per sample, and calls
#           the subtype with the highest mean probability. Stability = fraction of seeds
#           whose call equals the consensus call.
# INPUTS    args[1] pseudobulk_counts.tsv.gz (genes x samples, raw summed counts;
#                   written by fig2f_consensusov.py)
#           args[2] output csv; args[3] number of seeds (default 25; seeds 1..N)
# OUTPUTS   sample_id, consensus subtype, mean probabilities, margin (top1 - top2 of the
#           mean probabilities), stability, seed-42 call; and <out>_per_seed.csv.
# MANUSCRIPT PANEL(S): Fig 2f (input).   RUNTIME TIER: fast (~1 min per 10 seeds).
suppressPackageStartupMessages({
  library(edgeR); library(org.Hs.eg.db); library(consensusOV)
})
args <- commandArgs(trailingOnly = TRUE)
in_tsv <- args[1]; out_csv <- args[2]
n_seeds <- if (length(args) >= 3) as.integer(args[3]) else 25L

cnt <- as.matrix(read.delim(gzfile(in_tsv), row.names = 1, check.names = FALSE))
dge <- calcNormFactors(DGEList(counts = cnt), method = "TMM")
logcpm <- log2(cpm(dge, normalized.lib.sizes = TRUE, log = FALSE) + 1)
entrez <- suppressMessages(mapIds(org.Hs.eg.db, keys = rownames(logcpm), keytype = "SYMBOL",
                                  column = "ENTREZID", multiVals = "first"))
keep <- !is.na(entrez) & !duplicated(entrez)
expr <- logcpm[keep, , drop = FALSE]

call_once <- function(seed) {
  set.seed(seed)
  res <- get.subtypes(expr, entrez.ids = entrez[keep], method = "consensusOV")
  p <- as.data.frame(res$rf.probs)
  colnames(p) <- sub("_consensus$", "", colnames(p))
  list(probs = p, call = sub("_consensus$", "", as.character(res$consensusOV.subtypes)))
}
seeds <- unique(c(seq_len(n_seeds), 42L))
runs <- lapply(seeds, call_once)
names(runs) <- seeds
subtypes <- colnames(runs[[1]]$probs)
mean_p <- Reduce(`+`, lapply(runs[as.character(seq_len(n_seeds))], function(r) as.matrix(r$probs))) / n_seeds
consensus <- subtypes[max.col(mean_p, ties.method = "first")]
calls <- sapply(runs[as.character(seq_len(n_seeds))], function(r) r$call)
stability <- rowMeans(calls == consensus)
srt <- t(apply(mean_p, 1, sort, decreasing = TRUE))
out <- data.frame(sample_id = colnames(expr), subtype = consensus, mean_p,
                  margin = round(srt[, 1] - srt[, 2], 3), stability = round(stability, 3),
                  seed42 = runs[["42"]]$call, check.names = FALSE)
write.csv(out, out_csv, row.names = FALSE)
per_seed <- data.frame(sample_id = colnames(expr), calls, check.names = FALSE)
colnames(per_seed)[-1] <- paste0("seed_", seq_len(n_seeds))
write.csv(per_seed, sub("\\.csv$", "_per_seed.csv", out_csv), row.names = FALSE)
cat(sprintf("[R] %d seeds; consensus calls:\n", n_seeds)); print(table(consensus))
cat(sprintf("[R] stability: median %.2f; samples < 0.8: %d\n", median(stability), sum(stability < 0.8)))
