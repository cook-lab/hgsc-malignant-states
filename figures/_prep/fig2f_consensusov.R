#!/usr/bin/env Rscript
# Fig 2f prep — consensusOV calls on whole-sample pseudobulk (called by fig2f_consensusov.py)
# =========================================================================================
# PURPOSE   edgeR TMM normalisation -> log2(CPM + 1) -> Entrez IDs (org.Hs.eg.db)
#           -> consensusOV::get.subtypes(method = "consensusOV").
#           Genes without an Entrez ID are dropped; for duplicated Entrez IDs the
#           first symbol (atlas var order) is kept.
# INPUTS    args[1] pseudobulk_counts.tsv.gz (genes x samples, raw summed counts)
# OUTPUTS   args[2] consensusov_calls.csv: sample_id, subtype (DIF/IMR/MES/PRO),
#           the four random-forest probabilities, margin_top1_top2;
#           <dir>/gene_id_mapping.csv (symbol -> entrez used)
# MANUSCRIPT PANEL(S): Fig 2f (input; single seed 42). The Fig 2f subtypes are the 25-seed consensus of
#           fig2f_consensusov_multiseed.R (same pipeline); these calls fill the single-seed column of
#           fig2f_per_sample_joined.csv.   RUNTIME TIER: fast (< 1 min).
suppressPackageStartupMessages({
  library(edgeR); library(org.Hs.eg.db); library(consensusOV)
})
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
in_tsv <- args[1]; out_csv <- args[2]

cnt <- read.delim(gzfile(in_tsv), row.names = 1, check.names = FALSE)
cnt <- as.matrix(cnt)
cat(sprintf("[R] counts: %d genes x %d samples\n", nrow(cnt), ncol(cnt)))

dge <- calcNormFactors(DGEList(counts = cnt), method = "TMM")
logcpm <- log2(cpm(dge, normalized.lib.sizes = TRUE, log = FALSE) + 1)

entrez <- suppressMessages(mapIds(org.Hs.eg.db, keys = rownames(logcpm), keytype = "SYMBOL",
                                  column = "ENTREZID", multiVals = "first"))
map <- data.frame(symbol = rownames(logcpm), entrez = unname(entrez))
write.csv(map, file.path(dirname(out_csv), "gene_id_mapping.csv"), row.names = FALSE)
keep <- !is.na(entrez) & !duplicated(entrez)
expr <- logcpm[keep, , drop = FALSE]
cat(sprintf("[R] genes with unique Entrez IDs: %d / %d\n", nrow(expr), nrow(logcpm)))

res <- get.subtypes(expr, entrez.ids = entrez[keep], method = "consensusOV")
probs <- as.data.frame(res$rf.probs)
srt <- t(apply(as.matrix(probs), 1, sort, decreasing = TRUE))
out <- data.frame(sample_id = colnames(expr),
                  subtype = sub("_consensus$", "", as.character(res$consensusOV.subtypes)),
                  probs, margin_top1_top2 = round(srt[, 1] - srt[, 2], 3),
                  check.names = FALSE)
write.csv(out, out_csv, row.names = FALSE)
print(table(out$subtype))
