# Recovered copy (2026-10-22). This code ran inline (an `Rscript -` heredoc) from the mouse project folder
# 202608_kpca_secb_mouse/ on 2026-09-04 21:22 UTC and was never saved as a file; it is recovered from the analysis
# session log. Paths are relative to the mouse project folder, so run it from there.
# Writes output/18_infercnv/gene_positions.tsv, the gene order file of 18_infercnv.R. The shell command created the
# output folder first (mkdir -p output/18_infercnv); the dir.create line does that here. The rest is unchanged.
# The run printed: "wrote 19123 gene positions".

dir.create("output/18_infercnv", showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({library(org.Mm.eg.db); library(AnnotationDbi)})
# map Ensembl gene IDs -> chr, start, end (offline, from org.Mm.eg.db)
library(zellkonverter); library(SingleCellExperiment)
sce <- readH5AD("output/15_doublet_cleanup/kpca_all_clean.h5ad", verbose=FALSE)
ens <- as.character(rowData(sce)$gene_ids); sym <- rownames(sce)
names(sym) <- ens
chr <- mapIds(org.Mm.eg.db, ens, "CHR", "ENSEMBL")
st  <- suppressWarnings(mapIds(org.Mm.eg.db, ens, "CHRLOC", "ENSEMBL"))
en  <- suppressWarnings(mapIds(org.Mm.eg.db, ens, "CHRLOCEND", "ENSEMBL"))
pos <- data.frame(gene=sym, chr=paste0("chr",chr), start=abs(as.numeric(st)), stop=abs(as.numeric(en)),
                  stringsAsFactors=FALSE)
pos <- pos[!is.na(pos$chr) & !is.na(pos$start) & pos$chr %in% paste0("chr",c(1:19,"X","Y")),]
pos <- pos[!duplicated(pos$gene),]
pos <- pos[order(match(pos$chr,paste0("chr",c(1:19,"X","Y"))), pos$start),]
write.table(pos[,c("gene","chr","start","stop")], "output/18_infercnv/gene_positions.tsv",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("wrote", nrow(pos), "gene positions\n")
cat("chr coverage:", paste(names(table(pos$chr))[1:5],collapse=","), "...\n")
