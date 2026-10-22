# Recovered copy (2026-10-22). Figure S9g input: the per-cell consensusOV subtype and cell type of the five tumours
# of Geistlinger et al. (their published SingleCellExperiment objects sample<59|76|77|89|90>_sce.rds), with a join key
# tumour_barcode. It ran inline (`Rscript -e`) on 2026-06-12 04:52 UTC from
# 202601_final_atlas/output/spectrum_epithelial_crosstab/geistlinger_sce/ and was never saved as a file; it is
# recovered from the analysis session log. Run it from that folder. The code is unchanged.
# Output: ../geistlinger_cell_annotations.csv (read by figures/supplementary/15_figS9_composition_by_metadata.py).
# The run printed: "Total cells: 42253".

suppressMessages({library(SingleCellExperiment); library(SummarizedExperiment)})
tumors <- c(59,76,77,89,90)
out <- list()
for (t in tumors) {
  sce <- readRDS(paste0("sample",t,"_sce.rds"))
  cd <- as.data.frame(colData(sce))
  cat("T",t," cols: ", paste(colnames(cd),collapse=", "), "\n", sep="")
  bc <- if ("Barcode" %in% colnames(cd)) cd$Barcode else colnames(sce)
  df <- data.frame(tumor=paste0("T",t), Barcode=as.character(bc),
                   celltype=as.character(cd$celltype),
                   subtype=as.character(cd$subtype), stringsAsFactors=FALSE)
  out[[as.character(t)]] <- df
}
res <- do.call(rbind, out)
res$join_key <- paste0(res$tumor, "_", res$Barcode)
write.csv(res, "../geistlinger_cell_annotations.csv", row.names=FALSE)
cat("\nTotal cells:", nrow(res), "\n")
cat("celltype:\n"); print(table(res$celltype, useNA="ifany"))
cat("subtype:\n"); print(table(res$subtype, useNA="ifany"))
cat("subtype within EPI:\n"); print(table(res$subtype[res$celltype=="EPI"], useNA="ifany"))
