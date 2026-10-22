# 18 — inferCNV: confirm clusters 2/3/7/8 are tumour (aneuploid) vs stroma (diploid).
# Reference = immune cells (definitionally diploid). Query = the epithelial/EMT
# candidates + anchor clusters. Subsampled per cluster to keep it tractable.
# Tumour cells carry the KPCA.B clonal copy-number profile; host stroma does not,
# so a clear CNV signal in 2/3/7/8 confirms the EMT-tumour call from step 17.
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(infercnv)
})
set.seed(42)
PROJ <- "/Volumes/CookLab/Sarah/202608_kpca_secb_mouse"
out  <- file.path(PROJ, "output/18_infercnv")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

sce <- readH5AD(file.path(PROJ, "output/15_doublet_cleanup/kpca_all_clean.h5ad"),
                verbose = FALSE)
counts <- as(assay(sce, "counts"), "CsparseMatrix")
rownames(counts) <- rownames(sce)
leiden <- as.character(sce$leiden)

# ---- roles ---------------------------------------------------------------
TUMOUR_CAND <- c("2","3","7","8")           # in question (EMT candidates)
EPI_ANCHOR  <- c("10","13")                  # confident epithelial tumour core
CAF_ANCHOR  <- c("9","12")                   # confident host CAF
IMMUNE_REF  <- c("4","5","16","17")          # macrophage / T / plasma / monocyte = diploid reference

role <- rep(NA_character_, length(leiden))
for (cl in TUMOUR_CAND) role[leiden==cl] <- paste0("cand_",cl)
for (cl in EPI_ANCHOR)  role[leiden==cl] <- paste0("epi_",cl)
for (cl in CAF_ANCHOR)  role[leiden==cl] <- paste0("caf_",cl)
for (cl in IMMUNE_REF)  role[leiden==cl] <- "reference_immune"

# subsample each group to <=400 cells for tractability
keep <- unlist(lapply(unique(na.omit(role)), function(r){
  idx <- which(role==r); if(length(idx)>400) sample(idx,400) else idx
}))
counts <- counts[, keep]; role <- role[keep]
cat("cells:", ncol(counts), " groups:\n"); print(table(role))

ann <- data.frame(role, row.names = colnames(counts))
write.table(ann, file.path(out,"annotations.tsv"), sep="\t", quote=FALSE, col.names=FALSE)

icnv <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = file.path(out,"annotations.tsv"),
  gene_order_file   = file.path(out,"gene_positions.tsv"),
  ref_group_names   = "reference_immune")

icnv <- infercnv::run(
  icnv, cutoff = 0.1, out_dir = file.path(out,"run"),
  cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE,
  num_threads = 1, no_plot = FALSE, output_format = "png")

# ---- per-cell CNV signal: variance of the residual expr profile ----------
expr <- icnv@expr.data                       # genes x cells, reference-centred
cnv_signal <- apply(expr, 2, function(v) sum((v-1)^2)/length(v))  # deviation from ref
res <- data.frame(cell=colnames(expr), role=ann[colnames(expr),"role"],
                  cnv_signal=cnv_signal)
res$group <- sub("_[0-9]+$","",res$role)
write.csv(res, file.path(out,"18_cnv_signal_per_cell.csv"), row.names=FALSE)
cat("\nmean CNV signal (higher = more aneuploid = tumour):\n")
agg <- tapply(res$cnv_signal, res$role, mean)
for(r in names(sort(agg,decreasing=TRUE))) cat(sprintf("  %-18s %.5f (n=%d)\n", r, agg[r], sum(res$role==r)))
cat("\nwrote", file.path(out,"18_cnv_signal_per_cell.csv"), "\n")
