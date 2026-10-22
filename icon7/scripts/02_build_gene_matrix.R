# Trimmed copy (2026-10-21) of 2026_final_xenium_analysis/davids side quests/ICON7/scripts/02_build_gene_matrix.R,
# the script that ran. Kept: expr_gene_symbol.rds (read by 03; gene presence in Table S6) and the probe mapping it was
# built from (probe_annotation.tsv). Removed: the signature coverage report (01_signature_coverage.tsv).
# The untrimmed script is in the lab project folder.
# ============================================================================
# 02_build_gene_matrix.R
# Build a gene-symbol-level expression matrix from the GSE140082 normalized
# probe-level expression file. Strips out the interleaved detection p-value
# columns and collapses multiple probes per gene by max-mean expression.
#
# Outputs:
#   data/processed/expr_gene_symbol.rds  (gene x sample matrix, log2-normalized)
#   data/processed/probe_annotation.tsv  (probe -> symbol mapping used)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(illuminaHumanv4.db)
})

ROOT <- "/Users/dpcook/My Drive/Lab/Manuscripts/hgsc_epitypes/supplemental_analysis/ICON7"
EXPR_FILE <- file.path(ROOT, "data/raw/GSE140082_geo.normdata.gz")
META_FILE <- file.path(ROOT, "data/processed/sample_metadata.tsv")

OUT_EXPR <- file.path(ROOT, "data/processed/expr_gene_symbol.rds")
OUT_PROBE_ANN <- file.path(ROOT, "data/processed/probe_annotation.tsv")
OUT_COV <- file.path(ROOT, "results/tables/01_signature_coverage.tsv")

dir.create(dirname(OUT_COV), showWarnings = FALSE, recursive = TRUE)

# ---- Read expression file -------------------------------------------------
# File is CSV. First column = probe ID. Subsequent columns interleave
# expression and detection p-values per sample. fread handles .gz natively.
message("Reading expression matrix (this may take ~1 minute)...")
expr_raw <- fread(EXPR_FILE)

n_col <- ncol(expr_raw)
message(sprintf("Raw expression matrix: %d rows x %d cols", nrow(expr_raw), n_col))
message(sprintf("First 6 column names: %s",
                paste(names(expr_raw)[1:6], collapse = " | ")))
message(sprintf("Last 4 column names: %s",
                paste(names(expr_raw)[(n_col-3):n_col], collapse = " | ")))

# Identify detection-pval columns. The convention used by Winterhoff's deposit
# is that detection p-value columns include the substring "Detection" or
# "_DP" or share the sample name with a ".1" suffix. We try multiple
# detection rules and report which one matched.

cn <- names(expr_raw)
probe_col <- cn[1]

# Heuristic 1: explicit "Detection" or "Pval"
det_idx_h1 <- grep("Detection|Pval|pval|\\.detection", cn, ignore.case = TRUE)
# Heuristic 2: columns ending in ".1" (data.table renames duplicates) usually
# from a sample-name pair (X1234_A, X1234_A.1 = detection pval)
det_idx_h2 <- grep("\\.1$", cn)
# Heuristic 3: every second column starting from column 3 (interleaved)
det_idx_h3 <- if (n_col %% 2 == 1) seq(3, n_col, by = 2) else integer(0)

# Choose the heuristic that gives n_samples = 380 expression columns
candidate_sets <- list(h1 = det_idx_h1, h2 = det_idx_h2, h3 = det_idx_h3)
for (nm in names(candidate_sets)) {
  s <- candidate_sets[[nm]]
  n_expr_cols <- n_col - 1 - length(s)
  message(sprintf("Heuristic %s: %d detection-pval columns -> %d expression columns",
                  nm, length(s), n_expr_cols))
}

# Use h2 if it gives n_expr_cols == 380, else h3, else h1
n_target <- 380
chosen <- NULL
for (nm in c("h2", "h3", "h1")) {
  s <- candidate_sets[[nm]]
  if ((n_col - 1 - length(s)) == n_target && length(s) > 0) {
    chosen <- nm
    break
  }
}
if (is.null(chosen)) {
  message("WARNING: no heuristic produced 380 expression columns. Inspect manually.")
  # Default to h2 (rename-by-data.table convention)
  chosen <- "h2"
}
message(sprintf("Using heuristic: %s", chosen))
det_idx <- candidate_sets[[chosen]]
expr_idx <- setdiff(seq_len(n_col), c(1, det_idx))

probe_ids <- expr_raw[[probe_col]]
expr_mat <- as.matrix(expr_raw[, ..expr_idx])
rownames(expr_mat) <- probe_ids
# Clean column names — strip any leading "X" that data.table adds for headers
# starting with a digit, and any trailing ".1"/".2"
colnames(expr_mat) <- gsub("\\.1$", "", colnames(expr_mat))
# Keep an unmodified version too in case sample_title has leading X
message(sprintf("Expression matrix: %d probes x %d samples", nrow(expr_mat), ncol(expr_mat)))
message(sprintf("Expression value range: [%.2f, %.2f]; median = %.2f",
                min(expr_mat, na.rm=TRUE),
                max(expr_mat, na.rm=TRUE),
                median(expr_mat, na.rm=TRUE)))

# ---- Verify column alignment to metadata ----------------------------------
meta <- fread(META_FILE)

# Try aligning by sample_title (BeadChip well IDs). The expression columns
# may have a leading "X" added if sample_title starts with a digit.
expr_cols <- colnames(expr_mat)
candidates <- c(
  list(raw = expr_cols),
  list(stripX = sub("^X", "", expr_cols))
)
matched <- NULL
for (nm in names(candidates)) {
  hits <- sum(candidates[[nm]] %in% meta$sample_title)
  message(sprintf("Match expr cols (%s) to meta$sample_title: %d / %d",
                  nm, hits, length(expr_cols)))
  if (hits > length(expr_cols) * 0.9 && is.null(matched)) {
    matched <- nm
  }
}
if (!is.null(matched)) {
  colnames(expr_mat) <- candidates[[matched]]
  message(sprintf("Renamed expression columns using: %s", matched))
}

# ---- Probe -> gene symbol mapping -----------------------------------------
message("Annotating probes with illuminaHumanv4.db ...")
ann <- AnnotationDbi::select(illuminaHumanv4.db,
                             keys = probe_ids,
                             columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                             keytype = "PROBEID")
# Some probes map to multiple entries; keep first non-NA per probe
ann <- as.data.table(ann)
ann <- ann[!is.na(SYMBOL)]
setkey(ann, PROBEID)
ann_first <- ann[, .(SYMBOL = SYMBOL[1], ENTREZID = ENTREZID[1]), by = PROBEID]
message(sprintf("Annotated %d / %d probes to a gene symbol",
                nrow(ann_first), length(probe_ids)))

fwrite(ann_first, OUT_PROBE_ANN, sep = "\t")

# ---- Collapse to gene-level by max-mean -----------------------------------
# Subset to annotated probes; compute row means; pick the probe with the
# largest mean expression per gene.
keep_probes <- intersect(rownames(expr_mat), ann_first$PROBEID)
expr_sub <- expr_mat[keep_probes, , drop = FALSE]
probe2sym <- setNames(ann_first$SYMBOL, ann_first$PROBEID)[keep_probes]

row_means <- rowMeans(expr_sub, na.rm = TRUE)
order_idx <- order(probe2sym, -row_means)
expr_ord <- expr_sub[order_idx, , drop = FALSE]
sym_ord <- probe2sym[order_idx]
keep <- !duplicated(sym_ord)
gene_mat <- expr_ord[keep, , drop = FALSE]
rownames(gene_mat) <- sym_ord[keep]
message(sprintf("Gene-level matrix: %d genes x %d samples",
                nrow(gene_mat), ncol(gene_mat)))

# ---- Save ----------------------------------------------------------------
saveRDS(gene_mat, OUT_EXPR)
message(sprintf("Wrote: %s", OUT_EXPR))
