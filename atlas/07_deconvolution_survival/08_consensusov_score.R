# 08_consensusov_score.R — consensusOV TCGA HGSC subtype calls on pseudobulk
# ====================================================================
# HGSC malignant-states atlas backend.
#
# Purpose: call TCGA HGSC molecular subtypes (Differentiated / Immunoreactive /
# Mesenchymal / Proliferative) on the per-sample pseudobulk from
# 07_consensusov_export.py, using consensusOV (which aggregates Verhaak /
# Bentink / Helland / Konecny via an internal random forest). Runs both the
# bulk-tumor and epithelial-only pseudobulks.
#
# INPUTS (output_root/07_deconvolution_survival/consensusov/):
#   - pseudobulk_bulk_counts.tsv.gz, pseudobulk_epi_counts.tsv.gz
#   - pseudobulk_metadata.csv, gene_list.tsv
#
# OUTPUTS (output_root/07_deconvolution_survival/consensusov/):
#   - consensusov_calls_bulk.csv, consensusov_calls_epi.csv
#   - consensusov_classifiers_bulk.csv, consensusov_classifiers_epi.csv
#   - gene_id_mapping.csv
#
# MANUSCRIPT PANELS: upstream of Fig 3H (epitype x TCGA subtype).
#
# RUNTIME TIER: moderate (consensusOV RF on two pseudobulk matrices).
#
# SEEDING: set.seed(CFG$seed) before get.subtypes() so the consensusOV random
#   forest is deterministic (fixes the audit's non-determinism finding).
#
# Usage:
#   Rscript 08_consensusov_score.R

# config is 3 levels up: atlas/07_deconvolution_survival/08_*.R -> repo root
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."))
source(file.path(repo_root, "config", "config.R"))
set.seed(CFG$seed)

# Idempotent install. consensusOV was dropped from Bioconductor 3.20 (Apr 2026),
# so install from the canonical GitHub source repo (bhklab/consensusOV).
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")
for (pkg in c("org.Hs.eg.db", "edgeR", "AnnotationDbi")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
if (!requireNamespace("consensusOV", quietly = TRUE)) {
  remotes::install_github("bhklab/consensusOV", upgrade = "never")
}

suppressWarnings(suppressMessages({
  library(consensusOV)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(edgeR)
  library(data.table)
}))

# ── Paths ──────────────────────────────────────────────────────
in_dir  <- cfg_path("output_root", "07_deconvolution_survival", "consensusov")
out_dir <- in_dir
bulk_path <- file.path(in_dir, "pseudobulk_bulk_counts.tsv.gz")
epi_path  <- file.path(in_dir, "pseudobulk_epi_counts.tsv.gz")

cat(strrep("=", 70), "\n", sep = "")
cat("  consensusOV TCGA subtype calls on pseudobulk\n")
cat(strrep("=", 70), "\n", sep = "")

# ── Helpers ────────────────────────────────────────────────────
read_counts <- function(path) {
  cat(sprintf("  Reading %s ...\n", path))
  dt <- fread(path, sep = "\t", header = TRUE)
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- dt[[1]]
  m
}

normalize_logcpm <- function(counts) {
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = TRUE, prior.count = 1)
}

map_symbols_to_entrez <- function(symbols) {
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
                               columns = c("SYMBOL", "ENTREZID"))
  map <- map[!is.na(map$ENTREZID), , drop = FALSE]
  map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
  rownames(map) <- map$SYMBOL
  map
}

run_consensusov <- function(counts, label) {
  cat(sprintf("\n[%s] Pseudobulk: %d genes x %d samples\n", label, nrow(counts), ncol(counts)))
  keep_samples <- colSums(counts) > 0
  if (any(!keep_samples)) {
    cat(sprintf("  Dropping %d empty sample(s)\n", sum(!keep_samples)))
    counts <- counts[, keep_samples, drop = FALSE]
  }
  keep_genes <- rowSums(counts) >= 10
  cat(sprintf("  Keeping %d / %d genes (rowsum >= 10)\n", sum(keep_genes), length(keep_genes)))
  counts <- counts[keep_genes, , drop = FALSE]

  expr <- normalize_logcpm(counts)
  cat(sprintf("  log2(CPM+1) matrix: %d x %d\n", nrow(expr), ncol(expr)))

  sym_map <- map_symbols_to_entrez(rownames(expr))
  shared <- intersect(rownames(expr), sym_map$SYMBOL)
  expr <- expr[shared, , drop = FALSE]
  entrez <- sym_map[shared, "ENTREZID"]
  if (any(duplicated(entrez))) {
    ord <- order(rowMeans(expr), decreasing = TRUE)
    expr <- expr[ord, , drop = FALSE]; entrez <- entrez[ord]
    keep <- !duplicated(entrez)
    expr <- expr[keep, , drop = FALSE]; entrez <- entrez[keep]
  }
  rownames(expr) <- entrez
  cat(sprintf("  Mapped to %d Entrez genes\n", nrow(expr)))

  cat("  Running consensusOV (consensus method)...\n")
  res_consensus <- get.subtypes(expr, entrez.ids = entrez, method = "consensusOV")
  margins <- as.data.frame(res_consensus$rf.probs)
  margins$sample_id   <- colnames(expr)
  margins$consensusOV <- as.character(res_consensus$consensusOV.subtypes)
  prob_cols <- setdiff(colnames(margins), c("sample_id", "consensusOV"))
  margins$margin_top1_top2 <- apply(margins[, prob_cols, drop = FALSE], 1, function(p) {
    s <- sort(as.numeric(p), decreasing = TRUE)
    if (length(s) >= 2) s[1] - s[2] else NA_real_
  })

  classifiers <- c("Verhaak", "Bentink", "Helland", "Konecny")
  per_clf <- list()
  for (m in classifiers) {
    cat(sprintf("  Running %s ...\n", m))
    out <- tryCatch(get.subtypes(expr, entrez.ids = entrez, method = m),
                    error = function(e) {
                      cat(sprintf("    [warn] %s failed: %s\n", m, conditionMessage(e))); NULL })
    if (!is.null(out)) {
      lbl_field <- paste0(tolower(m), ".subtypes")
      per_clf[[m]] <- as.character(if (lbl_field %in% names(out)) out[[lbl_field]] else out[[1]])
    } else {
      per_clf[[m]] <- rep(NA_character_, ncol(expr))
    }
  }
  per_clf_df <- as.data.frame(per_clf, stringsAsFactors = FALSE)
  per_clf_df$sample_id   <- colnames(expr)
  per_clf_df$consensusOV <- margins$consensusOV
  per_clf_df <- per_clf_df[, c("sample_id", "consensusOV", classifiers)]

  list(margins = margins, per_clf = per_clf_df,
       gene_map = data.frame(symbol = shared, entrez = entrez, stringsAsFactors = FALSE))
}

# ── Run ────────────────────────────────────────────────────────
bulk_res <- run_consensusov(read_counts(bulk_path), "BULK")
epi_res  <- run_consensusov(read_counts(epi_path),  "EPI")

# ── Write ──────────────────────────────────────────────────────
cat("\n[writing outputs]\n")
fwrite(bulk_res$margins, file.path(out_dir, "consensusov_calls_bulk.csv"))
fwrite(epi_res$margins,  file.path(out_dir, "consensusov_calls_epi.csv"))
fwrite(bulk_res$per_clf, file.path(out_dir, "consensusov_classifiers_bulk.csv"))
fwrite(epi_res$per_clf,  file.path(out_dir, "consensusov_classifiers_epi.csv"))
fwrite(bulk_res$gene_map, file.path(out_dir, "gene_id_mapping.csv"))
cat("  consensusov_calls_{bulk,epi}.csv, consensusov_classifiers_{bulk,epi}.csv, gene_id_mapping.csv\n")
cat("\n[done]\n")
