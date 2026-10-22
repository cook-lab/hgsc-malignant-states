#!/usr/bin/env Rscript
# Step 19_cnv — Run CopyKAT on a single sample
# =============================================
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
#
# Reads per-sample count matrix + reference barcodes extracted by
# 19_cnv_extract.py, runs CopyKAT, clusters aneuploid cells into
# subclones (1-Pearson distance, Ward.D2, silhouette-picked k in 2..5,
# monoclonal fallback when max silhouette < SILHOUETTE_MIN).
#
# Progress-safe: writes DONE.txt on success and FAILED.txt on failure.
# Appends one line to $BATCH_LOG (set by caller) for live monitoring.
#
# Usage:
#   Rscript scripts/19_cnv_copykat.R <sample_dir> [n_cores]
#
# Input:
#   <sample_dir>/counts.mtx.gz, genes.txt, barcodes.csv, ref_barcodes.txt
#
# Output (all in <sample_dir>/):
#   copykat_prediction.csv        per-cell aneuploid/diploid
#   copykat_CNA_results.txt       bin x cell CNA matrix
#   copykat_subclones.csv         barcode, subclone (clone_k | diploid)
#   copykat_subclone_qc.csv       k, silhouette, n_aneuploid, monoclonal
#   {sample_id}_copykat_*         CopyKAT's native files (heatmap jpeg, etc.)
#   DONE.txt  or  FAILED.txt

SILHOUETTE_MIN <- 0.15
K_RANGE        <- 2:5
WIN_SIZE       <- 25

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 19_cnv_copykat.R <sample_dir> [n_cores]")
}
sample_dir <- sub("/$", "", args[1])
n_cores    <- ifelse(length(args) >= 2, as.integer(args[2]), 4L)
sample_id  <- basename(sample_dir)

cat(strrep("=", 60), "\n", sep = "")
cat("  CopyKAT - ", sample_id, "  |  cores: ", n_cores, "\n", sep = "")
cat(strrep("=", 60), "\n", sep = "")

if (file.exists(file.path(sample_dir, "DONE.txt"))) {
  cat("  Already done, skipping.\n"); quit(status = 0)
}
if (file.exists(file.path(sample_dir, "FAILED.txt"))) {
  file.remove(file.path(sample_dir, "FAILED.txt"))
}

batch_log <- Sys.getenv("BATCH_LOG", "")
t0 <- proc.time()

append_log <- function(status, extras = "") {
  if (nchar(batch_log) == 0) return(invisible())
  line <- sprintf("%s\t%s\t%s\t%.1fmin\t%s\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  sample_id, status,
                  (proc.time() - t0)[["elapsed"]] / 60, extras)
  cat(line, file = batch_log, append = TRUE)
}

tryCatch({

  suppressPackageStartupMessages({
    library(Matrix); library(copykat); library(cluster)
  })

  pred_path <- file.path(sample_dir, "copykat_prediction.csv")
  cna_path  <- file.path(sample_dir, "copykat_CNA_results.txt")

  if (file.exists(pred_path) && file.exists(cna_path)) {
    cat("  Re-using existing CopyKAT outputs (skipping copykat() call)\n")
    pred <- read.csv(pred_path, stringsAsFactors = FALSE)
    cna_mat <- as.matrix(read.table(cna_path, sep = "\t", header = TRUE,
                                    check.names = TRUE, row.names = NULL))
  } else {
    cat("  Reading inputs...\n")
    counts       <- readMM(file.path(sample_dir, "counts.mtx.gz"))
    genes        <- readLines(file.path(sample_dir, "genes.txt"))
    barcodes     <- read.csv(file.path(sample_dir, "barcodes.csv"),
                             stringsAsFactors = FALSE)
    ref_barcodes <- readLines(file.path(sample_dir, "ref_barcodes.txt"))

    cat("    dims: ", dim(counts)[1], " genes x ",
        dim(counts)[2], " cells   ref: ", length(ref_barcodes),
        "   epi: ", sum(barcodes$is_reference == "False"), "\n", sep = "")

    counts_dense <- as.matrix(counts)
    rownames(counts_dense) <- genes
    colnames(counts_dense) <- barcodes$barcode
    rm(counts); gc()

    cat("\n  Running CopyKAT (", n_cores, " cores)...\n", sep = "")
    old_wd <- getwd()
    setwd(sample_dir)

    copykat_result <- copykat(
      rawmat          = counts_dense,
      id.type         = "S",
      cell.line       = "no",
      ngene.chr       = 5,
      LOW.DR          = 0.05,
      UP.DR           = 0.1,
      win.size        = WIN_SIZE,
      norm.cell.names = ref_barcodes,
      KS.cut          = 0.1,
      sam.name        = sample_id,
      distance        = "euclidean",
      n.cores         = n_cores,
      output.seg      = FALSE
    )

    setwd(old_wd)
    rm(counts_dense); gc()

    pred <- copykat_result$prediction
    write.csv(pred, pred_path, row.names = FALSE)

    cna_mat <- copykat_result$CNAmat
    if (!is.null(cna_mat)) {
      write.table(cna_mat, cna_path, sep = "\t", quote = FALSE)
    }
  }

  n_aneu <- sum(pred$copykat.pred == "aneuploid",   na.rm = TRUE)
  n_dip  <- sum(pred$copykat.pred == "diploid",     na.rm = TRUE)
  n_nd   <- sum(pred$copykat.pred == "not.defined", na.rm = TRUE)
  cat("  Predictions -- aneu: ", n_aneu, "  dip: ", n_dip,
      "  nd: ", n_nd, "\n", sep = "")

  qc_k <- NA_integer_; qc_sil <- NA_real_; monoclonal <- FALSE

  if (!is.null(cna_mat) && ncol(cna_mat) > 1 && n_aneu >= 20) {
    # CopyKAT's CNAmat colnames are make.names()-mangled cell barcodes
    # (e.g. "49564-1" -> "X49564.1"). Also the first 3 columns are
    # genomic coords (chrom, chrompos, abspos) -- not cells. Map raw
    # barcodes through make.names() so the intersect hits, and keep
    # raw names for the output CSV (so Python can join by barcode).
    aneu_raw <- pred$cell.names[pred$copykat.pred == "aneuploid"]
    aneu_mangled <- make.names(aneu_raw)
    cell_cols <- setdiff(colnames(cna_mat),
                         c("chrom", "chrompos", "abspos"))
    keep <- aneu_mangled %in% cell_cols
    aneu_raw     <- aneu_raw[keep]
    aneu_mangled <- aneu_mangled[keep]
    cat("  Matched ", length(aneu_mangled), " / ", length(keep),
        " aneuploid cells to CNA matrix\n", sep = "")
    if (length(aneu_mangled) < 20) {
      stop(sprintf("Too few aneuploid cells after barcode match: %d",
                   length(aneu_mangled)))
    }
    cna_cells <- t(cna_mat[, aneu_mangled, drop = FALSE])
    # cna_mat bins may be character if read from disk; coerce to numeric
    storage.mode(cna_cells) <- "numeric"
    rownames(cna_cells) <- aneu_raw  # use raw names downstream

    cor_mat <- suppressWarnings(
      cor(t(cna_cells), use = "pairwise.complete.obs"))
    cor_mat[is.na(cor_mat)] <- 0
    d <- as.dist(1 - cor_mat)
    hc <- hclust(d, method = "ward.D2")

    best_k <- NA_integer_; best_sil <- -Inf
    for (k in K_RANGE) {
      if (k >= length(aneu_raw)) next
      cl <- cutree(hc, k = k)
      if (length(unique(cl)) < 2) next
      s <- tryCatch(mean(silhouette(cl, d)[, 3]),
                    error = function(e) -Inf)
      if (s > best_sil) { best_sil <- s; best_k <- k }
    }

    if (is.na(best_k) || best_sil < SILHOUETTE_MIN) {
      monoclonal <- TRUE
      subclone_labels <- setNames(rep(1L, length(aneu_raw)), aneu_raw)
      qc_k <- 1L
      qc_sil <- ifelse(is.finite(best_sil), best_sil, NA_real_)
      cat("  Monoclonal (silhouette < ", SILHOUETTE_MIN,
          "): k=1\n", sep = "")
    } else {
      subclone_labels <- cutree(hc, k = best_k)
      qc_k <- best_k; qc_sil <- best_sil
      cat("  Subclones: k=", best_k, "  silhouette=",
          round(best_sil, 3), "\n", sep = "")
    }

    # subclone_labels is named by raw barcode (from hclust labels
    # propagated through cor/dist). Align names explicitly to aneu_raw.
    subclone_df <- data.frame(
      barcode  = aneu_raw,
      subclone = paste0("clone_", subclone_labels[aneu_raw]),
      stringsAsFactors = FALSE)
    dip_cells <- pred$cell.names[pred$copykat.pred == "diploid"]
    if (length(dip_cells) > 0) {
      subclone_df <- rbind(subclone_df, data.frame(
        barcode = dip_cells, subclone = "diploid",
        stringsAsFactors = FALSE))
    }
    write.csv(subclone_df, file.path(sample_dir, "copykat_subclones.csv"),
              row.names = FALSE)

    for (sc in unique(subclone_df$subclone)) {
      cat("    ", sc, ": ", sum(subclone_df$subclone == sc),
          " cells\n", sep = "")
    }
  } else {
    cat("  Skipping subclone clustering (n_aneu<20 or no CNA mat)\n")
    subclone_df <- data.frame(
      barcode  = pred$cell.names,
      subclone = ifelse(pred$copykat.pred == "aneuploid",
                        "clone_1", "diploid"),
      stringsAsFactors = FALSE)
    write.csv(subclone_df, file.path(sample_dir, "copykat_subclones.csv"),
              row.names = FALSE)
    qc_k <- ifelse(n_aneu > 0, 1L, 0L); monoclonal <- TRUE
  }

  qc_df <- data.frame(
    sample_id   = sample_id,
    n_aneuploid = n_aneu, n_diploid = n_dip, n_nd = n_nd,
    k = qc_k, silhouette = qc_sil, monoclonal = monoclonal,
    elapsed_min = round((proc.time() - t0)[["elapsed"]] / 60, 2),
    stringsAsFactors = FALSE)
  write.csv(qc_df, file.path(sample_dir, "copykat_subclone_qc.csv"),
            row.names = FALSE)

  elapsed <- (proc.time() - t0)[["elapsed"]]
  writeLines(paste0("Completed in ", round(elapsed / 60, 1), " min"),
             file.path(sample_dir, "DONE.txt"))
  cat("\n  DONE in ", round(elapsed / 60, 1), " min\n", sep = "")

  append_log("DONE",
             sprintf("aneu=%d dip=%d k=%s sil=%s",
                     n_aneu, n_dip,
                     ifelse(is.na(qc_k), "NA", as.character(qc_k)),
                     ifelse(is.na(qc_sil), "NA", sprintf("%.3f", qc_sil))))

}, error = function(e) {
  msg <- conditionMessage(e)
  cat("\n  FAILED: ", msg, "\n")
  writeLines(paste0("Failed: ", msg),
             file.path(sample_dir, "FAILED.txt"))
  append_log("FAILED", substr(msg, 1, 200))
  quit(status = 1)
})
