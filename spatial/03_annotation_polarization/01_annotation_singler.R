# ============================================================================
# 01_annotation_singler.R — Cell type annotation via SingleR
# ============================================================================
# PURPOSE: Annotate every Xenium cell with SingleR against the atlas reference
#   (16 cell types, 1,000 cells per type). The reference is restricted to the
#   genes it shares with the Xenium panel, minus the 28 probe-QC exclusions
#   (429 genes), and SingleR runs once per SFE. Writes singler_label,
#   singler_pruned and singler_score into each SFE.
#
# SAMPLES: every SFE in sfe_names_all (TMA A, the arrays of cohort.tma_new and
#   WT01–WT08). An SFE that already carries singler_label is skipped.
#
# INPUTS:
#   - <sfe_dir>/sfe_*                                    (normalised SFEs)
#   - obj("xenium_ref") = xenium_celltype_downsampled.h5ad  (SingleR reference)
#   - <output_root>/05_probe_qc/genes_exclude_singler.txt  (28 excluded genes)
#
# OUTPUTS:
#   - SFEs updated with singler_label / singler_pruned / singler_score
#     (singler_score = best minus median reference score per cell)
#   - log messages "Shared genes (panel ∩ reference): N" and "After removing N
#     excluded genes: N genes for SingleR": tables/07_tableS7_xenium_gene_panel.py
#     parses them from this script's run log and checks them against Table S7,
#     so keep their wording
#
# MANUSCRIPT PANEL(S): the 16 SingleR types behind every Xenium cell type
#   (Figs 4–7, Figures S11–S18, Tables S7–S11), through 06_extract_singler.R.
#
# RUNTIME TIER: heavy (SingleR over millions of cells; SINGLER_WORKERS threads)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(SingleR)
library(zellkonverter)

message("=== Cell Type Annotation (SingleR) ===")

# ── 1. Setup & gene preparation ──────────────────────────────────────────────

exclude_genes <- readLines(file.path(out_dir, "05_probe_qc",
                                     "genes_exclude_singler.txt"))
exclude_genes <- trimws(exclude_genes[nchar(trimws(exclude_genes)) > 0])
message("Genes to exclude from SingleR: ", length(exclude_genes))

sfe_names <- sfe_names_all   # every SFE: TMA A, the cohort.tma_new arrays, WT01–WT08

# ── 2. Prepare reference ─────────────────────────────────────────────────────

message("\n--- Preparing scRNA-seq reference ---")

ref_path <- cfg_obj("xenium_ref")
message("Loading reference: ", ref_path)
ref_sce <- readH5AD(ref_path, reader = "R")

message(sprintf("  Reference: %s genes x %s cells",
                format(nrow(ref_sce), big.mark = ","),
                format(ncol(ref_sce), big.mark = ",")))

ref_labels <- ref_sce$xenium_celltype
message("  Cell types (", length(unique(ref_labels)), "):")
ref_tab <- sort(table(ref_labels), decreasing = TRUE)
for (ct in names(ref_tab)) {
  message(sprintf("    %-20s %s", ct, format(ref_tab[ct], big.mark = ",")))
}

first_sfe <- loadHDF5SummarizedExperiment(
  dir = file.path(sfe_dir, sfe_names[1])
)
panel_genes <- rownames(first_sfe)
rm(first_sfe); gc(verbose = FALSE)
message("  Xenium panel genes: ", length(panel_genes))

shared_genes <- intersect(rownames(ref_sce), panel_genes)
message("  Shared genes (panel ∩ reference): ", length(shared_genes))

shared_genes <- setdiff(shared_genes, exclude_genes)
n_excluded <- length(intersect(rownames(ref_sce), exclude_genes))
message("  After removing ", n_excluded, " excluded genes: ",
        length(shared_genes), " genes for SingleR")

ref_sce <- ref_sce[shared_genes, ]

if (!"logcounts" %in% assayNames(ref_sce)) {
  message("  Log-normalizing reference...")
  ref_sce <- logNormCounts(ref_sce)
} else {
  message("  Reference already has logcounts assay")
}

message(sprintf("  Final reference: %d genes x %s cells, %d cell types",
                nrow(ref_sce),
                format(ncol(ref_sce), big.mark = ","),
                length(unique(ref_labels))))

# ── 3. Run SingleR per-sample ────────────────────────────────────────────────

message("\n--- Running SingleR (per-sample) ---")

for (sfe_name in sfe_names) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Annotating: ", sfe_name)
  message(paste(rep("=", 60), collapse = ""))

  sfe_path <- file.path(sfe_dir, sfe_name)
  sfe <- loadHDF5SummarizedExperiment(dir = sfe_path)

  if ("singler_label" %in% colnames(colData(sfe))) {
    message("  Already annotated (singler_label found). Skipping.")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  n_cells <- ncol(sfe)
  message(sprintf("  Loaded: %s cells x %d genes",
                  format(n_cells, big.mark = ","), nrow(sfe)))

  genes_in_sfe <- intersect(shared_genes, rownames(sfe))
  if (length(genes_in_sfe) != length(shared_genes)) {
    warning("  Gene mismatch: ", length(genes_in_sfe), " / ",
            length(shared_genes), " shared genes found in SFE")
  }
  sfe_sub <- sfe[genes_in_sfe, ]

  message(sprintf("  Running SingleR on %s cells x %d genes...",
                  format(n_cells, big.mark = ","), nrow(sfe_sub)))
  t0 <- Sys.time()
  # Realize the (subsetted, 429-gene) logcounts into memory so SingleR's C++
  # threads (num.threads) work — the HDF5-backed DelayedArray breaks the parallel
  # path; dense is small here (~cells x 429). Threads via SINGLER_WORKERS (12).
  .nw <- as.integer(Sys.getenv("SINGLER_WORKERS", "12"))
  test_mat <- as.matrix(SummarizedExperiment::assay(sfe_sub, "logcounts"))
  singler_res <- SingleR(
    test           = test_mat,
    ref            = ref_sce,
    labels         = ref_labels,
    assay.type.ref  = "logcounts",
    num.threads    = .nw
  )
  rm(test_mat); gc(verbose = FALSE)
  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  message(sprintf("  SingleR completed in %s min", elapsed))

  sfe$singler_label  <- singler_res$labels
  sfe$singler_pruned <- singler_res$pruned.labels
  scores_mat <- singler_res$scores
  best_score <- apply(scores_mat, 1, max)
  median_score <- apply(scores_mat, 1, median)
  sfe$singler_score <- best_score - median_score

  label_tab <- sort(table(sfe$singler_label), decreasing = TRUE)
  n_pruned <- sum(is.na(sfe$singler_pruned))
  message(sprintf("  Assigned %d labels (%d pruned, %.1f%%)",
                  n_cells, n_pruned, 100 * n_pruned / n_cells))
  for (ct in names(label_tab)) {
    message(sprintf("    %-20s %6s (%5.1f%%)",
                    ct, format(label_tab[ct], big.mark = ","),
                    100 * label_tab[ct] / n_cells))
  }

  message("  Saving annotated SFE...")
  for (a in assayNames(sfe)) {
    assay(sfe, a) <- as(assay(sfe, a), "dgCMatrix")
  }

  save_sfe(sfe, sfe_name)

  rm(sfe, sfe_sub, singler_res, scores_mat, best_score, median_score)
  gc(verbose = FALSE)
}

# ── 4. Final summary ─────────────────────────────────────────────────────────

message("\n", paste(rep("=", 60), collapse = ""))
message("ANNOTATION SUMMARY")
message(paste(rep("=", 60), collapse = ""))

message(sprintf("Genes used for SingleR: %d (after excluding %d)",
                length(shared_genes), n_excluded))

message("\nDone. Review outputs:")
message("  Annotated SFEs: ", sfe_dir, "/sfe_*")
log_session()
