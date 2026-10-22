#!/usr/bin/env Rscript
# ============================================================================
# 05_clean_split_rctd.R — RCTD/SPLIT tables and the last write of the SFE cell_label
# ============================================================================
# PURPOSE: Run RCTD (spacexr, doublet mode) with SPLIT post-processing and
#   SPLIT::purify() on each SFE against the 16-type SingleR reference, and
#   write per-SFE tables to 06g_clean_split/. The --all / --writeback pass then
#   writes the result into each SFE: cell_label becomes the RCTD first type for
#   non-secretory cells and, for RCTD "Secretory epithelium" cells, SecA /
#   Intermediate / SecB from polarization_UCell and the thresholds in
#   threshold_summary.csv (cells without an RCTD call, and RCTD secretory cells
#   without a score, keep their previous label); the RCTD columns and a
#   "purified" assay are added.
#
#   No manuscript panel uses these labels or the purified assay: cell types
#   come from the SingleR label tables (06_extract_singler.R). The step is kept
#   because the Xenium raw-count export (spatial/05_count_export/as_run/) reads
#   <sfe>_rctd_meta.csv and <sfe>_purified_counts.rds, and takes its cell set
#   from the cells with a non-empty SFE cell_label after this writeback.
#
# SAMPLES: --all runs every SFE in sfe_names_all; a comma-separated list runs
#   those SFEs; no argument runs the first whole tissue only (pilot).
#   RERUN_WRITEBACK_ONLY="sfe_a,sfe_b" limits the writeback to the listed SFEs.
#
# DETERMINISM: set.seed(CFG$seed) before each SFE's RCTD run.
#
# USAGE:
#   Rscript spatial/03_annotation_polarization/05_clean_split_rctd.R            # pilot
#   Rscript .../05_clean_split_rctd.R --all                         # all SFEs, then writeback
#   Rscript .../05_clean_split_rctd.R sfe_<sample1>,sfe_<sample2> --suffix chunkA
#   Rscript .../05_clean_split_rctd.R --writeback                              # writeback only
#
# Reference: Marconato et al. 2026 Nat Methods; SPLIT @ bdsc-tds/SPLIT.
#
# INPUTS:
#   - <sfe_dir>/sfe_*  (SFEs with cell_label and polarization_UCell, steps 01–04)
#   - obj("xenium_ref") = xenium_celltype_downsampled.h5ad  (scRNA-seq reference)
#   - <output_root>/05_probe_qc/genes_exclude_singler.txt   (genes left out of RCTD)
#   - <output_root>/06f_reclassification_polarization/threshold_summary.csv
#
# OUTPUTS:
#   - <output_root>/06g_clean_split/<sfe>_rctd_meta.csv  (per cell: RCTD first and
#       second type, spot class, first-type weight; empty = RCTD reject)
#   - <output_root>/06g_clean_split/<sfe>_purified_counts.rds  (SPLIT-purified counts)
#   - <output_root>/06g_clean_split/06g_summary<suffix>.csv  (per-SFE RCTD counts of
#       the run: genes used, cells assigned / rejected, % secretory)
#   - SFEs (writeback): cell_label overwritten, the previous label kept as
#       cell_label_06f; rctd_first_type, rctd_second_type, rctd_spot_class,
#       rctd_first_weight; "purified" assay
#   - <output_root>/06g_clean_split/.sfe_update_complete  (checkpoint: when
#       present, the writeback is skipped)
#
# MANUSCRIPT PANEL(S): none directly. Its outputs define the cell set of the
#   Xenium raw-count export read by Figure 6f–i and Tables S9–S10.
#
# RUNTIME TIER: heavy (~21 GB peak per whole tissue; ~30-40 GB for TMA)
# ============================================================================

source("spatial/00_setup/00_setup.R")

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Matrix)
  library(data.table)
})

for (pkg in c("spacexr", "SPLIT")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Install from GitHub: bdsc-tds/SPLIT", pkg))
  }
}
library(spacexr)
library(SPLIT)

message("\n=== SPLIT/RCTD Annotation Cleaning ===")
message("[", Sys.time(), "] Starting...")

# ── CLI argument parsing ─────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

chunk_suffix <- ""
if ("--suffix" %in% args) {
  i <- which(args == "--suffix")
  if (i + 1 <= length(args)) {
    chunk_suffix <- paste0("_", args[i + 1])
    args <- args[-c(i, i + 1)]
  }
}

do_writeback_only <- "--writeback" %in% args
if (do_writeback_only) args <- args[args != "--writeback"]

do_run_all <- "--all" %in% args

# All SFEs for --all and for the writeback existence check.
ALL_SFES <- sfe_names_all

# Optional writeback scoping: RERUN_WRITEBACK_ONLY="sfe_a,sfe_b" restricts PART 2
# (the existence check and the SFE writeback) to exactly those SFEs, leaving
# every other SFE untouched. It also replaces the --all target list.
.wb_scope <- {
  .v <- Sys.getenv("RERUN_WRITEBACK_ONLY", "")
  if (nzchar(.v)) trimws(strsplit(.v, ",")[[1]]) else NULL
}
if (!is.null(.wb_scope)) ALL_SFES <- .wb_scope

PILOT_SAMPLES <- sfe_names_wt[1]

target_sfes <- if (do_writeback_only) {
  character(0)
} else if (length(args) == 0) {
  PILOT_SAMPLES
} else if (do_run_all) {
  ALL_SFES
} else {
  strsplit(args[1], ",")[[1]]
}

# SFEs to write back: sfe_names_all, or the RERUN_WRITEBACK_ONLY list.
ALL_SFES_WRITEBACK <- sfe_names_all
if (!is.null(.wb_scope)) ALL_SFES_WRITEBACK <- .wb_scope

# ── Paths & constants ────────────────────────────────────────────────────────

step_dir  <- file.path(out_dir, "06g_clean_split")
fig_out   <- file.path(step_dir, "figures")
flag_file <- file.path(step_dir, ".sfe_update_complete")

for (d in c(step_dir, fig_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sec_classes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

valid_cell_labels_extra <- c(
  "Ciliated epithelium", "Mesothelial", "Fibroblast", "Smooth muscle",
  "Pericyte", "Endothelial", "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

# The 16-type reference label for secretory cells (unsplit)
SECRETORY_REF_LABEL <- "Secretory epithelium"

# ── Load the frozen thresholds (04_reclassification_polarization.R) ──────────

thresh_path <- file.path(out_dir, "06f_reclassification_polarization", "threshold_summary.csv")
if (!file.exists(thresh_path)) {
  stop("Frozen thresholds not found: ", thresh_path, "\n  Run 04_reclassification_polarization.R first.")
}
thresh_dt <- fread(thresh_path)
t_low  <- as.numeric(thresh_dt[metric == "t_low_SecA_p75",  value])
t_high <- as.numeric(thresh_dt[metric == "t_high_SecB_p25", value])
message(sprintf("  Frozen thresholds: t_low = %.6f, t_high = %.6f", t_low, t_high))

classify_polarization <- function(pol, t_low, t_high) {
  fifelse(
    pol < t_low, "SecA epithelium",
    fifelse(pol < t_high, "Intermediate epithelium", "SecB epithelium")
  )
}

# ── Helpers ──────────────────────────────────────────────────────────────────

# spacexr disallows '/' in cell type names
sanitize_celltype <- function(x) gsub("/", "_SLASH_", as.character(x))
unsanitize_celltype <- function(x) gsub("_SLASH_", "/", as.character(x))


# ============================================================================
# PART 0: Build spacexr Reference
# ============================================================================

message("\n[", Sys.time(), "] === PART 0: Building spacexr Reference ===")

ref_path <- cfg_obj("xenium_ref")
message("  Loading reference: ", ref_path)
ref_sce <- readH5AD(ref_path, reader = "R")
message(sprintf("  Reference: %s genes x %s cells",
                format(nrow(ref_sce), big.mark = ","),
                format(ncol(ref_sce), big.mark = ",")))

ref_counts_raw <- as(assay(ref_sce, "X"), "CsparseMatrix")
if (any(ref_counts_raw@x %% 1 != 0)) {
  ref_counts_raw@x <- round(ref_counts_raw@x)
}
storage.mode(ref_counts_raw@x) <- "double"

ref_labels <- sanitize_celltype(ref_sce$xenium_celltype)
names(ref_labels) <- colnames(ref_counts_raw)
ref_labels <- factor(ref_labels)

# Drop rare types (< 25 cells)
ct_tab <- table(ref_labels)
keep_ct <- names(ct_tab)[ct_tab >= 25]
keep_cells <- names(ref_labels)[ref_labels %in% keep_ct]
ref_counts_raw <- ref_counts_raw[, keep_cells]
ref_labels <- droplevels(ref_labels[keep_cells])

message(sprintf("  Reference after filtering: %s cells, %d types",
                format(ncol(ref_counts_raw), big.mark = ","), nlevels(ref_labels)))
message("  Types: ", paste(levels(ref_labels), collapse = ", "))

spacexr_ref_full <- Reference(ref_counts_raw, cell_types = ref_labels,
                               min_UMI = 10, require_int = TRUE)

exclude_genes <- readLines(file.path(out_dir, "05_probe_qc", "genes_exclude_singler.txt"))
exclude_genes <- trimws(exclude_genes[nchar(trimws(exclude_genes)) > 0])
message("  Probe-QC excluded genes: ", length(exclude_genes))

rm(ref_sce); gc(verbose = FALSE)


# ============================================================================
# PART 1: RCTD + SPLIT Per Sample
# ============================================================================

if (length(target_sfes) > 0) {
  message("\n[", Sys.time(), "] === PART 1: RCTD + SPLIT ===")
  message("  Targets: ", paste(target_sfes, collapse = ", "))

  results_summary <- list()

  for (i in seq_along(target_sfes)) {
    nm <- target_sfes[i]
    message(sprintf("\n[%d/%d] %s", i, length(target_sfes), nm))
    t_total <- Sys.time()

    # Seed before the stochastic RCTD/SPLIT optimization for reproducible labels.
    set.seed(CFG$seed)

    meta_csv <- file.path(step_dir, paste0(nm, "_rctd_meta.csv"))
    if (file.exists(meta_csv)) {
      message("  Already processed, skipping (delete _rctd_meta.csv to re-run)")
      next
    }

    sfe <- load_sfe(nm)
    n_cells <- ncol(sfe)
    message(sprintf("  SFE: %s cells x %d genes", format(n_cells, big.mark = ","), nrow(sfe)))

    common <- intersect(rownames(sfe), rownames(slot(spacexr_ref_full, "counts")))
    use_genes <- setdiff(common, exclude_genes)
    message(sprintf("  Genes: shared=%d, usable (excl. probe-FAIL)=%d",
                    length(common), length(use_genes)))

    message("  Building SpatialRNA test object...")
    xen_counts <- as(assay(sfe, "counts"), "CsparseMatrix")[use_genes, , drop = FALSE]
    storage.mode(xen_counts@x) <- "double"

    coords_mat <- spatialCoords(sfe)
    coords_df <- data.frame(x = coords_mat[, 1], y = coords_mat[, 2],
                            row.names = colnames(sfe))

    test_rna <- SpatialRNA(coords = coords_df, counts = xen_counts,
                           require_int = TRUE)

    ref_panel_counts <- slot(spacexr_ref_full, "counts")[use_genes, , drop = FALSE]
    ref_panel <- Reference(ref_panel_counts,
                           cell_types = slot(spacexr_ref_full, "cell_types"),
                           min_UMI = 10, require_int = TRUE)

    rm(xen_counts); gc(verbose = FALSE)

    message("  Creating RCTD object...")
    n_cores <- max(1, parallel::detectCores() - 2)
    RCTD <- create.RCTD(
      test_rna, ref_panel,
      UMI_min = 10,
      counts_MIN = 10,
      UMI_min_sigma = 100,
      max_cores = n_cores,
      CELL_MIN_INSTANCE = 25
    )

    message("  Running RCTD doublet mode (this is the slow step)...")
    t0 <- Sys.time()
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    message(sprintf("  RCTD done in %.1f min",
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))

    message("  SPLIT post-processing RCTD...")
    RCTD <- SPLIT::run_post_process_RCTD(RCTD)

    message("  Running SPLIT::purify()...")
    t0 <- Sys.time()
    res_split <- SPLIT::purify(
      counts = assay(sfe, "counts")[use_genes, , drop = FALSE],
      rctd   = RCTD,
      DO_purify_singlets = TRUE
    )
    message(sprintf("  SPLIT purify done in %.1f min",
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))

    purified_rds <- file.path(step_dir, paste0(nm, "_purified_counts.rds"))
    saveRDS(res_split$purified_counts, purified_rds)
    message("  Saved: ", basename(purified_rds))

    all_barcodes <- colnames(sfe)
    cell_meta <- res_split$cell_meta
    if (is.null(cell_meta)) {
      cell_meta <- data.frame(row.names = colnames(res_split$purified_counts))
    }

    rctd_meta <- data.table(
      barcode = all_barcodes,
      rctd_first_type  = NA_character_,
      rctd_second_type = NA_character_,
      rctd_spot_class  = NA_character_,
      rctd_first_weight = NA_real_
    )

    meta_barcodes <- rownames(cell_meta)
    matched <- intersect(all_barcodes, meta_barcodes)
    idx <- match(matched, rctd_meta$barcode)

    if ("first_type" %in% colnames(cell_meta)) {
      rctd_meta[idx, rctd_first_type := unsanitize_celltype(cell_meta[matched, "first_type"])]
    }
    if ("second_type" %in% colnames(cell_meta)) {
      rctd_meta[idx, rctd_second_type := unsanitize_celltype(cell_meta[matched, "second_type"])]
    }
    if ("spot_class" %in% colnames(cell_meta)) {
      rctd_meta[idx, rctd_spot_class := as.character(cell_meta[matched, "spot_class"])]
    }
    if ("weight_first_type" %in% colnames(cell_meta)) {
      rctd_meta[idx, rctd_first_weight := as.numeric(cell_meta[matched, "weight_first_type"])]
    }

    fwrite(rctd_meta, meta_csv)
    message("  Saved: ", basename(meta_csv))

    elapsed <- round(as.numeric(difftime(Sys.time(), t_total, units = "mins")), 1)
    results_summary[[nm]] <- data.table(
      sfe_name     = nm,
      n_cells      = n_cells,
      n_genes_used = length(use_genes),
      n_rctd_assigned = sum(!is.na(rctd_meta$rctd_first_type)),
      n_rctd_reject   = sum(is.na(rctd_meta$rctd_first_type)),
      pct_secretory   = round(100 * sum(rctd_meta$rctd_first_type == SECRETORY_REF_LABEL,
                                         na.rm = TRUE) / n_cells, 2),
      elapsed_min  = elapsed
    )

    rm(sfe, RCTD, test_rna, ref_panel, res_split); gc(verbose = FALSE)
    message(sprintf("  Total elapsed: %.1f min", elapsed))
  }

  if (length(results_summary) > 0) {
    summary_dt <- rbindlist(results_summary)
    summary_csv <- file.path(step_dir, paste0("06g_summary", chunk_suffix, ".csv"))
    fwrite(summary_dt, summary_csv)
    message("\n  Summary saved: ", basename(summary_csv))
    print(summary_dt)
  }
}


# ============================================================================
# PART 2: Label Assignment & SFE Writeback
# ============================================================================

if (do_writeback_only || do_run_all) {

  message("\n[", Sys.time(), "] === PART 2: SFE Writeback ===")

  if (file.exists(flag_file)) {
    message("  Checkpoint flag found — SFEs already updated. Skipping writeback.")
  } else {

    # All RCTD meta CSVs must exist.
    missing <- character(0)
    for (nm in ALL_SFES) {
      csv <- file.path(step_dir, paste0(nm, "_rctd_meta.csv"))
      if (!file.exists(csv)) missing <- c(missing, nm)
    }
    if (length(missing) > 0) {
      stop("RCTD meta CSVs missing for: ", paste(missing, collapse = ", "),
           "\n  Run RCTD on these samples first before --writeback.")
    }

    for (nm in ALL_SFES_WRITEBACK) {
      t0 <- Sys.time()
      message("\n[", Sys.time(), "] Updating ", nm, "...")

      sfe <- load_sfe(nm)
      cd  <- colData(sfe)

      source_nm <- if (nm == "sfe_tma_filtered") "sfe_tma" else nm
      meta_csv <- file.path(step_dir, paste0(source_nm, "_rctd_meta.csv"))
      rctd_meta <- fread(meta_csv)

      sfe_barcodes <- colnames(sfe)
      rctd_meta <- rctd_meta[match(sfe_barcodes, barcode)]

      n_matched <- sum(!is.na(rctd_meta$barcode))
      message(sprintf("  Matched %s / %s barcodes from RCTD meta",
                      format(n_matched, big.mark = ","),
                      format(length(sfe_barcodes), big.mark = ",")))

      if (!"polarization_UCell" %in% colnames(cd)) {
        rm(sfe); gc(verbose = FALSE)
        stop(sprintf("SFE %s lacks polarization_UCell — rerun 03_ucell/04_reclassification first.", nm))
      }

      # Keep the current cell_label (from 04_reclassification_polarization.R) as cell_label_06f
      cd$cell_label_06f <- as.character(cd$cell_label)
      # Read a legacy "Transitioning epithelium" label as "Intermediate epithelium"
      # before the vector is used as the fallback label (no-op if absent).
      cd$cell_label_06f[cd$cell_label_06f == "Transitioning epithelium"] <- "Intermediate epithelium"

      cd$rctd_first_type   <- rctd_meta$rctd_first_type
      cd$rctd_second_type  <- rctd_meta$rctd_second_type
      cd$rctd_spot_class   <- rctd_meta$rctd_spot_class
      cd$rctd_first_weight <- rctd_meta$rctd_first_weight

      rctd_type <- as.character(rctd_meta$rctd_first_type)
      pol       <- as.numeric(cd$polarization_UCell)
      old_label <- as.character(cd$cell_label_06f)

      new_label <- character(length(sfe_barcodes))

      has_rctd <- !is.na(rctd_type)

      # Secretory: apply polarization thresholds
      is_sec_rctd <- has_rctd & rctd_type == SECRETORY_REF_LABEL
      has_pol <- !is.na(pol)

      new_label[is_sec_rctd & has_pol] <- classify_polarization(
        pol[is_sec_rctd & has_pol], t_low, t_high
      )
      # Secretory but no polarization score: fall back to old label
      new_label[is_sec_rctd & !has_pol] <- old_label[is_sec_rctd & !has_pol]

      # Non-secretory: use RCTD first_type directly
      is_nonsec_rctd <- has_rctd & rctd_type != SECRETORY_REF_LABEL
      new_label[is_nonsec_rctd] <- rctd_type[is_nonsec_rctd]

      # RCTD reject/NA: fall back to old label
      new_label[!has_rctd] <- old_label[!has_rctd]

      cd$cell_label <- new_label
      colData(sfe) <- cd

      allowed <- c(sec_classes, valid_cell_labels_extra)
      bad <- setdiff(unique(new_label), allowed)
      if (length(bad) > 0) {
        warning(sprintf("SFE %s has unexpected cell_label values: %s",
                        nm, paste(bad, collapse = ", ")))
      }

      # Add purified counts as new assay (expand to full gene set,
      # overwriting RCTD-modeled genes with purified values).
      purified_rds <- file.path(step_dir, paste0(source_nm, "_purified_counts.rds"))
      if (file.exists(purified_rds)) {
        message("  Adding purified counts assay...")
        purified_sub <- readRDS(purified_rds)

        if (ncol(purified_sub) == ncol(sfe)) {
          full_purified <- assay(sfe, "counts")
          purified_genes <- rownames(purified_sub)
          shared_genes <- intersect(purified_genes, rownames(full_purified))
          full_purified[shared_genes, ] <- purified_sub[shared_genes, ]
          assay(sfe, "purified") <- full_purified
          message(sprintf("  Purified assay: %d genes (%d SPLIT-corrected) x %s cells",
                          nrow(full_purified), length(shared_genes),
                          format(ncol(full_purified), big.mark = ",")))
          rm(full_purified)
        } else {
          message(sprintf("  Skipping purified assay (dimension mismatch: %d vs %d cells)",
                          ncol(purified_sub), ncol(sfe)))
        }
        rm(purified_sub); gc(verbose = FALSE)
      } else {
        message("  No purified counts RDS found, skipping assay.")
      }

      n_changed <- sum(new_label != old_label, na.rm = TRUE)
      n_sec_new <- sum(new_label %in% sec_classes)
      n_sec_old <- sum(old_label %in% sec_classes)
      message(sprintf("  %s: %s cells changed label (%s -> %s secretory cells)",
                      nm,
                      format(n_changed, big.mark = ","),
                      format(n_sec_old, big.mark = ","),
                      format(n_sec_new, big.mark = ",")))

      save_sfe(sfe, nm)
      rm(sfe); gc(verbose = FALSE)
      message(sprintf("  [%s] %s update complete in %.1f s",
                      Sys.time(), nm,
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }

    writeLines(as.character(Sys.time()), flag_file)
    message("\n[", Sys.time(), "] All SFEs updated. Checkpoint flag written.")
  }
}

# --- Session info -----------------------------------------------------------
log_session()
message("\n[", Sys.time(), "] === 06g clean split complete ===")
