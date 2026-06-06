#!/usr/bin/env Rscript
# ============================================================================
# 05_clean_split_rctd.R — SPLIT/RCTD annotation cleaning (final cell_label writer)
# ============================================================================
# PURPOSE: Apply SPLIT/RCTD deconvolution to correct transcript spillover in
#   Xenium FFPE data, producing decontaminated cell-type labels and purified
#   counts. RCTD first_type replaces SingleR labels for non-secretory cells;
#   for RCTD-confirmed secretory cells, the FROZEN 06f polarization thresholds
#   assign SecA / Intermediate / SecB from existing UCell scores. This is the
#   FINAL writer of cell_label.
#
# COHORT PIN: RCTD is run on sfe_tma + the 8 published whole tissues; writeback
#   also updates sfe_tma_filtered. FTE whole-tissue samples are excluded.
#
# DETERMINISM (audit fix): RCTD/SPLIT involve stochastic optimization. The
#   global seed (set.seed(CFG$seed)) is set before each sample's RCTD run so the
#   cell labels are reproducible.
#
# NAMING: secretory subtypes use "Intermediate epithelium" (was "Transitioning").
#
# USAGE:
#   Rscript spatial/03_annotation_polarization/05_clean_split_rctd.R            # pilot
#   Rscript .../05_clean_split_rctd.R --all                                     # all SFEs
#   Rscript .../05_clean_split_rctd.R sfe_OTB_2384,sfe_OTB_2417 --suffix chunkA
#   Rscript .../05_clean_split_rctd.R --writeback                              # Part 2+3 only
#
# Reference: Marconato et al. 2026 Nat Methods; SPLIT @ bdsc-tds/SPLIT.
#
# INPUTS:
#   - <sfe_dir>/sfe_*                                  (annotated/scored SFEs)
#   - obj("xenium_ref") = xenium_celltype_downsampled.h5ad  (scRNA-seq reference)
#   - <output_root>/05_probe_qc/genes_exclude_singler.txt
#   - <output_root>/06f_reclassification_polarization/threshold_summary.csv (frozen)
#
# OUTPUTS:
#   - <output_root>/06g_clean_split/{name}_rctd_results.rds   (final RCTD per SFE)
#   - <output_root>/06g_clean_split/{name}_purified_counts.rds, _rctd_meta.csv
#   - <output_root>/06g_clean_split/06g_summary*.csv, report HTML
#   - Updated SFEs: cell_label overwritten (FINAL); purified assay added
#
# MANUSCRIPT PANEL(S): final cell_label feeds all Fig 4–7 spatial panels.
#
# RUNTIME TIER: heavy (~21 GB peak per whole tissue; ~30-40 GB for TMA)
# ============================================================================

source("spatial/00_setup/00_setup.R")

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Matrix)
  library(data.table)
  library(base64enc)
  library(ragg)
  library(knitr)
  library(patchwork)
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

# All processable SFEs (published cohort; sfe_tma_filtered handled via barcode
# join in writeback, not RCTD).
ALL_SFES <- c(sfe_names_wt, "sfe_tma")

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

# All SFEs for writeback (includes sfe_tma_filtered).
ALL_SFES_WRITEBACK <- c("sfe_tma", "sfe_tma_filtered", sfe_names_wt)

# ── Paths & constants ────────────────────────────────────────────────────────

step_dir  <- file.path(out_dir, "06g_clean_split")
fig_out   <- file.path(step_dir, "figures")
html_out  <- file.path(step_dir, "06g_clean_split_report.html")
flag_file <- file.path(step_dir, ".sfe_update_complete")

for (d in c(step_dir, fig_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sec_classes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
class_pal   <- ref_palette[sec_classes]

valid_cell_labels_extra <- c(
  "Ciliated epithelium", "Mesothelial", "Fibroblast", "Smooth muscle",
  "Pericyte", "Endothelial", "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

# The 16-type reference label for secretory cells (unsplit)
SECRETORY_REF_LABEL <- "Secretory epithelium"

# ── Load frozen polarization thresholds from 06f ─────────────────────────────

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

save_and_embed <- function(p, name, width = 1200, height = 600, res = 150) {
  pdf_path <- file.path(fig_out, paste0(name, ".pdf"))
  ggsave(pdf_path, plot = p, width = width / res, height = height / res,
         device = cairo_pdf)
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = width, height = height, res = res)
  print(p)
  dev.off()
  b64 <- base64enc::base64encode(tmp); unlink(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;">', b64)
}

make_html_table <- function(df, digits = 3) {
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) df[[j]] <- round(df[[j]], digits)
  }
  knitr::kable(df, format = "html", row.names = FALSE,
               table.attr = 'class="styled-table"')
}


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

    rctd_rds <- file.path(step_dir, paste0(nm, "_rctd_results.rds"))
    saveRDS(RCTD, rctd_rds)
    message("  Saved: ", basename(rctd_rds))

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

    # All RCTD meta CSVs must exist (sfe_tma_filtered joins from sfe_tma).
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

      # Archive current cell_label (06f)
      cd$cell_label_06f <- as.character(cd$cell_label)

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

        if (nm == "sfe_tma_filtered") {
          shared_bc <- intersect(sfe_barcodes, colnames(purified_sub))
          purified_sub <- purified_sub[, match(sfe_barcodes[sfe_barcodes %in% shared_bc],
                                                colnames(purified_sub)), drop = FALSE]
          if (ncol(purified_sub) != length(sfe_barcodes)) {
            message(sprintf("  WARNING: %d / %d barcodes matched for purified counts",
                            ncol(purified_sub), length(sfe_barcodes)))
          }
        }

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


  # ==========================================================================
  # PART 3: QC Report
  # ==========================================================================

  message("\n[", Sys.time(), "] === PART 3: QC Report ===")

  summary_files <- list.files(step_dir, pattern = "^06g_summary.*\\.csv$",
                               full.names = TRUE)
  if (length(summary_files) > 0) {
    all_summaries <- rbindlist(lapply(summary_files, fread))
    all_summaries <- unique(all_summaries, by = "sfe_name")
    fwrite(all_summaries, file.path(step_dir, "06g_summary_combined.csv"))
    message("  Combined summary: ", nrow(all_summaries), " samples")
  }

  comp_list <- list()
  for (nm in ALL_SFES_WRITEBACK) {
    sfe <- load_sfe(nm)
    cd <- as.data.table(as.data.frame(colData(sfe)[, c("cell_label", "cell_label_06f",
                                                         "rctd_first_type", "rctd_spot_class",
                                                         "rctd_first_weight",
                                                         "singler_label",
                                                         "polarization_UCell")]))
    cd[, sfe_name := nm]
    comp_list[[nm]] <- cd
    rm(sfe); gc(verbose = FALSE)
  }
  comp_dt <- rbindlist(comp_list)

  # --- Fig A: RCTD spot_class distribution per sample -----------------------

  spot_class_dt <- comp_dt[, .N, by = .(sfe_name, rctd_spot_class)]
  spot_class_dt[, pct := 100 * N / sum(N), by = sfe_name]
  spot_class_dt[is.na(rctd_spot_class), rctd_spot_class := "NA/unassigned"]

  p_spot <- ggplot(spot_class_dt, aes(x = sfe_name, y = pct, fill = rctd_spot_class)) +
    geom_col(width = 0.7) +
    scale_fill_brewer(palette = "Set2", name = "RCTD spot class") +
    labs(x = NULL, y = "% of cells",
         title = "RCTD spot class distribution per SFE") +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  img_spot <- save_and_embed(p_spot, "rctd_spot_class_per_sfe",
                              width = 1400, height = 600, res = 150)

  # --- Fig B: Top label transitions (06f -> 06g) -----------------------------

  conf_dt <- comp_dt[, .N, by = .(cell_label_06f, cell_label)]
  conf_dt[, pct := 100 * N / sum(N)]

  n_total   <- nrow(comp_dt)
  n_changed <- sum(comp_dt$cell_label != comp_dt$cell_label_06f, na.rm = TRUE)
  pct_changed <- round(100 * n_changed / n_total, 1)

  top_changes <- conf_dt[cell_label_06f != cell_label][order(-N)][1:min(.N, 15)]

  p_conf_top <- ggplot(top_changes, aes(x = reorder(paste0(cell_label_06f, " -> ", cell_label), N),
                                         y = N)) +
    geom_col(fill = "#cc4444", width = 0.7) +
    geom_text(aes(label = format(N, big.mark = ",")), hjust = -0.1, size = 2.5) +
    coord_flip() +
    labs(x = NULL, y = "Number of cells",
         title = "Top label transitions (06f -> 06g)",
         subtitle = sprintf("%s / %s cells changed (%.1f%%)",
                            format(n_changed, big.mark = ","),
                            format(n_total, big.mark = ","), pct_changed)) +
    theme_lab(base_size = 8)

  img_conf_top <- save_and_embed(p_conf_top, "top_label_transitions",
                                   width = 1200, height = 700, res = 150)

  # --- Fig C: SecA/Intermediate/SecB composition before vs after -------------

  sec_old <- comp_dt[cell_label_06f %in% sec_classes,
                      .N, by = .(sfe_name, cell_label_06f)]
  sec_old[, pct := 100 * N / sum(N), by = sfe_name]
  sec_old[, scheme := "06f"]
  setnames(sec_old, "cell_label_06f", "class")

  sec_new <- comp_dt[cell_label %in% sec_classes,
                      .N, by = .(sfe_name, cell_label)]
  sec_new[, pct := 100 * N / sum(N), by = sfe_name]
  sec_new[, scheme := "06g"]
  setnames(sec_new, "cell_label", "class")

  sec_comp <- rbind(sec_old, sec_new)
  sec_comp[, class := factor(class, levels = sec_classes)]
  sec_comp[, scheme := factor(scheme, levels = c("06f", "06g"))]

  p_sec_comp <- ggplot(sec_comp, aes(x = scheme, y = pct, fill = class)) +
    geom_col(width = 0.75, position = "stack") +
    scale_fill_manual(values = class_pal, name = "Class") +
    facet_wrap(~ sfe_name, nrow = 1) +
    labs(x = NULL, y = "% of secretory cells",
         title = "Secretory composition per SFE: 06f vs 06g (SPLIT-cleaned)") +
    theme_lab(base_size = 7) +
    theme(axis.text.x = element_text(angle = 0, size = rel(0.9)),
          strip.text = element_text(size = rel(0.65)),
          legend.position = "bottom")

  img_sec_comp <- save_and_embed(p_sec_comp, "secretory_composition_06f_vs_06g",
                                   width = 2000, height = 550, res = 150)

  # --- Fig D: Per-sample % changed ------------------------------------------

  change_dt <- comp_dt[, .(
    stayed  = sum(cell_label == cell_label_06f, na.rm = TRUE),
    changed = sum(cell_label != cell_label_06f, na.rm = TRUE)
  ), by = sfe_name]
  change_long <- melt(change_dt, id.vars = "sfe_name",
                       variable.name = "status", value.name = "n")
  change_long[, pct := 100 * n / sum(n), by = sfe_name]

  p_change <- ggplot(change_long, aes(x = sfe_name, y = pct, fill = status)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(stayed = "grey70", changed = "#cc4444"),
                       name = "Status") +
    labs(x = NULL, y = "% of all cells",
         title = "Proportion of cells that changed label (06f -> 06g SPLIT)") +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  img_change <- save_and_embed(p_change, "per_sfe_change_bar",
                                width = 1200, height = 500, res = 150)

  # --- Fig E: RCTD first_type vs singler_label concordance ------------------

  concord_dt <- comp_dt[!is.na(rctd_first_type),
                         .(n_agree = sum(rctd_first_type == singler_label),
                           n_total = .N),
                         by = sfe_name]
  concord_dt[, pct_agree := round(100 * n_agree / n_total, 1)]

  p_concord <- ggplot(concord_dt, aes(x = sfe_name, y = pct_agree)) +
    geom_col(fill = "#228B22", width = 0.6) +
    geom_text(aes(label = paste0(pct_agree, "%")), vjust = -0.3, size = 2.8) +
    labs(x = NULL, y = "% agreement",
         title = "RCTD first_type vs original SingleR label concordance") +
    ylim(0, 105) +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  img_concord <- save_and_embed(p_concord, "rctd_vs_singler_concordance",
                                  width = 1200, height = 500, res = 150)

  # --- Biological plausibility checks ---------------------------------------

  biop <- comp_dt[cell_label %in% sec_classes, .(
    n_secretory      = .N,
    n_SecA           = sum(cell_label == "SecA epithelium"),
    n_Intermediate   = sum(cell_label == "Intermediate epithelium"),
    n_SecB           = sum(cell_label == "SecB epithelium"),
    pct_SecA_of_sec  = round(100 * sum(cell_label == "SecA epithelium") / .N, 2),
    pct_Int_of_sec   = round(100 * sum(cell_label == "Intermediate epithelium") / .N, 2),
    pct_SecB_of_sec  = round(100 * sum(cell_label == "SecB epithelium") / .N, 2),
    median_pol_SecA  = round(median(polarization_UCell[cell_label == "SecA epithelium"], na.rm = TRUE), 4),
    median_pol_Int   = round(median(polarization_UCell[cell_label == "Intermediate epithelium"], na.rm = TRUE), 4),
    median_pol_SecB  = round(median(polarization_UCell[cell_label == "SecB epithelium"], na.rm = TRUE), 4)
  ), by = sfe_name]
  biop[, flag_low_SecB   := pct_SecB_of_sec < 5]
  biop[, flag_high_SecB  := pct_SecB_of_sec > 30]
  biop[, flag_pol_SecA_sign := median_pol_SecA >= 0]
  biop[, flag_pol_SecB_sign := median_pol_SecB <= 0]
  fwrite(biop, file.path(step_dir, "biological_plausibility_checks_06g.csv"))

  # --- HTML tables ----------------------------------------------------------

  tbl_spot    <- make_html_table(as.data.frame(spot_class_dt), digits = 1)
  tbl_top_ch  <- make_html_table(as.data.frame(top_changes),   digits = 1)
  tbl_concord <- make_html_table(as.data.frame(concord_dt),    digits = 1)
  tbl_biop    <- make_html_table(as.data.frame(biop),          digits = 2)

  if (exists("all_summaries")) {
    tbl_summary <- make_html_table(as.data.frame(all_summaries), digits = 1)
  } else {
    tbl_summary <- "<p>No summary CSV found.</p>"
  }

  # --- HTML report ----------------------------------------------------------

  html_css <- '
<style>
  body {font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
    max-width:1200px;margin:0 auto;padding:30px 20px;background:white;color:#333;line-height:1.6;}
  h1 {border-bottom:2px solid #333;padding-bottom:10px;font-size:1.8em;}
  h2 {color:#444;border-bottom:1px solid #ddd;padding-bottom:5px;margin-top:40px;}
  h3 {color:#555;margin-top:25px;}
  .date {color:#888;font-size:0.9em;margin-bottom:20px;}
  .meta {background:#f8f8f8;padding:12px 18px;border-radius:5px;margin:15px 0;font-size:0.9em;}
  .toc  {background:#f8f8f8;padding:15px 25px;border-radius:5px;margin:20px 0;}
  .toc a{color:#0066cc;text-decoration:none;}
  .toc ul{list-style:none;padding-left:0;} .toc li{margin:6px 0;}
  img {border:1px solid #eee;border-radius:3px;margin:10px 0;}
  .note {background:#f0f7ff;border-left:3px solid #0066cc;padding:10px 15px;margin:15px 0;font-size:0.9em;}
  .result {background:#f0fff0;border-left:3px solid #228B22;padding:10px 15px;margin:15px 0;font-size:0.9em;}
  .styled-table{border-collapse:collapse;width:100%%;font-size:0.85em;margin:15px 0;}
  .styled-table th{background:#f4f4f4;text-align:left;padding:8px 12px;border-bottom:2px solid #ddd;}
  .styled-table td{padding:6px 12px;border-bottom:1px solid #eee;}
</style>
'

  html_toc <- '
<div class="toc"><strong>Table of Contents</strong>
<ul>
  <li><a href="#s1">1. RCTD/SPLIT Summary</a></li>
  <li><a href="#s2">2. RCTD Spot Class Distribution</a></li>
  <li><a href="#s3">3. Label Transitions (06f -> 06g)</a></li>
  <li><a href="#s4">4. RCTD vs SingleR Concordance</a></li>
  <li><a href="#s5">5. Secretory Composition Comparison</a></li>
  <li><a href="#s6">6. Per-SFE Change Rates</a></li>
  <li><a href="#s7">7. Biological Plausibility</a></li>
</ul></div>
'

  meta_block <- paste0(
    '<div class="meta">',
    '<strong>Method:</strong> RCTD doublet mode (spacexr) + SPLIT purification<br>',
    '<strong>Reference:</strong> ', basename(ref_path), ' (16 types)<br>',
    '<strong>Seed:</strong> ', CFG$seed, ' (set before each RCTD run)<br>',
    '<strong>Polarization thresholds (frozen from 06f):</strong> t_low = ',
    sprintf("%.6f", t_low), ', t_high = ', sprintf("%.6f", t_high), '<br>',
    '<strong>Secretory classification:</strong> RCTD first_type == "Secretory epithelium" -> ',
    'classify_polarization(polarization_UCell, t_low, t_high)<br>',
    '<strong>Total cells changed:</strong> ',
    format(n_changed, big.mark = ","), ' / ', format(n_total, big.mark = ","),
    ' (', pct_changed, '%)<br>',
    '<strong>SFEs updated:</strong> ', length(ALL_SFES_WRITEBACK),
    ' (cell_label overwritten; cell_label_06f preserves prior; purified assay added)',
    '</div>'
  )

  html_body <- paste0(
    '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="utf-8">\n',
    '<title>06g: SPLIT/RCTD Annotation Cleaning</title>\n',
    html_css,
    '\n</head>\n<body>\n',
    '<h1>06g: SPLIT/RCTD Annotation Cleaning</h1>\n',
    '<p class="date">Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>\n',
    '<p>SPLIT/RCTD deconvolution corrects transcript spillover in Xenium FFPE data. ',
    'Each cell is modeled as a mixture of primary + secondary types; contaminant transcripts ',
    'are removed. RCTD first_type replaces SingleR labels for non-secretory cells. ',
    'For RCTD-confirmed secretory cells, existing UCell polarization scores and frozen ',
    '06f thresholds assign SecA/Intermediate/SecB subtypes.</p>\n',
    meta_block, '\n', html_toc,

    '\n<h2 id="s1">1. RCTD/SPLIT Summary</h2>\n',
    '<p>Per-sample timing, gene counts, and assignment rates.</p>\n',
    tbl_summary, '\n',

    '\n<h2 id="s2">2. RCTD Spot Class Distribution</h2>\n',
    '<p>RCTD doublet mode classifies cells as singlet, doublet_certain, ',
    'doublet_uncertain, or reject.</p>\n',
    img_spot, '\n', tbl_spot, '\n',

    '\n<h2 id="s3">3. Label Transitions (06f -> 06g)</h2>\n',
    '<p>Top cell type transitions after SPLIT cleaning.</p>\n',
    img_conf_top, '\n',
    '<h3>Top transitions</h3>\n', tbl_top_ch, '\n',

    '\n<h2 id="s4">4. RCTD vs SingleR Concordance</h2>\n',
    '<p>Agreement between RCTD first_type and the original SingleR call ',
    '(both using the 16-type reference system).</p>\n',
    img_concord, '\n', tbl_concord, '\n',

    '\n<h2 id="s5">5. Secretory Composition Comparison</h2>\n',
    '<p>SecA/Intermediate/SecB proportions before (06f) and after (06g SPLIT). ',
    'Note: the total secretory pool may change if RCTD reclassifies cells.</p>\n',
    img_sec_comp, '\n',

    '\n<h2 id="s6">6. Per-SFE Change Rates</h2>\n',
    img_change, '\n',

    '\n<h2 id="s7">7. Biological Plausibility</h2>\n',
    '<p>Flags: <code>flag_low_SecB</code> if SecB%% < 5; ',
    '<code>flag_high_SecB</code> if > 30; polarization sign flags check ',
    'SecA/SecB medians sit on expected sides of 0.</p>\n',
    tbl_biop, '\n',

    '\n<hr>\n<p style="color:#888;font-size:0.8em;">',
    'Report generated by spatial/03_annotation_polarization/05_clean_split_rctd.R on ',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>\n',
    '</body>\n</html>'
  )

  writeLines(html_body, html_out)
  message("[", Sys.time(), "] HTML report written: ", html_out)

}

# --- Session info -----------------------------------------------------------
log_session()
message("\n[", Sys.time(), "] === 06g clean split complete ===")
