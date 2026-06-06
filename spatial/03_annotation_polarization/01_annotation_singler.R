# ============================================================================
# 01_annotation_singler.R — Cell type annotation via SingleR
# ============================================================================
# PURPOSE: Annotate all Xenium cells with SingleR using the matched scRNA-seq
#   reference (16K cells, 16 cell types). Prepares the reference on the
#   panel-gene intersection (excluding probe-QC-failed genes), runs SingleR
#   per-sample, validates against prior TMA annotations, and summarizes
#   composition.
#
# COHORT PIN: operates over sfe_tma + the 8 published whole tissues
#   (CFG$cohort$whole_tissue); FTE whole-tissue samples excluded.
#
# INPUTS:
#   - <sfe_dir>/sfe_tma, <sfe_dir>/sfe_<wt>            (normalized SFEs)
#   - obj("xenium_ref")  = xenium_celltype_downsampled.h5ad  (SingleR reference)
#   - <output_root>/05_probe_qc/genes_exclude_singler.txt  (28 excluded genes)
#   - <output_root>/05_probe_qc/genes_monitor_singler.txt
#   - <data_root>/.../previously processed/xenium_tma_final_sfe_hdf5  (old labels)
#
# OUTPUTS:
#   - SFEs updated with singler_label / singler_pruned / singler_score
#   - <output_root>/06_annotation/*.csv   (composition, validation, summary)
#   - <output_root>/figures/annotation/*.pdf
#
# MANUSCRIPT PANEL(S): annotation backend for Fig 4–7 composition/ROI panels;
#   first step of the noBCAM polarization chain.
#
# RUNTIME TIER: heavy (SingleR over millions of cells)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(SingleR)
library(zellkonverter)

message("=== Cell Type Annotation (SingleR) ===")

# ── 1. Setup & gene preparation ──────────────────────────────────────────────

out_path <- file.path(out_dir, "06_annotation")
fig_path <- file.path(fig_dir, "annotation")
for (d in c(out_path, fig_path)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

exclude_genes <- readLines(file.path(out_dir, "05_probe_qc",
                                     "genes_exclude_singler.txt"))
exclude_genes <- trimws(exclude_genes[nchar(trimws(exclude_genes)) > 0])
message("Genes to exclude from SingleR: ", length(exclude_genes))

monitor_genes <- readLines(file.path(out_dir, "05_probe_qc",
                                     "genes_monitor_singler.txt"))
monitor_genes <- trimws(monitor_genes[nchar(trimws(monitor_genes)) > 0])
message("Genes to monitor: ", length(monitor_genes))

sfe_names <- sfe_names_all   # sfe_tma + published 8 whole tissues (cohort PIN)

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

results_list <- list()

for (sfe_name in sfe_names) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Annotating: ", sfe_name)
  message(paste(rep("=", 60), collapse = ""))

  sfe_path <- file.path(sfe_dir, sfe_name)
  sfe <- loadHDF5SummarizedExperiment(dir = sfe_path)

  if ("singler_label" %in% colnames(colData(sfe))) {
    message("  Already annotated (singler_label found). Skipping.")
    results_list[[sfe_name]] <- data.table(
      sfe_name = sfe_name,
      n_cells  = ncol(sfe),
      status   = "skipped"
    )
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
  singler_res <- SingleR(
    test           = sfe_sub,
    ref            = ref_sce,
    labels         = ref_labels,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts"
  )
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

  results_list[[sfe_name]] <- data.table(
    sfe_name  = sfe_name,
    n_cells   = n_cells,
    n_pruned  = n_pruned,
    pct_pruned = round(100 * n_pruned / n_cells, 2),
    elapsed_min = as.numeric(elapsed),
    status    = "annotated"
  )

  rm(sfe, sfe_sub, singler_res, scores_mat, best_score, median_score)
  gc(verbose = FALSE)
}

# ── 4. Collect labels & composition summary ─────────────────────────────────

message("\n--- Collecting labels across all samples ---")

all_labels <- list()
for (sfe_name in sfe_names) {
  sfe <- loadHDF5SummarizedExperiment(dir = file.path(sfe_dir, sfe_name))
  cd <- as.data.table(colData(sfe)[, c("sample_id", "singler_label",
                                        "singler_pruned", "singler_score")])
  cd$sfe_name <- sfe_name
  all_labels[[sfe_name]] <- cd
  rm(sfe); gc(verbose = FALSE)
}
all_labels_dt <- rbindlist(all_labels)

comp <- all_labels_dt[, .N, by = .(sfe_name, singler_label)]
comp[, total := sum(N), by = sfe_name]
comp[, pct := round(100 * N / total, 2)]
comp <- comp[order(sfe_name, -N)]

fwrite(comp, file.path(out_path, "composition_per_sample.csv"))
message("Saved: composition_per_sample.csv")

overall <- all_labels_dt[, .N, by = singler_label][order(-N)]
overall[, pct := round(100 * N / sum(N), 2)]
message("\nOverall cell type composition:")
for (i in seq_len(nrow(overall))) {
  message(sprintf("  %-20s %8s (%5.1f%%)",
                  overall$singler_label[i],
                  format(overall$N[i], big.mark = ","),
                  overall$pct[i]))
}

summary_dt <- rbindlist(results_list, fill = TRUE)
fwrite(summary_dt, file.path(out_path, "annotation_summary.csv"))
message("Saved: annotation_summary.csv")

fwrite(all_labels_dt, file.path(out_path, "singler_results_per_sample.csv"))
message("Saved: singler_results_per_sample.csv (",
        format(nrow(all_labels_dt), big.mark = ","), " cells)")

# ── 5. TMA validation against prior annotations ──────────────────────────────

message("\n--- TMA validation against prior annotations ---")

old_tma_path <- file.path(data_dir, "previously processed",
                          "xenium_tma_final_sfe_hdf5")

if (dir.exists(old_tma_path)) {

  old_tma <- loadHDF5SummarizedExperiment(dir = old_tma_path)
  new_tma <- loadHDF5SummarizedExperiment(dir = file.path(sfe_dir, "sfe_tma"))

  message(sprintf("  Old TMA: %s cells, New TMA: %s cells",
                  format(ncol(old_tma), big.mark = ","),
                  format(ncol(new_tma), big.mark = ",")))

  shared_cells <- intersect(colnames(old_tma), colnames(new_tma))
  message(sprintf("  Matched cells: %s / %s (%.1f%%)",
                  format(length(shared_cells), big.mark = ","),
                  format(ncol(new_tma), big.mark = ","),
                  100 * length(shared_cells) / ncol(new_tma)))

  if (length(shared_cells) > 0) {

    old_labels <- as.character(old_tma[, shared_cells]$final_annotation)
    new_labels <- new_tma[, shared_cells]$singler_label

    message("  Old label levels: ", paste(sort(unique(old_labels)), collapse = ", "))
    message("  New label levels: ", paste(sort(unique(new_labels)), collapse = ", "))

    # Map old Level-2 TMA names -> reference Level-1 names for concordance.
    # Epithelial polarization classes (incl. Intermediate, formerly
    # "Transitioning") all collapse to "Secretory epithelium" at Level 1.
    old_to_ref <- c(
      "Epi_SecA"         = "Secretory epithelium",
      "Epi_SecB"         = "Secretory epithelium",
      "Epi_Intermediate" = "Secretory epithelium",
      "Epi_Cil"          = "Ciliated epithelium",
      "Epi"              = "Secretory epithelium",
      "Mesenchymal"      = "Fibroblast",
      "Pericyte_SM"      = "Pericyte",
      "Endothelial"      = "Endothelial",
      "T cell"           = "T cell",
      "NK cell"          = "NK cell",
      "B cell"           = "B cell",
      "Plasma cell"      = "Plasma cell",
      "Macrophage"       = "Macrophage",
      "DC"               = "Conventional dendritic cell",
      "Neutrophil"       = "Neutrophil",
      "Mast"             = "Mast cell"
    )

    old_labels_mapped <- ifelse(old_labels %in% names(old_to_ref),
                                old_to_ref[old_labels],
                                old_labels)

    conf_dt <- data.table(old_label = old_labels_mapped, new_label = new_labels)
    conf_mat <- conf_dt[, .N, by = .(old_label, new_label)]
    conf_wide <- dcast(conf_mat, old_label ~ new_label, value.var = "N",
                       fill = 0)

    fwrite(conf_wide, file.path(out_path, "tma_confusion_matrix.csv"))
    message("Saved: tma_confusion_matrix.csv")

    all_types <- union(unique(old_labels_mapped), unique(new_labels))
    concordance <- data.table(cell_type = character(),
                              n_old = integer(),
                              n_new = integer(),
                              n_agree = integer(),
                              precision = numeric(),
                              recall = numeric(),
                              f1 = numeric())

    for (ct in sort(all_types)) {
      is_old <- old_labels_mapped == ct
      is_new <- new_labels == ct
      n_old  <- sum(is_old, na.rm = TRUE)
      n_new  <- sum(is_new, na.rm = TRUE)
      n_agree <- sum(is_old & is_new, na.rm = TRUE)
      prec <- if (n_new > 0) n_agree / n_new else NA_real_
      rec  <- if (n_old > 0) n_agree / n_old else NA_real_
      f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
        2 * prec * rec / (prec + rec)
      } else NA_real_
      concordance <- rbind(concordance, data.table(
        cell_type = ct, n_old = n_old, n_new = n_new, n_agree = n_agree,
        precision = round(prec, 4), recall = round(rec, 4),
        f1 = round(f1, 4)
      ))
    }

    concordance <- concordance[order(-f1)]
    fwrite(concordance, file.path(out_path, "tma_validation.csv"))
    message("Saved: tma_validation.csv")

    overall_agree <- sum(old_labels_mapped == new_labels, na.rm = TRUE)
    overall_pct <- 100 * overall_agree / length(shared_cells)
    message(sprintf("\n  Overall concordance: %s / %s (%.1f%%)",
                    format(overall_agree, big.mark = ","),
                    format(length(shared_cells), big.mark = ","),
                    overall_pct))

    message("\n  Per-type concordance:")
    for (i in seq_len(nrow(concordance))) {
      r <- concordance[i]
      flag <- if (!is.na(r$f1) && r$f1 < 0.80) " *** LOW" else ""
      message(sprintf("    %-20s F1=%.3f  prec=%.3f  rec=%.3f  (n=%s)%s",
                      r$cell_type, r$f1, r$precision, r$recall,
                      format(r$n_old, big.mark = ","), flag))
    }

    # ── Validation figures ──────────────────────────────────────────────────

    conf_mat_m <- as.matrix(conf_wide[, -1])
    rownames(conf_mat_m) <- conf_wide$old_label
    conf_prop <- conf_mat_m / rowSums(conf_mat_m)

    pdf(file.path(fig_path, "tma_confusion_heatmap.pdf"),
        width = 10, height = 8)
    ht <- Heatmap(
      conf_prop,
      name = "Proportion",
      col  = circlize::colorRamp2(c(0, 0.5, 1),
                                  c("white", "#FDE725", "#440154")),
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      row_title       = "Old annotation (mapped)",
      column_title    = "New SingleR label",
      cell_fun = function(j, i, x, y, width, height, fill) {
        count <- conf_mat_m[i, j]
        if (count > 0) {
          grid.text(format(count, big.mark = ","),
                    x, y, gp = gpar(fontsize = 7,
                    col = ifelse(conf_prop[i, j] > 0.5, "white", "black")))
        }
      },
      row_names_gp    = gpar(fontsize = 9),
      column_names_gp = gpar(fontsize = 9),
      column_names_rot = 45,
      heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                  labels_gp = gpar(fontsize = 8))
    )
    draw(ht)
    dev.off()
    message("Saved: tma_confusion_heatmap.pdf")

    conc_plot <- concordance[!is.na(f1)] |>
      ggplot(aes(x = reorder(cell_type, f1), y = f1)) +
      geom_col(fill = "steelblue", width = 0.7) +
      geom_hline(yintercept = 0.80, linetype = "dashed", color = "red",
                 linewidth = 0.5) +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "F1 Score",
           title = "TMA annotation concordance (old vs. new)",
           subtitle = "Dashed line = 0.80 threshold") +
      theme_lab()

    ggsave(file.path(fig_path, "tma_concordance_barplot.pdf"),
           conc_plot, width = 6, height = 5)
    message("Saved: tma_concordance_barplot.pdf")

    rm(old_tma, new_tma, conf_dt, conf_mat, conf_wide, concordance)
    gc(verbose = FALSE)

  } else {
    warning("  No matching cell IDs between old and new TMA SFEs!")
  }

} else {
  message("  Old TMA SFE not found at: ", old_tma_path)
  message("  Skipping TMA validation.")
}

# ── 6. Composition visualization ─────────────────────────────────────────────

message("\n--- Composition visualization ---")

ref_types <- sort(unique(all_labels_dt$singler_label))
# Stacked bar plot uses the shared ref_palette (defined in 00_setup.R).

comp_plot_dt <- comp[, .(sfe_name, singler_label, pct)]

p_comp <- ggplot(comp_plot_dt,
                 aes(x = sfe_name, y = pct, fill = singler_label)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = ref_palette, name = "Cell type") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion (%)",
       title = "Cell type composition per sample",
       subtitle = "SingleR annotation with scRNA-seq reference") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(fig_path, "composition_barplot.pdf"),
       p_comp, width = 10, height = 6)
message("Saved: composition_barplot.pdf")

# ── 7. Final summary ─────────────────────────────────────────────────────────

message("\n", paste(rep("=", 60), collapse = ""))
message("ANNOTATION SUMMARY")
message(paste(rep("=", 60), collapse = ""))

total_cells <- sum(overall$N)
total_pruned <- sum(summary_dt$n_pruned, na.rm = TRUE)
message(sprintf("Total cells annotated: %s", format(total_cells, big.mark = ",")))
message(sprintf("Total pruned (low-confidence): %s (%.1f%%)",
                format(total_pruned, big.mark = ","),
                100 * total_pruned / total_cells))
message(sprintf("Genes used for SingleR: %d (after excluding %d)",
                length(shared_genes), n_excluded))
message(sprintf("Cell types assigned: %d", length(ref_types)))

message("\nDone. Review outputs:")
message("  Tables: ", out_path)
message("  Figures: ", fig_path)
message("  Annotated SFEs: ", sfe_dir, "/sfe_*")
log_session()
