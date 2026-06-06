# ============================================================================
# 04_pathway_scoring.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-cell pathway activity scoring with UCell (37 modules, rank-based, maxRank = n_genes to match 06d).
#
# INPUTS:
#   - SFEs (load_sfe): sfe_tma_filtered + 8 whole-tissue (cell_label)
#
# OUTPUTS:
#   - SFEs updated with pathway_* colData
#   - output/9b_scoring/pathway_gene_sets_v2.csv (Supp Data 6)
#   - output/9b_scoring/scoring_summary_v2.csv
#
# MANUSCRIPT PANEL(S): Supp Data 6 (UCell pathway gene sets); feeds Fig 5A/6B/6G/6H GAMs.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

library(UCell)
library(Seurat)

out_9b <- file.path(out_dir, "9b_scoring")
if (!dir.exists(out_9b)) dir.create(out_9b, recursive = TRUE)

message("\n=== 9b: Per-Cell Pathway Scoring (UCell v2) ===")

# ============================================================================
# PATHWAY GENE SETS — split into pro-pathway vs regulators where needed
# ============================================================================

pathway_sets <- list(

  # --- Unchanged (all genes same direction) ---
  proliferation = c("CCNB1", "CDC20", "CDK1", "MKI67", "PCNA", "STMN1", "TOP2A", "TUBB"),

  hypoxia = c("CAT", "ENO1", "EPAS1", "HIF1A", "LDHA", "NFE2L2", "PDK1", "PGK1",
              "PRDX1", "SLC2A1", "TXN", "VEGFA"),

  rtk_ras = c("ABL1", "BRAF", "CBL", "EGFR", "ERBB3", "ERBB4", "ERRFI1", "FGFR1",
              "FLT3", "HRAS", "JAK2", "KIT", "KRAS", "MAP2K1", "MAP2K2", "MAPK1",
              "MAPK3", "MET", "NF1", "NRAS", "NTRK2", "PDGFRA", "PTPN11", "RAC1",
              "RAF1", "RASA1", "RIT1", "SOS1"),

  tgfb = c("INHBA", "INHBB", "SMAD2", "SMAD3", "SMAD4", "TGFB1", "TGFB3",
           "TGFBR1", "TGFBR2"),

  hippo = c("CSNK1D", "CSNK1E", "FAT1", "LATS1", "LATS2", "MOB1A", "MOB1B",
            "SAV1", "STK4", "TAOK1", "TAOK3", "WWC1", "YAP1"),

  jak_stat = c("CSF3R", "EGFR", "IFNAR1", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2",
               "IL10", "IL10RA", "IL13RA1", "IL15", "IL2RA", "IL2RG", "IL4R",
               "IL6R", "IL7R", "JAK1", "JAK2", "JAK3", "LIFR", "OSM", "PDGFA",
               "PDGFB", "PDGFRA", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A",
               "STAT6", "TYK2"),

  type_i_ifn = c("IFI44L", "IFI6", "IFIT1", "IFIT3", "IFNAR1", "IFNAR2", "IRF9",
                 "ISG20", "MX1", "OAS1", "RSAD2", "SOCS1", "STAT1", "STAT2"),

  type_ii_ifn = c("CIITA", "CXCL9", "CXCL10", "CXCL11", "IDO1", "IFNG", "IFNGR1",
                  "IFNGR2", "IRF1", "STAT1"),

  cytotoxicity = c("CST7", "FASLG", "FGFBP2", "GZMH", "NCR1", "PRF1"),

  antigen_presentation = c("CALR", "CANX", "CIITA", "NLRC5", "PDIA3", "PSME1",
                           "PSME2", "PSME3", "TAP1", "TAP2", "TAPBP"),

  emt = c("ACTA2", "CDH2", "CTHRC1", "DCN", "INHBA", "MMP11", "POSTN", "SERPINE1",
          "SFRP2", "SFRP4", "TAGLN", "TGFBI", "VCAN"),

  angiogenesis = c("DLL4", "ENG", "F3", "GNG11", "HSPG2", "NRP2", "PDGFA", "PDGFB",
                   "PECAM1", "RAMP2", "VEGFA", "VWF"),

  chemokine = c("CCR1", "CCR5", "CCR7", "CCRL2", "CX3CR1", "CXCL3", "CXCL9",
                "CXCL10", "CXCL11", "CXCL17", "CXCR4", "CXCR6"),

  nfkb = c("BIRC3", "CD40", "FAS", "ICAM1", "NFKB1", "NFKB2", "NFKBIA", "RELA",
           "TNF", "TRAF6", "VCAM1"),

  # --- SPLIT: Wnt ---
  wnt_activating = c("WNT7A", "FZD1", "FZD3", "FZD4", "CTNNB1", "LRP6", "TCF7", "TCF7L2"),
  wnt_inhibitory = c("APC", "GSK3B", "DKK3", "SFRP2", "SFRP4", "TLE1", "TLE3",
                     "TLE4", "TCF7L1"),

  # --- SPLIT: Notch ---
  notch_activating = c("JAG1", "JAG2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "EP300"),
  notch_inhibitory = c("FBXW7", "KDM5A", "NCOR1", "NCOR2", "SPEN"),

  # --- SPLIT: PI3K/AKT/mTOR ---
  pi3k_activating = c("AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB", "PIK3R1",
                      "PIK3R3", "MTOR", "RHEB", "RICTOR"),
  pi3k_inhibitory = c("PTEN", "INPP4B", "TSC2", "STK11", "PPP2R1A"),

  # --- SPLIT: p53 ---
  p53_activating = c("TP53", "ATM", "CHEK2", "CDKN2A"),
  p53_inhibitory = c("MDM2", "MDM4"),

  # --- SPLIT: MYC ---
  myc_activating = c("MYC", "MYCL", "MAX", "MLX", "MLXIP"),
  myc_inhibitory = c("MGA", "MNT", "MXD1", "MXD4", "MXI1"),

  # --- SPLIT: Complement ---
  complement_activation = c("C1QC", "C1S", "C3", "C7"),
  complement_regulation = c("CD55", "CLU"),  # only 2 — will flag if UCell fails

  # --- SPLIT: Immune checkpoint ---
  checkpoint_inhibitory = c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "LGALS9",
                            "TIGIT", "BTLA", "TNFSF14"),
  checkpoint_costimulatory = c("CD80", "CD86", "CD226", "CD27", "CD4", "SH2D1A"),

  # --- SPLIT: Apoptosis ---
  apoptosis_pro = c("CASP3", "CASP8", "CASP10", "FAS", "FASLG", "TP53", "BAX"),
  apoptosis_anti = c("BCL2", "BCL2L1"),  # only 2 — will flag

  # --- NEW pathways ---
  glycolysis = c("ENO1", "LDHA", "PDK1", "PGK1", "SLC2A1", "SLC16A3"),

  lactate_metabolism = c("LDHA", "SLC16A3", "SLC2A1", "PDK1"),

  ecm_adhesion = c("CD44", "CDH2", "CDH6", "CLDN3", "DCN", "EPCAM", "ITGAV",
                   "ITGB5", "ITGB6", "POSTN", "VCAN"),

  matrix_remodeling = c("ADAM10", "ADAM17", "F3", "MMP7", "MMP11", "SERPINE1"),

  oxidative_stress = c("CAT", "GPX3", "NFE2L2", "PRDX1", "TXN"),

  survival_antiapoptotic = c("BCL2", "BCL2L1", "BIRC3", "CD55", "CLU"),

  adc_targets = c("TACSTD2", "NECTIN4", "MUC16", "MET", "EGFR", "ERBB2", "ERBB3", "VTCN1")
)

message(sprintf("Defined %d pathway modules (%d unique genes)",
                length(pathway_sets), length(unique(unlist(pathway_sets)))))

# --- Export gene set definitions ----
gene_set_df <- data.frame(
  pathway = rep(names(pathway_sets), lengths(pathway_sets)),
  gene    = unlist(pathway_sets, use.names = FALSE),
  stringsAsFactors = FALSE
)
write.csv(gene_set_df, file.path(out_9b, "pathway_gene_sets_v2.csv"), row.names = FALSE)
message("Saved pathway_gene_sets_v2.csv")

# --- SFE names ----
sfe_names <- c("sfe_tma_filtered",
               "sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
               "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
               "sfe_SP24_24824", "sfe_SP24_25573")

# ============================================================================
# SCORING LOOP — UCell on each SFE
# ============================================================================

summary_list <- list()

for (sname in sfe_names) {

  message("\n", strrep("-", 60))
  message("Processing ", sname, " ...")
  t0 <- Sys.time()

  sfe <- load_sfe(sname)
  n_total <- ncol(sfe)
  n_genes <- nrow(sfe)
  rn <- rownames(sfe)

  message(sprintf("  %s cells, %d genes", format(n_total, big.mark = ","), n_genes))

  # --- Remove old pathway_ columns if they exist ---
  old_pw <- grep("^pathway_", colnames(colData(sfe)), value = TRUE)
  if (length(old_pw) > 0) {
    message(sprintf("  Removing %d old pathway columns", length(old_pw)))
    for (col in old_pw) colData(sfe)[[col]] <- NULL
  }

  # --- Filter gene sets to what's on the panel ---
  filtered_sets <- lapply(pathway_sets, function(genes) intersect(genes, rn))
  n_per_set <- sapply(filtered_sets, length)

  # Report coverage
  for (nm in names(filtered_sets)) {
    ng <- n_per_set[nm]
    nt <- length(pathway_sets[[nm]])
    if (ng < 2) {
      message(sprintf("  WARNING: %s — only %d/%d genes, will score as NA", nm, ng, nt))
    } else if (ng < nt) {
      message(sprintf("  %s: %d/%d genes on panel", nm, ng, nt))
    }
  }

  # Drop sets with < 2 genes (UCell minimum)
  scoreable <- filtered_sets[n_per_set >= 2]
  not_scoreable <- names(filtered_sets)[n_per_set < 2]

  if (length(not_scoreable) > 0) {
    message(sprintf("  Skipping %d sets (< 2 genes): %s",
                    length(not_scoreable), paste(not_scoreable, collapse = ", ")))
  }

  message(sprintf("  Scoring %d pathway modules with UCell (maxRank = %d) ...",
                  length(scoreable), n_genes))

  # --- Convert to Seurat for UCell ---
  # Use counts slot (UCell ranks on raw counts, consistent with 06d)
  seurat <- CreateSeuratObject(
    counts = counts(sfe),
    meta.data = as.data.frame(colData(sfe))
  )

  # --- Run UCell ---
  seurat <- AddModuleScore_UCell(
    seurat,
    features = scoreable,
    maxRank = n_genes,  # all genes — same as 06d
    name = NULL         # no suffix added
  )

  # --- Transfer scores back to SFE colData ---
  # UCell with name=NULL uses exact list names (no _UCell suffix)
  n_scored <- 0

  for (i in seq_along(scoreable)) {
    set_name <- names(scoreable)[i]
    sfe_col <- paste0("pathway_", set_name)

    # Try exact name first, then with _UCell suffix as fallback
    if (set_name %in% colnames(seurat@meta.data)) {
      colData(sfe)[[sfe_col]] <- seurat@meta.data[[set_name]]
      n_scored <- n_scored + 1
    } else if (paste0(set_name, "_UCell") %in% colnames(seurat@meta.data)) {
      colData(sfe)[[sfe_col]] <- seurat@meta.data[[paste0(set_name, "_UCell")]]
      n_scored <- n_scored + 1
    } else {
      message(sprintf("  WARNING: %s not found in UCell output", set_name))
      colData(sfe)[[sfe_col]] <- NA_real_
    }
  }

  # Set NA for non-scoreable sets
  for (nm in not_scoreable) {
    colData(sfe)[[paste0("pathway_", nm)]] <- NA_real_
  }

  message(sprintf("  Scored %d/%d modules successfully", n_scored, length(pathway_sets)))

  # --- Quick summary stats ---
  pw_cols <- grep("^pathway_", colnames(colData(sfe)), value = TRUE)
  for (pc in pw_cols) {
    vals <- colData(sfe)[[pc]]
    if (all(is.na(vals))) {
      message(sprintf("  %-35s: ALL NA", pc))
    } else {
      message(sprintf("  %-35s: mean=%.4f  range=[%.4f, %.4f]",
                      pc, mean(vals, na.rm = TRUE),
                      min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
    }
  }

  # --- Realize and save ---
  message("  Realizing assays...")
  for (a in assayNames(sfe)) {
    assay(sfe, a) <- as(assay(sfe, a), "dgCMatrix")
  }

  message("  Saving SFE...")
  save_sfe(sfe, sname)

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  message(sprintf("  Done in %.1f min", elapsed))

  summary_list[[sname]] <- data.frame(
    sample_id = sname,
    n_cells = n_total,
    n_pathways_scored = n_scored,
    n_pathways_na = length(not_scoreable),
    elapsed_min = elapsed,
    stringsAsFactors = FALSE
  )

  rm(sfe, seurat)
  gc(verbose = FALSE)
}

# ============================================================================
# Summary
# ============================================================================

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_9b, "scoring_summary_v2.csv"), row.names = FALSE)

message("\n", strrep("=", 60))
message("=== 9b Complete (UCell v2) ===")
message(sprintf("Scored %s cells across %d samples with %d pathway modules",
                format(sum(summary_df$n_cells), big.mark = ","),
                nrow(summary_df),
                length(pathway_sets)))
print(summary_df)
message(strrep("=", 60))
