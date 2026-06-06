#!/usr/bin/env Rscript
# ============================================================================
# Supp Data 6 — UCell pathway gene-set definitions (37 modules)
# ============================================================================
# PURPOSE
#   Emit the canonical definitions of the 37 UCell pathway modules used for
#   per-cell pathway scoring in the spatial pipeline (canonical step 9b). This
#   table generator owns ONLY the gene-set definitions and their tidy export; the
#   heavy per-SFE UCell scoring lives in spatial/05_gradients_gams (9b proper).
#   The pathway_sets list is byte-identical to the canonical 9b script so the
#   exported definitions match the scored columns exactly.
#
# INPUTS
#   - none (gene-set definitions are inlined here, as in canonical 9b)
#
# OUTPUTS
#   - output_root/9b_scoring/pathway_gene_sets_v2.csv           (pathway, gene) long form
#   - supplemental/Supplemental_Table_6_UCell_pathway_gene_sets.csv
#         tidy: pathway, module_size, gene, gene_index
#
# MANUSCRIPT PANEL(S)
#   Supp Data 6 (UCell pathway gene sets). Underlies Fig 5A / 6B / 6G and SF12-14
#   (all GAMs scored on these modules).
#
# RUNTIME TIER
#   fast (definition export only).
# ============================================================================

# --- central config (tables/ is 1 level below repo root) ---
`%||%` <- function(a, b) if (is.null(a)) b else a
.this_file <- tryCatch(
  normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  error = function(e) NA_character_)
.script_dir <- if (!is.na(.this_file)) dirname(.this_file) else "tables"
source(file.path(.script_dir, "..", "config", "config.R"))

# ============================================================================
# PATHWAY GENE SETS — 37 modules (identical to canonical step 9b)
# Split into pro-pathway vs regulator modules where signalling is bidirectional.
# ============================================================================
pathway_sets <- list(

  # --- single-direction modules ---
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
  complement_regulation = c("CD55", "CLU"),

  # --- SPLIT: Immune checkpoint ---
  checkpoint_inhibitory = c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "LGALS9",
                            "TIGIT", "BTLA", "TNFSF14"),
  checkpoint_costimulatory = c("CD80", "CD86", "CD226", "CD27", "CD4", "SH2D1A"),

  # --- SPLIT: Apoptosis ---
  apoptosis_pro = c("CASP3", "CASP8", "CASP10", "FAS", "FASLG", "TP53", "BAX"),
  apoptosis_anti = c("BCL2", "BCL2L1"),

  # --- metabolic / adhesion / target modules ---
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

# ============================================================================
# Long-form gene-set definitions (matches canonical 9b export)
# ============================================================================
gene_set_df <- data.frame(
  pathway = rep(names(pathway_sets), lengths(pathway_sets)),
  gene    = unlist(pathway_sets, use.names = FALSE),
  stringsAsFactors = FALSE
)
out_long <- cfg_path("output_root", "9b_scoring", "pathway_gene_sets_v2.csv")
write.csv(gene_set_df, out_long, row.names = FALSE)
message("Saved: ", out_long)

# ============================================================================
# Supp Data 6 — tidy table with pinned module_size + gene index
# ============================================================================
supp6 <- data.frame(
  pathway     = rep(names(pathway_sets), lengths(pathway_sets)),
  module_size = rep(lengths(pathway_sets), lengths(pathway_sets)),
  gene        = unlist(pathway_sets, use.names = FALSE),
  stringsAsFactors = FALSE
)
supp6$gene_index <- ave(seq_len(nrow(supp6)), supp6$pathway, FUN = seq_along)
supp6 <- supp6[order(supp6$pathway, supp6$gene_index), ]

out_supp <- cfg_path("output_root", "supplemental",
                     "Supplemental_Table_6_UCell_pathway_gene_sets.csv")
write.csv(supp6, out_supp, row.names = FALSE)
message("Saved: ", out_supp, "  (", nrow(supp6), " rows)")
