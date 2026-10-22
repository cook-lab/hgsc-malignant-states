# ============================================================================
# 04_pathway_scoring.R
# ----------------------------------------------------------------------------
# PURPOSE: Defines the 37 UCell pathway modules (Table S9) and writes them as
#   9b_scoring/pathway_gene_sets_v2.csv, the gene-set file that
#   04c_pathway_scoring_allcells.R scores on every cell. This script wrote the
#   file in the final run (2026-09-15); 04b_pathway_gene_sets.R writes an
#   identical file.
#
# INPUTS:
#   - none (gene sets defined inline)
#
# OUTPUTS:
#   - <output_root>/9b_scoring/pathway_gene_sets_v2.csv  (pathway, gene)
#
# MANUSCRIPT PANEL(S): Tables S7 and S9 (module membership and definitions); the
#   gene sets scored by 04c (Figs 4f, 4j and 5; Figures S16–S18) and read by
#   figures/figure6/05_fig6fh_xenium_reciprocal_effects.R (Fig 6f,h).
# RUNTIME TIER: fast
# ============================================================================

# --- Config + shared setup ---------------------------------------------------
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # out_dir
set.seed(CFG$seed)

out_9b <- file.path(out_dir, "9b_scoring")
if (!dir.exists(out_9b)) dir.create(out_9b, recursive = TRUE)

# ============================================================================
# PATHWAY GENE SETS — split into pro-pathway vs regulators where needed
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
  complement_regulation = c("CD55", "CLU"),  # 2 genes

  # --- SPLIT: Immune checkpoint ---
  checkpoint_inhibitory = c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "LGALS9",
                            "TIGIT", "BTLA", "TNFSF14"),
  checkpoint_costimulatory = c("CD80", "CD86", "CD226", "CD27", "CD4", "SH2D1A"),

  # --- SPLIT: Apoptosis ---
  apoptosis_pro = c("CASP3", "CASP8", "CASP10", "FAS", "FASLG", "TP53", "BAX"),
  apoptosis_anti = c("BCL2", "BCL2L1"),  # 2 genes

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

# --- Export gene set definitions ----
gene_set_df <- data.frame(
  pathway = rep(names(pathway_sets), lengths(pathway_sets)),
  gene    = unlist(pathway_sets, use.names = FALSE),
  stringsAsFactors = FALSE
)
write.csv(gene_set_df, file.path(out_9b, "pathway_gene_sets_v2.csv"), row.names = FALSE)
message("Saved pathway_gene_sets_v2.csv")
