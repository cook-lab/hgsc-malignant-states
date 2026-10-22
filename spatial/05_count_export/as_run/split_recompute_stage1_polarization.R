# ============================================================================
# SPLIT-purified recompute — STAGE 1: SecA/SecB UCell polarization
# ----------------------------------------------------------------------------
# Recompute SecA_UCell / SecB_UCell on the SPLIT *purified counts*, for the
# SPLIT-defined secretory cells (rctd_first_type == "Secretory epithelium"),
# faithful to 03_ucell_scoring_noBCAM.R: UCell on the counts slot over the
# ~457-gene Xenium∩Atlas shared space, maxRank = n shared genes,
# polarization = SecB - SecA, then FROZEN thresholds -> SecA/Intermediate/SecB.
# Output: per-cell purified polarization scores + relabel for all 13 samples.
# ============================================================================
suppressWarnings(suppressMessages({
  library(Matrix); library(UCell); library(Seurat); library(data.table)
}))
G <- "/Users/sarahnersesian/202609_spatial_rerun/06g_clean_split"
OUT <- "/Users/sarahnersesian/202609_spatial_rerun/split_recompute"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
shared <- readLines("/Users/sarahnersesian/202609_spatial_rerun/06d_annotation_noBCAM/shared_gene_list_xenium_atlas.txt")
SECA <- c("MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2")
SECB <- c("KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2")
T_LOW <- 0.109429; T_HIGH <- 0.210501     # frozen (SecA_p75 / SecB_p25 on polarization)

samples <- sub("_purified_counts.rds$", "", basename(Sys.glob(file.path(G, "*_purified_counts.rds"))))
cat("Stage1: ", length(samples), " samples\n")
out_csv <- file.path(OUT, "purified_polarization_scores.csv")
if (file.exists(out_csv)) file.remove(out_csv)
first <- TRUE

for (nm in samples) {
  message("[", Sys.time(), "] ", nm)
  meta <- fread(file.path(G, paste0(nm, "_rctd_meta.csv")), select = c("barcode","rctd_first_type"))
  sec_bc <- meta[rctd_first_type == "Secretory epithelium", barcode]
  if (length(sec_bc) < 10) { message("  <10 secretory, skip"); next }
  pc <- readRDS(file.path(G, paste0(nm, "_purified_counts.rds")))
  sec_bc <- intersect(sec_bc, colnames(pc))
  gkeep <- intersect(shared, rownames(pc))
  m <- pc[gkeep, sec_bc, drop = FALSE]
  m <- as(m, "dgCMatrix")
  rm(pc); gc(verbose = FALSE)
  so <- CreateSeuratObject(counts = m)
  ng <- nrow(so)
  so <- AddModuleScore_UCell(so, features = list(SecA = intersect(SECA, gkeep),
                                                 SecB = intersect(SECB, gkeep)),
                             maxRank = ng, name = "")
  pol <- so$SecB - so$SecA
  lab <- ifelse(pol < T_LOW, "SecA epithelium",
                ifelse(pol < T_HIGH, "Intermediate epithelium", "SecB epithelium"))
  dt <- data.table(sample = nm, barcode = colnames(so),
                   SecA_UCell = so$SecA, SecB_UCell = so$SecB,
                   polarization_UCell = pol, cell_label_split = lab)
  fwrite(dt, out_csv, append = !first)
  first <- FALSE
  message("  ", nrow(dt), " secretory cells scored | SecA/Int/SecB = ",
          paste(table(factor(lab, c("SecA epithelium","Intermediate epithelium","SecB epithelium"))), collapse = "/"))
  rm(m, so); gc(verbose = FALSE)
}
message("\n[", Sys.time(), "] STAGE1_DONE -> ", out_csv)
tot <- fread(out_csv)
cat("TOTAL secretory cells:", nrow(tot), "\n")
print(tot[, .N, by = cell_label_split])
