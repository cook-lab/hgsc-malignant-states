# ============================================================================
# ICON7 / GSE140082 — Strategy 2 supplemental analysis
# 00_setup.R : ensure required packages are installed
# ============================================================================
# Run order:
#   00_setup.R   -> install / load packages
#   01_parse_metadata.R
#   02_build_gene_matrix.R
#   03_score_signatures.R
#   04_descriptive.R
#   05_survival_primary.R
#   06_survival_sensitivity.R
#   07_figures.R
# ============================================================================

req_cran <- c("data.table", "survival", "survminer", "ggplot2",
              "dplyr", "tidyr", "patchwork", "readr", "stringr")
req_bioc <- c("GEOquery", "illuminaHumanv4.db", "UCell",
              "AnnotationDbi", "Biobase")

ip <- rownames(installed.packages())

to_install_cran <- setdiff(req_cran, ip)
if (length(to_install_cran) > 0) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
}

if (!"BiocManager" %in% ip) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
to_install_bioc <- setdiff(req_bioc, rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
}

ip <- rownames(installed.packages())
for (p in c(req_cran, req_bioc, "BiocManager")) {
  message(sprintf("%-25s %s", p, p %in% ip))
}
