# Central config loader for R (spatial) scripts.
#
# Usage:
#   source(here::here("config", "config.R"))   # or relative path to this file
#   sfe <- loadHDF5SummarizedExperiment(cfg_obj("sfe_tma_filtered"))
#   out <- cfg_path("output_root", "06_spatial_stats")
#   set.seed(CFG$seed)
#
# Resolves ${VAR:-default} against the environment so DATA_ROOT / OUTPUT_ROOT
# can be overridden without editing config.yml.

if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")

.cfg_expand <- function(x) {
  if (is.character(x)) {
    m <- regmatches(x, gregexpr("\\$\\{([A-Za-z_][A-Za-z0-9_]*)(:-[^}]*)?\\}", x))[[1]]
    for (tok in m) {
      inner <- sub("^\\$\\{", "", sub("\\}$", "", tok))
      parts <- strsplit(inner, ":-", fixed = TRUE)[[1]]
      val <- Sys.getenv(parts[1], unset = if (length(parts) > 1) parts[2] else "")
      x <- sub(tok, val, x, fixed = TRUE)
    }
    x
  } else if (is.list(x)) {
    lapply(x, .cfg_expand)
  } else x
}

.config_file <- file.path(dirname(sys.frame(1)$ofile %||% "config/config.R"), "config.yml")
if (!file.exists(.config_file)) .config_file <- "config/config.yml"
`%||%` <- function(a, b) if (is.null(a)) b else a

CFG <- .cfg_expand(yaml::read_yaml(.config_file))
CFG$seed <- as.integer(CFG$seed %||% 42L)

cfg_obj  <- function(key) file.path(path.expand(CFG$paths$data_root), CFG$objects[[key]])
cfg_path <- function(root_key, ...) {
  base <- path.expand(CFG$paths[[root_key]])
  p <- file.path(base, ...)
  # auto-create dirs for any OUTPUT-side root (output_root, figures_dir, …); never for
  # data_root (read-only input). Directory-like path (no file extension) -> create it;
  # file-like path -> create its parent.
  if (root_key != "data_root") {
    target <- if (grepl("\\.[A-Za-z0-9]+$", basename(p))) dirname(p) else p
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
  }
  p
}
