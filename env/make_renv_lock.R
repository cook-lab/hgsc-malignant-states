#!/usr/bin/env Rscript
# ============================================================================
# renv.lock from the R library that ran the analyses (environment E1, env/README.md)
# ----------------------------------------------------------------------------
# PURPOSE
#   Records the R packages the repository's code loads, plus everything they need (recursive
#   Depends / Imports / LinkingTo), at the versions installed in the library that ran every
#   final-run R step (R 4.4.3, Bioconductor 3.20). renv is not installed in that library, so the
#   lockfile is written straight from each package's installed DESCRIPTION, in renv's lockfile
#   format (renv::restore() reads it). Each record has Package, Version, Source (CRAN,
#   Bioconductor, or GitHub with the installed commit) and Requirements; no Hash (renv adds it).
#   Packages are found by parsing the tracked R code (library / require / requireNamespace /
#   loadNamespace, pkg::fn, and `for (p in c(...)) requireNamespace(p)` loops), by a regex over
#   the R code that Python scripts run through Rscript, and by adding packages that other
#   packages load on demand (svglite for ggsave(.svg), R.utils for fread(.gz)), BiocManager,
#   BiocVersion, and the packages of upstream steps that ran in this library but whose code is
#   outside the repo (OUTSIDE_REPO below).
#   Also rewrites docs/r_key_package_versions.csv (the packages the Methods name).
#
# INPUTS
#   the R library: .Library of the R running this script (override with env var R_LOCK_LIB)
#   the repo's tracked *.R, *.py and *.sh files (git ls-files), minus renv_bootstrap.R and env/
#
# OUTPUTS
#   renv.lock
#   docs/r_key_package_versions.csv
#   stdout: package counts by source; packages the code names that are not installed;
#           packages with no recorded source
#
# MANUSCRIPT PANEL(S): none (environment record)
# RUNTIME TIER: fast (seconds)
#
# USAGE (from the repo root, with the R that ran the analyses)
#   Rscript env/make_renv_lock.R
# ============================================================================

suppressPackageStartupMessages(library(jsonlite))

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo <- normalizePath(file.path(dirname(self), ".."))
lib  <- Sys.getenv("R_LOCK_LIB", .Library)

# Packages other packages load on demand, so the code never names them
ON_DEMAND <- c("svglite",      # ggplot2::ggsave(*.svg)
               "R.utils",      # data.table::fread(*.gz)
               "BiocManager", "BiocVersion")
# Upstream steps that ran in this library, code outside the repo (cited in the Methods)
OUTSIDE_REPO <- c("scDblFinder",  # PDO object build (doublets), KPCA.B doublets
                  "infercnv")     # KPCA.B tumour-cell confirmation
# Packages the Methods name (docs/r_key_package_versions.csv); only those in the lock are written
KEY <- c("SpatialFeatureExperiment", "Voyager", "SingleCellExperiment", "SingleR", "celldex",
         "UCell", "scuttle", "scran", "scater", "edgeR", "limma", "mgcv", "spdep", "sf", "dbscan",
         "RANN", "survival", "survRM2", "data.table", "arrow", "ggplot2", "patchwork",
         "ComplexHeatmap", "Seurat", "copykat", "BayesPrism", "consensusOV", "scDblFinder",
         "infercnv", "spacexr", "SPLIT")
# Installed without source metadata but not from Bioconductor (their DESCRIPTION has biocViews)
NOT_BIOC <- c("SecAct", "RidgeFast")

# ---- installed library -------------------------------------------------------
ip <- installed.packages(lib.loc = lib, noCache = TRUE)
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
rownames(ip) <- ip[, "Package"]
base_pkgs <- rownames(ip)[!is.na(ip[, "Priority"]) & ip[, "Priority"] == "base"]

# ---- packages the code names ------------------------------------------------
files <- system2("git", c("-C", shQuote(repo), "ls-files", "--", "'*.R'", "'*.py'", "'*.sh'"),
                 stdout = TRUE)
files <- files[!grepl("^env/", files) & files != "renv_bootstrap.R"]
LOADERS <- c("library", "require", "requireNamespace", "loadNamespace", "attachNamespace")

pkgs_in_r <- function(exprs) {
  found <- character()
  walk <- function(x) {
    if (is.call(x)) {
      fn <- if (is.symbol(x[[1]])) as.character(x[[1]]) else ""
      if (fn %in% c("::", ":::")) { found <<- c(found, as.character(x[[2]])); return(invisible()) }
      if (fn %in% LOADERS && length(x) >= 2) {
        a  <- as.list(x)[-1]
        nm <- names(a); if (is.null(nm)) nm <- rep("", length(a))
        i  <- if ("package" %in% nm) which(nm == "package") else which(nm == "")[1]
        co <- "character.only" %in% nm && isTRUE(a[[which(nm == "character.only")]])
        if (!is.na(i)) {
          a1 <- a[[i]]
          if (is.character(a1)) found <<- c(found, a1)
          else if (is.symbol(a1) && !co && fn %in% c("library", "require")) found <<- c(found, as.character(a1))
        }
      }
      if (fn == "for" && length(x) == 4 && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("c"))) {
        strs <- unlist(Filter(is.character, as.list(x[[3]])[-1]))
        body <- paste(deparse(x[[4]]), collapse = " ")
        pat  <- paste0("(", paste(LOADERS, collapse = "|"), ")\\(\\s*", as.character(x[[2]]), "\\b")
        if (length(strs) && grepl(pat, body)) found <<- c(found, strs)
      }
    }
    if (is.call(x) || is.pairlist(x) || is.expression(x))
      for (i in seq_along(x)) if (!identical(x[[i]], quote(expr = ))) walk(x[[i]])   # skip empty arguments
  }
  walk(exprs)
  found
}
pkgs_in_text <- function(txt) {
  txt <- paste(txt, collapse = "\n")
  ld <- regmatches(txt, gregexpr(paste0("(?:", paste(LOADERS, collapse = "|"),
                                        ")\\(\\s*[\"']?[A-Za-z][A-Za-z0-9.]*"), txt, perl = TRUE))[[1]]
  ns <- regmatches(txt, gregexpr("\\b[A-Za-z][A-Za-z0-9.]*:::?(?=[A-Za-z._])", txt, perl = TRUE))[[1]]
  c(sub("^.*\\(\\s*[\"']?", "", ld), sub(":::?$", "", ns))
}

named <- list()
for (f in files) {
  path <- file.path(repo, f)
  txt  <- readLines(path, warn = FALSE, encoding = "UTF-8")
  p <- if (grepl("\\.R$", f)) {
    tryCatch(pkgs_in_r(parse(text = txt, keep.source = FALSE)), error = function(e) {
      message("  regex fallback for ", f, ": ", conditionMessage(e)); pkgs_in_text(txt) })
  } else pkgs_in_text(txt)
  for (x in unique(p)) named[[x]] <- c(named[[x]], f)
}
named_installed <- intersect(names(named), rownames(ip))
missing_pkgs    <- setdiff(names(named), c(rownames(ip), "R", "base"))
# Text matches in Python/shell files are noisy: report a missing package only if an R file names it
missing_pkgs <- missing_pkgs[vapply(missing_pkgs, function(p) any(grepl("\\.R$", named[[p]])), logical(1))]

# ---- closure -------------------------------------------------------------------
roots <- setdiff(union(named_installed, intersect(c(ON_DEMAND, OUTSIDE_REPO), rownames(ip))), base_pkgs)
deps  <- tools::package_dependencies(roots, db = ip, which = c("Depends", "Imports", "LinkingTo"),
                                     recursive = TRUE)
pkgs  <- sort(setdiff(intersect(unique(c(roots, unlist(deps))), rownames(ip)), base_pkgs),
              method = "radix")

# ---- records -------------------------------------------------------------------
desc <- function(p) { d <- read.dcf(file.path(lib, p, "DESCRIPTION")); setNames(as.list(d[1, ]), colnames(d)) }
reqs <- function(d) {
  f <- unlist(d[intersect(c("Depends", "Imports", "LinkingTo"), names(d))])
  if (!length(f)) return(character())
  x <- trimws(sub("\\(.*$", "", trimws(unlist(strsplit(f, ",")))))
  sort(unique(x[nzchar(x)]), method = "radix")
}
source_of <- function(p, d) {
  if (identical(d$RemoteType, "github")) return("GitHub")
  if (identical(d$Repository, "CRAN")) return("Repository")
  if (!is.null(d$Repository) && grepl("^Bioconductor", d$Repository)) return("Bioconductor")
  if (!(p %in% NOT_BIOC) && (!is.null(d$git_url) || !is.null(d$biocViews))) return("Bioconductor")
  if (!is.null(d$Repository)) return("Repository")
  "unknown"
}
records <- list()
for (p in pkgs) {
  d   <- desc(p)
  src <- source_of(p, d)
  r   <- list(Package = p, Version = d$Version, Source = src)
  if (src %in% c("Repository", "Bioconductor") && !is.null(d$Repository)) r$Repository <- d$Repository
  if (src == "GitHub")
    for (f in c("RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteSubdir", "RemoteRef", "RemoteSha"))
      if (!is.null(d[[f]])) r[[f]] <- d[[f]]
  r$Requirements <- I(reqs(d))
  records[[p]] <- r
}

bioc <- sub("\\.[0-9]+$", "", as.character(packageVersion("BiocVersion", lib.loc = lib)))
lock <- list(
  R = list(Version = paste(R.version$major, R.version$minor, sep = "."),
           Repositories = list(list(Name = "CRAN", URL = "https://cloud.r-project.org"))),
  Bioconductor = list(Version = bioc),
  Packages = records)
writeLines(toJSON(lock, pretty = TRUE, auto_unbox = TRUE), file.path(repo, "renv.lock"))

src_label <- function(r) switch(r$Source,
  Repository   = r$Repository,
  Bioconductor = paste("Bioconductor", bioc),
  GitHub       = sprintf("GitHub %s/%s@%s", r$RemoteUsername, r$RemoteRepo, substr(r$RemoteSha, 1, 7)),
  r$Source)
key <- intersect(KEY, pkgs)
write.csv(data.frame(Package = key,
                     Version = vapply(key, function(p) records[[p]]$Version, ""),
                     Source  = vapply(key, function(p) src_label(records[[p]]), "")),
          file.path(repo, "docs", "r_key_package_versions.csv"), row.names = FALSE)

# ---- report --------------------------------------------------------------------
cat(sprintf("R %s, Bioconductor %s, library %s\n", lock$R$Version, bioc, lib))
cat(sprintf("code names %d installed packages; %d with dependencies -> renv.lock\n",
            length(named_installed), length(pkgs)))
print(table(vapply(records, function(r) r$Source, "")))
cat("\nGitHub:\n")
for (r in Filter(function(r) r$Source == "GitHub", records))
  cat(sprintf("  %-14s %-12s %s\n", r$Package, r$Version, src_label(r)))
unk <- names(Filter(function(r) r$Source == "unknown", records))
cat("\nno recorded source:", if (length(unk)) paste(unk, collapse = ", ") else "none", "\n")
cat("\nnamed in R code but not installed in this library:\n")
for (p in sort(missing_pkgs)) cat(sprintf("  %-22s %s\n", p, paste(head(named[[p]], 3), collapse = ", ")))
