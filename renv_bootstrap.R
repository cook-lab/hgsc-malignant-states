#!/usr/bin/env Rscript
# Bootstrap the R environment without renv (R 4.4.x).
# Installs every package recorded in renv.lock (the R 4.4.3 / Bioconductor 3.20 library that ran the
# analyses; env/README.md), preferring binary packages. CRAN and Bioconductor packages come at the
# CURRENT CRAN / Bioconductor 3.20 versions, so some can differ from renv.lock (the differences are
# listed at the end); GitHub packages are installed at the commit recorded in renv.lock.
# For the exact versions use renv instead:  Rscript -e 'renv::restore(lockfile = "renv.lock")'
# Run from the repo root:  Rscript renv_bootstrap.R

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
for (p in c("jsonlite", "remotes", "BiocManager"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
if (!startsWith(paste(R.version$major, R.version$minor, sep = "."), sub("\\.[0-9]+$", "", lock$R$Version)))
  warning("renv.lock was written with R ", lock$R$Version, "; this is R ", getRversion())
BiocManager::install(version = lock$Bioconductor$Version, ask = FALSE, update = FALSE)

recs <- lock$Packages
src  <- vapply(recs, function(r) r$Source, "")
have <- rownames(installed.packages())
cran <- setdiff(names(recs)[src == "Repository"], have)
bioc <- setdiff(names(recs)[src == "Bioconductor"], have)
if (length(cran)) install.packages(cran)
if (length(bioc)) BiocManager::install(bioc, ask = FALSE, update = FALSE)
for (r in recs[src == "GitHub"])
  if (!requireNamespace(r$Package, quietly = TRUE))
    remotes::install_github(sprintf("%s/%s@%s", r$RemoteUsername, r$RemoteRepo, r$RemoteSha),
                            subdir = r$RemoteSubdir, upgrade = "never")

cat("\n--- installed vs renv.lock ---\n")
ip   <- installed.packages()
miss <- setdiff(names(recs), rownames(ip))
diff <- Filter(function(r) r$Package %in% rownames(ip) && ip[r$Package, "Version"] != r$Version, recs)
cat(sprintf("%d of %d packages installed; %d at a different version\n",
            length(recs) - length(miss), length(recs), length(diff)))
if (length(miss)) cat("missing:", paste(miss, collapse = ", "), "\n")
for (r in diff) cat(sprintf("  %-26s lock %-12s installed %s\n", r$Package, r$Version, ip[r$Package, "Version"]))
