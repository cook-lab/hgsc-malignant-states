# env/: the environments that ran the analyses

The two lock files at the repo root are written from the environments that produced the results
(regenerated 2026-10-22 from the code in this repository; the package versions are those installed on
2026-10-17):

- `renv.lock` (R)
- `environment.yml` (Python)

This folder records the other environments involved and holds the scripts that write the locks.

| ID | Environment | Ran | Record |
|---|---|---|---|
| E1 | Analysis Mac: R 4.4.3 (2025-02-28), Bioconductor 3.20, one system library | See [E1](#e1-r-443-bioconductor-320) below | `../renv.lock`: 442 packages, i.e. those the code loads plus their dependencies |
| E2 | Analysis Mac: conda base env, Python 3.12.9, scientific stack installed with pip | See [E2](#e2-python-3129) below | `../environment.yml` (92 pins: the imports plus their requirements); `E2_base_full.yml` (whole env) |
| E3 | Analysis Mac: conda env `ovarian_atlas`, Python 3.12.3 | Atlas data prep, QC, and Scrublet doublet scores (2026-02-10 to 02-17). Scrublet = scanpy 1.11.4 built-in. mygene 3.2.2 | `E3_ovarian_atlas_full.yml` |
| E4 | GPU cluster | seurat_v3 HVGs, CellAssign, scVI, scANVI (2026-02-13/14) | scvi-tools **1.4.1**, from the model metadata. Python, CUDA and scikit-misc versions were not recorded |
| E5 | Second Mac: R 4.5.2 | ICON7 cohort build and UCell scores (Fig 7k–n input); copies of that code, trimmed to what the paper uses, are in `icon7/scripts/` | `renv_R4.5.2_bioc3.21_previous.lock`: the repo's renv.lock before 2026-10-17, which matches that R 4.5.2 library and not E1. The UCell version used for ICON7 is not confirmed |
| E6 | Vendor pipelines | Xenium Onboard Analysis 3.2.0.7 / 3.3.0.1 / 4.0.2.2, Space Ranger 4.0.1, Cell Ranger 9.0.0, CellBender 0.3.2, QuPath 0.7.0 | the manuscript's Methods |

### E1 (R 4.4.3, Bioconductor 3.20)

All of the final-run R:

- Xenium SFE build, QC, SingleR and UCell
- Figs 4–7, Figures S11–S18 and the inputs to Tables S10–S12

It also ran these upstream steps:

- CopyKAT, BayesPrism, atlas UCell, and consensusOV (Fig 2f)
- PDO object build (scDblFinder)
- KPCA.B (scDblFinder, inferCNV)

### E2 (Python 3.12.9)

- All of the final-run Python.
- Atlas steps from 2026-02-18 onward: refilter, UMAP/Leiden, NMF, decoupler, TCGA scoring and the SingleR reference export.
- The scFEA all-cell re-run (PyTorch) and PyDESeq2.

**Versions.** The locks pin the versions installed on 2026-10-17. Neither environment changed during
the final run: the newest R package is from 2026-08-10, the newest E2 package from 2026-09-11, and
the run started on 2026-09-13. Earlier atlas steps (Feb–May 2026) may have run on older versions of
a few packages that were updated later (in E2):

- anndata, installed 2026-05-23
- scipy, installed 2026-08-31
- igraph, installed 2026-09-11

**Restore.** The commands are the same as in the main README:

```bash
conda env create -f environment.yml && conda activate epitype-py
Rscript -e 'install.packages("renv"); renv::restore(lockfile = "renv.lock")'   # R 4.4.x
```

Notes on what to expect:

- **Pins.** The Python pins were taken on macOS arm64. GitHub packages restore at their recorded
  commits.
- **renv.lock format.** renv was not installed in E1, so `make_renv_lock.R` writes the lockfile in
  renv's format straight from the installed DESCRIPTION files. The records have no `Hash`; renv adds it.
- **Without renv.** `Rscript renv_bootstrap.R` installs the same packages at current CRAN /
  Bioconductor 3.20 versions and lists any version that differs.
- **Missing from `environment.yml`.** One repo import is not in E2: mygene, used for atlas data prep (E3).

**Regenerate.** Run these from the repo root, with the environments above:

| Command | Writes |
|---|---|
| `Rscript env/make_renv_lock.R` | `renv.lock` and `docs/r_key_package_versions.csv` |
| `bash env/export_conda_envs.sh` | `E2_base_full.yml` and `E3_ovarian_atlas_full.yml`, then runs `make_environment_yml.py` → `environment.yml` |
