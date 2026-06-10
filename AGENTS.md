# AGENTS.md — working in this repository

Orientation for coding agents (and people) working in this repo. Read this first, then
`README.md` for the full reproduction walkthrough. This is a **code** repository — the large data
objects are deposited separately (see *Data*, below).

## What this is

Reproducible analysis code for the manuscript *"Tumour architecture shapes polarized epithelial
states that predict survival in high-grade serous ovarian cancer"* (Nersesian, Abou-Hamad, …,
Cook, 2026). It integrates HGSC scRNA-seq (13 studies, ~2M cells) into an atlas, resolves a
secretory-epithelial **polarization axis** — **SecA** (proliferative, progenitor-like) ↔ **SecB**
(quiescent, mucosal-injury / keratin program) — and validates it with Xenium spatial
transcriptomics, patient-derived organoids, and TCGA-OV.

The repo is **object-based**: the heavy `atlas/` and `spatial/` backends produce standardized
objects + per-stage caches; the lightweight `figures/` and `tables/` scripts load those and render
one panel/table each.

## Layout

```
config/   Single source of truth: paths, seed, cohort, signatures, polarization schema
          (config.yml + config.py / config.R loaders)
shared/   Canonical SecA/SecB signatures (signatures.yml) + CellAssign markers
atlas/    scRNA-seq backend (Python / scanpy)            01_preprocess_qc … 08_xenium_reference
spatial/  Xenium/Visium backend (R / SpatialFeatureExperiment)  00_setup … 09_external_validation
figures/  ONE script per manuscript panel (load cache → plot → export); figureN/ + supplementary/
tables/   Supplemental-data table generators
docs/     METHODS.md (expanded methods), package versions
data/     Pointers only — the large objects are NOT in the repo (see data/README.md)
output/   Everything is written here (OUTPUT_ROOT); gitignored
```

Every stage directory has a `README.md` with an ordered *script → inputs→outputs → panel* table.
`figures/README.md` is the full **panel → script** map (Fig 1A…7G, SF1…SF14, Supp Data 1–7).

## Conventions — follow these

- **Paths: never hardcode.** Resolve everything through `config/`.
  - Python: `from config.config import CFG, obj, path, SEED` → `obj("atlas_final")`,
    `path("output_root", "03_epithelial_nmf", "x.csv")`.
  - R: `source(".../config/config.R")` → `cfg_obj("sfe_tma_filtered")`,
    `cfg_path("output_root", "06_spatial_stats")`.
- **Data lives outside the repo.** Inputs resolve under `DATA_ROOT`; set it (env var, or edit
  `config.yml`) to the deposited bundle / a mounted drive. Outputs go under `OUTPUT_ROOT`
  (default `./output`). `ORGANOIDS_ROOT` points at the PDO substrate. Note the cache contract:
  **figure/table scripts read their input caches from `DATA_ROOT`** (the deposited bundle preserves
  the original directory layout the scripts expect); **backend stages write fresh results to
  `OUTPUT_ROOT`**.
- **Naming:** polarization labels are **SecA / Intermediate / SecB / Ciliated**. "Intermediate"
  was "Transitioning" in some original scripts — always use **Intermediate**.
- **Signatures:** load SecA/SecB from `shared/signatures.yml` (the noBCAM 7-gene sets). Do **not**
  inline divergent gene lists.
- **Cohort:** the published cohort lives only in `config.yml` (`cohort:` — 8 whole-tissue samples,
  FTE whole-tissue excluded, 97-patient TMA). Don't re-enumerate samples inline.
- **Determinism:** seed stochastic steps from `config.seed` / `CFG$seed` (= 42).
- **Polarization schema:** NMF **Factor_2** defines SecB; partition = SecA `<p50`,
  Intermediate `p50–p75`, SecB `>=p75` (see `config.yml`).
- **Script headers:** every backend/figure script opens with a
  `PURPOSE / INPUTS / OUTPUTS / MANUSCRIPT PANEL(S) / RUNTIME TIER` docstring naming its caches and
  producer — read it before editing, and keep the format if you add scripts.

## Integration (trust boundary)

The full pipeline from raw data **including multi-study integration** is in
`atlas/01_preprocess_qc/` (`03_preprocess_hvg.py` … `07_finalize.py`; five methods benchmarked,
**scANVI selected**). Re-running it needs a GPU and is **not required** to reproduce downstream
results — the integrated atlas (`hgsc_atlas_scanvi.h5ad`) is **deposited as the entry object**.
**Start downstream analysis from the deposited objects** (`config.objects`, e.g. `atlas_final`,
`sfe_tma_filtered`). GPU training is not bit-identical across hardware; the deposited object is the
canonical reference.

## Reproduce one panel

1. Environments: `conda env create -f environment.yml && conda activate epitype-py`;
   `Rscript -e 'renv::restore(lockfile="renv.lock")'`.
2. `export DATA_ROOT=/path/to/deposited/bundle` (see `data/README.md`).
3. Look the panel up in `figures/README.md`, then run that one script:
   `PYTHONPATH=. python figures/figure2/03_atlas_volcano_secA_secB.py` (or `Rscript figures/…R`).
   `run_panel.sh <relpath-from-repo-root>` wraps this with logging and reports new outputs.

## Where to look next

`README.md` (walkthrough) · `figures/README.md` (panel→script map) · each stage `README.md` ·
`docs/METHODS.md` (expanded methods) · `config/config.yml` (all paths, cohort, seed, signatures).
