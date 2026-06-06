# Refactor conventions (READ BEFORE MIGRATING ANY SCRIPT)

This repo is the cleaned, reproducible code for the HGSC epitype manuscript. Migration =
take a CANONICAL source script (see `_repro_refactor/reports/CANONICAL_SET.md`), move it
here, and clean it to these conventions. Preserve the validated analytical logic; do not
re-derive results.

## Source → destination
- Canonical script list + their roles: `_repro_refactor/reports/CANONICAL_SET.md`.
- Panel→generator map: `_repro_refactor/reports/LINEAGE.md`.
- Atlas backend (Python) → `atlas/<NN_stage>/`; Xenium backend (R) → `spatial/<NN_stage>/`.
- One reproduction script per manuscript figure/panel → `figures/figureN/` (or `supplementary/`).
- Supp-data table generators → `tables/`.

## Mandatory cleanups for every migrated script
1. **Paths**: remove ALL hardcoded `/Volumes/CookLab/Sarah/...` (and `/Volumes/Nosepass/...`).
   - Python: `from config.config import obj, path, SEED` → `obj("atlas_final")`, `path("output_root","03_epithelial_nmf","x.csv")`.
   - R: `source(<rel>/config/config.R)` → `cfg_obj("sfe_tma_filtered")`, `cfg_path("output_root","06_spatial_stats")`.
2. **Header docstring** (every script) stating: purpose; INPUTS (objects/caches); OUTPUTS (tables/figures); which MANUSCRIPT PANEL(S) it supports; runtime tier (fast/moderate/heavy).
3. **Seeds**: set the global seed from config at the top of any stochastic step
   (`np.random.seed(SEED)` / `set.seed(CFG$seed)`). REQUIRED for: xenium 06f/06g cell-labeling,
   `20_consensusov_score.R`, NMF inits, any subsampling. (Fixes the audit's non-determinism finding.)
4. **Naming**: standardize the epithelial label **"Transitioning" → "Intermediate"** everywhere
   (matches the manuscript). Keep SecA/SecB/Ciliated as-is.
5. **Cohort**: where whole-tissue samples are enumerated, use `CFG.cohort.whole_tissue`
   (the published 8) and EXCLUDE `cohort.fte_exclude_wt` (FT1-1, EAOC-1-FTE) from the
   whole-tissue arm. (Fixes the cohort-drift finding.) FTE TMA cores (n=15) stay.
6. **Signatures**: load SecA/SecB from `shared/signatures.yml` (the noBCAM 7-gene set) —
   do not inline divergent gene lists. (Atlas: use the noBCAM-matching scoring, not `18b_v2`.)
7. **Dead code**: delete commented-out exploration, superseded branches, and `_v2/_noBCAM`
   suffixes — keep only the canonical path. Keep analytical parameters identical.
8. **Documented quirks to fix**: scFEA flux label-map bug (int `Module_id` vs `M_x` keys —
   make labels human-readable); add `torch` usage note where step 21/scFEA runs.

## Integration (the official cluster scripts, included as runnable)
The integration **is included** in `atlas/01_preprocess_qc/` as the **official cluster scripts**
(`01_preprocess` … `05_finalize`), migrated to steps `03_preprocess_hvg.py` → `04_cellassign.py`
→ `{05a_integrate_scvi, 05b_integrate_scanvi, 05c_integrate_mrvi, 05d_integrate_sysvi,
05e_integrate_harmony}.py` → `06_benchmark.py` → `07_finalize.py`. These are the authoritative
originals, path-centralised — **not reconstructions**; analytical parameters were not altered.
Per-cell QC + per-sample **Scrublet < 0.3** are applied upstream in `02_concat_qc_doublets.py`,
whose output `atlas_concatenated_filtered.h5ad` is the raw input to the chain. CellAssign
(`04_cellassign.py`) uses `shared/cellassign_markers_v3.csv` (**81-gene × 16-type**, the single
canonical marker matrix). **Five integration methods are benchmarked** (scib-metrics, step 06)
and **scANVI is SELECTED** for the manuscript (`07_finalize.py --method scanvi` → `integrated_scanvi.h5ad`,
obsm `X_scanvi`). Re-running the integration needs a GPU and is computationally expensive, so it
is optional: the deposited output (`hgsc_atlas_scanvi.h5ad`) is a config entry-point object, and
downstream analysis can start from it. The code is included so every step from raw data is
scrutable.

## Per-directory README
Each stage dir gets a short README: ordered scripts, inputs→outputs, which figures depend on it.

## Known discrepancies (record, don't silently "fix" the science)
Fig 6J circularity (published 8/8 p=0.008 vs reproducible 6/8 p=0.042), Fig 5F LIANA autocrine
(under investigation), SF7 caption "3735"→375, BayesPrism multivariate covariate set. These go in
`docs/REPRODUCIBILITY.md` as known issues — the code is migrated faithfully; corrections are the
authors' call.
