# Reproducibility

This document summarizes the reproducibility audit for the HGSC epitype manuscript and tracks the known issues carried into the refactored monorepo. The full audit record — scorecard, per-area evidence, and environment spec — is at:

**`_repro_refactor/reports/REPRODUCIBILITY_REPORT.md`** (Phase 3 deliverable, compiled 2026-06-05)

## Verdict: strongly reproducible

The manuscript is **strongly reproducible**. The defining result — the secretory polarization axis — reproduces **exactly**: re-running the core NMF (`atlas/03_epithelial_nmf/01_epithelial_nmf.py`) recovered SecA = Factor 3 and SecB = Factor 2 with identical assignments, gene-loading correlation **1.000** ref↔new for both factors, and 15/15 overlap in each factor's top-15 genes. The survival story is bit-reproducible: Xenium polarization Cox HRs (Fig 7A OS HR=1.45, 7B PFS HR=1.42) and TCGA-OV validation (Fig 7G OS HR=2.00 p=0.004, PFS HR=1.79 p=0.009; KM p=0.010) match the reference to printed precision.

Quantitatively:
- **Figures:** 42/67 panel rows reproduced first-pass → **54/67** after one path-typo symlink fix (the 12 unblocked panels had already matched the ledger in statistics before a render crash).
- **Backends:** **21/34 exact** (several byte-identical), **11 within-tolerance**, **1 divergent** (Fig 6J), 0 unresolved.
- **Spatial backend re-ran clean post-synthesis:** `44` autocorrelation regenerated all 7 tables (WT Lee's L p=0.001, TMA per-core significant → Fig 4E/F); GAM steps `16b`/`19d`/`19e` all exit 0.

**Integration:** the integration that produces the 1,980,703-cell atlas embedding **is included in the repository** (the official cluster scripts `atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`). Per-cell QC + per-sample Scrublet (< 0.3) are applied upstream in `02_concat_qc_doublets.py` (→ `atlas_concatenated_filtered.h5ad`, the chain's raw input); CellAssign annotations (`cellassign_markers_v3.csv`, 81×16) seed semi-supervised integration; **five methods are benchmarked** (scVI, scANVI, MrVI, SysVI, Harmony; scib-metrics) and **scANVI is selected** (`07_finalize --method scanvi` → `integrated_scanvi.h5ad`, obsm `X_scanvi` → `08_refilter_umap.py` → deposited `hgsc_atlas_scanvi.h5ad`). We did **not re-run** the GPU integration during validation — a validation-time choice, since it is expensive — and used the deposited integrated object as input; but the integration is fully included and scrutable, and is reproducible from raw data via the included code. The four stochastic tools (CopyKAT, BayesPrism, LIANA, scFEA) were flanking-audited (cached outputs checked for consistency; deterministic input-prep + consumers re-run), not re-executed at their stochastic cores. The `2026_organoids/` directory is external (Fig 3 TRUST-EXISTING per author decision).

## Known issues

### Material — require author action before publication
- **Fig 6J — macrophage circularity (circularity 8/8 vs 6/8).** Published: 8/8 samples, p=0.008. Reproducible value, and the reference cache itself: **6/8, p=0.042**. Direction (circularity ↑ toward SecB-rich) holds; magnitude/significance overstated. Code migrated faithfully (`figures/figure6/08_fig_macrophage_apoptosis_prolif_morphology.R`), not changed. **Action:** correct the statistic to 6/8 p=0.042, or trace and justify the 8/8 source.
- **Fig 5F — autocrine LR shift (RESOLVED; statistically unsupported).** Full investigation: `_repro_refactor/reports/FIG5F_INVESTIGATION.md`. The generator **is** in the repo (`figures/figure5/04_atlas_seca_secb_autocrine_shift.py`) and reproduces the panel **exactly**. The issue is the threshold: Fig 5F filters on **`lrscore > 0.5`** — a raw expression-magnitude floor that passes **67.7% of all interactions** — **not** a significance/specificity test. (The manuscript's "iarscore > 0.5" is a mis-transcription of `lrscore`; `iarscore` exists nowhere in the pipeline.) At that magnitude level the asymmetry is real and the category deltas reproduce exactly (DAMP/TLR +12, ECM −10); the named DAMP/TLR pairs (HMGB1→TLR4/TLR2, S100A8/9→TLR4) clear lrscore>0.5 in SecB and are absent in SecA. **But under proper significance** (`magnitude_rank ≤ 0.05`), the shift collapses — zero TLR autocrine pairs are significant in either pole, and the only robust ECM autocrine signal (collagen→CD44/CD93/ITGA5) is shared by both. **Action (author):** either reframe as a non-inferential expression-magnitude observation (and fix iarscore→lrscore), or soften/drop the "DAMP replaces ECM autocrine signaling" mechanistic narrative. The analysis is **retained** in the refactor as-is.

### Clarifications / corrections
- **SF7 caption (3735 → 375).** Caption states "n = 3735 clones across 248 samples"; the actual informative-clone count is **375** (factor-of-10 digit-insertion typo). All other SF7/SF4C numbers verified exact. **Action:** correct "3735" → "375".
- **BayesPrism multivariate covariate set.** The exported "full multivariate" model defaults to including platinum sensitivity, which dominates and renders polarization non-significant. The manuscript-stated model (epithelial fraction + stage + age), re-fit here, **retains significance** (OS p=0.039, PFS p=0.041). **Action:** clarify in Methods that the multivariate adjustment is epi fraction + stage + age (not platinum sensitivity).

### Reproducibility-hardening (applied / documented in the refactor)
- **Cohort drift handling.** The live cache had 2 extra FTE whole-tissue samples (`FT1-1`, `EAOC-1-FTE`) added post-preprint, inflating FT/WT/TMA n. The refactor **pins the published cohort** via `CFG$cohort$whole_tissue` (8 WT; FTE whole-tissue excluded, 15 FTE TMA cores retained), resolving 4A/4B/4E n inflation and the SF11 assertion. See `_repro_refactor/reports/COHORT_DRIFT_NOTE.md`.
- **Non-determinism now seeded.** `set.seed(CFG$seed)` added to Xenium cell-labeling `spatial/03_annotation_polarization/04_reclassification_polarization.R` and `05_clean_split_rctd.R` (source of TMA n/p and per-sample GAM consistency-count drift) and to `atlas/07_deconvolution_survival/08_consensusov_score.R` (RF flips 2 samples IMR↔MES). Python steps seed from config `SEED`.
- **HR reciprocal convention.** Scripts compute polarization as SecA-high (HR<1, protective); the manuscript reports the SecB-high reciprocal. 5-yr-clipped HRs invert exactly to the published values. The convention is documented so the mapping is explicit (`figures/figure7/`, `atlas/07_deconvolution_survival/05_signature_survival.py`).
- Cosmetic: scFEA flux-label-map bug fixed (`atlas/04_functional/01_epitype_functional_characterization.py`); `torch` dependency needed only for the scFEA step; "Transitioning" → "Intermediate" naming standardized across the polarization display vocabulary.

## Environment
- **Python (`epitype-py`):** 3.12.3, scanpy 1.11.4, anndata 0.12.2, numpy 1.26.4, numba 0.59 (+ lifelines, decoupler, liana).
- **R 4.5.2 stack:** SingleR 2.10.0, UCell 2.12.0, mgcv 1.9-4, spdep 1.4-2, spatstat 3.6-1, survival 3.8-6, BayesPrism 2.2.3, consensusOV 1.30.0, copykat 1.1.0, SpatialFeatureExperiment 1.10.1, Voyager 1.10.0 (see `docs/r_key_package_versions.csv`).
