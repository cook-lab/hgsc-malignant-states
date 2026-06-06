# Migration Summary — Phase 4

Finalizes the Phase-4 migration of the HGSC epitype codebase into the refactored monorepo. Migration record: the manifests in `_repro_refactor/content_map/migration_shards/`. Canonical migrate list: `_repro_refactor/reports/CANONICAL_SET.md`. Compiled 2026-06-05.

## Scripts migrated per area

| Area | Migrated scripts | Notes |
|---|---|---|
| `atlas/` (Python/R backend) | **41** | 01_preprocess_qc 5, 02_annotation 12, 03_epithelial_nmf 5, 04_functional 1, 05_cnv 4, 06_cellcomm 3, 07_deconvolution_survival 9, 08_xenium_reference 2 |
| `spatial/` (R/SFE backend) | **50** | 00_setup 1, 01_build_sfe 2, 02_qc 4, 03_annotation_polarization 5, 04_neighborhood 11, 05_gradients_gams 5, 06_spatial_stats 1, 07_morphometry 10, 08_clinical_survival 8, 09_external_validation 3 |
| `figures/` (per-panel generators + prep) | **69** | _prep 3, figure1 6, figure2 5, figure3 7, figure4 10, figure5 6, figure6 9, figure7 3, supplementary 20 |
| `tables/` (Supp Data generators) | **7** | Supp Data 1-7 (4 committed-as-is + 3 newly authored stubs) |
| **Total** | **167** | + `config/` (config.py, config.R, config.yml) and `shared/signatures.yml` |

Counts exclude `._*` AppleDouble sidecars and `__pycache__`. The `figures/` count includes shared prep helpers and the supplementary-side cache copies (`00_extract_atlas_obs.py`, `00b_extract_integration_umaps.py`), which is why it exceeds the ~50 unique per-panel generators in CANONICAL_SET (c).

## Conventions applied across the migration
- Hardcoded `/Volumes/CookLab/Sarah/...` paths → central config (`config/config.py`, `config/config.R`, object keys + `path()`/`cfg_path()`).
- Header docstrings (purpose / INPUTS / OUTPUTS / PANELS / runtime tier) on every script.
- RNG seeded from config `SEED` / `set.seed(CFG$seed)` (distinct local null seeds preserved where exact results depend on them, e.g. `29c` set.seed(29), `40b` bootstrap 20260508).
- SecA/SecB signatures loaded from `shared/signatures.yml` (canonical noBCAM 7-gene set), replacing divergent inline lists.
- "Transitioning" → "Intermediate" for the polarization **display** vocabulary only; validated level-2 annotation label strings preserved verbatim.
- Notebooks (`.ipynb`) converted to `.py` preserving logic exactly.

## Completeness check

### CANONICAL_SET scripts NOT migrated (with reason)
- **scVI→scANVI integration** (CANONICAL (a) step 4): the **original cluster scripts have now been located and migrated** (superseding the earlier interim reconstruction). The authoritative `atlas_00_aggregate` … `atlas_05_process` chain is migrated into `atlas/01_preprocess_qc/` as steps `02_aggregate.py` → `03_preprocess_hvg.py` → `04_scvi.py` → `05_cellassign.py` → `06_scanvi.py` → `07_process.py` (plus `04b_harmony_comparison.py`, an optional comparison from `atlas_02_harmony`). These are the real path-centralised originals (analytical parameters unchanged), not reconstructions; the interim `03_integration_scvi_scanvi.py` was removed. The original `atlas_04a` CellAssign step uses `shared/cellassign_markers.csv` (53-gene × 11-type), which differs from the 81-gene × 16-type `cellassign_markers_v3.csv` used by the later deposited integration — flagged in `atlas/01_preprocess_qc/README.md` along with the QC/Scrublet-placement and object-naming reconciliation points. The integrated atlas remains deposited as a convenience entry object; the GPU steps are optional to re-run. Superseded originals not migrated: `atlas_04_annotation.py` (CellAssign+scANVI monolith, replaced by the 04a/04b split) and `atlas_integration.py` / `atlas_integration2.py` (older all-in-one versions).
- **`42a..42f` FTE UCell/GAM** (CANONICAL (b) line 78, "pending manual resolution"): retained as backend only *if* the SF10C FTE gallery is reproduced; left unresolved at synthesis and **not migrated**. SF10C currently runs in TMA-cores mode (documented low-confidence flag).
- **Probe-QC heavy Rmd** (referenced by `spatial/02_qc/05_probe_qc.R`): the compute-heavy notebook lived in the excluded `scripts/sandbox/` and was **not migrated**; `05_probe_qc.R` invokes it via the `PROBE_QC_RMD` env var. This is the one convention partially deferred.
- **`19_cnv_*.sh` SLURM/batch wrappers** (drivers for `atlas/05_cnv/02_cnv_copykat.R`): documented batch wrappers, intentionally not migrated as numbered pipeline scripts.
- **To-author tidy-reshape stubs** for Supp Data 1/2/3: CANONICAL flagged these as "to author in refactor" (no committed exact-schema generator). They are now present as newly authored committed stubs (`tables/05`, `06`, `07`), so this gap is **closed**.

All other CANONICAL_SET entries — atlas backend (a), Xenium backend (b) including `13/13b/13c/13d/13f/13g/13m`, per-figure generators (c), and Supp-data generators (d) — have a migrated destination in the manifests and a corresponding file on disk.

### Panels without an in-repo generator (external / known gaps, not migration failures)
- **Fig 3A/B/D/E/F/G** — external organoid data (TRUST-EXISTING per author decision; scripts present but depend on `organoids_root`).
- **Fig 3C** — external flow cytometry; no in-repo generator (documented).
- **Fig 5H** — cell-type map absent from the canonical 5I gene-panel generator (not invented).

### Repo stage dirs lacking a README
**None.** All numbered pipeline stage directories have a `README.md`:
- `atlas/01_preprocess_qc` … `atlas/08_xenium_reference` (8/8)
- `spatial/00_setup` … `spatial/09_external_validation` (10/10)
- `figures/_prep`, `figures/figure1` … `figures/figure7`, `figures/supplementary` (9/9)
- `tables/` (present)

Repo root `README.md` and `data/README.md` are present. Top-level *area* dirs (`atlas/`, `spatial/`, `figures/`) and infra dirs (`config/`, `shared/`, `output/`, `docs/`) do not have an area-level README, but every pipeline **stage** dir does, so the stage-README check passes.

## Deliverables produced this phase
- `figures/README.md` — full panel→script map (Fig 1A..7G, SF1..SF14, Supp Data 1-7) with external/known-gap flags.
- `docs/REPRODUCIBILITY.md` — audit verdict + known-issues list (pointer to the full report).
- `docs/MIGRATION_SUMMARY.md` — this file.
