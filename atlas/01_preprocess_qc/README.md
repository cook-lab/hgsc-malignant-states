# atlas/01_preprocess_qc — preprocessing, QC, integration (raw → integrated atlas)

The **full preprocessing chain from raw per-study data through the integrated atlas object**
that is the entry point for all downstream analysis. The integration sub-sequence (steps 03–07)
consists of the **official cluster scripts**, path-centralised into this repo. The
scVI/scANVI/MrVI/SysVI/Harmony steps are GPU-heavy and **optional to re-run**: the integrated
atlas object is deposited as a config entry object, so downstream work can start from it.

## Ordered scripts

| # | Script | Official source | Inputs | Outputs |
|---|--------|-----------------|--------|---------|
| 01 | `01_dataprep.py` | (repo) | per-study raw matrices (`DATA_ROOT/2026_final_atlas/raw/<study>/`), `atlas_metadata.xlsx` | harmonised per-study h5ad → `processed/<study>.h5ad` |
| 02 | `02_concat_qc_doublets.py` | (repo) | 13 per-study h5ad | `processed/atlas_concatenated_filtered.h5ad` (per-cell QC: total_counts≥500, n_genes≥300, sample≥500 cells; **Scrublet < 0.3**) + QC violins. **This IS the QC step.** |
| 03 | `03_preprocess_hvg.py` | `01_preprocess.py` | `atlas_concatenated_filtered.h5ad` | `…/integration/anndata/preprocessed.h5ad` (HVG 4000 seurat_v3 batch=sample_id, subset; X=lognorm, `layers["counts"]`=raw) |
| 04 | `04_cellassign.py` | `02_cellassign.py` | `atlas_concatenated_filtered.h5ad`, `shared/cellassign_markers_v3.csv` (81×16) | `…/integration/cellassign/{predictions,probabilities}.csv` + QC plots |
| 05a | `05a_integrate_scvi.py` | `03a_integrate_scvi.py` | `preprocessed.h5ad`, predictions | `…/integration/models/scvi/`, `…/embeddings/scvi/embedding.npz`. **GPU; comparison + base for scANVI** |
| 05b | `05b_integrate_scanvi.py` | `03b_integrate_scanvi.py` | `preprocessed.h5ad`, predictions, scvi model | `…/integration/models/scanvi/`, `…/embeddings/scanvi/embedding.npz` (obsm `X_scanvi`). **GPU; ★ SELECTED METHOD** |
| 05c | `05c_integrate_mrvi.py` | `03c_integrate_mrvi.py` | `preprocessed.h5ad`, predictions | `…/embeddings/mrvi/embedding.npz` (u + z). **GPU; comparison** |
| 05d | `05d_integrate_sysvi.py` | `03d_integrate_sysvi.py` | `preprocessed.h5ad`, predictions | `…/embeddings/sysvi/embedding.npz`. **GPU; comparison** |
| 05e | `05e_integrate_harmony.py` | `03e_integrate_harmony.py` | `preprocessed.h5ad`, predictions | `…/embeddings/harmony/embedding.npz`. **CPU; comparison** |
| 06 | `06_benchmark.py` | `04_benchmark.py` | `preprocessed.h5ad`, predictions, all `embeddings/<m>/embedding.npz` | `…/integration/benchmark/metrics_{raw,scaled}.csv` + summary plot (scib-metrics; basis for selecting scANVI) |
| 07 | `07_finalize.py` | `05_finalize.py` | `atlas_concatenated_filtered.h5ad`, predictions, `embeddings/<method>/embedding.npz` | `…/integration/anndata/integrated_<method>.h5ad` (default `--method scanvi` → `integrated_scanvi.h5ad`, obsm `X_scanvi`, normalize/log1p base2, UMAP min_dist=0.3, Leiden res=0.2) |
| 08 | `08_refilter_umap.py` | (repo; was `03b`) | `integrated_scanvi.h5ad` (step 07 / deposited 20260213 run) | `obj("atlas_scanvi")` = `hgsc_atlas_scanvi.h5ad` (post-integration Scrublet<0.25, neighbours on `X_scanvi`, UMAP min_dist=0.2) + QC UMAPs |
| 09 | `09_atlas_qc_review.py` | (repo; was `04`) | `obj("atlas_scanvi")` | per-study QC summary CSVs + QC histograms/UMAPs |
| 10 | `10_umap_suite.py` | (repo; was `05`) | `obj("atlas_scanvi")` | 9 publication metadata UMAP SVGs |

The integration sub-sequence is **02 → 03 → 04 → {05a–05e} → 06 → 07**. **scANVI (05b) is the
SELECTED method** carried into the manuscript; scVI/MrVI/SysVI/Harmony + the `06_benchmark`
table are the **method comparison** (SF2 / Methods). All five integration methods read the same
`preprocessed.h5ad` and CellAssign `predictions.csv`. Integration-chain intermediates
(`preprocessed.h5ad`, CellAssign CSVs, per-method `models/` + `embeddings/`, `benchmark/`,
`integrated_<method>.h5ad`) are **regenerated** under `output_root/01_preprocess_qc/integration/`
— not deposited primary inputs. The QC/post-integration steps (`02`, `08`–`10`) write under
`output_root/01_preprocess_qc/`.

## These are the official scripts (optional GPU re-run)

Steps 03–07 are the **official cluster scripts** (`01_preprocess.py` … `05_finalize.py`),
migrated verbatim with only the hardcoded cluster paths
(`/home/snersesi/projects/def-dcook/active/hgsc_atlas/...`) replaced by central config.
**Analytical parameters were not altered.** Re-running the scVI/scANVI/MrVI/SysVI/CellAssign
steps needs a GPU and is computationally expensive; it is **not required** to reproduce
downstream results because the integrated atlas object is **deposited as a config entry object**
(`obj("atlas_scanvi")`). The full code (all five methods + benchmark) is included so every step
from raw data is open to scrutiny and can be independently re-executed. (Minor cosmetic
cleanups: `07_finalize` had an accidental duplicate `sc.pp.neighbors` call — idempotent,
collapsed to one — and corrupted emoji in some print strings, cleaned; no logic changed.)

## Resolved provenance (formerly open reconciliation flags)

The earlier reconciliation flags are now **RESOLVED** by the official scripts:

- **QC / Scrublet placement — RESOLVED.** The official `01_preprocess.py` and `02_cellassign.py`
  both read `atlas_concatenated_filtered.h5ad` directly, i.e. the output of
  `02_concat_qc_doublets.py`. So per-cell QC + per-sample **Scrublet < 0.3** (the Methods'
  2,398,571 → 2,326,532 filter) IS the QC step of the integration chain — it sits between
  concatenation (02) and HVG/integration (03+). No missing step.
- **Object-naming / pipeline-version bridge — RESOLVED.** `07_finalize.py --method scanvi` writes
  `integrated_scanvi.h5ad` with obsm key **`X_scanvi` (lowercase)** — exactly what
  `08_refilter_umap.py` reads (`sc.pp.neighbors(..., use_rep="X_scanvi")`). Step 08 then applies
  the stricter Scrublet < 0.25 + UMAP and writes the deposited `hgsc_atlas_scanvi.h5ad`
  (`obj("atlas_scanvi")`). The chain is continuous: 07 (`integrated_scanvi.h5ad`) → 08
  (`hgsc_atlas_scanvi.h5ad`).
- **CellAssign marker matrix — RESOLVED.** The official `02_cellassign.py` uses
  `shared/cellassign_markers_v3.csv` (**81 genes × 16 cell types**, identical to the staged
  official copy). This is the single canonical marker matrix; the older 53×11
  `cellassign_markers.csv` was an earlier vintage and has been removed.
- **Integration method — RESOLVED.** Five methods are run and benchmarked (scib-metrics, step 06);
  **scANVI is selected** for the manuscript (`07_finalize` default `--method scanvi`).

## Conventions applied

Central config for all paths: the chain's raw input is
`path("data_root", "2026_final_atlas", "processed", "atlas_concatenated_filtered.h5ad")` (the
step-02 output); markers → `shared/cellassign_markers_v3.csv`; all regenerated
intermediates/models/embeddings → `output_root/01_preprocess_qc/integration/`. Header docstrings
(purpose / INPUTS / OUTPUTS / runtime tier / manuscript role) added to each migrated script.
Analytical parameters of the official scripts preserved EXACTLY. `matplotlib` forced to the Agg
backend in the migrated scripts (headless cluster QC plotting).

## Figures supported

QC/overview substrate only. SF1 (per-study QC violins), SF2 (integration benchmark + UMAPs),
Fig 1A (atlas UMAP) and SF3 (metadata UMAPs) are rendered by per-panel scripts in `figures/`
from `obj("atlas_scanvi")` / `obj("atlas_final")`.
