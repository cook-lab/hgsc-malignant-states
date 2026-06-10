# figure_icon7_bevacizumab — PLACEHOLDER (panel letter TBD)

ICON7 / GSE140082 bevacizumab external validation, slated as an **expanded
Figure 7** (sibling to Fig 7E–G TCGA-OV). The final manuscript panel is **not
yet assigned**, so this directory uses a placeholder name — rename it to
`figureN/` (or `supplementary/`) once decided.

## The result (canonical)

Baseline SecB–SecA polarization is **prognostic under chemotherapy**
(chemo-arm continuous OS HR 1.38, p=0.038) and **abolished by bevacizumab**
(bev arm flat) — a per-arm slope reversal demonstrating treatment-modifier
interaction; with a dual-timescale time-varying decomposition. Full narrative,
figures, and caveats: the source module
`2026_final_xenium_analysis/davids side quests/ICON7/report.html`.

## Scripts

| Order | Script | Panel | Inputs → Outputs |
|---|---|---|---|
| 01 | `01_icon7_bev_slope_reversal.R` | TBD | `cfg_obj("icon7_cohort")` (per-patient `cohort_filtered.tsv`) → `figures_dir/figure_icon7_bevacizumab/icon7_per_arm_slope_forest.{svg,pdf}` (+ `_data.csv`); per-arm Cox computed in-script (mirrors `figure7/03_tcga_km_forest.R`) |
| S1 | `S1_secb_bev_confound_robustness.R` | none (supplementary) | `cfg_obj("icon7_cohort")` + `cfg_obj("icon7_expr")` → `…/supplementary_novelty/{A1_*,A2_*}.tsv` + `suppl_signature_battery_forest.{pdf,svg}` |

`01_…` reproduces the single headline panel (the per-arm slope reversal as an
HR-per-1-SD forest). The other narrative panels (4-group KM extremes,
median-rescue bars, time-varying forest, composite) are produced by the source
backend `…/ICON7/scripts/04_figures.R` (Section A) and can be ported here as
`02_…`, `03_…` when the figure layout is finalized.

`S1_…` is **supplementary robustness** — it produces no main panel and does not
change `01`. It stress-tests two reviewer objections (is this just the known
clinical high-risk → bev-benefit effect? is it just a generic hypoxia signature?),
both answered No with the trial data. It is not auto-run with the figures; invoke
it explicitly. Part B (the hypoxia/angiogenesis battery) needs `UCell` + `msigdbr`
and is skipped (with a flag) if they're absent. The detailed narrative, literature
context, and positioning notes (incl. citing/distinguishing Kommoss 2017,
PMID 28159814, on the same dataset) are kept in the authors' internal manuscript notes.

## Status / TODO before this is final

- [ ] Assign the manuscript panel letter(s); rename this dir + update headers and `figures/README.md`.
- [ ] Decide which panels are main vs supplementary (slope reversal, KM extremes, time-varying, composite).
- [ ] Migrate the ICON7 **backend** (the `01_prepare` → `03_sensitivity` stage scripts) into the repo as a numbered stage and have this figure read its deposited cache (see `…/ICON7/MIGRATION.md`). Currently `icon7_cohort` points at the source module's `data/processed/`.
- [ ] Confirm `config.yml` `objects: icon7_cohort` resolves under the deposited `DATA_ROOT` bundle.

## Config

Reads `cfg_obj("icon7_cohort")` and (optionally) `cfg_obj("icon7_results_dir")`,
added to `config/config.yml`. No hardcoded paths. `S1_…` additionally reads
`cfg_obj("icon7_expr")` (the gene-symbol × sample matrix, for Part B signature scoring).
