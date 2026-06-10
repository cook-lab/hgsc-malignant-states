# figure_ccle_smac_mimetics — PLACEHOLDER (panel letter TBD)

CCLE / DepMap 24Q4 + PRISM Repurposing Harmonized 25Q2 side-quest, slated as a
**supplementary / mechanistic** panel (a pharmacological vulnerability of the
SecB epitype, sibling to the in vivo niche analyses). The final manuscript panel
is **not yet assigned**, so this directory uses a placeholder name — rename it to
`figureN/` (or move into `supplementary/`) once decided.

## The result (canonical)

**SecB-polarized HGSC cell lines are selectively sensitive to SMAC mimetics.**
In PRISM Repurposing Harmonized 25Q2, all three available SMAC mimetics have
**negative** Spearman correlations of viability AUC with SecB polarization
(lower AUC = more sensitive): **LCL-161** ρ = −1.00 (exact p = 0.017, n = 5),
**GDC-0152** ρ = −0.71 (exact p = 0.088, n = 7; broadest coverage), and
**birinapant** ρ = −0.50 (n = 5). Chemotherapy negative controls (carboplatin,
doxorubicin) do **not** show the pattern. By contrast, the pre-specified
single-gene **BIRC3 Chronos dependency is null / opposite** (ρ = +0.22
permissive; +0.37 strict) — interpreted as **IAP-family redundancy**
(cIAP1/cIAP2/XIAP) that pan-IAP SMAC mimetics bypass but a single-gene CRISPR
knockout does not reveal. The upstream, non-redundant TNF-receptor-signalosome
adapter **TRAF2 trends** in the predicted direction (ρ = −0.45, p = 0.08),
consistent with the SMAC-mimetic signal. Full narrative + caveats: the source
module `2026_final_xenium_analysis/davids side quests/ccle_depmap/report.html`
(and its `README.md`).

## Scripts

| Order | Script | Role | Inputs → Outputs |
|---|---|---|---|
| 01 | `01_score_hgsc_lines.R` | PREREQUISITE — regenerate per-line scores from raw | `cfg_obj("ccle_model_meta")` (Model.csv), `cfg_obj("ccle_expression")` (OmicsExpressionProteinCodingGenesTPMLogp1.csv) → `figures_dir/figure_ccle_smac_mimetics/hgsc_line_scores.tsv`, `hgsc_model_meta.tsv`, `01_polarization_ranking.pdf`, `01_polarization_ranking_UCell.pdf`. UCell + per-gene z-mean (z-mean primary); regenerates the deposited `ccle_line_scores` cache deterministically. |
| 02 | `02_prism_smac_mimetics.R` | **HEADLINE** — PRISM SMAC-mimetic AUC correlation | `cfg_obj("ccle_line_scores")`, `cfg_obj("ccle_prism_auc")` (REPURPOSINGAUCMatrix.csv), `cfg_obj("ccle_prism_conditions")` (REPURPOSINGCollapsedConditions.csv) → `smac_mimetic_correlations.tsv`, **`smac_combined_panel.pdf`** (the headline 3-SMAC facet figure), `smac_lead_scatter.pdf`, `scatter_<COMPOUND>.pdf` (e.g. `scatter_LCL_161.pdf`, `scatter_GDC_0152.pdf`, `scatter_BIRINAPANT.pdf`). |
| 03 | `03_crispr_targeted_nfkb.R` | SUPPORT — targeted NF-κB / BIRC3 / IAP Chronos | `cfg_obj("ccle_line_scores")`, `cfg_obj("ccle_crispr")` (CRISPRGeneEffect.csv) → `nfkb_panel_correlations.tsv`, `forest_nfkb_panel.pdf`, `scatter_BIRC3_polarization.pdf`, `scatter_BIRC3_polarization_strict.pdf`, `scatter_<GENE>_polarization.pdf` (top gene by p_zmean → `scatter_TRAF2_polarization.pdf`). 27-gene panel; permissive + strict Domcke-9; Fisher-z CIs. |

Run order: **01 → (02, 03)**. Script 01 is the prerequisite that produces the
per-line scores; 02 (headline) and 03 (support) both read those scores. In the
deposited bundle 02/03 read the canonical cache `cfg_obj("ccle_line_scores")`
directly, so they can run without re-running 01.

## Streamlining note (what was dropped vs the source module)

The source module `…/ccle_depmap/scripts/` had more files; only the canonical
path was migrated. Dropped:

- **`03_unbiased_dependency_scan.R`** — the genome-wide (~17.5k-gene) Spearman
  scan. The report itself flags it as **n = 16 noise**: the only q < 0.1 hit
  (DEFB121, a defensin) is leave-one-out-unstable and has never been essential
  anywhere in CCLE (Chronos range −0.44 to +0.70); the top hits are dominated by
  defensins/olfactory receptors that are not real essentialities. Statistical
  artifacts at this n, so not migrated.
- **`00_find_smac_compounds.py`** — a Python discovery helper that located the
  SMAC-mimetic compound names in PRISM. Its result is already **baked into
  script 02's `smac_targets` list** (LCL-161, birinapant, GDC-0152, etc.), so the
  helper is redundant.
- **`01b_score_diagnostics.R`, `04b_prism_diag.R`, `05_inspect_top_hits.R`** —
  diagnostic / inspection scripts (UCell-saturation diagnostics, exact-p
  diagnostics, top-hit inspection); not part of the canonical figure path.
- **Unused IC50 matrix load** — source `04_prism_smac_mimetics.R` read
  `REPURPOSINGLog2IC50Matrix.csv` into `ic50` but never used it. Removed from the
  migrated script 02 (so `ccle_prism_ic50` is **not** a config key).

All analytical logic that was kept (the exact-Spearman settings, BH adjustment,
27-gene panel, negative controls, Fisher-z CIs, plot geoms, `theme_classic(base_size = 9/8)`)
is preserved byte-for-byte from the source; signatures now load from
`shared/signatures.yml` (identical 7-gene noBCAM sets) and outputs/inputs route
through the repo config.

## Status / TODO before this is final

- [ ] Assign the manuscript panel letter(s); rename this dir + update script
      headers and `figures/README.md` (move out of "Pending / placeholder").
- [ ] Decide main vs supplementary (likely supplementary / mechanistic).
- [ ] Deposit under the `DATA_ROOT` bundle at full deposit time: the per-line
      scores cache (`hgsc_line_scores.tsv` → `cfg_obj("ccle_line_scores")`) and
      the raw matrices (Model.csv, the TPM matrix, CRISPRGeneEffect.csv, the two
      PRISM files) so `cfg_obj("ccle_*")` resolve. Confirm the six keys are added
      to `config.yml` (owned by the config-wiring agent).
- [ ] Verify script 01 regenerates the deposited `ccle_line_scores` identically
      (provenance check).

## Config (6 cfg keys read)

| Key | Source file | Used by |
|---|---|---|
| `ccle_line_scores` | `hgsc_line_scores.tsv` (deposited; regenerated by 01) | 02, 03 |
| `ccle_model_meta` | `Model.csv` | 01 |
| `ccle_expression` | `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 01 |
| `ccle_crispr` | `CRISPRGeneEffect.csv` | 03 |
| `ccle_prism_auc` | `REPURPOSINGAUCMatrix.csv` | 02 |
| `ccle_prism_conditions` | `REPURPOSINGCollapsedConditions.csv` | 02 |

No hardcoded paths; inputs via `cfg_obj(...)`, all outputs via
`cfg_path("figures_dir", "figure_ccle_smac_mimetics", ...)`.

## Data sources (public; documented for reproducibility)

DepMap 24Q4 Public (Figshare; downloaded 2026-06-04):

- `Model.csv` — https://ndownloader.figshare.com/files/51065297
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` — https://ndownloader.figshare.com/files/51065489
- `CRISPRGeneEffect.csv` — https://ndownloader.figshare.com/files/51064667

PRISM Repurposing Harmonized Secondary Screen 25Q2 (DepMap portal):

- `REPURPOSINGAUCMatrix.csv` — https://depmap.org/portal/download/api/download?file_name=downloads-by-canonical-id%2Fprocessed-files-for-ctd2-gdsc1-gdsc2-and-repurposing-secondary-4e65.3%2FREPURPOSINGAUCMatrix.csv
- `REPURPOSINGLog2ViabilityCollapsedConditions.csv` (saved as `REPURPOSINGCollapsedConditions.csv`) — same path with that filename

DepMap file manifest API used to look up the URLs:
`https://depmap.org/portal/download/api/downloads`

Derived files (`hgsc_line_scores.tsv`, `hgsc_model_meta.tsv`, all figures/tables)
are reproducible from the scripts above.
