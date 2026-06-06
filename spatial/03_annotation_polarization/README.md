# spatial/03_annotation_polarization — cell typing + secretory polarization

The canonical **noBCAM** annotation/polarization chain. Produces the final
`cell_label` (SecA / Intermediate / SecB / Ciliated + 15 non-epithelial types)
written into every SFE, used by all Fig 4–7 spatial panels.

## Ordered scripts
1. `01_annotation_singler.R` — SingleR 16-type annotation against the matched
   scRNA-seq reference (`obj("xenium_ref")`), excluding probe-QC-failed genes;
   TMA validation vs prior labels.
   - IN: normalized SFEs + xenium reference + `genes_exclude_singler.txt`
   - OUT: `singler_label` in SFEs; `06_annotation/*.csv`
2. `02_adaptive_secretory_noBCAM.R` — initial SecA/Intermediate/SecB subtype
   from the SecB/SecA logcounts ratio (thresholds 1.0/2.0) using the shared
   noBCAM signatures.
   - IN: annotated SFEs
   - OUT: `secretory_subtype`/`cell_label`/scores in SFEs;
     `06b_adaptive_secretory_noBCAM/*.csv`
3. `03_ucell_scoring_noBCAM.R` — UCell SecA/SecB scoring on the Xenium∩atlas
   shared gene space; writes `polarization_UCell` back to SFEs.
   - IN: annotated SFEs + atlas `18_ucell_atlas/` exports (path via
     `ATLAS_UCELL_DIR`)
   - OUT: `06d_annotation_noBCAM/{xenium,atlas}_ucell_scores.csv`; SFE scores
4. `04_reclassification_polarization.R` — derive atlas-calibrated polarization
   thresholds (SecA p75 / SecB p25), classify, write `cell_label`, freeze
   thresholds.
   - IN: 06d UCell score CSVs + SFEs
   - OUT: `06f_reclassification_polarization/threshold_summary.csv` (FROZEN),
     `reclassified_xenium_scores.csv` (override consumed by most Fig 4–6 panels)
5. `05_clean_split_rctd.R` — RCTD doublet-mode + SPLIT purification; final
   `cell_label` writer. Non-secretory cells take RCTD first_type; RCTD-confirmed
   secretory cells get SecA/Intermediate/SecB from the frozen 06f thresholds.
   - IN: SFEs + xenium reference + frozen 06f thresholds
   - OUT: `06g_clean_split/{name}_rctd_results.rds` + purified counts/labels;
     final `cell_label` + `purified` assay in SFEs

> Renumbered from source names `06_annotation`, `06b_adaptive_secretory_production_noBCAM`,
> `06d_annotation_noBCAM`, `06f_reclassification_polarization`, `06g_clean_split`.

## Conventions applied (chain-specific)
- **Determinism fix:** `04` seeds the atlas sub-sample and `05` calls
  `set.seed(CFG$seed)` before each RCTD/SPLIT run (the audit's non-determinism
  finding) so cell labels reproduce.
- **noBCAM signatures:** SecA/SecB come from `shared/signatures.yml` via
  `00_setup.R` (no inlined gene lists).
- **Label rename:** "Transitioning" -> "Intermediate" throughout (label
  strings, palette keys, atlas NMF labels remapped on read, summary columns).
- **Cohort PIN:** RCTD/scoring run on sfe_tma + the 8 published whole tissues;
  writeback also covers sfe_tma_filtered. FTE whole-tissue samples excluded.

## Figures supported
Final `cell_label` + polarization feed Fig 4 (composition/ROI), Fig 5
(polarization gradients), Fig 6 (niche/macrophage), Fig 7 (clinical), and
SF10–SF14.
