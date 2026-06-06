# spatial/02_qc — QC, filtering, normalization

Computes QC metrics, filters cells, removes control probes, normalizes, runs
cross-platform probe QC, and produces the canonical filtered TMA entry-point.

## Ordered scripts
1. `03_qc_metrics.R` — per-cell QC metrics across all samples (no filtering).
   - IN: `<sfe_dir>/sfe_tma`, `<sfe_dir>/sfe_<wt>`
   - OUT: `03_04_qc/qc_metrics_all_samples.csv`, `qc_summary_per_sample.csv`
2. `04_qc_filter_normalize.R` — cell filters (counts>=10, genes>=5,
   neg<=5%, area 10–500 um^2), remove 64 control probes, add `logcounts` +
   `logcounts_area`.
   - IN: raw SFEs
   - OUT: SFEs updated in place; `03_04_qc/filtering_summary.csv`
3. `05_probe_qc.R` — cross-platform probe QC (Supp Data 5). Renders the probe-QC
   notebook (path via `PROBE_QC_RMD` env var) and verifies outputs.
   - IN: SFEs + xenium reference
   - OUT: `05_probe_qc/{probe_qc_full,flag_qc,gene_exclusion_decisions,celltype_coverage}.csv`,
     `genes_exclude_singler.txt`, `genes_monitor_singler.txt`
4. `06_filter_tma.R` — apply core-level QC exclusions; write the canonical TMA.
   - IN: `<sfe_dir>/sfe_tma`, `07_core_qc/core_qc_summary.csv`,
     `clinical_data_clean.csv`
   - OUT: `<sfe_dir>/sfe_tma_filtered`, `07_core_qc/excluded_cores_documentation.csv`

> `06_filter_tma.R` was `07b_filter_tma.R` in the source tree (renumbered for
> clean ordering within this stage).

## Cohort PIN
Scripts 03/04 operate over `sfe_names_all` (sfe_tma + the published 8 whole
tissues); the FTE whole-tissue samples are excluded.

## Figures supported
QC backend. `05_probe_qc.R` -> Supp Data 5. `06_filter_tma.R` produces
`sfe_tma_filtered`, the canonical TMA object feeding Fig 4–7 TMA panels,
SF10B, SF12.
