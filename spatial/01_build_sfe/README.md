# spatial/01_build_sfe — build raw SpatialFeatureExperiment objects

Constructs the raw, unfiltered SFEs from Xenium data. All cells kept; only raw
counts stored (normalization happens in `spatial/02_qc/`).

## Ordered scripts
1. `01_build_sfe_tma.R` — load fresh TMA_1 + TMA_2, transfer core/patient
   metadata from the previously processed SFE (cell-ID match + spatial-NN
   fallback), shift TMA_2 coordinates, merge.
   - IN: raw TMA slides + `previously processed/.../se.rds`
   - OUT: `<sfe_dir>/sfe_tma`
2. `02_build_sfe_whole_tissue.R` — build one SFE per published whole tissue.
   OTB_2457_2384 is DBSCAN-split into OTB_2457 + OTB_2384.
   - IN: raw whole-tissue runs under `data/xenium/whole_tissue/`
   - OUT: `<sfe_dir>/sfe_<sample>` for the 8 published whole tissues

## Cohort PIN
`02_build_sfe_whole_tissue.R` enumerates exactly the published 8 whole tissues
(`CFG$cohort$whole_tissue`) and **excludes** FT1-1 and EAOC-1-FTE
(`CFG$cohort$fte_exclude_wt`) from the whole-tissue arm. FTE TMA cores remain in
`sfe_tma`.

## Figures supported
Backend for all TMA + whole-tissue panels (Fig 4–7, SF10–SF14). Upstream of
`sfe_tma_filtered` (the canonical TMA entry-point, built in `spatial/02_qc/`).
