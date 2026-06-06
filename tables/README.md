# tables/ — Supplemental Data table generators

Generators for the manuscript's Supplemental Data items (Supp Data 1-7). Each
script resolves paths/objects through the central config (`config/config.py` for
Python, `config/config.R` for R) and writes its output under
`${OUTPUT_ROOT}/supplemental/`. SecA/SecB signatures come from
`shared/signatures.yml` (noBCAM 7-gene sets); the epithelial "Transitioning"
label is standardized to "Intermediate" throughout.

Run any generator from the repo root, e.g.:
```
python tables/01_supp_data_4_nmf_loadings.py            # or: --figures-only
Rscript tables/02_supp_data_5_xenium_gene_panel.R
```

## Ordered scripts (inputs -> outputs)

| # | Script | Input (config-resolved) | Output | Supp Data |
|---|--------|-------------------------|--------|-----------|
| 01 | `01_supp_data_4_nmf_loadings.py` | `atlas_celltype_dir/hgsc_atlas_final_epithelial.h5ad`; `shared/signatures.yml` | `output/11d_epithelial_nmf/11d_nmf_{usage,loadings,...}.csv`; `output/supplemental/supplementary_table_NMF_factor_genes.csv` | SD4 |
| 02 | `02_supp_data_5_xenium_gene_panel.R` | `output/05_probe_qc/{gene_exclusion_decisions,probe_qc_full}.csv` | `output/supplemental/Supplemental_Table_5_Xenium_Gene_Panel.csv` | SD5 |
| 03 | `03_supp_data_6_pathway_gene_sets.R` | (inlined 37-module definitions) | `output/9b_scoring/pathway_gene_sets_v2.csv`; `output/supplemental/Supplemental_Table_6_UCell_pathway_gene_sets.csv` | SD6 |
| 04 | `04_supp_data_7_autocrine_lr.py` | `output/17_cellcomm_nmf/tables/17b_liana_global.csv` | `output/17c_secA_secB_communication_nmf/tables/*.csv`; `output/supplemental/Supplemental_Table_7_autocrine_LR_pairs.csv` | SD7 |
| 05 | `05_supp_data_1_cell_metadata.py` | `atlas_final` (obs) | `output/supplemental/T1_atlas_metadata.csv` | SD1 |
| 06 | `06_supp_data_2_level1_de.py` | `atlas_final` (X + `celltype_level1`) | `output/supplemental/T2_celltypelevel1_markers.csv` | SD2 |
| 07 | `07_supp_data_3_ciliated_vs_secretory_de.py` | `atlas_celltype_dir/hgsc_atlas_final_epithelial.h5ad` | `output/supplemental/T3_epithelial_markers.csv` | SD3 |

## Figures supported
These are data tables, not figure panels. The shared caches they produce/consume
also underpin: SD4's `11d_nmf_usage.csv` -> Fig 1G/1H, SF5/SF6; SD6's pathway
modules -> Fig 5A/6B/6G and SF12-14; SD7's `17b_liana_global.csv` -> Fig 5F.

## Notes
- Filename typos fixed vs the canonical outputs: `T2_*.csv.csv` -> `T2_*.csv`;
  `T3_*markeres.csv` -> `T3_*markers.csv`. SD5/SD6/SD7 filenames are pinned with
  their `_5/_6/_7` suffixes.
- SD1 must run after the canonical Zhang-2022 `treatment_status` fix (atlas step
  12a) is baked into `hgsc_atlas_final.h5ad`.
- SD1/SD2/SD3 are newly authored committed export stubs (the canonical analysis
  exported these interactively); SD4/SD5/SD6/SD7 are migrated from committed
  canonical generators with logic preserved.
- SD6 owns only the gene-set *definitions* + tidy export; the heavy per-SFE UCell
  scoring from canonical 9b lives in the spatial pipeline shard.
