# atlas/06_cellcomm

LIANA+ cell-cell communication at NMF epitype resolution, and the SecA-vs-SecB
differential-communication / mechanism analyses.

## Ordered scripts (inputs → outputs)

| # | Script | Inputs | Key outputs |
|---|--------|--------|-------------|
| 01 | `01_cellcomm_nmf.py` | `atlas_final` h5ad (obs + `layers/counts`); `03_epithelial_nmf/celltype_nmf_mapping.csv` | `output/06_cellcomm/tables/17b_liana_global.csv` (+ top50, interaction counts, subsampling); LIANA figures |
| 02 | `02_secA_secB_communication_nmf.py` | `tables/17b_liana_global.csv` | `17c_secA_secB_communication_nmf/tables/17c_*.csv`; differential L-R + TME-partner figures; HTML |
| 03 | `03_secA_secB_mechanisms_nmf.py` | `tables/17b_liana_global.csv`; `fig_secretory_polarization/data/panel_i_deg_results.csv` | `17d_secA_secB_mechanisms_nmf/tables/17d_*.csv`; pathway/level2/DEG-crossref figures; HTML |

## Figures / panels supported

- `17b_liana_global.csv` → **Fig 5F** (autocrine shift), **Supp Data 7** (autocrine L-R pairs)
- 02/03 differential + mechanism tables → supporting Fig 5F / Supp Data 7

## Conventions applied

Central config paths; header docstrings; subsample seeded from config
(`RANDOM_SEED = SEED`; the LIANA permutation seed is a fixed analytical parameter,
preserved as published); epithelial label **Transitioning → Intermediate**
throughout the NMF/level1 maps.

### Migration note (03)
The original `17d` was a runtime wrapper that string-patched and `exec()`-ed the
(Leiden-based, non-canonical) `16c` script. Script 03 bakes those exact patches
in as a standalone script — NMF label schema, the `17b` LIANA input, the `17d`
output prefixes — so there is no fragile runtime exec and no dependency on an
un-migrated `16c`. Analytical logic is preserved verbatim.
