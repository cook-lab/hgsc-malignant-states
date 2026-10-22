# Figure 6: SecB-niche macrophages

SecB epithelium co-opts a hypoxic macrophage state into a pro-dissemination niche. Panels a–e use Visium HD of whole
tissue WT01. Panel f compares the atlas macrophage subtypes with Xenium SecB-niche macrophages. Panel g measures the
SecB neighbours of hypoxic versus resident-like macrophages in Xenium (whole tissues and TMA patients). Panel h is a
cross-platform gene matrix (Xenium, Visium HD, atlas), and panel i a Xenium ROI of WT02.

## How this folder runs

- Steps 01–11 write tables to `$OUTPUT_ROOT/fig6_singler/`. Step 13 draws **every** panel in one run, so every input
  of every panel must exist before it starts.
- `13_fig6_draw_panels.py` imports `fig6_roi_picker/roi_server.py` as a library: palettes, the Xenium ROI loader and the
  cell polygons. Run on its own, the same file is the interactive ROI picker that writes `fig6_rois.json`.
- The panel keys inside `13_fig6_draw_panels.py` are not the published letters. They map to the published panels (and
  to the file names in its `NAMES`) as follows:

  | Key in 13 | a | b | c | m | d | e | f | g | h | k |
  |---|---|---|---|---|---|---|---|---|---|---|
  | Panel | 6a | 6b | 6c | 6d | 6e | 6f | 6g | 6h | 6i, top row | 6i, bottom row |

  The bottom row of 6i is written as `Fig6j_ROI_SecB_CD55_IL32_MMP7_MCT4` (tiles `Fig6j1`–`Fig6j4`). It is panel i;
  there is no panel 6j.

## Requires

These steps live in other folders and must run first. Run them with the same `OUTPUT_ROOT` as this folder.

| Upstream step | Writes | Used by |
|---|---|---|
| Xenium backend: `spatial/01_build_sfe/` → `spatial/02_qc/04_qc_filter_normalize.R` → `spatial/03_annotation_polarization/01_annotation_singler.R` … `04_reclassification_polarization.R` → `06_extract_singler.R sfe_<sample>`. SingleR reads the probe-QC exclusion list `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt` (from `spatial/02_qc/as_run/explore_probe_qc.Rmd`, copied by hand; see [`spatial/02_qc/README.md`](../../spatial/02_qc/README.md)) | `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` (SingleR cell type, x, y, `polarization_UCell`, TMA core and patient) | 04–07, picker |
| Xenium raw-count export: `spatial/03_annotation_polarization/05_clean_split_rctd.R` → `spatial/05_count_export/as_run/split_recompute_stage1_polarization.R` → `spatial/05_count_export/as_run/02_deg_precompute.R sfe_<sample>` | `$OUTPUT_ROOT/fig6_viewer/degdata/<sample>_{counts.mtx.gz, cells.csv.gz, genes.txt}` (raw counts, 477 genes) for WT01–WT08, `tma` (TMA A), `COEUR_TMA_G`, `COEUR_TMA_H` (TMA B) | 04–07, picker |
| Final TMA core set: `spatial/02_qc/05b_core_qc.R` (once per array) → `spatial/02_qc/06b_final_tma_cores.py`; loader `figures/_prep/final_tma_cores.py` | `$OUTPUT_ROOT/07_core_qc/final_tma_{cores,patients,counts}.csv` | 06 (TMA arm of 6g) |
| Pathway gene sets: `spatial/04_pathway_scoring/04_pathway_scoring.R` (or `04b_pathway_gene_sets.R`, same definitions) | `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv`, copied by hand to `$OUTPUT_ROOT/fig6_singler/inputs/` | 05 |
| Nath 2021 relabelling: `figures/_prep/nath2021_fix.py` | `$OUTPUT_ROOT/_prep_caches/nath2021_cell_correction.parquet` | 08 |
| Fig 2g pseudobulk DE: `figures/figure2/04_fig2g_pseudobulk_volcano.py` | `$OUTPUT_ROOT/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv` | 11, 13 (6h atlas SecB row) |
| Atlas lineage subsets: `atlas/02_annotation/as_run/9a_atlas_celltypelevel1_v2.py` | deposited `$DATA_ROOT/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_{macrophage,epithelial}.h5ad` | 08 (macrophage, see Environment), 10 (epithelial) |
| Atlas epithelial states: `atlas/03_epithelial_nmf/as_run/fig_secretory_polarization_00_prepare_data.py` | deposited `$DATA_ROOT/2026_final_atlas/output/fig_secretory_polarization/data/meta.parquet` (`schema_nmf`) | 10 |
| Visium HD Space Ranger outputs, staged by hand (see Not in this folder) | `$OUTPUT_ROOT/fig6_singler/visium_local/{WT01B1, WT03B1, WT01B1_8um}/` | 01–03 |
| Raw Xenium cell boundaries of WT01 and WT02 (deposit) | `$DATA_ROOT/202609_final_hgsc_spatial/data/HGSC xenium /<run>/cell_boundaries.parquet` (the folder name ends with a space; runs `WT05_WT01` and `WT02`) | picker, 13 (6i polygons) |

## Run order

Paths are under `$OUTPUT_ROOT/fig6_singler/` unless written otherwise.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 01 | `01_fig6abcd_visium_hd_bins.py` | `visium_local/WT01B1/`, `visium_local/WT03B1/` (16 µm `filtered_feature_bc_matrix.h5`, `tissue_positions.parquet`) | `visium_n2_bins_WT01.csv.gz` (class, polarization and macrophage program per bin; CD47, SPP1, MIF, CD55, SLC16A3, ICAM1); `visium_n2_bins_WT03.csv.gz` (same for the second section; not read by any panel) | a–d |
| 01b | `01b_fig6c_visium_hd_bins_extra_genes.py` | `visium_n2_bins_WT01.csv.gz`; `visium_local/WT01B1/` (16 µm) | `visium_n2_bins_WT01.csv.gz`, in place, with log1p CP10k of ITGB6, ITGAV, CA9, CTSL and TIMP1 | c (CA9) |
| 02 | `02_fig6e_visium_hd_niche_volcano.py 16`, then `… 8` | `visium_local/WT01B1/` (16 µm) or `visium_local/WT01B1_8um/` | `visiumWT01_niche_vs_distal_volcano.csv`, `visiumWT01_niche_vs_distal_volcano_8um.csv` | e (16 µm); Table S10 (both) |
| 03 | `03_fig6h_visium_hd_candidates.py` | `visium_local/WT01B1/`, `visium_local/WT03B1/` | `visium_candidates_2sections.csv` | h (Visium rows, WT01) |
| 04 | `04_fig6fhi_xenium_niche_macrophages.R` | cell tables, count export (8 whole tissues) | `mac_niche_WT.csv.gz` (niche / control / other, per macrophage), `mac_counts_WT.csv`, `mac_signature_deg_WT.csv` (pseudobulk limma-voom) | f, h, i; Table S10 |
| – | `fig6_roi_picker/roi_server.py` (+ `index.html`), interactive, once | count export, cell tables, `mac_niche_WT.csv.gz`, `visium_n2_bins_WT01.csv.gz`, raw Xenium `cell_boundaries.parquet` | `fig6_rois.json` (appended on each save); polygon cache `$OUTPUT_ROOT/fig6_roi_picker/cache/` | ROI of a–d and i |
| 05 | `05_fig6fh_xenium_reciprocal_effects.R` | cell tables, count export, `mac_niche_WT.csv.gz`, `inputs/pathway_gene_sets_v2.csv` | `reciprocal_persample.csv`, `reciprocal_summary.csv` (comparisons `mac_niche_vs_ctrl`, `secb_contact`, `secb_vs_seca`) | f (Xenium row), h (Xenium rows) |
| 06 | `06_fig6g_macrophage_secB_neighbours.py` | cell tables and count export (WT01–WT08, `tma`, `COEUR_TMA_G`, `COEUR_TMA_H`); final TMA cores | `fig6_colocal_percell.csv.gz` (per macrophage), `fig6_colocal_units.csv` (one row per whole tissue or TMA patient) | g |
| 07 | `07_fig6h_xenium_pct_expressing.py` | cell tables, count export, `mac_niche_WT.csv.gz` | `xenium_pct_persample.csv`, `xenium_pct_summary.csv` | h (dot size, Xenium rows) |
| 08 | `08_fig6fh_atlas_macrophage_pseudobulk.py [macrophage h5ad]` | atlas macrophage object (309,111 cells), Nath 2021 cache | `atlas_pb/{counts_sum,lognorm_sum,nnz}.npy`, `atlas_pb/groups.csv`, `atlas_pb/genes.txt` (sample × level-2 type sums, 2,108 groups) | f, h |
| 08b | `08b_fig6f_atlas_subtype_dotmatrix.py` | `atlas_pb/` | `atlas_subtype_dotmatrix.csv` (mean log1p CP10k and % expressing of 15 genes in 6 macrophage / monocyte types) | f (dot plot) |
| 09 | `09_fig6h_atlas_hypoxic_vs_other.py` | `atlas_pb/` | `atlas_hypoxic_vs_other_paired.csv` | h (atlas macrophage row) |
| 10 | `10_fig6h_atlas_epithelial_pct.py` | atlas epithelial object (counts), `meta.parquet` (`schema_nmf`) | `atlas_epi_pct.csv` | h (dot size, atlas SecB row) |
| 11 | `11_fig6h_concordance.py` | `reciprocal_summary.csv`, `visium_candidates_2sections.csv`, `atlas_pb/`, `$OUTPUT_ROOT/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv` | `fig6_concordance_table.csv` (call per platform and a verdict per gene) | h |
| 13 | `13_fig6_draw_panels.py ["example 2 (700 um)"]` | the outputs of 01–11; `atlas_subtype_dotmatrix.csv`; `fig6_rois.json`; the Fig 2g DE table; WT02 cell polygons through `roi_server.py` | `$OUTPUT_ROOT/final_panels/figure6/`, one file per panel (`.pdf`, `.png`, `.svg`): `Fig6a_visium_section_SecA_SecB`, `Fig6b_visium_ROI_SecA_SecB`, `Fig6c_visium_ROI_CA9`, `Fig6d_visium_ROI_macrophage_niche`, `Fig6e_visium_niche_macrophage_volcano`, `Fig6f_atlas_hypoxic_macrophage_dotplot`, `Fig6g_hypoxic_macrophage_SecB_neighbours`, `Fig6h_crossplatform_matrix`, `Fig6i_ROI_macrophage_VEGFA_CTSL_CD44_ICAM1` (+ tiles `Fig6i1`–`Fig6i5`), `Fig6j_ROI_SecB_CD55_IL32_MMP7_MCT4` (+ tiles `Fig6j1`–`Fig6j4`); `fig6h_matrix_standardized.csv` | a–i |

`13` writes the panel files to `$OUTPUT_ROOT/final_panels/figure6/` and does not create that folder (run
`mkdir -p $OUTPUT_ROOT/final_panels/figure6` first). It draws the full page once, without saving it, to record each
panel's crop region, then draws every panel and tile alone on the same page grid. `fig6h_matrix_standardized.csv` holds
the plotted (standardized) values of the 6h matrix.

## Environment

| Setting | Step | Value in the final run |
|---|---|---|
| `OUTPUT_ROOT` | all | Output root (default `./output`): `$OUTPUT_ROOT/singler_recompute` (cell tables), `$OUTPUT_ROOT/fig6_viewer/degdata` (count export), `$OUTPUT_ROOT/fig6_singler` (this folder's tables and `visium_local/`), `$OUTPUT_ROOT/final_panels/figure6` (panel files of 13) |
| `DATA_ROOT` | 10, picker, 13 | Deposit root: the atlas epithelial object and `meta.parquet` (10); the raw Xenium cell boundaries (picker; 13 through `roi_server.py`) |
| first argument | 08 | Path to the atlas macrophage object. Without it, 08 reads `$OUTPUT_ROOT/atlas_local/hgsc_atlas_final_macrophage.h5ad`, the local copy of the deposited `$DATA_ROOT/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_macrophage.h5ad` that the final run used. 08 writes to `$OUTPUT_ROOT/fig6_singler/atlas_pb/` and reads the Nath 2021 cache from `$OUTPUT_ROOT/_prep_caches/` |
| `FINAL_TMA_CORES` | 06 | default, `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` |
| first argument | 02 | `16` (Fig 6e, Table S10), then `8` (Table S10) |
| first argument | 13 | ROI name in `fig6_rois.json`; default `example 2 (700 um)` (the 6i ROI). The Visium window of a–d is the **last** saved entry with `kind: visium` |
| `PORT` | picker | 8793 (default) |
| `LANG`, `LC_ALL` | upstream R steps that read `config.yml` | `en_US.UTF-8` |
| Interpreters | all | Python 3.12.9 (E2), R 4.4.3 (E1); see `env/README.md` |

Run `figures/_prep/nath2021_fix.py` before 08: if the Nath 2021 cache is missing, 08 builds it from the Nath cells of
the macrophage object only.

## Panel → script

| Panel | Drawn by | Numbers computed by |
|---|---|---|
| 6a | `13` (key `a`) | `01`: bin classes. Pure epithelial bin = ≥ 5 epithelial-marker UMIs and no macrophage UMIs; SecA-high / SecB-high = bottom / top tertile of z(SecB) − z(SecA) among pure epithelial bins (16 µm bins). The box is the Visium window from `fig6_rois.json` |
| 6b | `13` (key `b`) | `01`; the SecA and SecB territory contours are drawn in `13` |
| 6c | `13` (key `c`) | CA9 column of `visium_n2_bins_WT01.csv.gz` (`01b`); the 3 × 3-bin mean is taken in `13` |
| 6d | `13` (key `m`) | `01`: macrophage bin = ≥ 3 macrophage UMIs and ≥ epithelial UMIs; SecB-niche = ≥ 2 SecB-high bins within 3 bins (about 48 µm) |
| 6e | `13` (key `d`) | `02 16`: niche vs distal macrophage bins, log2 fold change of mean CP10k, bin-level Mann–Whitney with Benjamini–Hochberg FDR, carry-over rule (`mac_sourced`). The macrophage-content adjustment (median log2FC of 7 pan-macrophage genes) and the stromal class are applied in `13`; Table S10 applies the same rules |
| 6f | `13` (key `e`) | Dot plot: `08` → `08b` (`atlas_subtype_dotmatrix.csv`). Bottom row: `04` → `05` (`mac_niche_vs_ctrl`: per-section log2FC, two-sided Wilcoxon signed-rank over n = 8 whole tissues, ≥ 20 cells per group) |
| 6g | `13` (key `f`) | `06`: program = z(9 glycolytic genes) − z(4 resident genes), z-scored over all macrophages of each whole tissue or TMA array; per macrophage, SecB share of the secretory cells (SecA, Intermediate, SecB) within 50 µm (≥ 3 neighbours); tertiles within each unit; unit value = mean per tertile. A whole tissue or TMA core enters only with ≥ 30 macrophages, ≥ 50 secretory cells and ≥ 5 endothelial or pericyte cells; TMA cores are pooled per patient (WT n = 8, TMA n = 311 patients). The paired Wilcoxon tests are computed in `13` |
| 6h | `13` (key `g`) | Macrophage rows: Xenium `04` → `05` (`mac_niche_vs_ctrl`, n = 8) with % expressing from `07`; Visium `03` (WT01 column); atlas `08` → `09` (hypoxic vs other macrophages, paired per sample, n = 181). Middle row, Xenium SecB cells within 30 µm of a SecB-niche macrophage vs SecB cells with no macrophage within 30 µm: `05` (`secb_contact`, n = 8). SecB rows: Xenium `05` (`secb_vs_seca`) with `07`; Visium `03`; atlas Fig 2g DESeq2 (n = 261 tumours, `figures/figure2/04_fig2g_pseudobulk_volcano.py`) with % expressing from `10`. `11` supplies atlas macrophage values for genes outside `09`'s list. The gene list and groups are fixed in `13` (`BEATS`); effects are standardized per row in `13` |
| 6i | `13` (keys `h`, `k`) | ROI `example 2 (700 um)` from the picker; niche labels from `04`; colours = log1p CP10k from the count export, 2nd–98th percentile |

Two Xenium niche-vs-control analyses exist and both are used. `05` gives per-section effects with a Wilcoxon test
across the 8 whole tissues (6f bottom row, 6h). `04` gives a pseudobulk limma-voom model paired by section, on the 7
tissues with ≥ 50 macrophages per group (Table S10).

### Numbers in the text

- "twice the proportion of SecB neighbors": `06` prints the unit medians per tertile (whole tissues 7.9% → 15.2%, TMA
  4.0% → 7.7%; 1.92-fold).
- SPP1 up in SecB-niche macrophage bins (Fig 6e) and MMP9 (Table S10): `02` (16 µm and 8 µm).
- CA9 along the SecB territories (Fig 6b–d): the CA9 column of `visium_n2_bins_WT01.csv.gz` (`01b`).
- Gene statements of the two niche paragraphs (VEGFA, CD44, CTSL, TGM2 in macrophages; BIRC3, CD47, CD55, IL32,
  LGALS9, SLC16A3, ICAM1, MMP7, TIMP1, TGM2 in SecB cells): `05`, `03`, `09`, Fig 2g, summarized by `11`.
- Statements at macrophage contact (BIRC3, LGALS9, ICAM1; 30 µm): `05` (`secb_contact`).
- Methods "Visium HD" rules: `01`, `02`. Methods "Macrophage niches" rules: `04`, `05`, `06`.

## Not in this folder

- `01b` and `08b` ran inline in the final run's output folder (2026-09-25) and are recovered from the analysis session
  log. Run on the final run's inputs, `01b` gives the deposited `visium_n2_bins_WT01.csv.gz` and `08b` the deposited
  `atlas_subtype_dotmatrix.csv`, value for value (checked 2026-10-22). Re-running `01` drops the five columns, so `01b`
  must follow it. No script prints the 6f legend's cell count (n = 295,296); it is the sum of `n_cells` over the six
  plotted subtypes in `atlas_pb/groups.csv`.
- **Staging of the Visium HD Space Ranger outputs** (by hand) into `$OUTPUT_ROOT/fig6_singler/visium_local/`, one
  folder per section named by study code: `WT01B1/` and `WT03B1/` (16 µm `filtered_feature_bc_matrix.h5` and
  `tissue_positions.parquet` in the folder itself) and `WT01B1_8um/` (8 µm `filtered_feature_bc_matrix.h5` and
  `spatial/tissue_positions.parquet`).
- **Hand copy:** `inputs/pathway_gene_sets_v2.csv` (read by `05`) is a copy of
  `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv`:
  `mkdir -p $OUTPUT_ROOT/fig6_singler/inputs && cp $OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv $OUTPUT_ROOT/fig6_singler/inputs/`.
- **Second Visium HD section.** `01` and `03` also process WT03. The per-gene calls in `fig6_concordance_table.csv`
  (`11`) use both sections; the 6h Visium rows show WT01, and no panel reads `visium_n2_bins_WT03.csv.gz`.
- **Interactive ROI picker.** `python figures/figure6/fig6_roi_picker/roi_server.py` serves http://127.0.0.1:8793
  (WT01 and WT02 Xenium, and the WT01 Visium bins). **Save** appends the ROI to
  `$OUTPUT_ROOT/fig6_singler/fig6_rois.json`, a list of entries with `name`, `sample`, `x`, `y`, `size` (µm) and `rows`
  (plus `kind: visium` and `aspect` for a Visium window). The final run used:
  - `visium roi`: `sample` `VISIUM`, a 1.6 mm square centred on bin column 226, row 156 (panels a–d);
  - `example 2 (700 um)`: `sample` `WT02`, 700 µm, centred on x = 6549 µm, y = 2314 µm (panel i).

  The file is append-only and `13` takes the last `kind: visium` entry, so a new Visium save moves panels a–d.
  `sample` must be the study code that names the cell tables (`WT01`, `WT02`).
- **Cell set of the count export.** Every Xenium script in this folder keeps only the cells present in the count export
  (`spatial/05_count_export/as_run/02_deg_precompute.R`, the script that ran, unchanged). That export keeps cells with an
  SFE `cell_label`. For TMA B it holds fewer cells than the SingleR tables, so the TMA arm of 6g runs on that subset.
- **Assembly.** The manuscript page was assembled by hand in Affinity Designer from the per-panel files.
