# Figure 3: secretory polarization in patient-derived organoids and the KPCA.B mouse model

Figure 3 asks whether in vitro and in vivo models of HGSC reproduce SecA ↔ SecB polarization. Panels a–e show eight
patient-derived organoid (PDO) models and OPTO98 cultured over time, with cytokines and in different media. Panels f–h
score public organoid and spheroid datasets for SecB markers and SecB function (HIF induction, quiescence). Panels i–l
show one KPCA.B orthotopic syngeneic mouse tumour: scRNA-seq cell types, epithelial polarization, Krt7, and KRT7 IHC.

Only a, b, c–e, j, k and l are drawn by scripts in this folder. Panels f–h and i are drawn by scripts in `models/`:
trimmed copies of the scripts that ran in the organoid and mouse project folders (`202603_organoids`,
`202608_kpca_secb_mouse`), which keep the untrimmed scripts (see [Not in this folder](#not-in-this-folder)).

## Run order

**Requires** (run once, in other folders):
- `atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py` → `atlas/03_epithelial_nmf/as_run/17_prepare_nmf_labels.py`
  → `atlas/03_epithelial_nmf/as_run/18a_ucell_atlas_export.py` → `atlas/03_epithelial_nmf/as_run/18b_ucell_atlas_scoring.R`:
  the deposited atlas UCell export `$DATA_ROOT/2026_final_atlas/output/18_ucell_atlas/` (`atlas_ucell_scores.csv`,
  `atlas_secretory_metadata.csv`, `shared_gene_list.txt` with 14,271 genes). It gives the atlas reference lines of
  b–e and the gene space in which PDOs and OPTO98 are scored. `18b` intersects the atlas genes with those of the PDO
  object, so `models/organoids/scripts/03_create_seurat_untreated_baseline.R` runs before it.
- `models/organoids/scripts/`: `03_create_seurat_untreated_baseline.R` → `04_secretory_polarization_v4_UCell_aligned.R`
  → `04v5_secB_only_comparison.py` (PDO scores for a, b) and `08_opto98_ucell_aligned.R` (OPTO98 scores for c–e).
  Outputs are under `$ORGANOIDS_ROOT/output/`. See [`models/organoids/README.md`](../../models/organoids/README.md).
- `models/kpca_mouse/scripts/`: `01_load_qc_kpca.py` … `33_consolidate_handoff.py` → the KPCA.B handoff object
  `kpca_final_annotated.h5ad` (20,598 cells), and `02_orthologs.py` → `ortholog_map.csv`. See
  [`models/kpca_mouse/README.md`](../../models/kpca_mouse/README.md).

Run from the repository root. Renders, and the two CSVs of step 2, go to `$OUTPUT_ROOT/figures/figure3/`; the tables
of step 4 go to `$OUTPUT_ROOT/figure3/`.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `01_fig3a_pdo_models_umap.R` | `$ORGANOIDS_ROOT/output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds` (UMAP, `pdo_model`) | `organoids_models_umap_letter.{png,svg}` | a |
| 2 | `02_fig3b_pdo_secB_by_model.py` | `$ORGANOIDS_ROOT/output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv`; atlas `atlas_ucell_scores.csv` + `atlas_secretory_metadata.csv` | `secB_by_model_violin_letter.{png,svg}`, `secB_by_model_violin_letter_models.csv` (per-model n, median, mean, % zero), `secB_by_model_violin_letter_stats.csv` (the ascites vs primary test) | b |
| 3 | `03_fig3cde_opto98_secB_violins.py` | `$ORGANOIDS_ROOT/output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv`; the same atlas files | `secB_DEF_combined_violin.{png,svg}`, which is panels c–e (D/E/F in the file name and in the code stand for c/d/e) | c, d, e |
| 4 | `04_fig3j_kpca_polarization_7plus7.py` (uses `_kpca_umap_panel.py`) | `$KPCA_H5AD`, `$KPCA_ORTHO`, `shared/signatures.yml` | `$OUTPUT_ROOT/figure3/kpca_polarization_7plus7_cells.csv` (per-cell SecA / SecB UCell and polarization), `$OUTPUT_ROOT/figure3/kpca_polarization_7plus7_stats.csv` (gene detection; agreement with the stored 5 + 4-gene score, quoted in Table S6); `Fig3_j_mouse_polarization_7plus7.{svg,png}` | j |
| 5 | `05_fig3k_kpca_krt7_umap.py` (uses `_kpca_umap_panel.py`) | `$KPCA_H5AD` | `Fig3_k_mouse_KRT7.{svg,png}` | k |
| 6 | `06_fig3l_krt7_ihc_scale_bar.py --downsample 2 --um 100 <snapshot>/krt7_d26_t8.svg` | the QuPath snapshot `krt7_d26_t8.svg` (path given on the command line); for a plausibility check only, the scan sizes in `$DATA_ROOT/20260606_murine_secb_ihc/output/qupath_k7/data/*/server.json` | `Fig3_l_KRT7_IHC_krt7_d26_t8.{png,tif,svg}` (the SVG holds the image plus a vector scale bar) | l |

`_kpca_umap_panel.py` is a helper, not run on its own. It holds the shared 1.8 in UMAP square of j and k, with the
colour bar underneath.

Notes on the steps:
- **Step 2.** The atlas metadata call the middle class "Transitioning epithelium". The script renames it to
  Intermediate epithelium on read. The dashed atlas lines are medians of a random 20,000-cell sample per class
  (seed 42 from `config.yml`). The same holds for step 3.
- **Step 3.** Panel e shows the eight media BM, BM1, BM2, BM02, BM03, CGM, CGM + ascites and HPLM; PBS is drawn in d.
- **Step 4.** It rescores the 12,429 tumour cells (`celltype == "Tumour"`, the cells with `epi_umap_1/2`) with UCell
  on the mouse orthologues of the 7 + 7 signature genes in `shared/signatures.yml`. This uses maxRank 1500 on the
  log-normalised matrix, the implementation of `models/kpca_mouse/scripts/24_secb_whole_umap.py`. Polarization is
  SecB − SecA. The script asserts the cell count, and that its UCell reproduces the 5 + 4-gene scores stored in the
  handoff object. Runs in about a minute.
- **Step 6.** At downsample 2 the snapshot has 0.4402 µm/px (scan resolution 0.2200804 µm/px, the `FULL_UM_PER_PX`
  constant), so 100 µm = 227.2 px. The script prints, for downsample 2 and 10, the implied field size and whether it
  fits the QuPath scan sizes. If the QuPath project is not available, the check has no scans to compare with (it
  then prints "LARGER than every scan") and the panel is still written.

## Environment

| Variable | Steps | Value in the final run (written generically) |
|---|---|---|
| `OUTPUT_ROOT` | all | `$OUTPUT_ROOT`; renders go to `$OUTPUT_ROOT/figures/figure3/` |
| `ORGANOIDS_ROOT` | 1–3 | the organoid project folder `202603_organoids`, i.e. `ORGANOIDS_ROOT=$DATA_ROOT/202603_organoids` |
| `DATA_ROOT` | 2–6 | `$DATA_ROOT`, holding `2026_final_atlas/output/18_ucell_atlas/` (2, 3), the mouse project folder `202608_kpca_secb_mouse` (4, 5) and the IHC project folder `20260606_murine_secb_ihc` (6) |
| `KPCA_H5AD` | 4, 5 | not set; the default is the mouse project's handoff object, `$DATA_ROOT/202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad` |
| `KPCA_ORTHO` | 4 | not set; the default is `$DATA_ROOT/202608_kpca_secb_mouse/output/02_orthologs/ortholog_map.csv` |
| `LANG`, `LC_ALL` | 1 (R) | `en_US.UTF-8`, so that R reads `config.yml` in full |

- `ORGANOIDS_ROOT` does not follow `DATA_ROOT`: its default in `config.yml` is the placeholder
  `/path/to/deposit/202603_organoids`, so set it. With the placeholders, steps 1–5 stop with an error; they do not
  fall back to other files (step 6 only loses its plausibility check).
- Steps 2 and 3 read the deposited atlas UCell export in `$DATA_ROOT/2026_final_atlas/output/18_ucell_atlas/`; no
  script in this repository rewrites it.
- Seeds: steps 1–3 use `seed: 42` from `config.yml`; steps 4 and 5 use seed 0 for the point order, as in the mouse
  project's panels.

## Panel → script

- **3a** PDO UMAP (8 models, 10,455 cells). Drawn by `01_fig3a_pdo_models_umap.R`. The object, QC and UMAP come from
  `models/organoids/scripts/03_create_seurat_untreated_baseline.R`.
- **3b** SecB UCell per PDO model. Drawn by `02_fig3b_pdo_secB_by_model.py`, which also runs the test: a two-sided
  exact Mann–Whitney U test on the model medians, with the model as the unit (4 ascites vs 4 primary; U = 16,
  p = 0.029). Scores come from `models/organoids/scripts/04_secretory_polarization_v4_UCell_aligned.R` →
  `models/organoids/scripts/04v5_secB_only_comparison.py`. The "fewer than 1% of cells reaching the SecB signature
  threshold" in the text comes from the `atlas_class` column written by
  `models/organoids/scripts/04v4_comparison_report.py` (not drawn).
- **3c** OPTO98 over days 2–12, **3d** cytokines (PBS, WNT7A, TGFβ, TNFα, IFNγ) and **3e** media. Drawn together by
  `03_fig3cde_opto98_secB_violins.py`. Scores come from `models/organoids/scripts/08_opto98_ucell_aligned.R` (11,543
  cells). The panels show no tests. `models/organoids/scripts/08b_opto98_condition_comparison.py` writes two-sided
  Mann–Whitney tests of each day vs day 2, each cytokine vs PBS and each medium vs BM to
  `opto98_statistical_tests.csv`, and prints a Kruskal–Wallis test across days.
- **3f** SecB marker vs function. Drawn by `models/organoids/scripts/12h_panels_fgh_layout.py`
  (`fig3F_marker_vs_function`). Scores come from `models/organoids/scripts/12_gse284507_spheroid_organoid_secAB.R`,
  which needs `models/organoids/scripts/11c_atlas_pseudobulk.py` and `models/organoids/scripts/10a_geo_cluster_annotate.R`.
- **3g, 3h** HIF induction (g) and quiescence (h, negative proliferation), with the atlas SecB reference line. Drawn
  together by `models/organoids/scripts/12h_panels_fgh_layout.py` (`fig3GH_shared`, g on top). Scores come from the
  same chain as f.
- **3i** KPCA.B UMAP by cell type (20,598 cells). Drawn by `models/kpca_mouse/scripts/41b_export_1p8in.py`. It writes
  `41_mouse_umap_compartment_1p8in.{svg,png}`, which is panel i (coloured by `celltype`, despite the file name). Cell
  types come from `models/kpca_mouse/scripts/27_final_annotation_split.py`, carried into the handoff object by
  `models/kpca_mouse/scripts/33_consolidate_handoff.py`.
- **3j** Epithelial polarization (n = 12,429). Computed and drawn by `04_fig3j_kpca_polarization_7plus7.py`, with
  orthologues from `models/kpca_mouse/scripts/02_orthologs.py`. The colour scale is symmetric about 0, with limits
  at the 99th percentile of |polarization|.
- **3k** Krt7 on the epithelial UMAP. Drawn by `05_fig3k_kpca_krt7_umap.py`. Colour limits run from 0 to the 99th
  percentile.
- **3l** KRT7 IHC, day-26 tumour, 100 µm scale bar. Drawn by `06_fig3l_krt7_ihc_scale_bar.py` from a manual QuPath
  export.

Figure S10 uses the same KPCA.B cells. It is drawn by `figures/supplementary/16_figS10_kpca_mouse_vs_human_programs.py`,
from the handoff object and the outputs of mouse steps `02`, `31` and `32`.

## Not in this folder

- **Panels f–h:** `models/organoids/scripts/12h_panels_fgh_layout.py`. Its input and output folders are written in
  the script (`SRC`, `OUT`); the final run wrote to `$OUTPUT_ROOT/figures/figure3/`.
- **Panel i:** `models/kpca_mouse/scripts/41b_export_1p8in.py`. It writes into the mouse project's
  `output/41_fig3_umaps/`. Panels j and k are drawn by steps 4 and 5 here.
- **Manual steps:**
  - The Fig 3l image `krt7_d26_t8.svg` is a QuPath viewer snapshot of a Zeiss Z1 20× scan, exported by hand. No code
    exists for this export.
  - Each render was copied by hand under a `Fig3_<letter>_*` name, and the page was assembled in Affinity Designer.
- **Downloaded data (no script):** the DepMap next-generation organoid data and GSE299200 / GSE284507 (f–h), in the
  organoid project's `data/` folder.
- **Missing code** (details in the two `models/` READMEs):
  - **[code missing]** The Parse split-pipe run for the PDOs. It made the `DGE_filtered` matrices that
    `models/organoids/scripts/03_create_seurat_untreated_baseline.R` reads; they are not on disk. Its product survives
    as the PDO object `seurat_untreated_baseline.rds` (a).
  - **[code missing]** Cell Ranger, MULTI-seq demultiplexing and QC for OPTO98. These wrote the OPTO98 Seurat object
    `seurat_opto98.rds` in the organoid project's `data/` folder (c–e).
  - KPCA.B (i–k): nothing is missing. The steps that ran as inline code (the inferCNV gene positions, steps 21 and
    22, the epithelial split) are recovered in `models/kpca_mouse/scripts/`, and the Cell Ranger multi configuration
    and CellBender command are in `models/kpca_mouse/upstream/`.
