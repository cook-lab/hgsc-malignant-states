# models/organoids: patient-derived organoids (PDOs) and OPTO98

Trimmed copies of the scripts that ran, from `scripts/` of the organoid project folder `202603_organoids` (see
[`../README.md`](../README.md)). The code that wrote the outputs the paper uses is kept as it ran; plots and analyses
not in the paper were removed. The untrimmed scripts stay in the project folder. `03`, `10_markers_config.R` and
`11c` had nothing to remove and are unchanged copies.

The scripts build the PDO snRNA-seq object, score PDOs and the OPTO98 perturbations with UCell in the atlas gene
space, and score PDOs, atlas pseudobulks and three public datasets for SecB markers and function. Figure 3a–e reads
their outputs through `$ORGANOIDS_ROOT` (= the organoid project folder). Fig 3f–h is drawn here by
`12h_panels_fgh_layout.py`.

Paths below are relative to the organoid project folder. `<atlas>` is the atlas project folder `202601_final_atlas`,
which is `$DATA_ROOT/2026_final_atlas` in the repository config.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| – | `00_setup.R` | – | (sourced) Seurat, `project_dir`, `output_dir`; gene lists `secA_genes` (8, with BCAM, which 04 and 08 drop) and `secB_genes` (7) | 04, 08, 10a; `tables/06_tableS6_secA_secB_signatures.py` (Table S6) parses the gene lists of the copy in the deposited project folder, `$DATA_ROOT/202603_organoids/scripts/00_setup.R` |
| 1 | `03_create_seurat_untreated_baseline.R` | `data/single_cell/all-sample/DGE_filtered/{count_matrix.mtx, all_genes.csv, cell_metadata_updated.csv}` **[code missing]**: split-pipe output, not on disk | `output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds` (8 models, 10,455 cells, `pdo_model`, UMAP); an identical copy sits at the project root as `seurat_untreated_baseline.rds` | Fig 3a (`figures/figure3/01_…`); atlas 18b; 04; 12 |
| 2 | `atlas/03_epithelial_nmf/as_run/18b_ucell_atlas_scoring.R` (atlas, listed for order) | the PDO object (gene names only) and the atlas export of 18a | `<atlas>/output/18_ucell_atlas/shared_gene_list.txt` (14,271 genes = atlas genes ∩ PDO genes), `atlas_ucell_scores.csv` | 04, 04v4, 04v5, 08, 12 |
| 3 | `04_secretory_polarization_v4_UCell_aligned.R` | PDO object; `<atlas>/output/18_ucell_atlas/shared_gene_list.txt` | `output/02_Secretory_Polarization_v4_UCell_aligned/{organoid_ucell_scores_v4_aligned.csv, seurat_secretory_scored_v4_aligned.rds}` | 04v4, 04v5 |
| 4 | `04v4_comparison_report.py` | v4 scores; atlas UCell scores and metadata | `output/02_Secretory_Polarization_v4_UCell_aligned/organoid_ucell_scores_v4_with_atlas_class.csv` | text: "fewer than 1% of cells reaching the SecB signature threshold" |
| 5 | `04v5_secB_only_comparison.py` | v4 scores; atlas UCell scores and metadata | `output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv` | Fig 3b (`figures/figure3/02_…`); Table S6 parses `SECB_GENES` of the copy in `$DATA_ROOT/202603_organoids/scripts/` |
| 6 | `08_opto98_ucell_aligned.R` | `data/seurat_opto98.rds` **[code missing]**: writer not on disk; `<atlas>/output/18_ucell_atlas/shared_gene_list.txt` | `output/08_OPTO98_SecB_Conditions/{opto98_ucell_scores_aligned.csv, seurat_opto98_scored.rds}` (11,543 cells) | Fig 3c–e (`figures/figure3/03_…`); 08b |
| 7 | `08b_opto98_condition_comparison.py` | OPTO98 scores | `output/08_OPTO98_SecB_Conditions/opto98_statistical_tests.csv` | text on Fig 3c–e |
| 8 | `11c_atlas_pseudobulk.py` | `<atlas>/output/18_ucell_atlas/{atlas_secretory_counts.mtx.gz, atlas_secretory_genes.tsv, atlas_secretory_barcodes.tsv, atlas_secretory_metadata.csv}` | `data/11_nextgen/atlas_pseudobulk_{counts,meta}.csv` (369 groups) | 12 |
| 9 | `10a_geo_cluster_annotate.R` (sources `00_setup.R`, `10_markers_config.R`) | `data/GSE299200/filtered/{PBMCs_control, Tumor_PBMC_hypoxia, Tumor_PBMC_physoxia}/` (10x matrices) | `output/10_HGSC_public_SecAB_scoring/{10a_seurat_clustered_UNLABELED.rds, 10a_cluster_module_scores.csv, 10a_seurat_epithelial.rds}` | 12 |
| 10 | `12_gse284507_spheroid_organoid_secAB.R` | `seurat_untreated_baseline.rds` (project root); `data/11_nextgen/{atlas_pseudobulk_counts.csv, atlas_pseudobulk_meta.csv, next_gen_expression.csv, depmap_model_metadata.csv}`; `output/10_HGSC_public_SecAB_scoring/10a_seurat_epithelial.rds`; `data/GSE284507/GSE284507_normalized_counts.xlsx`; `<atlas>/output/18_ucell_atlas/shared_gene_list.txt` | `output/12_GSE284507_spheroid_organoid/12_scores_all_sources.csv` | 12h |
| 11 | `12h_panels_fgh_layout.py` | `output/12_GSE284507_spheroid_organoid/12_scores_all_sources.csv` | `fig3F_marker_vs_function.{svg,png}` (f), `fig3GH_shared.{svg,png}` (g above h). Written to the folder in `OUT` (final run: `$OUTPUT_ROOT/figures/figure3/`) | Fig 3f–h |

`10_markers_config.R` is not run on its own. It holds the GSE299200 annotation markers sourced by `10a`.

Order constraints: step 1 comes before the atlas step 18b (step 2), which reads the PDO genes. Steps 3–6 need step 2;
step 7 needs step 6. Step 8 needs only the atlas export and step 9 only the downloaded GSE299200 data. Step 10 needs
steps 1, 2, 8 and 9 and the downloaded DepMap and GSE284507 files; step 11 needs step 10.

## What each step does

- **03 (PDO object).**
  - Keeps the 8 untreated models (`samples_keep`): OPTO98, PDO66, OPTO112 and OPTO129 from primary tumors, and
    OCAD93, OCAD96, OCAD97 and OCAD106 from ascites.
  - Removes genes without a symbol, lncRNA, miRNA, snoRNA and ribosomal-protein genes; mitochondrial genes are kept.
  - Calls doublets with scDblFinder per sample (seed 42), then keeps singlets with > 500 UMIs, > 300 genes and < 20%
    mitochondrial reads: 10,455 cells.
  - Normalizes, takes 2,000 variable genes, 30 PCs and a UMAP (Fig 3a).
  - Uses paths relative to the project folder, so run it from there.
- **04 (PDO UCell).** UCell (maxRank 1500, raw counts) on the 14,271 shared genes. SecA uses the 7 genes left after
  BCAM is removed and SecB the 7 SecB genes, the same as the atlas scores of 18b. Writes `SecA_UCell`, `SecB_UCell` and
  `sec_polarization` per cell.
- **04v4 (the "< 1%" statistic).** Classes each PDO cell by the atlas polarization distribution: SecA below the atlas
  50th percentile, SecB at or above the 75th percentile (`atlas_class`). The share of all PDO cells with
  `atlas_class == "SecB"` is the "fewer than 1%" of the Results: 11 of 10,455 cells (0.105%).
- **04v5 (SecB-only table).** Uses the SecB UCell score only, derives atlas thresholds (Youden and the atlas SecB 25th
  percentile), and writes the per-cell table with `model` and `SecB_UCell` that Fig 3b plots. The ascites vs primary
  test (p = 0.029) is computed by `figures/figure3/02_fig3b_pdo_secB_by_model.py`, with the model as the unit.
- **08 (OPTO98 UCell).**
  - Drops the cells that MULTI-seq demultiplexing called `multiplet` (4,017) or `negative` (43) in the `barcode`
    column of `seurat_opto98.rds`, leaving 11,543 cells.
  - Scores them as in 04 (maxRank 1500, raw counts, 14,271-gene space).
  - The condition of each cell is in `experiment` (Growth, Treatment, Media) and `sample_id`.
- **08b (OPTO98 condition statistics).** Two-sided Mann–Whitney tests on per-cell SecB UCell: each day vs day 2, each
  cytokine vs PBS, each medium vs BM. They are written to `opto98_statistical_tests.csv`; a Kruskal–Wallis test
  across days is printed. These back the text on Fig 3c–e; the panels show no tests.
- **11c (atlas pseudobulk).** Sums raw counts of atlas secretory cells per patient × `celltype_nmf`, keeping groups of
  ≥ 50 cells (369 groups). Heavy: it streams the 5.5 GB count matrix.
- **10a (GSE299200).**
  - Processes primary HGSC cells co-cultured with autologous PBMCs, in hypoxia and physoxia, with a PBMC-only
    control. QC is as in 03; the data are integrated and clustered at resolution 1.0.
  - Annotation follows the authors' cell-type scheme. The epithelial compartment is kept for scoring.
  - The cluster → cell-type map was reviewed by hand against the per-cluster module scores of Section 7
    (`10a_cluster_module_scores.csv`) and written into Section 8 of the script. It refers to the cluster numbers of
    the cached `10a_seurat_clustered_UNLABELED.rds`, so a run from scratch must check the map.
- **12 (one UCell space; Methods "Organoid polarization quantification").**
  - Builds one count matrix on the genes common to all sources and the 14,271-gene list:
    - 369 atlas pseudobulks;
    - 8 PDO pseudobulks (summed per model);
    - DepMap next-generation organoids (TPM, un-logged);
    - GSE299200 epithelial pseudobulks per condition (2);
    - GSE284507 organoids and spheroids (28 samples).
  - Converts to CPM and runs one UCell call (maxRank = all common genes) with four modules: SecB marker (7 genes), HIF
    (NDRG1, SLC16A3, VEGFA, CA9, BNIP3, P4HA1, PDK1, SLC2A1), proliferation (MKI67, TOP2A, PCNA, CCNB1, CDK1) and
    SecA (7).
  - Function = z(HIF) − z(proliferation), z-scored over all samples.
- **12h (Fig 3f–h).**
  - f plots SecB marker against function for the in vitro systems. The shaded "SecB niche" starts at the atlas
    SecB-pseudobulk medians of both scores.
  - g and h plot HIF and quiescence (− proliferation) per system, with the atlas SecB median as a reference line.

## Paths in the scripts

Edit these to run outside the lab file server:

| Script | Path constants |
|---|---|
| `00_setup.R` | `project_dir` (earlier folder name `organoids`) |
| `03_…` | relative paths `data/single_cell/…`, `output/…`; run from the project folder |
| `04_…` | sources `scripts/00_setup.R` (run from the project folder); `shared_genes_file` |
| `04v4_…`, `04v5_…` | `ORG_DIR` (earlier name `organoids`), `ATLAS_DIR` |
| `08_…` | sources `scripts/00_setup.R`, then sets `project_dir` (earlier name `2026_organoids`); `shared_genes_file` |
| `08b_…` | `PROJECT` (earlier name `2026_organoids`) |
| `10a_…` | the absolute path of `00_setup.R`; `project_dir` |
| `11c_…` | `AD` (atlas export), `OUT` |
| `12_…` | `proj`, `atlasdir` |
| `12h_…` | `SRC`, `OUT` |

## Panel and Methods numbers fed by this folder

| Item | Steps |
|---|---|
| Fig 3a; Methods: QC thresholds, 10,455 cells | 03 |
| Methods: organoids scored in the 14,271-gene space shared with the atlas | atlas 18b (`shared_gene_list.txt`) → 04, 08 |
| Methods: organoid panels show the SecB UCell score | 04 → 04v5 (PDOs), 08 (OPTO98) |
| Fig 3b; text: ascites vs primary, p = 0.029 | 03 → 04 → 04v5 → `figures/figure3/02_…` |
| Text: fewer than 1% of PDO cells reach the SecB threshold | 03 → 04 → 04v4 |
| Fig 3c–e; Methods: 4,017 doublets, 43 unassigned, 11,543 cells | `seurat_opto98.rds` → 08 → `figures/figure3/03_…` |
| Text on Fig 3c–e (decline over days, IFNγ, media) | 08 → 08b |
| Fig 3f–h; Methods: organoid polarization quantification | 11c, 10a → 12 → 12h |

## Missing code

- **[code missing]** The Parse split-pipe run (Parse Evercode WT v3). It wrote the `DGE_filtered` matrices
  (`count_matrix.mtx`, `all_genes.csv`, `cell_metadata_updated.csv`) that step 03 reads. The matrices are not on disk:
  `data/single_cell/` is empty. What survives is the PDO object
  `output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds`.
- **[code missing]** Cell Ranger, MULTI-seq demultiplexing and QC for the OPTO98 library. These wrote
  `data/seurat_opto98.rds` (15,603 cells; `barcode` holds the MULTI-seq call, plus `experiment`, `sample_id`, `Phase`
  and `seurat_clusters`). The object is in the project folder; the laboratory has had it since 2024-07-25, and no
  script on the lab drives writes it. It was built on 2024-07-25 with Seurat 5.0.2, and it
  records its later steps: `NormalizeData` (LogNormalize), `SCTransform` (v2, regressing `S.Score` and `G2M.Score`,
  3,000 variable features), `RunPCA`, `RunUMAP` and `FindNeighbors` on PCs 1–30, and `FindClusters` (resolution 0.2;
  a resolution-0.1 column is also present). The MULTI-seq calls and the QC filters are not recorded; the kept cells
  all have ≥ 1,002 genes and < 30% mitochondrial reads.
- **Downloads (no script).** The external data were downloaded by hand:
  - DepMap next-generation organoids: `data/11_nextgen/next_gen_expression.csv` and `depmap_model_metadata.csv`
  - GSE299200: `data/GSE299200/filtered/`
  - GSE284507: `data/GSE284507/GSE284507_normalized_counts.xlsx`

  The files are in the project folder.
