# models/: patient-derived organoids and the KPCA.B mouse model

These folders hold trimmed copies of the scripts that ran in two lab project folders and wrote the objects and
tables behind Figure 3 and Figure S10:

| Folder | Lab project folder | What it covers |
|---|---|---|
| [`organoids/scripts/`](organoids/README.md) | `202603_organoids` | PDO snRNA-seq object (8 models, 10,455 cells) and PDO UCell scores; OPTO98 conditions (11,543 cells); scoring of PDOs, atlas pseudobulks, DepMap organoids, GSE299200 and GSE284507 in one UCell space; the Fig 3f–h panels |
| [`kpca_mouse/scripts/`](kpca_mouse/README.md) | `202608_kpca_secb_mouse` | KPCA.B scRNA-seq: QC, doublets, inferCNV, cell types (20,598 cells, 12,429 tumour cells), orthologues, de novo NMF and its alignment to the human NMF factors, the handoff object, and the Fig 3i panel |

**What "trimmed" means here.**
- The code that wrote the outputs the paper uses is kept as it ran, line for line. Plots, reports and analyses that
  are not in the paper were removed, with the imports and settings that only they used. The untrimmed scripts stay
  in the lab project folders. Do not edit the copies otherwise; they record what was run.
- Each trimmed script opens with a note: the project file it was copied from, what it still writes and what was
  removed. The original header below the note is unchanged, so it still describes the untrimmed script. Scripts
  with nothing to remove are unchanged copies. Every copy keeps the file name of its original; the folder READMEs
  say where the originals are.
- Their paths are those of the original project folders: absolute paths on the lab file server, or paths relative to
  the project root. They do not use `config/` or `DATA_ROOT`. To rerun one, run it from a copy of its project folder,
  or change the path constants at the top of the script; each README lists them.
- Several organoid scripts name the organoid project by earlier folder names (`organoids`, `2026_organoids`), and some
  name the atlas project `2026_final_atlas` instead of `202601_final_atlas`.

**How the repository uses their outputs.** The figure scripts in `figures/figure3/` and
`figures/supplementary/16_figS10_kpca_mouse_vs_human_programs.py` do not run anything in `models/`. They read the
outputs:
- **Organoids.** `figures/figure3/01_…`, `02_…` and `03_…` (Fig 3a–e) read the PDO object and the PDO and OPTO98 score
  tables under `$ORGANOIDS_ROOT/output/`. `ORGANOIDS_ROOT` is the deposited organoid project folder,
  `$DATA_ROOT/202603_organoids` (the default of `run_panel.sh`; `config.yml` has a placeholder, so set it when running
  a script directly).
- **Mouse.** `figures/figure3/04_…` and `05_…` (Fig 3j, k) read the KPCA.B handoff object
  `$DATA_ROOT/202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad` (override: `$KPCA_H5AD`); `04_…`
  also reads `output/02_orthologs/ortholog_map.csv` (override: `$KPCA_ORTHO`). The Figure S10 script reads the same
  two files and `output/31_fulltumour_nmf/` and `output/32_secb_align/`, under `$DATA_ROOT/202608_kpca_secb_mouse/`.
- **Table S6.** `tables/06_tableS6_secA_secB_signatures.py` parses the SecA / SecB lists of `00_setup.R` and
  `04v5_secB_only_comparison.py` from the deposited organoid project folder (`$DATA_ROOT/202603_organoids/scripts/`).
  The copies here carry the same lists.
- **Panels drawn here.** Two panels are drawn by scripts in `models/` and not by a repository figure script:
  - Fig 3f–h: `organoids/scripts/12h_panels_fgh_layout.py`
  - Fig 3i: `kpca_mouse/scripts/41b_export_1p8in.py`

**Atlas inputs.** Both pipelines read deposited atlas outputs written by the scripts in
`atlas/03_epithelial_nmf/as_run/`:
- the atlas UCell export `output/18_ucell_atlas/` (organoids);
- the NMF loadings `output/11d_epithelial_nmf/11d_nmf_loadings.csv` and gene classes
  `output/11e_nmf_characterization/11e_gene_classification.csv` (mouse).

**Environments.** The R steps (PDO object with scDblFinder, UCell, GSE299200 processing, inferCNV) ran in E1 (R 4.4.3,
Bioconductor 3.20). The Python steps ran in E2 (Python 3.12.9). Cell Ranger, CellBender and QuPath are vendor tools
(E6). See [`env/README.md`](../env/README.md).

**Missing code.** Each README tags with **[code missing]** the steps whose code is not on disk. For the organoids
these are the Parse split-pipe run and the OPTO98 Cell Ranger / MULTI-seq / QC; the public datasets were downloaded
by hand. The mouse has no missing code: the steps that ran as inline code are recovered from the session log, and
`kpca_mouse/upstream/` records the Cell Ranger multi configuration and the CellBender command.
