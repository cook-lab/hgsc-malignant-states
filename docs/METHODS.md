# Methods

> **Provenance note.** This expanded Methods document reflects everything actually
> performed in the canonical analysis pipeline (`atlas/` and `spatial/` of the
> `hgsc-malignant-states` repository) and supersedes the lighter Methods in the
> 2026-04-27 manuscript draft. Where the published draft was silent or imprecise,
> the exact parameter, software version, threshold, or random seed has been
> inserted from the canonical code. Discrepancies between the draft and the code
> are flagged inline in **[NOTE]** call-outs and collated at the end. All stochastic
> steps in the pipeline are seeded from a single global constant, `SEED = 42`
> (`config/config.py`), unless a different fixed seed is stated explicitly.

---

## Single-cell RNA sequencing atlas construction

Publicly available single-cell RNA sequencing (scRNA-seq) datasets from 13
independent studies of high-grade serous ovarian carcinoma (HGSC) were compiled
(Denisenko 2024, Geistlinger 2020, Hornburg 2021, Loret 2022, Luo 2024,
Nath 2021, Olalekan 2021, Olbrecht 2021, Regner 2021, Vázquez-García 2022,
Xu 2022, Zhang 2022, Zheng 2023). Raw count matrices were obtained in their
original formats (10x MTX, H5, CSV, H5AD) and gene identifiers were harmonized
to HGNC symbols using the MyGene Python API (v3.2.2). A standardized 16-column
metadata schema (barcode, sample identifier, study, patient identifier, sample
number, treatment status, histological subtype, stage, anatomic site, metastatic
site, age, treatment response, BRCA status, HRD status, TP53 status, reference)
was enforced across all studies; missing fields were encoded as `NA`. Gene sets
were intersected across all 13 studies, retaining 14,470 protein-coding genes
common to all datasets. Count matrices were stored as integer CSR sparse
matrices and concatenated on the gene intersection (`scanpy.concat`,
`join="inner"`), yielding 2,731,632 cells.

Cells from non-serous histologies were removed, after which per-cell quality
control removed cells with fewer than 500 total UMI counts or fewer than 300
detected genes, and samples with fewer than 500 cells were excluded, leaving
2,398,571 cells across 148 patients. Doublets were identified using Scrublet
(v0.2.3) run independently per sample (samples with fewer than 1,000 cells were
skipped) with an adaptive expected doublet rate scaled as
0.008 × (n_cells / 1,000) and `random_state = 42`. Cells with Scrublet scores
≥ 0.30 were removed and the low-cell-sample filter re-applied, leaving
2,326,532 cells. **[NOTE: two distinct Scrublet thresholds are applied — a
pre-integration cut at 0.30 (here) and a stricter post-integration cut at 0.25
(below); supplementary QC plots display the 0.25 line.]**

### Integration (CellAssign → scVI → scANVI)

Integration was performed once on a GPU cluster; the full integration code is
included in the repository (`atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`).
Re-running it is computationally expensive and not required to reproduce downstream
results — the integrated atlas object is deposited as an entry object — but the step
is included so it can be independently re-executed. For integration, 4,000 highly variable genes were
selected per batch using the Seurat v3 variance-stabilizing transform
(scanpy v1.11.4, `flavor="seurat_v3"`, `batch_key="sample_id"`). Initial cell-type
priors were generated with CellAssign using a curated marker matrix of 81 genes
across 16 cell types. An scVI model (negative-binomial likelihood,
`n_layers = 2`) was trained, followed by scANVI semi-supervised integration that
used the CellAssign predictions as label priors (`batch_key="sample_id"`,
`unlabeled_category="Unknown"`); both models were trained for up to 800 epochs
with early stopping (`patience = 10`, `monitor="elbo_validation"`) on an NVIDIA
V100 GPU (CUDA 12.6.2). The resulting 10-dimensional scANVI latent
representation (`X_scanvi`) was used for all downstream analyses. A
post-integration doublet filter removed cells with Scrublet scores ≥ 0.25,
leaving 2,294,893 cells. UMAP embeddings were computed on the scANVI latent
space (`scanpy.pp.neighbors`, `use_rep="X_scanvi"`, `n_neighbors = 10`;
`min_dist = 0.2`).

## Cell-type annotation

Cell-type annotation followed a hierarchical scheme. For coarse compartment
assignment (`celltype_level1`), Leiden clustering (resolution = 0.3) on the
scANVI neighbor graph produced clusters that were mapped to 13 cell-type
compartments (Epithelial, Mesothelial, Fibroblast, Smooth muscle, Pericyte,
Endothelial, T/NK cell, B cell, Plasma cell, Macrophage, Dendritic cell,
Neutrophil, Mast cell) based on dominant CellAssign composition per cluster.
Two artifact clusters (38 cells) were removed.

For fine-grained annotation (`celltype_level2`), each compartment was
independently subclustered: neighbors were recomputed per compartment
(`n_neighbors = 15`, `use_rep="X_scanvi"`), followed by Leiden clustering at
compartment-specific resolutions (range 0.05–0.5). Cluster identities were
assigned from differentially expressed markers (Wilcoxon rank-sum test,
BH-adjusted p < 0.05, log2 fold change > 0) and canonical marker expression.
Residual doublet and artifact clusters were identified and excluded, yielding a
final atlas of **1,980,703 cells** across 51 cell-type subtypes, 13 compartments,
148 patients, 371 samples and 12 anatomic sites. scANVI annotation confidence was
high (99.0% of cells assigned with prediction probability > 0.5). The epithelial
compartment comprised 575,366 cells; ciliated epithelial cells (6,678 cells,
1.2% of epithelium) were identified by FOXJ1, CAPS, RSPH1 and TPPP3, and the
remaining secretory epithelium by PAX8, MUC16 and WFDC2.

## NMF-based epithelial subtyping

To identify coordinated transcriptional programs within the epithelial
compartment, non-negative matrix factorization (NMF;
`sklearn.decomposition.NMF`, `init="nndsvda"`, `max_iter = 500`,
`random_state = 42`) was applied with **k = 10** factors. The input consisted of
3,000 highly variable genes (scanpy, `flavor="seurat_v3"`) selected from
log-normalized counts (`normalize_total`, `target_sum = 10,000`, then `log1p`);
the curated SecA/SecB signature genes (below) were forced into the HVG set to
guarantee their representation. NMF was fit on a stratified, seeded subsample of
50,000 epithelial cells (per-subtype allocation proportional to abundance,
minimum 100 cells per subtype), and the learned gene-loading matrix (H) was used
to project all 575,366 epithelial cells; per-cell factor usage (W) was
row-normalized to sum to 1.

Factor 3 was identified as the **SecA** (progenitor-like) program based on
enrichment for SOX17, WT1, PBX1, MECOM, LGR5, LPAR3 and RCN2, and Factor 2 as
the **SecB** (differentiated/adaptive) program marked by KRT17, KRT19, KRT7,
TACSTD2, SLPI, LCN2 and PRSS22. Factor identity was assigned automatically by
ranking each factor's mean signature-gene loading (normalized to the factor
maximum). To classify genes along the polarization axis, the 90th percentile of
each factor's loading distribution was used as a threshold, yielding 177
SecA-specific, 177 SecB-specific and 124 shared genes.

To assess metaprogram stability, consensus NMF was additionally performed
following Gavish et al.: per-patient NMF was run across k = 5–10 with 30 random
initializations per k; consensus programs were identified by hierarchical
clustering (average linkage, cosine-similarity threshold 0.5) and metaprograms
defined by clustering all patient consensus programs (cosine-similarity
threshold 0.3; minimum 5 programs from ≥ 3 patients). Two metaprograms were
recovered (a general-secretory program and a small metallothionein program);
critically, SecA and SecB did **not** separate into distinct metaprograms,
confirming that secretory polarization represents a continuous gradient rather
than discrete cell states.

Epithelial cells were partitioned along the polarization axis using percentile
cut-offs computed on **Factor 2 usage among non-ciliated epithelial cells**:
**SecA (below the 50th percentile), Intermediate (50th–75th percentile), and SecB
(at or above the 75th percentile)**. Ciliated epithelial cells were assigned
independently by canonical marker expression. This canonical epitype schema
(`celltype_nmf`) was written back into the atlas object and is the substrate for
all downstream epitype comparisons. **[NOTE: the published draft uses the legacy
label "Transitioning" for the 50th–75th-percentile class; the canonical schema
renames this to "Intermediate".]**

## Diffusion pseudotime analysis

Diffusion pseudotime (DPT) was computed on a stratified subsample of 40,000
epithelial cells. A neighbor graph was constructed on the scANVI latent space
(`n_neighbors = 30`), followed by diffusion-map computation (`n_comps = 10`) and
PAGA graph abstraction. The root cell was selected as the medoid of the SecA
cluster in diffusion-map space, and DPT values were projected to all epithelial
cells by inverse-distance-weighted interpolation from the 15 nearest neighbors
in scANVI space.

## Differential gene expression analysis

Differentially expressed genes between cell types and between epithelial states
were identified with the Wilcoxon rank-sum test
(`scanpy.tl.rank_genes_groups`, one-versus-rest). Genes with BH-adjusted
p < 0.05 and log2 fold change > 0 were retained as positive markers (yielding
51,592 significant gene–cell-type pairs across the 13 compartments, and 10,029
significant DEGs for the ciliated-versus-secretory contrast; Supplemental Data
2–3). A curated, cross-platform-portable signature of seven SecA markers
(MECOM, FBXO21, LGR5, LPAR3, PBX1, SOX17, RCN2) and seven SecB markers (KRT17,
KRT19, KRT7, LCN2, PRSS22, SLPI, TACSTD2) was defined for downstream UCell
scoring; this "noBCAM" 7-gene formulation deliberately excludes BCAM to maintain
robust cross-platform behavior.

## Pathway activity and functional characterization

Pathway and functional activity was quantified across epithelial states using
four complementary approaches. PROGENy pathway activity (14 pathways) was
inferred with the multivariate linear model (MLM) method via decoupler.
Gene-set enrichment was scored with MSigDB Hallmark gene sets via
`scanpy.tl.score_genes`. Transcription-factor activity was estimated from
DoRothEA regulons (confidence levels A, B, C) using the MLM method. Metabolic
flux was estimated with scFEA (Chang et al., 2021) using the human metabolic
model (**168 modules, 70 compounds**), trained on a stratified subsample of
20,000 epithelial cells (epochs = 100, learning rate = 0.008); scFEA runs in a
dedicated PyTorch sub-environment.

Cell-cycle phase was scored with the Tirosh S-phase (43 genes) and G2M (54 genes)
gene sets, with assignment based on 75th-percentile thresholds. Extended
functional profiling additionally scored EMT, stemness, DNA-damage response,
chromosomal instability (CIN70) and senescence (SenMayo) signatures via
`scanpy.tl.score_genes`. Statistical comparisons across epithelial states used
the Kruskal-Wallis H-test with Bonferroni correction, followed by pairwise
Mann-Whitney U tests.

## Copy-number variation inference

Copy-number variation (CNV) profiles were inferred from scRNA-seq counts using
CopyKAT (v1.1.0) run independently for each of 251 samples, with non-epithelial
cells supplied as the diploid reference and parameters `id.type = "S"`,
`cell.line = "no"`, `ngene.chr = 5`, `win.size = 25`, `KS.cut = 0.1`,
`distance = "euclidean"` (`set.seed(42)` for determinism). For samples with ≥ 20
aneuploid cells, subclonal structure was resolved by Ward.D2 hierarchical
clustering on 1 − Pearson correlation distances across genomic bins, with the
clone number k selected over {2, 3, 4, 5} by maximum mean silhouette width;
samples whose maximum silhouette fell below 0.15 were designated monoclonal.
Independence of epithelial state from clonal identity was tested with chi-square
or Fisher's exact tests with BH-FDR correction, supplemented by multinomial
logistic regression (5-fold cross-validated AUROC). Within-clone coexistence of
epithelial states was quantified by Shannon entropy. Across 248 samples with
sufficient cells (3,735 clones), 140 (58%) were monoclonal, 54 (22%) clonally
driven and 54 (22%) mixed; all three secretory states co-existed within
individual clones.

## Cell–cell communication analysis

Ligand–receptor interactions were inferred using LIANA (≥ v1.6.0) with the
multi-method rank-aggregate algorithm (`li.mt.rank_aggregate`) and the consensus
ligand–receptor resource (`resource_name="consensus"`, `expr_prop = 0.1`,
`use_raw = False`). Analysis was performed on a stratified, seeded subsample of
the atlas (maximum cells per cell type with a per-type minimum;
`n_perms` as configured, `seed = 1337`). Interactions with `magnitude_rank` ≤ 0.05
were considered significant, and differential communication between SecA and SecB
epithelium was quantified as the log2 fold change of significant-interaction
counts per signalling category. **[NOTE: the LIANA permutation seed (1337) is a
fixed analytical parameter preserved as published and is intentionally distinct
from the global SEED = 42.]**

## TCGA molecular-subtype assignment

TCGA molecular subtypes were assigned to atlas samples using ConsensusOV
(v1.30.0). Pseudobulk count matrices were generated per sample by summing raw
counts across all cells, TMM-normalized (edgeR), converted to log2(CPM + 1), and
mapped to Entrez IDs (org.Hs.eg.db). ConsensusOV classification was run with
default parameters, producing per-subtype posterior probabilities (differentiated,
immunoreactive, mesenchymal, proliferative).

## Patient-derived organoid culture and scRNA-seq

Eight HGSC patient-derived organoid (PDO) models were cultured: OPTO98, OPTO112
and OPTO129 (primary tumour-derived); OCAD93, OCAD96, OCAD97 and OCAD106
(ascites-derived; Princess Margaret Cancer Centre); and PDO66 (primary
tumour-derived; Western University). PDOs were maintained in Advanced DMEM/F-12
supplemented with GlutaMAX (2 mM), HEPES (10 mM), antibiotic/antimycotic (1×),
B27 (1×), N-acetyl-L-cysteine (1.25 mM), EGF (20 ng/mL), b-FGF (100 ng/mL),
FGF-10 (100 ng/mL), Noggin (100 ng/mL), nicotinamide (1 mM), Y-27632 (10 µM),
forskolin (10 µM) and β-estradiol (200 nM); PDO66 medium was additionally
supplemented with N2 (1×), SB431542 (0.5 µM) and BMP2 (100 ng/mL). PDOs were
embedded in Matrigel (Corning CB-40230) at 20,000–50,000 cells per 50 µL droplet
and passaged every 1–4 weeks following TrypLE dissociation (~45 min, 37 °C).

For scRNA-seq, cells were fixed and barcoded with the Parse Biosciences Evercode
Cell Fixation and WT v3 kits (~100,000 cells per sample fixed, permeabilized,
40-µm-strained and stored at −80 °C). Library preparation comprised three rounds
of combinatorial barcoding, cDNA amplification, 0.8× SPRI cleanup,
fragmentation, A-tailing and index PCR, with sequencing on an Illumina NovaSeqX
(1.5B flow cell). Reads were aligned to GRCh38 with split-pipe (v1.6.4) and the
count matrix filtered to 34,526 protein-coding genes. Cells were retained with
> 500 UMI counts, > 300 detected genes and < 20% mitochondrial reads; doublets
(3.4%) were removed with scDblFinder (v1.20.2), retaining 95.9% of cells, which
were processed in Seurat v5 (log-normalization, 2,000 HVGs, 50 principal
components, 30-dimension UMAP). PDO cells were scored for secretory polarization
with UCell (`AddModuleScore_UCell`) using the atlas-defined 7-gene SecA and SecB
signatures, with polarization computed as SecB_UCell − SecA_UCell. Time-course
(OPTO98, day 2–12) and perturbation experiments (OPTO98 treated with PBS, IFNγ,
TNFα, TGFβ or WNT7A) were scored identically. KRT19 protein expression was
quantified by flow cytometry in a subset of models.

## Targeted spatial transcriptomics

A custom 10x Genomics Xenium Prime panel was designed comprising **477
biological genes** spanning cell-type identity, epithelial programs, signalling
pathways, immune function, therapeutic targets and molecular classifiers
relevant to HGSC (Supplemental Data 5). Panel quality was assessed by
cross-platform comparison with the scRNA-seq atlas (per-gene
correlation-of-correlations, top-20-neighbor Jaccard, and scRNA/Xenium detection
rates); **28 genes with aberrant cross-platform behavior (> 3 QC flags) were
excluded from annotation, yielding 449 genes for SingleR classification**, while
all panel genes were retained for downstream expression analysis (20 genes were
present on the Xenium panel but absent from the atlas).

The panel was applied to three spatial cohorts: **8 whole HGSC tissue sections**
(5 treatment-naïve, 3 post-chemotherapy; ~1.9 million cells; samples OTB_2384,
SP24_25573, OTB_2432, OTB_2454, OTB_2457, OTB_2461, OTB_2417, SP24_24824); a
**tissue microarray (TMA) of 97 primary treatment-naïve HGSC patients**
represented by duplicate cores (590,090 cells); and **fallopian-tube epithelium
(FTE) cores (34,281 cells)**. Spatial objects were assembled in the Bioconductor
SpatialFeatureExperiment framework (v1.10.1; Voyager v1.10.0). **[NOTE on the
FTE cohort: 15 FTE cores were profiled and are used for cohort-level composition
summaries, but the canonical re-runnable cohort excludes two FTE whole-tissue
samples added after the preprint (FT1-1 and EAOC-1-FTE) to avoid cohort drift;
consequently spatial-clustering/polarization panels display 13 FTE cores while
cohort-level text and donut summaries reference 15. This reconciles the n = 13
vs n = 15 discrepancy in the figure legends.]**

### Probe and per-cell quality control

Transcript- and probe-level QC was performed per sample (negative-control probe
and unassigned-codeword rates), and per-cell QC removed low-quality cells before
normalization and log-transformation within the SpatialFeatureExperiment objects
(scater/scran v1.36.0). TMA objects were additionally split and filtered per
core.

## Cell-type annotation and polarization in spatial data

Cell types were annotated with SingleR (v2.10.0) against a downsampled scRNA-seq
atlas reference of **16,000 cells (1,000 cells per type × 16 atlas-derived
`xenium_celltype` labels)**; reference markers were derived by Wilcoxon
`rank_genes_groups`. SingleR was applied independently per sample on
log-normalized counts using the shared quality-controlled panel genes
(`assay.type.test = "logcounts"`, `assay.type.ref = "logcounts"`), with
fine-tuning and label pruning enabled (default delta-based pruning); per-cell
confidence was recorded as the difference between the best and median assignment
scores.

Secretory epithelial cells were subclassified along the polarization axis with
UCell (v2.12.0; `AddModuleScore_UCell` / `ScoreSignatures_UCell`) using the
7-gene SecA and SecB signatures and all shared panel genes for ranking
(`maxRank = nrow`). A single bivariate polarization score was computed as
**polarization_UCell = SecB_UCell − SecA_UCell**, capturing both poles (and
matching the atlas NMF "Intermediate" definition, unlike a SecB-only cut).
Atlas-calibrated, frozen thresholds were applied: `t_low` = the 75th percentile
of atlas SecA-class polarization, and `t_high` = the 25th percentile of atlas
SecB-class polarization, giving **SecA (polarization < t_low), Intermediate
(t_low ≤ polarization < t_high), and SecB (polarization ≥ t_high)**. Combined
with SingleR labels this produced **18 annotated spatial cell types**. The
thresholds were validated against the atlas (binary SecA-vs-SecB AUC for the
polarization score ≈ 0.94) and written to a frozen threshold-summary table.

## Spatial neighborhood analysis

For each secretory epithelial cell, the local neighborhood was defined as all
cells within a 50-µm radius using fixed-radius nearest-neighbor search
(`dbscan::frNN`, dbscan v1.2.4); neighborhoods were computed per core for TMA
data and per sample for whole tissue. Neighborhood features comprised cell-type
proportions (18 types), cell density, and mean UCell pathway scores (over all
neighbors and stratified by neighbor cell type). Across the cohort,
**1,651,007 secretory cells** were characterized within their 50-µm
neighborhoods.

Neighborhoods were classified by k-means clustering (**k = 10**) on cell-type
composition, fit on a stratified subsample of ~100,000 cells
(`nstart = 50`, `iter.max = 200`, `set.seed(42)`); k was selected with
within-cluster sum of squares, the Davies-Bouldin index and silhouette width.
All remaining cells were assigned to the nearest centroid (`FNN::get.knnx`,
k = 1). **[NOTE: the published draft states `nstart = 25`, `iter.max = 100`; the
canonical production script uses `nstart = 50`, `iter.max = 200`.]** In the
whole-tissue cohort, SecB cells were strongly self-aggregating (87% resided in
SecB-dominated neighborhoods).

## Co-localization and spatial association analysis

Spatial co-localization was assessed at three tiers. First, k-nearest-neighbor
enrichment (k = 20; `FNN::get.knnx`) quantified observed-versus-expected neighbor
proportions under complete spatial randomness, reported as log2 enrichment with
consensus (median) across TMA cores or whole-tissue samples; tumour-versus-FTE
contrasts used Wilcoxon rank-sum tests with BH-FDR correction. Second, Ripley's
cross-K function (`spatstat.explore::envelope`/`Kcross`, border correction,
rmax = 200 µm; 199 permutations for TMA, 19 for whole tissue) was computed for
12 focal cell-type pairs. Third, cross-context consistency between TMA and
whole-tissue enrichment matrices was assessed by Spearman correlation.

Bivariate spatial autocorrelation of SecA and SecB UCell scores was quantified
using **Lee's bivariate L statistic** with row-standardized fixed-radius spatial
weights (50 µm; `spdep::dnearneigh`, `spdep::nb2listw`, `style = "W"`;
spdep v1.4-2). Statistical inference used **999 permutations** with a two-sided
empirical p-value. Univariate global Moran's I (`spdep::moran.mc`, 999
permutations) and local indicators (LISA via `localmoran_perm`, and bivariate
BiLISA) were computed per cell to classify spatial regimes
(SecA-dominant / SecB-dominant / both / neither) with BH-FDR correction. SecA and
SecB were spatially segregated (Lee's L < 0; whole tissue p < 0.001, TMA
p < 0.05).

## Generalized additive modeling of polarization gradients

Generalized additive models (GAMs; `mgcv::gam`, mgcv v1.9-4) mapped
transcriptional and spatial features as a function of the epithelial
polarization score. Pooled models were fit as
`y ~ s(polarization_UCell, k = 10)` with REML estimation on a stratified
subsample of 50,000 whole-tissue secretory cells, and per-sample models
(`k = 5`) provided directional-consistency assessment. Distributional families
were chosen by response type: **beta regression (`betar`, logit link) for
proportions, negative binomial (`nb`) for counts, and Gaussian for continuous
scores** (e.g., UCell pathway scores). TMA validation used per-core Spearman
correlations for features significant in whole tissue (BH-FDR < 0.05, minimum
100 secretory cells per core). **[NOTE: the published draft states pooled
`k = 20` and per-sample `k = 10`; the canonical niche-succession GAM script uses
pooled `k = 10` and per-sample `k = 5`.]**

## Distance to vasculature

For each secretory epithelial cell, the Euclidean distance to the nearest
vascular cell (pericyte or endothelial) was computed by exact nearest-neighbor
search (`RANN::nn2`, k = 1; RANN v2.6.2), per sample for whole tissue and per
core for TMA (summarized as the per-patient mean of per-core medians). SecA,
Intermediate and SecB cells were compared by paired Wilcoxon signed-rank tests
(whole tissue) and pooled Wilcoxon rank-sum tests (TMA); distance increased
monotonically SecA → SecB (whole tissue p = 0.016, TMA p < 0.0001), and local
cell density decreased correspondingly (whole tissue p = 0.008, TMA p < 0.0001).

## UCell pathway scoring in spatial data

Pathway activity was scored with UCell (`AddModuleScore_UCell`) using **37 custom
gene modules** (Supplemental Data 6) covering proliferation, hypoxia, glycolysis,
NF-κB, RTK/RAS, TGFβ, type I and type II interferon, EMT, angiogenesis,
complement, immune checkpoint and therapeutic-target programs. Scoring was
performed on raw counts with `maxRank` equal to the number of panel genes, and
modules required a minimum of 2 panel genes to be scored.

## Immune phenotyping by neighborhood

Immune-cell composition (macrophage, T, NK, B, plasma) was modeled as a function
of local epithelial hypoxia and glycolysis scores using GAMs as above. Niche
metabolic stress was defined as the mean of epithelial hypoxia and glycolysis
UCell scores within 50 µm of each tumour-microenvironment cell, z-scored within
sample (whole tissue) or within core (TMA) and binned into tertiles. Immune-cell
presence as a function of niche metabolic stress was additionally modeled with
generalized linear mixed models (`lme4::glmer`, binomial family, sample as a
random intercept; lme4). T-cell and NK-cell exhaustion were scored with
composite signatures (T cells: PDCD1, HAVCR2, LAG3, TIGIT, CTLA4; NK cells:
HAVCR2, TIGIT, LAG3) and compared between high-stress (top decile) and
low-stress (bottom decile) epithelial neighborhoods using paired Wilcoxon
signed-rank tests (T-cell exhaustion p = 0.039, n = 7; NK-cell exhaustion
p = 0.031, n = 7). Macrophage transcriptional programs were assessed by
niche-specific pseudobulk differential expression (Wilcoxon rank-sum on
log2(CPM + 1), minimum 20 macrophages per group, BH-FDR < 0.05,
|log2FC| > 0.25).

## Nuclear and cellular morphometry

Cell- and nucleus-segmentation polygons from Xenium were used to compute
morphometric features with the sf R package (v1.1-1): cell area, nuclear area,
nuclear-to-cytoplasmic ratio, perimeter, circularity (4π × area / perimeter²),
solidity (area / convex-hull area), eccentricity (from PCA of polygon vertex
coordinates, computed in chunks of 25,000 cells) and nuclear-centroid offset.
Geometries were validated (`st_make_valid`) before measurement. Quality filters
excluded cells with area outside 10–5,000 µm² (cell) or 5–1,000 µm² (nucleus).
Comparisons across SecA, Intermediate and SecB cells used pairwise Wilcoxon
tests with a minimum of 30 cells per epitype per sample. Nuclear area and N:C
ratio decreased from SecA to SecB (whole tissue n = 8, p = 0.008; TMA n = 62,
p < 0.0001), and macrophages in SecB-rich niches were larger and more circular
than those in SecA-rich niches (each p = 0.008, 8/8 samples concordant).

## Survival analysis

Clinical outcomes were assessed in the TMA cohort (97 patients). Cell-type
densities (cells/mm²) and the mean per-patient polarization UCell score were
computed per patient. Kaplan-Meier curves for 5-year overall survival (OS) and
progression-free survival (PFS) were compared by log-rank tests after
stratification at the median, and univariate Cox proportional-hazards regression
(survival v3.8-6; survminer v0.5.2) was fit for each of the 18 cell-type
densities and for continuous polarization. Key associations: SecA proportion
(OS HR = 0.71, 95% CI 0.55–0.92, p = 0.009; PFS HR = 0.74, 0.60–0.91,
p = 0.005), SecB proportion (OS HR = 1.31, 1.04–1.64, p = 0.023; PFS HR = 1.28,
1.06–1.56, p = 0.011) and polarization (OS HR = 1.45, 1.07–1.96, p = 0.018;
PFS HR = 1.42, 1.12–1.79, p = 0.004). Per-core Xenium transcript levels were
correlated with TMA immunofluorescence protein MFI (KRT7, KRT19, KRT18, VIM,
E-cadherin) by Spearman correlation.

### External validation (TCGA-OV)

The 7-gene SecA and SecB signatures were applied to TCGA-OV bulk RNA-seq
(n = 314 tumours) after a pre-registered cohort lock. UCell scoring was performed
on log2(TPM + 1) expression, and polarization was computed as SecB_UCell −
SecA_UCell (z-scored across samples; an SecA/SecB log2-ratio score was also
derived). Patients were stratified into tertiles by polarization score and
compared by Kaplan-Meier analysis with log-rank tests (high-versus-low bins:
5-year OS and PFS p = 0.01). Cox proportional-hazards regression was fit
univariately and multivariately, forcing age and stage (stage collapsed to
{II/III, IV}) and adding data-driven covariates by cross-validation; the
polarization score remained significant after adjustment (univariate OS
HR = 2.00, 1.25–3.20, p = 0.004; PFS HR = 1.79, 1.15–2.78, p = 0.009). Tumour
epithelial fraction was estimated by BayesPrism (v2.2.3) deconvolution from a
12-cell-type atlas reference (300 cells per group) and included as a covariate.
Robustness was assessed by a 10,000-shuffle label-permutation test of the Cox
model, non-parametric bootstrap 95% confidence intervals (1,000 resamples;
`set.seed(20260508)`), leave-one-patient-out influence analysis,
leave-one-gene-out signature stability, and alternative scoring methods
(singscore, GSVA/ssGSEA). **[NOTE: the canonical TCGA script references
CIBERSORTx-based scores alongside BayesPrism; the manuscript should confirm which
deconvolution estimate is reported in the final multivariate model — both appear
in the code path.]**

## Statistical analysis

All statistical tests were two-sided unless otherwise specified. Multiple-testing
correction used the Benjamini-Hochberg method for FDR control or Bonferroni
correction, as indicated. Spatial permutation inference used 999 permutations
(Lee's L, Moran's I, LISA/BiLISA). Non-parametric tests (Wilcoxon rank-sum,
Wilcoxon signed-rank, Kruskal-Wallis, Mann-Whitney U) were used throughout owing
to the non-normal distributions typical of single-cell and spatial data. Effect
sizes included Cramér's V for contingency tables, median differences for pairwise
comparisons, and hazard ratios with 95% confidence intervals for survival
analyses. Across the cohort, 93.8% of patients with ≥ 100 epithelial cells had
representation from SecA, Intermediate and SecB states.

### Reproducibility

All stochastic steps in the re-runnable portion of the pipeline (HVG/NMF
subsampling and initialization, Scrublet, CopyKAT, k-means, GAM subsampling,
bootstrap and permutation procedures) were seeded from a single global constant
(`SEED = 42`, `config/config.py`), except where a different fixed analytical seed
is preserved as published (LIANA `seed = 1337`; TCGA bootstrap
`set.seed(20260508)`). Scrublet, CopyKAT and k-means seeds were added during the
reproducibility-hardening refactor where the original notebooks left them
implicit. The scVI/scANVI integration code is included in the repository
(`atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`) but, being an expensive GPU
job, was not re-run as part of routine downstream reproduction; the canonical cohort comprises the 8 published
whole HGSC tissues, the 97-patient TMA, and the 15 FTE cores (13 of which enter
the spatial-clustering/polarization panels, per the FTE note above).

## Software and versions

Analyses were performed in **Python 3.12.3** and **R 4.5.2**.

**Python stack:** scanpy 1.11.4, anndata 0.12.2, numpy 1.26.4, numba 0.59,
scikit-learn, statsmodels, scikit-misc (Seurat-v3 HVG), leidenalg/python-igraph,
umap-learn, decoupler (PROGENy/Hallmark/DoRothEA, MLM), liana ≥ 1.6.0
(cell–cell communication), lifelines (survival), mygene 3.2.2 (gene-symbol
harmonization), scrublet 0.2.3, and PyTorch (scFEA metabolic flux, in a dedicated
sub-environment). numpy < 2 and numba 0.59 are pinned for compatibility with
decoupler's JIT code paths.

**R/Bioconductor stack:** SingleR 2.10.0, UCell 2.12.0, mgcv 1.9-4, spdep 1.4-2,
spatstat 3.6-1 (spatstat.explore 3.8-1), sf 1.1-1, dbscan 1.2.4, RANN 2.6.2,
FNN, lme4, survival 3.8-6, survminer 0.5.2, copykat 1.1.0, BayesPrism 2.2.3,
consensusOV 1.30.0, spacexr 2.2.1 (RCTD), SpatialFeatureExperiment 1.10.1,
Voyager 1.10.0, scater/scran 1.36.0, ComplexHeatmap 2.24.1, data.table 1.18.4,
arrow 24.0.0, tidyverse 2.0.0, ggplot2 4.0.3, singscore and GSVA (robustness).
Integration (scVI/scANVI) was run with scvi-tools on an NVIDIA
V100 GPU (CUDA 12.6.2); the integration code is included in the repository
(`atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`). Organoid scRNA-seq was processed with Parse Biosciences
split-pipe 1.6.4 and Seurat v5; PDO doublets were removed with scDblFinder
1.20.2.

## Data and code availability

All scRNA-seq datasets used in atlas construction are publicly available under
the accessions listed in Supplemental Data 1 (per-cell metadata) and the study
list (Figure 1B / Supplemental Figure 2). Xenium spatial transcriptomics data,
processed atlas objects, and analysis code (`hgsc-malignant-states`;
https://www.github.com/cook-lab/hgsc-malignant-states) will be deposited upon
publication.

---

## Appendix: discrepancies and items flagged for the author

1. **GAM smoothing dimension** — draft: pooled `k = 20`, per-sample `k = 10`;
   code: pooled `k = 10`, per-sample `k = 5` (`spatial/05_gradients_gams/01_niche_succession_gams.R`).
   The code values are used above. Confirm intended values.
2. **Neighborhood k-means tuning** — draft: `nstart = 25`, `iter.max = 100`;
   code: `nstart = 50`, `iter.max = 200`
   (`spatial/04_neighborhood/02_neighborhood_k10_production.R`). Code values used.
3. **FTE cohort n = 13 vs 15** — reconciled in code: 15 cores profiled; 2 FTE
   whole-tissue samples (FT1-1, EAOC-1-FTE) excluded from the canonical
   re-runnable cohort as post-preprint additions, so 13 enter clustering/
   polarization panels. Confirm legends are updated consistently.
4. **CopyKAT sample n** — 251 samples had CopyKAT run; 248 had sufficient cells
   for within-clone epitype analysis (3,735 clones). Both numbers are correct in
   context; ensure figures cite the right one.
5. **Deconvolution method in TCGA multivariate model** — both BayesPrism and
   CIBERSORTx scores appear in `spatial/09_external_validation/01_tcga_external_validation.R`;
   confirm which estimate is reported in the final adjusted Cox model. (Atlas
   side also includes `atlas/07_deconvolution_survival/01_cibersort_reference.py`.)
6. **TCGA high-polarization bin n** — manuscript text says n = 105 while the
   Fig 7E/F legend/panel says n = 106 (low = 104). The KM figure binning differs
   from the tertile description in the original Methods; the tertile procedure is
   described above, but the figure bins should be reconciled with the reported n.
7. **Could not fully determine from code (flag for author):**
   (a) exact `n_perms` value passed to LIANA at production runtime (it is a CLI
   argument; the draft cites `n_perms = 100`); (b) some scVI/scANVI integration
   hyperparameters in `atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py` are
   taken from the draft text / scvi-tools defaults (marked CONFIRM in that script)
   because the original cluster script was not in the working copy; (c) precise scFEA
   epochs/learning-rate (epochs = 100, lr = 0.008) are from the draft, run in a
   vendored sub-environment not in the main `environment.yml`.
