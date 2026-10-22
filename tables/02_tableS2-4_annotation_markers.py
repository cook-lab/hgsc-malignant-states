#!/usr/bin/env python3
"""Supplementary Tables S2-S4: the marker/DEG tables produced DURING atlas annotation (author decision 2026-09-28).

PURPOSE / INPUTS
  Table S2  level-1 lineage markers  : 202601_final_atlas/output/06_celltype_level1/markers/level1_<type>.csv
                                       (scanpy rank_genes_groups, Wilcoxon, each lineage vs all other cells; 51,592 pairs)
  Table S3  level-2 markers           : 202601_final_atlas/output/08_celltype_level2/markers/<compartment>_<type>.csv
                                       (scripts/08_finalize_celltype_level2.py: Wilcoxon within each level-1 compartment,
                                       each level-2 cluster vs the other clusters of that compartment; padj < 0.05, logFC > 0)
                                       -> restricted to the 51 final level-2 types (the 'Excluded' clusters are dropped)
  Table S4  ciliated vs secretory     : 202601_final_atlas/output/supplemental/T3_epithelial_markeres.csv (10,029 DEGs)
  hgsc_atlas_final.h5ad (obs celltype_level1 / celltype_level2, via h5py): the 51 final level-2 types and their cell counts
OUTPUTS (supp_xlsx.write_table)
  SuppTable2_atlas_level1_markers (Top50, All_markers); SuppTable3_atlas_level2_markers_within_level1 (Top50,
  All_markers, Cell_types); SuppTable4_atlas_ciliated_vs_secretory_markers (All_DEGs, Canonical_markers)
MANUSCRIPT ITEM(S)
  Tables S2, S3, S4; Results "51 discrete cell types and states" (asserted: 51 level-2 types).
RUNTIME
  fast (reads the marker CSVs and the h5ad obs).
"""
import os, re, sys, glob, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from supp_xlsx import write_table, OUT
import h5py

A = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas")
L1_ORDER = ["Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial", "T/NK cell",
            "B cell", "Plasma cell", "Macrophage", "DC", "Neutrophil", "Mast cell"]
COLS = {"names": "gene", "scores": "wilcoxon_score", "logfoldchanges": "log2FC", "pvals": "pval", "pvals_adj": "pval_adj"}

# final level-1 / level-2 labels and cell counts (final atlas, all cells)
h = h5py.File(f"{A}/hgsc_atlas_final.h5ad", "r")
def cat(k):
    g = h["obs"][k]; c = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
    return pd.Categorical.from_codes(g["codes"][:], c)
obs = pd.DataFrame({"celltype_level1": cat("celltype_level1"), "celltype_level2": cat("celltype_level2")})
L2 = obs.groupby(["celltype_level1", "celltype_level2"], observed=True).size().reset_index(name="n_cells_final_atlas")
assert len(L2) == 51 and L2.n_cells_final_atlas.sum() == 1980703

def top_wide(df, keys, n=50):
    cols = {}
    for k, g in df.groupby(keys, sort=False):
        cols[k] = g.sort_values("wilcoxon_score", ascending=False).gene.head(n).tolist()
    return cols

# ---------------- Table 2 ----------------
T2 = []
for L in L1_ORDER:
    d = pd.read_csv(f"{A}/output/06_celltype_level1/markers/level1_{L.replace(' ', '_').replace('/', '-')}.csv").rename(columns=COLS)
    d = d[(d.pval_adj < 0.05) & (d.log2FC > 0)].sort_values("wilcoxon_score", ascending=False)
    d.insert(0, "celltype_level1", L); d.insert(1, "rank", range(1, len(d) + 1)); T2.append(d)
T2 = pd.concat(T2, ignore_index=True)
assert len(T2) == 51592, len(T2)
May = pd.read_csv(f"{A}/output/supplemental/T2_celltypelevel1_markers.csv.csv")
assert len(May) == len(T2) and set(zip(May.celltype_level1, May.gene)) == set(zip(T2.celltype_level1, T2.gene))
w = top_wide(T2, "celltype_level1")
top2 = pd.DataFrame({L: pd.Series(w[L]) for L in L1_ORDER}); top2.insert(0, "rank", range(1, len(top2) + 1))
write_table(2, "atlas_level1_markers",
    "Supplementary Table 2. Differentially expressed markers for cell type level 1 annotation.",
    ["Wilcoxon rank-sum test (scanpy rank_genes_groups, method = 'wilcoxon') of each of the 13 level-1 cell type compartments versus all other cells, "
     "run at the level-1 annotation step of the atlas (Leiden clusters of the scANVI-integrated atlas). Upregulated genes with Benjamini-Hochberg adjusted p < 0.05 "
     f"and log2 fold change > 0 are reported: {len(T2):,} significant gene-cell type pairs.",
     "Columns: celltype_level1, rank (by Wilcoxon score within cell type), gene, wilcoxon_score (z statistic), log2FC, pval, pval_adj (Benjamini-Hochberg).",
     "Sheets: Top50 (row 1 = cell type; the 50 highest-scoring markers per cell type), All_markers (all significant pairs).",
     "Source: 202601_final_atlas/output/06_celltype_level1/markers/level1_<cell type>.csv (identical gene-cell type pairs to the May 2026 Supplemental Data 2).",
     "Build: tables/02_tableS2-4_annotation_markers.py"],
    [("Top50", top2), ("All_markers", T2)])

# ---------------- Table 3 ----------------
norm = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())
files = {norm(os.path.basename(f)[:-4].split("_", 1)[1]): f for f in glob.glob(f"{A}/output/08_celltype_level2/markers/*.csv")}
T3, single = [], []
for L in L1_ORDER:
    for _, r in L2[L2.celltype_level1 == L].sort_values("celltype_level2").iterrows():
        d = pd.read_csv(files[norm(r.celltype_level2)])
        if "scores" not in d:            # compartment with a single level-2 type: no within-compartment DE was possible
            single.append(r.celltype_level2); continue
        d = d.rename(columns=COLS)
        d = d[(d.pval_adj < 0.05) & (d.log2FC > 0)].sort_values("wilcoxon_score", ascending=False)
        d.insert(0, "celltype_level1", L); d.insert(1, "celltype_level2", r.celltype_level2); d.insert(2, "rank", range(1, len(d) + 1))
        T3.append(d)
T3 = pd.concat(T3, ignore_index=True)
order = [(L, t) for L in L1_ORDER for t in sorted(L2[L2.celltype_level1 == L].celltype_level2) if t not in single]
w = top_wide(T3, ["celltype_level1", "celltype_level2"])
top3 = pd.DataFrame({f"{L} | {t}": pd.Series(w[(L, t)]) for L, t in order}); top3.insert(0, "rank", range(1, len(top3) + 1))
write_table(3, "atlas_level2_markers_within_level1",
    "Supplementary Table 3. Differentially expressed markers for cell type level 2 annotation within each level 1 compartment.",
    ["Wilcoxon rank-sum test (scanpy rank_genes_groups) run within each level-1 compartment at the level-2 annotation step: each level-2 cluster versus the other "
     "clusters of the same compartment (scripts/08_finalize_celltype_level2.py). Upregulated genes with Benjamini-Hochberg adjusted p < 0.05 and log2 fold change > 0 "
     f"are reported for the {len(order)} final level-2 cell types with a within-compartment comparison: {len(T3):,} significant gene-cell type pairs.",
     f"Compartments with a single level-2 type ({', '.join(single)}) have no within-compartment comparison; their markers are the level-1 markers in Supplementary Table 2. "
     "Clusters excluded from the final atlas (doublet/low-quality 'Excluded' clusters) are not reported, although they formed part of the comparison background.",
     "Columns: celltype_level1, celltype_level2, rank (by Wilcoxon score), gene, wilcoxon_score, log2FC, pval, pval_adj (Benjamini-Hochberg).",
     "Sheets: Top50 (column header = 'level 1 | level 2'; the 50 highest-scoring markers per level-2 type), All_markers, Cell_types (final level-2 types with cell counts in the final atlas).",
     "Source: 202601_final_atlas/output/08_celltype_level2/markers/<compartment>_<cell type>.csv. Build: tables/02_tableS2-4_annotation_markers.py"],
    [("Top50", top3), ("All_markers", T3), ("Cell_types", L2)])

# ---------------- Table 4 ----------------
T4 = pd.read_csv(f"{A}/output/supplemental/T3_epithelial_markeres.csv")
assert len(T4) == 10029
T4 = T4.sort_values(["direction", "wilcoxon_score"], ascending=[True, False])
for c in ["pct_ciliated", "pct_secretory"]: T4[c] = (100 * T4[c]).round(2)
CAN = ["FOXJ1", "CAPS", "RSPH1", "TPPP3", "PAX8", "MUC16", "WFDC2"]
write_table(4, "atlas_ciliated_vs_secretory_markers",
    "Supplementary Table 4. Differentially expressed markers for ciliated and secretory clusters.",
    ["Ciliated and secretory epithelial cell clusters of the atlas epithelial compartment were compared by Wilcoxon rank-sum test. "
     f"Genes with Benjamini-Hochberg adjusted p < 0.05 are reported ({len(T4):,} DEGs: "
     f"{(T4.direction.str.startswith('Cil')).sum():,} ciliated-up, {(~T4.direction.str.startswith('Cil')).sum():,} secretory-up).",
     "Columns: gene, log2FC (ciliated vs secretory), wilcoxon_score, pval, pval_adj (Benjamini-Hochberg), direction of enrichment, pct_ciliated / pct_secretory "
     "(% of cells expressing the gene in each group).",
     "Sheets: All_DEGs, Canonical_markers (the ciliated and secretory markers named in the text).",
     "Source: 202601_final_atlas/output/supplemental/T3_epithelial_markeres.csv (May 2026 Supplemental Data 3). Build: tables/02_tableS2-4_annotation_markers.py"],
    [("All_DEGs", T4), ("Canonical_markers", T4[T4.gene.isin(CAN)].set_index("gene").loc[[g for g in CAN if g in set(T4.gene)]].reset_index())])
print("T2", len(T2), "T3", len(T3), "single", single, "T4", len(T4))
