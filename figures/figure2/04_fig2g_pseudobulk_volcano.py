#!/usr/bin/env python3
"""
Figure 2g — Pseudobulk volcano, SecA vs SecB (paired within tumours)
====================================================================
PURPOSE
    Wide horizontal volcano of the SecA-vs-SecB differential expression computed
    on PSEUDOBULK profiles (raw counts summed per sample_id x epitype), with
    tumour SAMPLES as the replicate unit (SecA and SecB paired within each tumour)
    and DESeq2 (pydeseq2) for testing. The 7 + 7 SecA/SecB scoring signature genes
    (shared/signatures.yml) are drawn large / outlined and labelled; top DEGs per
    side highlighted.

REPRODUCIBILITY (recorded from the deposited caches, 2026 final atlas):
    - Epitype = epithelial state (NMF Factor_2 percentile classes) and sample_id
      from figures/_prep/atlas_final_tables.load_epithelial() (final atlas, with
      the Nath 2021 sample mapping fix).
    - Inclusion: SecA/SecB epithelial cells; a (sample x epitype) group is used
      only with >= MIN_CELLS_PER_GROUP (=10) cells; a sample is used only if it has
      BOTH epitypes at that threshold (paired presence).
    - FINAL n: 261 tumours with both epitypes -> 522 pseudobulk profiles
      (261 SecA + 261 SecB). Gene filter: total pseudobulk count >= 10.
    - DESeq2 design: ~sample + epitype (paired: each tumour contributes one SecA and one SecB profile, so the tumour
      is a blocking factor and SecB is compared with SecA within tumours); contrast = (epitype, SecB, SecA) so
      log2FC > 0 is up in SecB. Seed from config. Final run: 3,601 DE genes at padj < 0.05 and |log2FC| > 0.5;
      all 14 signature genes significant in the expected direction.
    (Re-run prints the realized counts; assert guards the expected n.)

INPUTS
    - obj-style epithelial counts:
        data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
        (layer "counts" = raw integer counts ONLY; its obs metadata are not used)
    - state + sample_id per cell: figures/_prep/atlas_final_tables.py (load_epithelial)
    - SecA/SecB 7-gene signatures: shared/signatures.yml

OUTPUTS
    - figures_dir/figure2/atlas_volcano_secA_secB_pseudobulk.{svg,png}
    - output_root/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv (full DE table; also read by
      tables/06_tableS6_secA_secB_signatures.py and the Fig 6h atlas row)
    - printed: tumours / profiles, DE genes and signature genes significant

MANUSCRIPT PANEL(S): Fig 2g.

RUNTIME TIER: moderate (one chunked pass over the epithelial counts layer, then
    a DESeq2 fit on 522 x ~14.5k genes with 262 coefficients; about 15 minutes, no GPU).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ---------- Paths ----------
EPI_H5AD = path("data_root", "2026_final_atlas", "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
SIG_YML = Path(__file__).resolve().parents[2] / "shared" / "signatures.yml"
OUT_SVG = path("figures_dir", "figure2", "atlas_volcano_secA_secB_pseudobulk.svg")
OUT_PNG = path("figures_dir", "figure2", "atlas_volcano_secA_secB_pseudobulk.png")
OUT_DE = path("output_root", "figure2", "atlas_volcano_secA_secB_pseudobulk_deseq2.csv")

# ---------- Signatures ----------
SIG = yaml.safe_load(open(SIG_YML))
SECA_MARKERS, SECB_MARKERS = SIG["SecA"], SIG["SecB"]
ALL_MARKERS = set(SECA_MARKERS + SECB_MARKERS)

MIN_CELLS_PER_GROUP = 10
GENE_MIN_TOTAL = 10
EXPECTED_SAMPLES = 261   # recorded; asserted below (guards against silent cache drift)

# ---------- Style ----------
# Standardized fonts: base/axis-titles/legend/annotations = 8pt (FA/FN);
# tick (axis-type) labels = 6pt (FK).
FA, FK, FN = 8, 6, 8
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "legend.fontsize": FN, "figure.dpi": 450, "savefig.dpi": 450, "pdf.fonttype": 42,
    "svg.fonttype": "none", "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False,
})
SECA_COLOR, SECB_COLOR, NS_COLOR = "#E6A141", "#9A7D55", "#D9D9D9"
N_TOP = 10
LFC_THRESH, PADJ_THRESH = 0.5, 0.05

# ============================================================================
# 1. Build pseudobulk (sample_id x epitype) raw-count matrix
# ============================================================================
print("[1] Per-cell epitype (state) + sample labels (load_epithelial) ...")
meta = load_epithelial(columns=["sample_id", "state"])
adata = ad.read_h5ad(EPI_H5AD, backed="r")
obs = adata.obs
assert obs.index.isin(meta.index).all(), "epithelial h5ad cells missing from load_epithelial()"
epitype = meta["state"].astype(str).reindex(obs.index).values
sample = meta["sample_id"].astype(str).reindex(obs.index).values

df = pd.DataFrame({"sample": sample, "epitype": epitype}, index=np.arange(len(obs)))
df = df[df["epitype"].isin(["SecA", "SecB"])].copy()

gsize = df.groupby(["sample", "epitype"]).size()
ok = gsize[gsize >= MIN_CELLS_PER_GROUP].reset_index(name="n")
both = ok.groupby("sample")["epitype"].nunique()
paired_samples = both[both == 2].index
ok = ok[ok["sample"].isin(paired_samples)]
group_ids = (ok["sample"] + "|" + ok["epitype"]).tolist()
gidx = {g: i for i, g in enumerate(group_ids)}
n_groups = len(group_ids)
print(f"    tumours with both epitypes (>= {MIN_CELLS_PER_GROUP} cells each): {len(paired_samples)}")
print(f"    pseudobulk profiles (samples x epitype): {n_groups}")
assert len(paired_samples) == EXPECTED_SAMPLES, (
    f"pseudobulk sample count {len(paired_samples)} != recorded {EXPECTED_SAMPLES}; "
    "inputs changed — update the header/EXPECTED_SAMPLES after review.")

df["group"] = df["sample"] + "|" + df["epitype"]
df = df[df["group"].isin(gidx)]
row_group = np.full(len(obs), -1, dtype=np.int64)
row_group[df.index.values] = df["group"].map(gidx).values

print("[2] Summing raw counts per group (chunked backed read) ...")
counts = adata.layers["counts"]
n_cells, n_genes = adata.shape
acc = np.zeros((n_groups, n_genes), dtype=np.float64)
CHUNK = 50000
for i in range(0, n_cells, CHUNK):
    j = min(i + CHUNK, n_cells)
    codes = row_group[i:j]
    sel = np.where(codes >= 0)[0]
    if sel.size == 0:
        continue
    sub = counts[i:j]
    sub = sub.tocsr() if sp.issparse(sub) else sp.csr_matrix(sub)
    sub = sub[sel]
    ind = sp.csr_matrix((np.ones(sel.size), (codes[sel], np.arange(sel.size))),
                        shape=(n_groups, sel.size))
    acc += (ind @ sub).toarray()

genes = np.asarray(adata.var_names)
pb = pd.DataFrame(np.rint(acc).astype(int), index=group_ids, columns=genes)
pmeta = pd.DataFrame({"epitype": pd.Categorical([g.split("|")[1] for g in group_ids],
                                                categories=["SecA", "SecB"]),
                      "sample": ["s" + g.split("|")[0] for g in group_ids]}, index=group_ids)
pb = pb.loc[:, pb.sum(axis=0) >= GENE_MIN_TOTAL]
print(f"    pseudobulk matrix: {pb.shape[0]} profiles x {pb.shape[1]} genes")

# ============================================================================
# 2. DESeq2 (pydeseq2): design ~sample + epitype (paired), contrast SecB vs SecA
# ============================================================================
FORMULA = "~sample + epitype"   # paired within tumours
print(f"[3] DESeq2 ({FORMULA}; contrast SecB vs SecA) ...")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

np.random.seed(SEED)
try:
    dds = DeseqDataSet(counts=pb, metadata=pmeta, design=FORMULA, quiet=True)
except TypeError:
    dds = DeseqDataSet(counts=pb, metadata=pmeta, design_factors=["sample", "epitype"], quiet=True)
dds.deseq2()
res = DeseqStats(dds, contrast=["epitype", "SecB", "SecA"], quiet=True)
res.summary()
de = res.results_df.copy()
de.index.name = "gene"
de = de.reset_index()
os.makedirs(os.path.dirname(OUT_DE), exist_ok=True)
de.to_csv(OUT_DE, index=False)

# signature check (recorded in the run log)
sig_ok = 0
for lab, gl in [("SecA", SECA_MARKERS), ("SecB", SECB_MARKERS)]:
    for g in gl:
        r = de[de["gene"] == g]
        if len(r) and (r["padj"].iloc[0] < PADJ_THRESH) and \
           (np.sign(r["log2FoldChange"].iloc[0]) == (1 if lab == "SecB" else -1)):
            sig_ok += 1
print(f"    noBCAM signature genes significant & correct-direction: {sig_ok}/{len(ALL_MARKERS)}")

# ============================================================================
# 3. Volcano
# ============================================================================
d = de.dropna(subset=["padj", "log2FoldChange"]).copy()
d["neg_log10p"] = -np.log10(d["padj"].clip(lower=1e-300))
YMAX = float(np.ceil(np.nanpercentile(d["neg_log10p"], 99.5)))
XMIN, XMAX = -6, 6
d["lfc_plot"] = d["log2FoldChange"].clip(XMIN + 0.1, XMAX - 0.1)
d["p_plot"] = d["neg_log10p"].clip(upper=YMAX)
d["category"] = "ns"
d.loc[(d["log2FoldChange"] < -LFC_THRESH) & (d["padj"] < PADJ_THRESH), "category"] = "SecA"
d.loc[(d["log2FoldChange"] > LFC_THRESH) & (d["padj"] < PADJ_THRESH), "category"] = "SecB"
d["is_marker"] = d["gene"].isin(ALL_MARKERS)
n_sig = int((d["category"] != "ns").sum())
print(f"    DE genes (padj<{PADJ_THRESH}, |LFC|>{LFC_THRESH}): {n_sig}")

seca_sig = d[(d.category == "SecA") & ~d.is_marker]
secb_sig = d[(d.category == "SecB") & ~d.is_marker]
top_seca = seca_sig.nlargest(N_TOP, "neg_log10p")["gene"].tolist()
top_secb = secb_sig.nlargest(N_TOP, "neg_log10p")["gene"].tolist()
top_all = set(top_seca + top_secb)

fig, ax = plt.subplots(figsize=(3.4, 2.9))
ns = d[d.category == "ns"]
ax.scatter(ns.lfc_plot, ns.p_plot, s=1.5, c=NS_COLOR, alpha=0.3, linewidths=0, rasterized=True, zorder=1)
for cat, col in [("SecA", SECA_COLOR), ("SecB", SECB_COLOR)]:
    bg = d[(d.category == cat) & ~d.is_marker & ~d.gene.isin(top_all)]
    ax.scatter(bg.lfc_plot, bg.p_plot, s=3, c=col, alpha=0.4, linewidths=0, rasterized=True, zorder=2)
    tp = d[d.gene.isin(top_seca if cat == "SecA" else top_secb)]
    ax.scatter(tp.lfc_plot, tp.p_plot, s=8, c=col, alpha=0.6, linewidths=0, rasterized=True, zorder=3)
for gl, col in [(SECA_MARKERS, SECA_COLOR), (SECB_MARKERS, SECB_COLOR)]:
    mm = d[d.gene.isin(gl)]
    ax.scatter(mm.lfc_plot, mm.p_plot, s=22, c=col, edgecolors="black", linewidths=0.4, zorder=5)

ax.axhline(-np.log10(PADJ_THRESH), color="grey", lw=0.4, ls="--", zorder=0)
ax.axvline(-LFC_THRESH, color="grey", lw=0.4, ls="--", zorder=0)
ax.axvline(LFC_THRESH, color="grey", lw=0.4, ls="--", zorder=0)
ax.axhline(YMAX, color="grey", lw=0.3, ls=":", alpha=0.4, zorder=0)
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(-3, YMAX + 14)  # extra headroom so repelled marker labels clear the y-cap

# SecA marker labels: anchored in a clean staggered column in the LEFT margin,
# left-aligned, each tied to its point by a thin leader line (keeps the centre of
# the volcano clear). Ordered top->bottom by point height to avoid leader crossings.
seca_rows = []
for gene in SECA_MARKERS:
    r = d[d.gene == gene]
    if len(r):
        r = r.iloc[0]
        seca_rows.append((gene, float(r.lfc_plot), float(r.p_plot)))
seca_rows.sort(key=lambda t: t[2], reverse=True)
label_x = XMIN + 0.15  # just inside the left edge
col_ys = np.linspace(YMAX - 3, 5, len(seca_rows)) if len(seca_rows) > 1 \
    else [0.5 * YMAX]
for (gene, gx, gy), ly in zip(seca_rows, col_ys):
    ax.annotate(gene, xy=(gx, gy), xytext=(label_x, ly), fontsize=FA,
                fontweight="bold", fontstyle="italic", color="#8B6914",
                ha="left", va="center", zorder=7,
                arrowprops=dict(arrowstyle="-", color="#999999", linewidth=0.3,
                                shrinkA=1, shrinkB=3))

# SecB marker labels: repelled in place (upper-right), leader lines to points.
texts = []
for gene in SECB_MARKERS:
    r = d[d.gene == gene]
    if len(r):
        r = r.iloc[0]
        texts.append(ax.text(r.lfc_plot, r.p_plot, gene, fontsize=FA, fontweight="bold",
                             fontstyle="italic", color="#5C4A2A", ha="center", va="center", zorder=7))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#999999", linewidth=0.3, shrinkA=1, shrinkB=3),
            expand=(1.6, 1.9), force_text=(0.6, 1.1), force_static=(0.35, 0.6),
            force_pull=(0.02, 0.02), only_move="xy", ensure_inside_axes=True,
            min_arrow_len=2, iter_lim=1500, time_lim=20)

ax.set_xlabel("log$_2$ fold change (SecA ← → SecB)", fontsize=FA)
ax.set_ylabel("$-$log$_{10}$ adj. p-value", fontsize=FA)
ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(width=0.5, length=2)
# No in-plot SecA/SecB marker legend: the legend is added in the figure layout (Affinity).

fig.tight_layout()
fig.savefig(OUT_SVG, format="svg")
fig.savefig(OUT_PNG, format="png")
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}\n  Saved: {OUT_DE}")
print(f"  [pseudobulk] {len(paired_samples)} tumours -> {n_groups} profiles; "
      f"{n_sig} DE genes; signature {sig_ok}/{len(ALL_MARKERS)}")
