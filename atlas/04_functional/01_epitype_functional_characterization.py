#!/usr/bin/env python3
"""
Epitype functional characterization (PROGENy / Hallmark / DoRothEA / flux)
==========================================================================
HGSC malignant-states atlas backend.

Freshly computes PROGENy pathway activity, Hallmark gene-set scores, DoRothEA
TF activity, and scFEA metabolic flux for the 4 NMF-based epithelial epitypes
(SecA / Intermediate / SecB / Ciliated), then renders z-scored heatmaps,
violin grids, and radar plots.

Label assignment (NMF Factor_2 percentile schema; thresholds on non-ciliated epi):
  - SecA:         Factor_2 < p50
  - Intermediate: p50 <= Factor_2 < p75   (was "Transitioning" in originals)
  - SecB:         Factor_2 >= p75
  - Ciliated:     celltype_level2 == "Ciliated epithelial cell"

INPUTS:
  - <data_root>/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
  - output_root/03_epithelial_nmf/11d_nmf_usage.csv   (Factor_2)
  - <data_root>/2026_final_atlas/tools/scFEA/data/*   (M168 metabolic model)

OUTPUTS (output_root/04_functional/):
  - 21_{progeny,hallmark,dorothea,flux}_zscored.csv  (group z-scores; KEY caches)
  - 21_*_means.csv / _radar.csv / _per_cell.parquet
  - 21_*_heatmap / _violins / _radar_* .svg/pdf

MANUSCRIPT PANELS: Fig 2A (radars), SF8 (functional heatmaps).

RUNTIME TIER: heavy. NOTE: the metabolic-flux step trains a PyTorch graph
  neural network (scFEA, Chang et al. 2021); requires `torch` (CPU is fine,
  ~20k cells x 168 modules). decoupler is required for PROGENy/DoRothEA.

SEEDING: np.random.seed(SEED) + RandomState(SEED) for the flux downsample;
  torch.manual_seed(SEED) before scFEA training (determinism).

scFEA FLUX LABEL-MAP FIX:
  The original keyed the readable label map by the integer `Module_id`
  column (1,2,...), but the flux matrix columns are the string module index
  (M_1, M_2, ...). The keys never matched, so heatmaps/radars showed raw
  "M_x" instead of "Compound_IN->Compound_OUT". Fixed here by reading the
  annotation with the module index (M_x) as the key so labels match.

Usage:
    python 01_epitype_functional_characterization.py
"""

import os
import sys
import time
import gc
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import seaborn as sns

import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

H5AD       = path("data_root", "2026_final_atlas", "celltype_h5ad",
                  "hgsc_atlas_final_epithelial.h5ad")
NMF_CSV    = path("output_root", "03_epithelial_nmf", "11d_nmf_usage.csv")
OUT_DIR    = path("output_root", "04_functional")
SCFEA_DATA = path("data_root", "2026_final_atlas", "tools", "scFEA", "data")
FLUX_ANNOT = os.path.join(SCFEA_DATA, "Human_M168_information.symbols.csv")
os.makedirs(OUT_DIR, exist_ok=True)

# scFEA hyperparameters (validated — do not change)
DOWNSAMPLE_N = 20_000
EPOCH        = 100
LEARN_RATE   = 0.008
LAMB_BA      = 1
LAMB_NG      = 1
LAMB_CELL    = 1
LAMB_MOD     = 1e-2
RNG          = np.random.RandomState(SEED)

# ============================================================================
# COOK LAB STYLE v1.2
# ============================================================================

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         8,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   7,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ============================================================================
# GROUPS & PALETTE
# ============================================================================

ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]
PALETTE = {
    "SecA":         "#E6A141",
    "Intermediate": "#C08E48",
    "SecB":         "#9A7D55",
    "Ciliated":     "#E07850",
}


# ============================================================================
# HELPERS
# ============================================================================

def zscore_df(df):
    mu = df.mean(axis=0)
    sd = df.std(axis=0).replace(0, 1)
    return (df - mu) / sd


def compute_group_means(scores_df, group_labels, feature_cols=None):
    if feature_cols is None:
        feature_cols = scores_df.columns.tolist()
    numeric_cols = [c for c in feature_cols
                    if scores_df[c].dtype in ("float64", "float32", "int64", "int32")]
    df = scores_df[numeric_cols].copy()
    df["group"] = group_labels
    return df.groupby("group").mean(numeric_only=True).reindex(ORDER)


def select_top_differential(group_means, n=8):
    ranges = group_means.max(axis=0) - group_means.min(axis=0)
    return ranges.nlargest(n).index.tolist()


def save_figs(fig, fname):
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(OUT_DIR, f"{fname}.{ext}"), format=ext)
    plt.close(fig)


# ============================================================================
# HEATMAP
# ============================================================================

def plot_full_heatmap(group_means_z, title, fname, figwidth=None):
    df = group_means_z.copy()
    n_features = df.shape[1]
    if n_features > 2:
        link = linkage(pdist(df.T.values, metric="euclidean"), method="ward")
        df = df.iloc[:, leaves_list(link)]
    if figwidth is None:
        figwidth = max(6, 0.35 * n_features + 1.5)
    figheight = max(2.5, 0.6 * len(ORDER) + 1.5)
    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    im = ax.imshow(df.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(n_features))
    ax.set_xticklabels(df.columns, rotation=90, ha="center", fontsize=6)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER, fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Z-score", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10)
    fig.tight_layout()
    save_figs(fig, fname)
    print(f"    Saved {fname}.svg/pdf ({n_features} features)")


# ============================================================================
# VIOLIN GRID
# ============================================================================

def plot_full_violins(per_cell_df, group_labels, features, title, fname,
                      max_features=30, ncols=5):
    if len(features) > max_features:
        tmp = per_cell_df[features].copy()
        tmp["group"] = group_labels
        gm = tmp.groupby("group").mean(numeric_only=True).reindex(ORDER)
        ranges = gm.max() - gm.min()
        features = ranges.nlargest(max_features).index.tolist()

    n = len(features)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.2, nrows * 1.8))
    axes = np.atleast_2d(axes)
    plot_df = per_cell_df[features].copy()
    plot_df["group"] = group_labels

    for i, feat in enumerate(features):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        sns.violinplot(data=plot_df, x="group", y=feat, order=ORDER,
                       palette=PALETTE, ax=ax, linewidth=0.5, cut=0,
                       inner="quartile", density_norm="width")
        ax.set_title(feat, fontsize=7, fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=5, rotation=45)
        ax.tick_params(axis="y", labelsize=5)
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].set_visible(False)

    fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)
    fig.tight_layout()
    save_figs(fig, fname)
    print(f"    Saved {fname}.svg/pdf ({len(features)} features)")


# ============================================================================
# RADAR
# ============================================================================

def radar_plot(ax, df_z, title):
    features = df_z.columns.tolist()
    n = len(features)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    for grp in ORDER:
        if grp not in df_z.index:
            continue
        vals = df_z.loc[grp].values.tolist()
        vals += vals[:1]
        ax.plot(angles, vals, linewidth=1.2, color=PALETTE[grp])
        ax.fill(angles, vals, alpha=0.12, color=PALETTE[grp])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(features, fontsize=6)
    ax.set_rlabel_position(30)
    ax.tick_params(axis="y", labelsize=6, pad=2)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.grid(linewidth=0.4, alpha=0.5)
    ax.spines["polar"].set_linewidth(0.4)
    ax.set_title(title, fontsize=9, fontweight="bold", pad=14)


# ============================================================================
# scFEA NETWORK (Chang et al. 2021, Genome Research)
# ============================================================================

class FLUX(nn.Module):
    """Graph neural network for metabolic flux estimation."""
    def __init__(self, matrix, n_modules, f_in=50, f_out=1):
        super(FLUX, self).__init__()
        self.inSize = f_in
        self.m_encoder = nn.ModuleList([
            nn.Sequential(
                nn.Linear(self.inSize, 8, bias=False), nn.Tanhshrink(),
                nn.Linear(8, f_out), nn.Tanhshrink(),
            )
            for _ in range(n_modules)
        ])

    def updateC(self, m, n_comps, cmMat):
        c = torch.zeros((m.shape[0], n_comps))
        for i in range(c.shape[1]):
            c[:, i] = torch.sum(m * cmMat[i, :], dim=1)
        return c

    def forward(self, x, n_modules, n_genes, n_comps, cmMat):
        for i in range(n_modules):
            x_block = x[:, i * n_genes: (i + 1) * n_genes]
            subnet = self.m_encoder[i]
            m = subnet(x_block) if i == 0 else torch.cat((m, subnet(x_block)), 1)
        c = self.updateC(m, n_comps, cmMat)
        return m, c


class FluxDataset(Dataset):
    def __init__(self, X, geneExprScale, module_scale):
        self.X = X
        self.geneExprScale = geneExprScale
        self.module_scale = module_scale

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx, :], self.geneExprScale[idx], self.module_scale[idx, :]


def pearsonr_torch(x, y):
    xm = x - torch.mean(x)
    ym = y - torch.mean(y)
    r_num = torch.sum(xm * ym)
    r_den = torch.sqrt(torch.sum(xm ** 2) * torch.sum(ym ** 2))
    return r_num / (r_den + 1e-8)


def scfea_loss(m, c, lamb1, lamb2, lamb3, lamb4, geneScale, moduleScale):
    total1 = torch.sum(torch.pow(c, 2), dim=1)
    error = torch.abs(m) - m
    total2 = torch.sum(error, dim=1)
    diff = torch.pow(torch.sum(m, dim=1) - geneScale, 2)
    if (diff > 0).sum() == m.shape[0]:
        total3 = torch.pow(diff, 0.5)
    else:
        total3 = diff
    if lamb4 > 0:
        corr = torch.FloatTensor(np.ones(m.shape[0]))
        for i in range(m.shape[0]):
            corr[i] = pearsonr_torch(m[i, :], moduleScale[i, :])
        total4 = torch.FloatTensor(np.ones(m.shape[0])) - torch.abs(corr)
    else:
        total4 = torch.FloatTensor(np.zeros(m.shape[0]))
    return (torch.sum(lamb1 * total1) + torch.sum(lamb2 * total2) +
            torch.sum(lamb3 * total3) + torch.sum(lamb4 * total4))


# ============================================================================
# AXIS SCORING + VISUALIZATION HELPER
# ============================================================================

def score_and_visualize(per_cell_scores, group_labels, axis_name, prefix):
    feats = per_cell_scores.columns.tolist()
    means = compute_group_means(per_cell_scores, group_labels)
    z = zscore_df(means)
    top8 = select_top_differential(means, n=8)
    radar_z = zscore_df(means[top8])

    means.to_csv(os.path.join(OUT_DIR, f"{prefix}_means.csv"))
    z.to_csv(os.path.join(OUT_DIR, f"{prefix}_zscored.csv"))
    radar_z.to_csv(os.path.join(OUT_DIR, f"{prefix}_radar.csv"))
    per_cell_scores.to_parquet(os.path.join(OUT_DIR, f"{prefix}_per_cell.parquet"))

    print(f"  Plotting {axis_name} heatmap ({len(feats)} features)...")
    plot_full_heatmap(z, f"{axis_name} (all {len(feats)} features)", f"{prefix}_heatmap")
    print(f"  Plotting {axis_name} violins...")
    plot_full_violins(per_cell_scores, group_labels, feats,
                      axis_name, f"{prefix}_violins", max_features=30)
    return radar_z


# ============================================================================
# MAIN
# ============================================================================

print("=" * 70)
print("Epitype Functional Characterization (Fresh Computation)")
print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

np.random.seed(SEED)

# ── STEP 1: Load h5ad & assign NMF-based epitype labels ──────
t0 = time.time()
print("\n[1/8] Loading epithelial h5ad (full, not backed)...")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes  ({time.time()-t0:.0f}s)")

xmax = adata.X.max() if not sparse.issparse(adata.X) else adata.X.data.max()
print(f"  X max = {xmax:.2f} (expected ~12.7 for log1p)")
if xmax > 50:
    print("  WARNING: X appears raw — normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

print("  Loading NMF Factor 2...")
nmf_usage = pd.read_csv(NMF_CSV, index_col=0)
f2 = nmf_usage["Factor_2"]
shared = adata.obs_names.intersection(f2.index)
print(f"  Shared cells: {len(shared):,}")

adata = adata[shared].copy()
f2 = f2.loc[shared]

is_ciliated = adata.obs["celltype_level2"] == "Ciliated epithelial cell"
non_cil_f2 = f2[~is_ciliated]
p50 = np.percentile(non_cil_f2, 50)
p75 = np.percentile(non_cil_f2, 75)
print(f"  Thresholds: p50={p50:.4f}, p75={p75:.4f}")

epitype = pd.Series("SecA", index=shared, name="epitype")
epitype[f2 >= p50] = "Intermediate"
epitype[f2 >= p75] = "SecB"
epitype[is_ciliated] = "Ciliated"
adata.obs["epitype"] = epitype.values
for g in ORDER:
    n = (epitype == g).sum()
    print(f"    {g}: {n:,} ({n/len(epitype)*100:.1f}%)")

group_labels = epitype.values
print(f"  Step 1 done ({time.time()-t0:.0f}s)\n")

# ── STEP 2: PROGENy ──────────────────────────────────────────
progeny_done = os.path.exists(os.path.join(OUT_DIR, "21_progeny_per_cell.parquet"))
if progeny_done:
    print("[2/8] PROGENy — loading existing results...")
    progeny = pd.read_parquet(os.path.join(OUT_DIR, "21_progeny_per_cell.parquet"))
    common_prog = progeny.index.intersection(pd.Index(adata.obs_names))
    progeny = progeny.loc[common_prog]
    prog_labels = adata.obs.loc[common_prog, "epitype"].values
    prog_means = compute_group_means(progeny, prog_labels)
    top8 = select_top_differential(prog_means, n=8)
    prog_radar_z = zscore_df(prog_means[top8])
    del progeny
    gc.collect()
else:
    print("[2/8] PROGENy pathway activity (decoupler MLM)...")
    t1 = time.time()
    progeny_net = dc.op.progeny(organism="human")
    print(f"  Network: {progeny_net.shape[0]:,} links, "
          f"{progeny_net['source'].nunique()} pathways")
    dc.mt.mlm(adata, net=progeny_net)
    progeny = adata.obsm["score_mlm"].copy()
    print(f"  Scored: {progeny.shape[1]} pathways x {progeny.shape[0]:,} cells  "
          f"({time.time()-t1:.0f}s)")
    prog_radar_z = score_and_visualize(progeny, group_labels,
                                       "PROGENy Pathway Activity", "21_progeny")
    del progeny
    gc.collect()

# ── STEP 3: Hallmark ─────────────────────────────────────────
print(f"\n[3/8] Hallmark gene sets (score_genes on {adata.n_obs:,} cells)...")
t2 = time.time()
hallmark_net = dc.op.hallmark(organism="human")
hall_sets = {}
for gs_name, grp in hallmark_net.groupby("source"):
    genes = [g for g in grp["target"].values if g in adata.var_names]
    if len(genes) >= 5:
        hall_sets[gs_name] = genes
print(f"  {len(hall_sets)} gene sets with >=5 overlapping genes")

hall_scores = {}
for gs_name, genes in hall_sets.items():
    key = f"hall_{gs_name}"
    sc.tl.score_genes(adata, gene_list=genes, score_name=key, use_raw=False)
    hall_scores[gs_name.replace("HALLMARK_", "")] = adata.obs[key].values

hallmark = pd.DataFrame(hall_scores, index=adata.obs_names)
print(f"  Scored: {hallmark.shape[1]} gene sets x {hallmark.shape[0]:,} cells  "
      f"({time.time()-t2:.0f}s)")
hall_radar_z = score_and_visualize(hallmark, group_labels, "Hallmark Gene Sets",
                                   "21_hallmark")
del hallmark
gc.collect()
for key in list(adata.obsm.keys()):
    if key.startswith("score_"):
        del adata.obsm[key]
gc.collect()

# ── STEP 4: DoRothEA ─────────────────────────────────────────
print(f"\n[4/8] DoRothEA TF activity (decoupler MLM, full {adata.n_obs:,} cells)...")
t3 = time.time()
dorothea_net = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
print(f"  Network: {dorothea_net.shape[0]:,} links, "
      f"{dorothea_net['source'].nunique()} TFs")
dc.mt.mlm(adata, net=dorothea_net)
dorothea = adata.obsm["score_mlm"].copy()
print(f"  Scored: {dorothea.shape[1]} TFs x {dorothea.shape[0]:,} cells  "
      f"({time.time()-t3:.0f}s)")
doro_radar_z = score_and_visualize(dorothea, group_labels, "DoRothEA TF Activity",
                                   "21_dorothea")
del dorothea
gc.collect()

# ── STEP 5: Metabolic flux via scFEA ─────────────────────────
print("\n[5/8] Metabolic flux (scFEA, fresh training)...")
t4 = time.time()

print(f"  Downsampling to {DOWNSAMPLE_N:,} cells (stratified)...")
labels_arr = adata.obs["epitype"].values
unique_cats, counts = np.unique(labels_arr, return_counts=True)
frac = DOWNSAMPLE_N / len(labels_arr)
indices = []
for cat, count in zip(unique_cats, counts):
    cat_idx = np.where(labels_arr == cat)[0]
    n_keep = max(1, int(round(count * frac)))
    indices.append(RNG.choice(cat_idx, size=min(n_keep, count), replace=False))
indices = np.concatenate(indices)
RNG.shuffle(indices)

adata_flux = adata[indices].copy()
flux_epitype = adata_flux.obs["epitype"].values
print(f"  Downsampled: {adata_flux.n_obs:,} cells")
for g in ORDER:
    print(f"    {g}: {(flux_epitype == g).sum():,}")

print("  Preparing expression matrix...")
X_flux = adata_flux.X
if sparse.issparse(X_flux):
    X_flux = X_flux.toarray()
geneExpr = pd.DataFrame(X_flux, index=adata_flux.obs_names, columns=adata_flux.var_names)
if geneExpr.max().max() > 50:
    print("  Values > 50 detected, applying log2(x+1)...")
    geneExpr = np.log2(geneExpr + 1)

del adata, adata_flux, X_flux
gc.collect()

print("  Loading scFEA metabolic model (168 modules)...")
moduleGene = pd.read_csv(os.path.join(SCFEA_DATA, "module_gene_m168.csv"), index_col=0)
cmMat = pd.read_csv(os.path.join(SCFEA_DATA, "cmMat_c70_m168.csv"), header=None).values
cName = pd.read_csv(os.path.join(SCFEA_DATA, "cName_c70_m168.csv")).columns

moduleLen = np.array([moduleGene.iloc[i, :].notna().sum()
                      for i in range(moduleGene.shape[0])])

module_gene_all = set()
for i in range(moduleGene.shape[0]):
    for j in range(moduleGene.shape[1]):
        if pd.notna(moduleGene.iloc[i, j]):
            module_gene_all.add(moduleGene.iloc[i, j])

gene_overlap = sorted(set(geneExpr.columns) & module_gene_all)
print(f"  Module genes: {len(module_gene_all)} | Data genes: {geneExpr.shape[1]}")
print(f"  Overlap: {len(gene_overlap)} genes")

geneExpr = geneExpr[gene_overlap]
gene_names = list(geneExpr.columns)
cell_names = list(geneExpr.index)
n_modules = moduleGene.shape[0]
n_genes = len(gene_names)
n_cells = len(cell_names)
n_comps = cmMat.shape[0]

geneExprSum = geneExpr.sum(axis=1)
stand = geneExprSum.mean()
geneExprScale = torch.FloatTensor((geneExprSum / stand).values)

print("  Building module-blocked matrix...")
blocks = []
for i in range(n_modules):
    genes = [g for g in moduleGene.iloc[i, :].values if pd.notna(g)]
    temp = geneExpr.copy()
    temp.loc[:, [g for g in gene_names if g not in genes]] = 0
    blocks.append(temp.values)
X_blocked = np.hstack(blocks).astype(np.float32)
X_tensor = torch.FloatTensor(X_blocked)
del blocks, X_blocked
gc.collect()

module_sums = []
for i in range(n_modules):
    genes = [g for g in moduleGene.iloc[i, :].values if pd.notna(g)]
    overlap = [g for g in genes if g in gene_names]
    module_sums.append(geneExpr[overlap].sum(axis=1).values if overlap
                       else np.zeros(n_cells))
module_scale = torch.FloatTensor(np.column_stack(module_sums) / moduleLen)
cmMat_tensor = torch.FloatTensor(cmMat)

print(f"  Matrix: {n_cells} cells x {n_modules} modules x {n_genes} genes")

print(f"  Training scFEA ({EPOCH} epochs)...")
t_train = time.time()
torch.manual_seed(SEED)
net = FLUX(X_tensor, n_modules, f_in=n_genes, f_out=1)
optimizer = torch.optim.Adam(net.parameters(), lr=LEARN_RATE)
dataset = FluxDataset(X_tensor, geneExprScale, module_scale)
train_loader = DataLoader(dataset, batch_size=n_cells, shuffle=False,
                          num_workers=0, pin_memory=False)
net.train()
for epoch in tqdm(range(EPOCH), desc="  Training"):
    for X_batch, scale_batch, mscale_batch in train_loader:
        X_batch = Variable(X_batch.float())
        scale_batch = Variable(scale_batch.float())
        mscale_batch = Variable(mscale_batch.float())
        out_m, out_c = net(X_batch, n_modules, n_genes, n_comps, cmMat_tensor)
        loss = scfea_loss(out_m, out_c, LAMB_BA, LAMB_NG, LAMB_CELL, LAMB_MOD,
                          scale_batch, mscale_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
print(f"  Training done ({time.time()-t_train:.0f}s)")

print("  Predicting per-cell flux...")
net.eval()
test_loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0)
flux_matrix = np.zeros((n_cells, n_modules), dtype=np.float32)
balance_matrix = np.zeros((n_cells, n_comps), dtype=np.float32)
with torch.no_grad():
    for i, (X_i, _, _) in enumerate(test_loader):
        out_m, out_c = net(X_i.float(), n_modules, n_genes, n_comps, cmMat_tensor)
        flux_matrix[i, :] = out_m.numpy()
        balance_matrix[i, :] = out_c.numpy()

flux_df = pd.DataFrame(flux_matrix, index=cell_names, columns=moduleGene.index)
balance_df = pd.DataFrame(balance_matrix, index=cell_names, columns=cName)
balance_df.to_parquet(os.path.join(OUT_DIR, "21_balance_per_cell.parquet"))
print(f"  Flux: {flux_df.shape[0]} cells x {n_modules} modules")

del net, optimizer, dataset, train_loader, test_loader
del X_tensor, geneExprScale, module_scale, cmMat_tensor
del flux_matrix, balance_matrix, balance_df, geneExpr
gc.collect()

# ── 5g. Score and visualize flux ─────────────────────────────
# FIX: key the readable label map by the module INDEX (M_x), which matches
# flux_df.columns (= moduleGene.index). The original keyed by the integer
# Module_id column (1,2,...) and never matched the M_x flux columns.
flux_label_map = {}
if os.path.exists(FLUX_ANNOT):
    annot = pd.read_csv(FLUX_ANNOT, index_col=0)  # index = M_1, M_2, ... (matches flux cols)
    if "Compound_IN_name" in annot.columns and "Compound_OUT_name" in annot.columns:
        for mid, row in annot.iterrows():
            cin = row.get("Compound_IN_name", "?")
            cout = row.get("Compound_OUT_name", "?")
            flux_label_map[mid] = f"{cin}->{cout}"

flux_feats = flux_df.columns.tolist()
flux_means = compute_group_means(flux_df, flux_epitype)
flux_z = zscore_df(flux_means)
flux_top = select_top_differential(flux_means, n=8)
flux_radar_z = zscore_df(flux_means[flux_top])

if flux_label_map:
    flux_radar_z.columns = [flux_label_map.get(c, c) for c in flux_radar_z.columns]

flux_means.to_csv(os.path.join(OUT_DIR, "21_flux_means.csv"))
flux_z.to_csv(os.path.join(OUT_DIR, "21_flux_zscored.csv"))
flux_radar_z.to_csv(os.path.join(OUT_DIR, "21_flux_radar.csv"))
flux_df.to_parquet(os.path.join(OUT_DIR, "21_flux_per_cell.parquet"))

flux_z_labeled = flux_z.copy()
if flux_label_map:
    flux_z_labeled.columns = [flux_label_map.get(c, c) for c in flux_z_labeled.columns]

print(f"  Plotting flux heatmap ({len(flux_feats)} modules)...")
plot_full_heatmap(flux_z_labeled, f"Metabolic Flux (all {len(flux_feats)} modules)",
                  "21_flux_heatmap")
print("  Plotting flux violins...")
plot_full_violins(flux_df, flux_epitype, flux_feats,
                  "Metabolic Flux", "21_flux_violins", max_features=30)
del flux_df
gc.collect()
print(f"  Flux done ({time.time()-t4:.0f}s)")

# ── STEP 6: Combined 2x2 radar panel ─────────────────────────
print("\n[6/8] Combined 2x2 radar panel...")
radar_data = [
    (prog_radar_z, "PROGENy"),
    (hall_radar_z, "Hallmark"),
    (doro_radar_z, "DoRothEA"),
    (flux_radar_z, "Metabolic flux"),
]
fig, axes = plt.subplots(2, 2, figsize=(9, 9), subplot_kw=dict(polar=True))
for idx, (df_z, title) in enumerate(radar_data):
    r, c = divmod(idx, 2)
    radar_plot(axes[r, c], df_z, title)
handles = [Line2D([0], [0], color=PALETTE[g], linewidth=1.5, label=g) for g in ORDER]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           fontsize=8, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.04, 1, 1])
save_figs(fig, "21_radar_combined")
print("  Saved 21_radar_combined.svg/pdf")

# ── STEP 7: Individual radar panels ──────────────────────────
print("\n[7/8] Individual radar panels...")
for df_z, title in radar_data:
    stem = f"21_radar_{title.lower().replace(' ', '_')}"
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))
    radar_plot(ax, df_z, title)
    handles = [Line2D([0], [0], color=PALETTE[g], linewidth=1.5, label=g) for g in ORDER]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.35, 1.1),
              frameon=False, fontsize=7)
    fig.tight_layout()
    save_figs(fig, stem)
    print(f"  Saved {stem}.svg/pdf")

# ── STEP 8: Summary statistics ───────────────────────────────
print("\n[8/8] Summary statistics...")
summary_rows = []
for name, means_path in [("PROGENy", "21_progeny_means.csv"),
                         ("Hallmark", "21_hallmark_means.csv"),
                         ("DoRothEA", "21_dorothea_means.csv"),
                         ("Flux", "21_flux_means.csv")]:
    fpath = os.path.join(OUT_DIR, means_path)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, index_col=0)
        summary_rows.append({
            "axis": name, "n_features": df.shape[1],
            "top_SecA": df.loc["SecA"].idxmax() if "SecA" in df.index else "",
            "top_Intermediate": df.loc["Intermediate"].idxmax() if "Intermediate" in df.index else "",
            "top_SecB": df.loc["SecB"].idxmax() if "SecB" in df.index else "",
            "top_Ciliated": df.loc["Ciliated"].idxmax() if "Ciliated" in df.index else "",
        })
if summary_rows:
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(OUT_DIR, "21_summary.csv"), index=False)
    print(summary.to_string(index=False))

print(f"\n{'='*70}")
print(f"DONE — Epitype Functional Characterization ({(time.time()-t0)/60:.1f} min)")
print(f"  Output directory: {OUT_DIR}")
print(f"{'='*70}")
