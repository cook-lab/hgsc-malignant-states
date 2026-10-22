#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/21_epitype_functional_characterization.py, the script that ran.
# Kept: the per-cell PROGENy, Hallmark and DoRothEA scores (21_*_per_cell.parquet; Fig 2a, Figure S8), the per-state z-scores
# (21_*_zscored.csv) and the 20,000-cell scFEA flux, which the Fig 2a prep compares against. Removed: heatmaps, violins, radars,
# the means, radar, balance and summary outputs and the closing report. The untrimmed script is in the lab project folder.
"""
Step 21 — Epitype Functional Characterization (Fresh Computation)
==================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Freshly computes PROGENy, Hallmark, DoRothEA, and metabolic flux scores
for all 4 NMF-based epithelial subtypes directly from the epithelial h5ad.

Label assignment (NMF Factor 2 percentile schema):
  - SecA:         non-ciliated, Factor_2 < p50
  - Intermediate: non-ciliated, p50 <= Factor_2 < p75
  - SecB:         non-ciliated, Factor_2 >= p75
  - Ciliated:     celltype_level2 == "Ciliated epithelial cell"

Computation methods:
  - PROGENy:  decoupler MLM (dc.op.progeny + dc.mt.mlm)
  - Hallmark: scanpy score_genes (gene sets from dc.op.hallmark)
  - DoRothEA: decoupler MLM (dc.op.dorothea A/B/C + dc.mt.mlm)
  - Flux:     scFEA graph neural network (embedded, 168 modules, 20k cells)

Visualisations (all features shown in their entirety):
  - Full z-scored heatmaps with hierarchical column clustering
  - Violin grids (top 30 most differential features per axis)
  - Radar plots (top 8 per axis, individual + 2x2 combined)

Outputs -> output/21_epitype_functional_characterization/

Usage:
    python scripts/21_epitype_functional_characterization.py
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
from scipy import sparse

import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

BASE       = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD       = os.path.join(BASE, "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
NMF_CSV    = os.path.join(BASE, "output", "11d_epithelial_nmf", "11d_nmf_usage.csv")
OUT_DIR    = os.path.join(BASE, "output", "21_epitype_functional_characterization")
SCFEA_DATA = os.path.join(BASE, "tools", "scFEA", "data")
os.makedirs(OUT_DIR, exist_ok=True)

# scFEA hyperparameters (matching 11v2_epithelial_flux.py)
DOWNSAMPLE_N = 20_000
EPOCH        = 100
LEARN_RATE   = 0.008
LAMB_BA      = 1
LAMB_NG      = 1
LAMB_CELL    = 1
LAMB_MOD     = 1e-2
RNG          = np.random.RandomState(42)

# ============================================================================
# GROUPS & PALETTE
# ============================================================================

ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]

# ============================================================================
# HELPERS
# ============================================================================

def zscore_df(df):
    """Z-score each column across groups."""
    mu = df.mean(axis=0)
    sd = df.std(axis=0).replace(0, 1)
    return (df - mu) / sd


def compute_group_means(scores_df, group_labels, feature_cols=None):
    """Per-group means for given features."""
    if feature_cols is None:
        feature_cols = scores_df.columns.tolist()
    numeric_cols = [c for c in feature_cols
                    if scores_df[c].dtype in ("float64", "float32", "int64", "int32")]
    df = scores_df[numeric_cols].copy()
    df["group"] = group_labels
    return df.groupby("group").mean(numeric_only=True).reindex(ORDER)


def select_top_differential(group_means, n=8):
    """Top n features by range across groups."""
    ranges = group_means.max(axis=0) - group_means.min(axis=0)
    return ranges.nlargest(n).index.tolist()


# ============================================================================
# scFEA NETWORK (from Chang et al. 2021, Genome Research)
# Copied verbatim from scripts/11v2_epithelial_flux.py
# ============================================================================

class FLUX(nn.Module):
    """Graph neural network for metabolic flux estimation."""
    def __init__(self, matrix, n_modules, f_in=50, f_out=1):
        super(FLUX, self).__init__()
        self.inSize = f_in
        self.m_encoder = nn.ModuleList([
            nn.Sequential(
                nn.Linear(self.inSize, 8, bias=False),
                nn.Tanhshrink(),
                nn.Linear(8, f_out),
                nn.Tanhshrink(),
            )
            for _ in range(n_modules)
        ])

    def updateC(self, m, n_comps, cmMat):
        c = torch.zeros((m.shape[0], n_comps))
        for i in range(c.shape[1]):
            tmp = m * cmMat[i, :]
            c[:, i] = torch.sum(tmp, dim=1)
        return c

    def forward(self, x, n_modules, n_genes, n_comps, cmMat):
        for i in range(n_modules):
            x_block = x[:, i * n_genes: (i + 1) * n_genes]
            subnet = self.m_encoder[i]
            if i == 0:
                m = subnet(x_block)
            else:
                m = torch.cat((m, subnet(x_block)), 1)
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
    mean_x = torch.mean(x)
    mean_y = torch.mean(y)
    xm = x - mean_x
    ym = y - mean_y
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
        corr = torch.abs(corr)
        total4 = torch.FloatTensor(np.ones(m.shape[0])) - corr
    else:
        total4 = torch.FloatTensor(np.zeros(m.shape[0]))
    loss = (torch.sum(lamb1 * total1) + torch.sum(lamb2 * total2) +
            torch.sum(lamb3 * total3) + torch.sum(lamb4 * total4))
    return loss


# ============================================================================
# AXIS SCORING + VISUALIZATION HELPER
# ============================================================================

def score_and_visualize(per_cell_scores, group_labels, axis_name, prefix):
    """
    Given per-cell scores DataFrame and group labels:
    - Compute group means, z-score, select top 8
    - Save CSVs, plot heatmap + violins
    - Return radar z-scored DataFrame
    """
    feats = per_cell_scores.columns.tolist()
    means = compute_group_means(per_cell_scores, group_labels)
    z = zscore_df(means)
    top8 = select_top_differential(means, n=8)
    radar_z = zscore_df(means[top8])

    # Save data
    z.to_csv(os.path.join(OUT_DIR, f"{prefix}_zscored.csv"))
    per_cell_scores.to_parquet(os.path.join(OUT_DIR, f"{prefix}_per_cell.parquet"))

    return radar_z


# ============================================================================
# MAIN
# ============================================================================

print("=" * 70)
print("Step 21 — Epitype Functional Characterization (Fresh Computation)")
print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# ---------------------------------------------------------------------------
# STEP 1: Load h5ad & assign NMF-based epitype labels
# ---------------------------------------------------------------------------
t0 = time.time()
print("\n[1/8] Loading epithelial h5ad (full, not backed)...")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes  ({time.time()-t0:.0f}s)")

# Verify log-normalized
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

# Subset to shared cells
adata = adata[shared].copy()
f2 = f2.loc[shared]

# Assign epitype labels
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

# ---------------------------------------------------------------------------
# STEP 2: PROGENy pathway activity (fresh)
# ---------------------------------------------------------------------------
progeny_done = os.path.exists(os.path.join(OUT_DIR, "21_progeny_per_cell.parquet"))
if progeny_done:
    print("[2/8] PROGENy — loading existing results (already computed)...")
    t1 = time.time()
    progeny = pd.read_parquet(os.path.join(OUT_DIR, "21_progeny_per_cell.parquet"))
    # Reindex to current shared cells (in case of slight differences)
    common_prog = progeny.index.intersection(pd.Index(adata.obs_names))
    progeny = progeny.loc[common_prog]
    prog_labels = adata.obs.loc[common_prog, "epitype"].values
    prog_means = compute_group_means(progeny, prog_labels)
    prog_z = zscore_df(prog_means)
    top8 = select_top_differential(prog_means, n=8)
    prog_radar_z = zscore_df(prog_means[top8])
    print(f"  Loaded: {progeny.shape[1]} pathways x {progeny.shape[0]:,} cells")
    del progeny
    gc.collect()
else:
    print("[2/8] PROGENy pathway activity (decoupler MLM)...")
    t1 = time.time()

    progeny_net = dc.op.progeny(organism="human")
    print(f"  Network: {progeny_net.shape[0]:,} links, "
          f"{progeny_net['source'].nunique()} pathways")
    dc.mt.mlm(adata, net=progeny_net)
    progeny = adata.obsm["score_mlm"].copy()  # copy before DoRothEA overwrites
    print(f"  Scored: {progeny.shape[1]} pathways x {progeny.shape[0]:,} cells  "
          f"({time.time()-t1:.0f}s)")

    prog_radar_z = score_and_visualize(progeny, group_labels, "PROGENy Pathway Activity",
                                        "21_progeny")
    del progeny
    gc.collect()

# ---------------------------------------------------------------------------
# STEP 3: Hallmark gene sets (fresh — sc.tl.score_genes, full data)
# ---------------------------------------------------------------------------
print(f"\n[3/8] Hallmark gene sets (score_genes on {adata.n_obs:,} cells)...")
t2 = time.time()

# Get hallmark gene sets from decoupler
hallmark_net = dc.op.hallmark(organism="human")
hall_sets = {}
for gs_name, grp in hallmark_net.groupby("source"):
    genes = [g for g in grp["target"].values if g in adata.var_names]
    if len(genes) >= 5:
        hall_sets[gs_name] = genes
print(f"  {len(hall_sets)} gene sets with ≥5 overlapping genes (of {hallmark_net['source'].nunique()} total)")

# Score each set with scanpy
hall_scores = {}
for gs_name, genes in hall_sets.items():
    key = f"hall_{gs_name}"
    sc.tl.score_genes(adata, gene_list=genes, score_name=key, use_raw=False)
    # Clean name for display: remove HALLMARK_ prefix
    display_name = gs_name.replace("HALLMARK_", "")
    hall_scores[display_name] = adata.obs[key].values

hallmark = pd.DataFrame(hall_scores, index=adata.obs_names)
print(f"  Scored: {hallmark.shape[1]} gene sets x {hallmark.shape[0]:,} cells  "
      f"({time.time()-t2:.0f}s)")

hall_radar_z = score_and_visualize(hallmark, group_labels, "Hallmark Gene Sets",
                                    "21_hallmark")
del hallmark
gc.collect()

# Clear any obsm from previous steps to free memory
for key in list(adata.obsm.keys()):
    if key.startswith("score_"):
        del adata.obsm[key]
gc.collect()

# ---------------------------------------------------------------------------
# STEP 4: DoRothEA TF activity (fresh, MLM on full data)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# STEP 5: Metabolic flux via scFEA (fresh)
# ---------------------------------------------------------------------------
print(f"\n[5/8] Metabolic flux (scFEA, fresh training)...")
t4 = time.time()

# 5a. Downsample to DOWNSAMPLE_N cells (stratified by epitype)
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

# 5b. Prepare expression matrix
print("  Preparing expression matrix...")
X_flux = adata_flux.X
if sparse.issparse(X_flux):
    X_flux = X_flux.toarray()
geneExpr = pd.DataFrame(X_flux, index=adata_flux.obs_names,
                         columns=adata_flux.var_names)
if geneExpr.max().max() > 50:
    print("  Values > 50 detected, applying log2(x+1)...")
    geneExpr = np.log2(geneExpr + 1)

# Free adata to reclaim memory before building blocked matrix
del adata, adata_flux, X_flux
gc.collect()

# 5c. Load scFEA metabolic model
print("  Loading scFEA metabolic model (168 modules)...")
moduleGene = pd.read_csv(os.path.join(SCFEA_DATA, "module_gene_m168.csv"),
                         index_col=0)
cmMat = pd.read_csv(os.path.join(SCFEA_DATA, "cmMat_c70_m168.csv"),
                     header=None).values
cName = pd.read_csv(os.path.join(SCFEA_DATA, "cName_c70_m168.csv")).columns

moduleLen = np.array([moduleGene.iloc[i, :].notna().sum()
                      for i in range(moduleGene.shape[0])])

# Find overlapping genes
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

# Gene scaling
geneExprSum = geneExpr.sum(axis=1)
stand = geneExprSum.mean()
geneExprScale = torch.FloatTensor((geneExprSum / stand).values)

# 5d. Build module-blocked expression matrix
print("  Building module-blocked matrix...")
blocks = []
for i in range(n_modules):
    genes = [g for g in moduleGene.iloc[i, :].values if pd.notna(g)]
    temp = geneExpr.copy()
    mask_cols = [g for g in gene_names if g not in genes]
    temp.loc[:, mask_cols] = 0
    blocks.append(temp.values)

X_blocked = np.hstack(blocks).astype(np.float32)
X_tensor = torch.FloatTensor(X_blocked)
del blocks, X_blocked
gc.collect()

# Module scale
module_sums = []
for i in range(n_modules):
    genes = [g for g in moduleGene.iloc[i, :].values if pd.notna(g)]
    overlap = [g for g in genes if g in gene_names]
    if overlap:
        module_sums.append(geneExpr[overlap].sum(axis=1).values)
    else:
        module_sums.append(np.zeros(n_cells))
module_scale = np.column_stack(module_sums) / moduleLen
module_scale = torch.FloatTensor(module_scale)
cmMat_tensor = torch.FloatTensor(cmMat)

print(f"  Matrix: {n_cells} cells x {n_modules} modules x {n_genes} genes")

# 5e. Train scFEA network
print(f"  Training scFEA ({EPOCH} epochs)...")
t_train = time.time()

torch.manual_seed(16)
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
        loss = scfea_loss(out_m, out_c,
                          LAMB_BA, LAMB_NG, LAMB_CELL, LAMB_MOD,
                          scale_batch, mscale_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

print(f"  Training done ({time.time()-t_train:.0f}s)")

# 5f. Predict per-cell flux
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

flux_df = pd.DataFrame(flux_matrix, index=cell_names,
                        columns=moduleGene.index)
balance_df = pd.DataFrame(balance_matrix, index=cell_names,
                           columns=cName)

print(f"  Flux: {flux_df.shape[0]} cells x {n_modules} modules")
print(f"  Balance: {balance_df.shape[0]} cells x {n_comps} compounds")

# Clean up scFEA objects
del net, optimizer, dataset, train_loader, test_loader
del X_tensor, geneExprScale, module_scale, cmMat_tensor
del flux_matrix, balance_matrix, balance_df, geneExpr
gc.collect()

# 5g. Score and visualize flux
flux_feats = flux_df.columns.tolist()
flux_means = compute_group_means(flux_df, flux_epitype)
flux_z = zscore_df(flux_means)
flux_top = select_top_differential(flux_means, n=8)
flux_radar_z = zscore_df(flux_means[flux_top])

# Save data
flux_z.to_csv(os.path.join(OUT_DIR, "21_flux_zscored.csv"))
flux_df.to_parquet(os.path.join(OUT_DIR, "21_flux_per_cell.parquet"))

del flux_df
gc.collect()
print(f"  Flux done ({time.time()-t4:.0f}s)")
