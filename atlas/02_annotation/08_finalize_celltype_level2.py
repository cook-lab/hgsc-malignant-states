#!/usr/bin/env python3
"""
Atlas 02 — Step 08 (finalize): per-compartment level-2 subclustering + labels

PURPOSE
    For each level-1 compartment, subcluster (neighbours/UMAP/Leiden on the scANVI
    latent space), compute markers, assign celltype_level2 labels, and write the
    barcode->level2 maps + per-label marker CSVs consumed by 08_merge / 08c.

INPUTS
    DATA_ROOT/2026_final_atlas/hgsc_atlas_celltype_level1.h5ad

OUTPUTS
    output_root/02_annotation/08_celltype_level2/{barcode_maps,markers}/*.csv
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad

MANUSCRIPT PANEL(S)
    Annotation backend; defines the level-2 schema for Fig 1B-E and SF4.

RUNTIME TIER
    heavy (per-compartment neighbours + UMAP + Leiden).
"""

import gc
import os
import sys
import time
import warnings
from collections import OrderedDict
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================================
# PATHS
# ============================================================================

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

ATLAS_H5AD  = path("data_root", "2026_final_atlas", "hgsc_atlas_celltype_level1.h5ad")
OUT_DIR     = path("output_root", "02_annotation", "08_celltype_level2")
MAP_DIR     = os.path.join(OUT_DIR, "barcode_maps")
MARKER_DIR  = os.path.join(OUT_DIR, "markers")

for d in [OUT_DIR, MAP_DIR, MARKER_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================================
# COMPARTMENT REGISTRY
#
# Format:
#   key: short name (command line arg)
#   level1: celltype_level1 value
#   resolution: Leiden resolution (None = no clustering, assign single label)
#   cluster_map: {leiden_cluster_str: celltype_level2_label}
#   no_leiden_label: for compartments with no subclustering, the single label
# ============================================================================

COMPARTMENTS = OrderedDict()

# ── Endothelial ─────────────────────────────────────────────
COMPARTMENTS["endothelial"] = {
    "level1": "Endothelial",
    "resolution": 0.2,
    "cluster_map": {
        "0": "Activated capillary endothelial cell",
        "1": "Angiogenic endothelial cell",
        "2": "Arterial endothelial cell",
        "3": "Excluded endothelial cell",
        "4": "Lymphatic endothelial cell",
    },
}

# ── Mast cell ───────────────────────────────────────────────
COMPARTMENTS["mastcell"] = {
    "level1": "Mast cell",
    "resolution": 0.05,
    "cluster_map": {
        "0": "Mast cell",
        "1": "Excluded mast cell",
        "2": "Cycling mast cell",
    },
}

# ── Neutrophil (no further subsetting) ──────────────────────
COMPARTMENTS["neutrophil"] = {
    "level1": "Neutrophil",
    "resolution": None,
    "no_leiden_label": "Neutrophil",
    "cluster_map": {},
}

# ── Pericyte (no further subsetting) ────────────────────────
COMPARTMENTS["pericyte"] = {
    "level1": "Pericyte",
    "resolution": None,
    "no_leiden_label": "Pericyte",
    "cluster_map": {},
}

# ── Plasma cell ─────────────────────────────────────────────
COMPARTMENTS["plasmacell"] = {
    "level1": "Plasma cell",
    "resolution": 0.075,
    "cluster_map": {
        "0": "Plasma cell",
        "1": "Excluded plasma cell",
        "2": "Cycling plasma cell",
    },
}

# ── Smooth muscle ───────────────────────────────────────────
COMPARTMENTS["smoothmuscle"] = {
    "level1": "Smooth muscle",
    "resolution": 0.2,
    "cluster_map": {
        "0": "Homeostatic smooth muscle",
        "1": "Contractile smooth muscle cell",
        "2": "ECM-producing smooth muscle cell",
        "3": "Excluded smooth muscle cell_1",
    },
}

# ── DC ──────────────────────────────────────────────────────
COMPARTMENTS["dc"] = {
    "level1": "DC",
    "resolution": 0.05,
    "cluster_map": {
        "0": "Plasmacytoid dendritic cell",
        "1": "Conventional dendritic cell",
        "2": "Mature dendritic cell",
        "3": "Excluded dendritic cell",
    },
}

# ── B cell ──────────────────────────────────────────────────
COMPARTMENTS["bcell"] = {
    "level1": "B cell",
    "resolution": 0.125,
    "cluster_map": {
        "0": "Activated B cell",
        "1": "Naive B cell",
        "2": "Excluded B cell",
        "3": "Germinal centre B cell",
    },
}

# ── Mesothelial ─────────────────────────────────────────────
COMPARTMENTS["mesothelial"] = {
    "level1": "Mesothelial",
    "resolution": 0.1,
    "cluster_map": {
        "0": "Mesothelial cell",
        "1": "Hypoxic mesothelial cell",
        "2": "Excluded mesothelial cell_1",
        "3": "Excluded mesothelial cell_2",
    },
}

# ── Fibroblast ──────────────────────────────────────────────
COMPARTMENTS["fibroblast"] = {
    "level1": "Fibroblast",
    "resolution": 0.3,
    "cluster_map": {
        "0": "Myo-fibroblastic cancer-associated fibroblast",
        "1": "Activated fibroblast",
        "2": "Fibroblast",
        "3": "Excluded fibroblast_1",
        "4": "Complement secreting fibroblast",
        "5": "Inflammatory cancer-associated fibroblast",
        "6": "Cycling fibroblast",
        "7": "Ovarian stromal cell",
        "8": "Neuronal cell",
    },
}

# ── Macrophage ──────────────────────────────────────────────
COMPARTMENTS["macrophage"] = {
    "level1": "Macrophage",
    "resolution": 0.3,
    "cluster_map": {
        "0": "C1Q tissue-resident macrophage",
        "1": "Monocyte-derived macrophage",
        "2": "Excluded macrophage_1",
        "3": "Hypoxic macrophage",
        "4": "Excluded macrophage_2",
        "5": "Cycling macrophage",
        "6": "Perivascular macrophage",
        "7": "Myeloid-derived dendritic cell",
        "8": "Classical monocyte",
        "9": "Excluded macrophage_3",
    },
}

# ── Epithelial ──────────────────────────────────────────────
COMPARTMENTS["epithelial"] = {
    "level1": "Epithelial",
    "resolution": 0.4,
    "cluster_map": {
        "0": "Transitioning secretory epithelial cell",
        "1": "Adaptive secretory epithelial cell",
        "2": "Epithelial cell_1",
        "3": "Excluded epithelial cell_1",
        "4": "Epithelial cell_2",
        "5": "Excluded epithelial cell_2",
        "6": "Ciliated epithelial cell",
        "7": "Proliferative epithelial cell",
        "8": "Excluded epithelial cell_3",
    },
}

# ── T/NK cell ───────────────────────────────────────────────
COMPARTMENTS["tnkcell"] = {
    "level1": "T/NK cell",
    "resolution": 0.5,
    "cluster_map": {
        "0":  "CD8 memory T cell",
        "1":  "CD4 naive T cell",
        "2":  "Quiescent T cell",
        "3":  "CD56bright NK cell",
        "4":  "CD4 Regulatory T cell",
        "5":  "Excluded T/NK cell_1",
        "6":  "CD56dim NK cell",
        "7":  "Excluded T/NK cell_2",
        "8":  "Cycling T/NK cell",
        "9":  "T/NK cell",
        "10": "Excluded T/NK cell_3",
        "11": "Innate lymphoid cell",
    },
}

# Run order: smallest → largest (to catch config errors fast)
RUN_ORDER = [
    "mastcell", "neutrophil", "pericyte", "dc", "plasmacell",
    "smoothmuscle", "bcell", "mesothelial", "endothelial",
    "fibroblast", "macrophage", "tnkcell", "epithelial",
]


# ============================================================================
# CORE PIPELINE: PROCESS ONE COMPARTMENT
# ============================================================================

def process_compartment(comp_key, config, atlas=None):
    """Run the full finalization pipeline for one compartment.

    If atlas is provided, uses it directly (avoids re-loading for --all).
    """
    t0 = time.time()
    level1_name    = config["level1"]
    resolution     = config["resolution"]
    cluster_map    = config["cluster_map"]
    no_leiden      = resolution is None
    no_leiden_label = config.get("no_leiden_label", None)

    print(f"\n{'='*70}")
    if no_leiden:
        print(f"  Finalizing: {level1_name} ({comp_key}) — no subclustering")
    else:
        print(f"  Finalizing: {level1_name} ({comp_key}) at r={resolution}")
    print(f"{'='*70}\n")

    # ------------------------------------------------------------------
    # 1. Load and subset
    # ------------------------------------------------------------------
    print("1. Loading atlas and subsetting...", flush=True)
    free_atlas = False
    if atlas is None:
        atlas = ad.read_h5ad(ATLAS_H5AD)
        free_atlas = True
        print(f"   Atlas: {atlas.shape[0]:,} x {atlas.shape[1]:,}")

    mask = atlas.obs["celltype_level1"] == level1_name
    adata = atlas[mask].copy()
    if free_atlas:
        del atlas
        gc.collect()

    n_cells = len(adata)
    print(f"   {level1_name}: {n_cells:,} cells")

    # ------------------------------------------------------------------
    # 2. Recompute neighbors + UMAP on X_scanvi
    #    (matches 07_resolution_explorer.py exactly)
    # ------------------------------------------------------------------
    print("\n2. Computing neighbors (n=15, X_scanvi) → UMAP...", flush=True)
    sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15)
    sc.tl.umap(adata, random_state=SEED)
    print("   Done")

    # ------------------------------------------------------------------
    # 3. Normalize for marker detection
    # ------------------------------------------------------------------
    print("\n3. Normalizing for marker detection...", flush=True)
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("   Done (counts → normalize_total 1e4 → log1p)")

    # ------------------------------------------------------------------
    # 4. Leiden (or skip for no-subcluster compartments)
    # ------------------------------------------------------------------
    if no_leiden:
        print(f"\n4. No subclustering — assigning all cells → '{no_leiden_label}'")
        adata.obs["celltype_level2"] = no_leiden_label
        adata.obs["celltype_level2"] = pd.Categorical(adata.obs["celltype_level2"])
        leiden_key = None
    else:
        leiden_key = f"leiden_r{resolution}"
        print(f"\n4. Leiden clustering at r={resolution}...", flush=True)
        # Match 07_resolution_explorer.py: default flavor/iterations/seed
        sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key, random_state=SEED)
        n_cl = adata.obs[leiden_key].nunique()
        print(f"   {n_cl} clusters")

        # Show cluster sizes
        clusters = sorted(adata.obs[leiden_key].unique(), key=lambda x: int(x))
        print(f"\n   Cluster sizes:")
        for cl in clusters:
            n = (adata.obs[leiden_key] == cl).sum()
            label = cluster_map.get(cl, "??? UNMAPPED")
            print(f"     {cl:>3}  {n:>9,}  → {label}")

        # Validate
        unmapped = [cl for cl in clusters if cl not in cluster_map]
        if unmapped:
            print(f"\n   ERROR: Unmapped clusters: {unmapped}")
            print(f"   Leiden produced: {clusters}")
            print(f"   cluster_map has: {list(cluster_map.keys())}")
            raise ValueError(f"Cluster map incomplete — missing: {unmapped}")

        extra = [k for k in cluster_map if k not in clusters]
        if extra:
            print(f"\n   WARNING: cluster_map has extra keys not in Leiden: {extra}")

        # Assign
        print("\n5. Assigning celltype_level2 labels...", flush=True)
        adata.obs["celltype_level2"] = adata.obs[leiden_key].map(cluster_map)
        adata.obs["celltype_level2"] = pd.Categorical(adata.obs["celltype_level2"])

    # Distribution
    print(f"\n   celltype_level2 distribution:")
    for label, count in adata.obs["celltype_level2"].value_counts().items():
        pct = 100 * count / n_cells
        print(f"     {label:50s} {count:>8,}  ({pct:5.1f}%)")

    # ------------------------------------------------------------------
    # 6. Wilcoxon DEG — all significant markers per level2 cluster
    # ------------------------------------------------------------------
    step = "5" if no_leiden else "6"
    print(f"\n{step}. Computing DEGs per celltype_level2 (Wilcoxon)...", flush=True)

    level2_labels = sorted(adata.obs["celltype_level2"].unique())
    n_groups = len(level2_labels)

    if n_groups == 1:
        # Single group — no DEG possible, save all expressed genes as markers
        print(f"   Single group — saving expressed gene list instead of DEGs")
        label = level2_labels[0]
        # Get mean expression per gene
        mean_expr = np.asarray(adata.X.mean(axis=0)).ravel()
        gene_df = pd.DataFrame({
            "names": adata.var_names,
            "mean_expression": mean_expr,
        }).sort_values("mean_expression", ascending=False)
        gene_df = gene_df[gene_df["mean_expression"] > 0]

        safe_label = label.replace("/", "-").replace(" ", "_")
        fname = f"{comp_key}_{safe_label}.csv"
        gene_df.to_csv(os.path.join(MARKER_DIR, fname), index=False)
        print(f"   {label:50s} → {len(gene_df):>5} expressed genes → {fname}")
    else:
        sc.tl.rank_genes_groups(adata, groupby="celltype_level2", method="wilcoxon",
                                 use_raw=False)
        for label in level2_labels:
            result = sc.get.rank_genes_groups_df(adata, group=label)

            # Filter to significant markers (padj < 0.05 and logFC > 0)
            sig = result[(result["pvals_adj"] < 0.05) & (result["logfoldchanges"] > 0)].copy()
            sig = sig.sort_values("scores", ascending=False)

            safe_label = label.replace("/", "-").replace(" ", "_")
            fname = f"{comp_key}_{safe_label}.csv"
            sig.to_csv(os.path.join(MARKER_DIR, fname), index=False)

            print(f"   {label:50s} → {len(sig):>5} sig markers → {fname}")

    # ------------------------------------------------------------------
    # 7. Save barcode map
    # ------------------------------------------------------------------
    step_n = "6" if no_leiden else "7"
    print(f"\n{step_n}. Saving barcode → celltype_level2 map...", flush=True)
    barcode_data = {
        "barcode": adata.obs.index,
        "celltype_level2": adata.obs["celltype_level2"].values,
    }
    if leiden_key is not None:
        barcode_data["leiden_cluster"] = adata.obs[leiden_key].values
    barcode_map = pd.DataFrame(barcode_data)
    map_path = os.path.join(MAP_DIR, f"{comp_key}_barcode_to_level2.csv")
    barcode_map.to_csv(map_path, index=False)
    print(f"   Saved: {map_path}")
    print(f"   {len(barcode_map):,} barcodes")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  DONE — {level1_name} finalized in {elapsed/60:.1f} min")
    print(f"{'='*70}")
    if resolution is not None:
        print(f"  Resolution:    r={resolution}")
    else:
        print(f"  Resolution:    (none — single label)")
    print(f"  Cells:         {n_cells:,}")
    print(f"  Level2 labels: {n_groups}")
    print(f"  Barcode map:   {map_path}")
    print(f"  Marker CSVs:   {MARKER_DIR}/{comp_key}_*.csv")
    print()

    del adata
    gc.collect()

    return True


# ============================================================================
# STATUS: Show which compartments have been completed
# ============================================================================

def show_status():
    """Show progress across all compartments."""
    print(f"\n{'='*70}")
    print(f"  Status — celltype_level2 finalization")
    print(f"{'='*70}\n")
    done = 0
    for comp_key in RUN_ORDER:
        map_path = os.path.join(MAP_DIR, f"{comp_key}_barcode_to_level2.csv")
        if os.path.exists(map_path):
            df = pd.read_csv(map_path)
            n = len(df)
            n_labels = df["celltype_level2"].nunique()
            print(f"  [DONE]    {comp_key:15s}  {n:>9,} barcodes  {n_labels} labels")
            done += 1
        else:
            print(f"  [PENDING] {comp_key:15s}")
    print(f"\n  {done}/{len(RUN_ORDER)} compartments complete")


# ============================================================================
# MERGE: Combine all barcode maps into atlas
# ============================================================================

def merge_into_atlas():
    """Read all barcode maps and merge celltype_level2 into the atlas."""
    print(f"\n{'='*70}")
    print(f"  Merging celltype_level2 into atlas")
    print(f"{'='*70}\n")

    # 1. Load all barcode maps
    print("1. Loading barcode maps...", flush=True)
    map_files = sorted([
        f for f in os.listdir(MAP_DIR)
        if f.endswith("_barcode_to_level2.csv") and not f.startswith(".")
    ])
    print(f"   Found {len(map_files)} map files")

    dfs = []
    for f in map_files:
        df = pd.read_csv(os.path.join(MAP_DIR, f))
        compartment = f.replace("_barcode_to_level2.csv", "")
        n_unique = df["celltype_level2"].nunique()
        print(f"   {compartment:20s}: {len(df):>8,} barcodes, {n_unique} level2 labels")
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    print(f"\n   Total barcodes: {len(merged):,}")

    # Check for duplicates
    n_dup = merged["barcode"].duplicated().sum()
    if n_dup > 0:
        print(f"   WARNING: {n_dup} duplicate barcodes found!")
        merged = merged.drop_duplicates(subset="barcode", keep="last")
        print(f"   After dedup: {len(merged):,}")
    else:
        print("   No duplicate barcodes — OK")

    barcode_to_level2 = merged.set_index("barcode")["celltype_level2"]

    # 2. Load atlas
    print("\n2. Loading celltype_level1 atlas...", flush=True)
    t0 = time.time()
    adata = ad.read_h5ad(ATLAS_H5AD)
    print(f"   {adata.shape[0]:,} cells x {adata.shape[1]:,} genes ({time.time()-t0:.0f}s)")

    # 3. Map
    print("\n3. Mapping celltype_level2...", flush=True)
    adata.obs["celltype_level2"] = barcode_to_level2.reindex(adata.obs.index).values

    n_mapped = adata.obs["celltype_level2"].notna().sum()
    n_missing = adata.obs["celltype_level2"].isna().sum()
    print(f"   Mapped:  {n_mapped:,} ({100 * n_mapped / adata.n_obs:.1f}%)")
    print(f"   Missing: {n_missing:,} ({100 * n_missing / adata.n_obs:.1f}%)")

    if n_missing > 0:
        missing_mask = adata.obs["celltype_level2"].isna()
        print("\n   Missing by celltype_level1:")
        print(adata.obs.loc[missing_mask, "celltype_level1"].value_counts().to_string())
        adata.obs["celltype_level2"] = adata.obs["celltype_level2"].fillna("Unmapped")

    adata.obs["celltype_level2"] = pd.Categorical(adata.obs["celltype_level2"])

    # 4. Summary
    print("\n4. celltype_level2 distribution:")
    vc = adata.obs["celltype_level2"].value_counts()
    for label, count in vc.items():
        pct = 100 * count / adata.n_obs
        print(f"   {label:50s} {count:>8,}  ({pct:5.1f}%)")
    print(f"\n   Total unique level2 labels: {vc.shape[0]}")

    # 5. Save
    output_h5ad = obj("atlas_celltype_l2")
    print(f"\n5. Saving: {output_h5ad}", flush=True)
    t0 = time.time()
    adata.write_h5ad(output_h5ad)
    fsize = os.path.getsize(output_h5ad) / 1e9
    print(f"   File size: {fsize:.1f} GB ({time.time()-t0:.0f}s)")

    print(f"\n{'='*70}")
    print(f"  DONE — celltype_level2 atlas saved")
    print(f"{'='*70}")
    print(f"  Cells:  {adata.n_obs:,}")
    print(f"  Level2: {vc.shape[0]} unique labels")
    print(f"  Mapped: {n_mapped:,} | Missing: {n_missing:,}")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print(f"{'='*70}")
    print(f"  08_finalize_celltype_level2.py — HGSC Atlas")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")

    if len(sys.argv) < 2:
        print("\nUsage:")
        print("  python scripts/08_finalize_celltype_level2.py <compartment>")
        print("  python scripts/08_finalize_celltype_level2.py --all")
        print("  python scripts/08_finalize_celltype_level2.py --merge")
        print("  python scripts/08_finalize_celltype_level2.py --status")
        print(f"\nRegistered compartments ({len(COMPARTMENTS)}):")
        for k, v in COMPARTMENTS.items():
            r = v["resolution"] if v["resolution"] is not None else "none"
            print(f"  {k:15s}  r={r}")
        sys.exit(0)

    arg = sys.argv[1].lower().strip()

    if arg == "--merge":
        merge_into_atlas()
    elif arg == "--status":
        show_status()
    elif arg == "--all":
        print(f"\nRunning ALL {len(RUN_ORDER)} compartments sequentially...")
        print(f"Order: {' → '.join(RUN_ORDER)}\n")

        # Load atlas ONCE for all compartments
        print("Loading atlas (shared across all compartments)...", flush=True)
        t_load = time.time()
        atlas = ad.read_h5ad(ATLAS_H5AD)
        print(f"Atlas: {atlas.shape[0]:,} x {atlas.shape[1]:,} ({time.time()-t_load:.0f}s)\n")

        results = {}
        t_total = time.time()
        for i, comp_key in enumerate(RUN_ORDER, 1):
            config = COMPARTMENTS[comp_key]
            print(f"\n[{i}/{len(RUN_ORDER)}] {comp_key}")
            try:
                process_compartment(comp_key, config, atlas=atlas)
                results[comp_key] = "SUCCESS"
            except Exception as e:
                import traceback
                traceback.print_exc()
                results[comp_key] = f"FAILED: {e}"
            gc.collect()

        # Summary
        elapsed = time.time() - t_total
        print(f"\n{'='*70}")
        print(f"  ALL COMPARTMENTS — {elapsed/60:.1f} min total")
        print(f"{'='*70}")
        for comp, status in results.items():
            print(f"  {comp:15s}: {status}")
        n_ok = sum(1 for s in results.values() if s == "SUCCESS")
        n_fail = len(results) - n_ok
        print(f"\n  {n_ok}/{len(results)} succeeded, {n_fail} failed")

        del atlas
        gc.collect()

    elif arg in COMPARTMENTS:
        process_compartment(arg, COMPARTMENTS[arg])
    else:
        print(f"\nERROR: Unknown argument '{arg}'")
        print(f"Registered: {', '.join(COMPARTMENTS.keys())}")
        print(f"Or: --all, --merge, --status")
        sys.exit(1)
