#!/usr/bin/env python3
"""
Step 19_cnv — Extract per-sample count matrices for CopyKAT
============================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Extracts per-sample raw count matrices (epithelial + reference cells)
for CopyKAT CNV inference. Reference cells = fibroblast + T/NK + endothelial
from the same sample, subsampled to MAX_REF_CELLS per sample so CopyKAT
runtime stays bounded (diploid baseline saturates well below 500 cells).

Epithelial cells labeled with NMF percentile SecA/Transitioning/SecB/Ciliated.

By default extracts ALL eligible samples (>=MIN_EPI_CELLS epi +
>=MIN_REF_CELLS ref). Use --n-samples N to limit for a pilot.

Input:
    celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
    celltype_h5ad/hgsc_atlas_final_fibroblast.h5ad
    celltype_h5ad/hgsc_atlas_final_t_nk_cell.h5ad
    celltype_h5ad/hgsc_atlas_final_endothelial.h5ad
    output/11d_epithelial_nmf/11d_nmf_usage.csv

Output:
    output/19_cnv/per_sample/{sample_id}/
        counts.mtx.gz, genes.txt, barcodes.csv, ref_barcodes.txt
    output/19_cnv/tables/sample_manifest.csv

Usage:
    python scripts/19_cnv_extract.py                     # all eligible
    python scripts/19_cnv_extract.py --n-samples 10      # pilot
"""

import argparse
import gzip
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

# ── Paths ────────────────────────────────────────────────────
BASE    = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_EPI = os.path.join(BASE, "celltype_h5ad",
                        "hgsc_atlas_final_epithelial.h5ad")
H5AD_REF = {
    "fibroblast":  os.path.join(BASE, "celltype_h5ad",
                                "hgsc_atlas_final_fibroblast.h5ad"),
    "t_nk_cell":   os.path.join(BASE, "celltype_h5ad",
                                "hgsc_atlas_final_t_nk_cell.h5ad"),
    "endothelial": os.path.join(BASE, "celltype_h5ad",
                                "hgsc_atlas_final_endothelial.h5ad"),
}
NMF_CSV = os.path.join(BASE, "output", "11d_epithelial_nmf",
                       "11d_nmf_usage.csv")
OUT_DIR = os.path.join(BASE, "output", "19_cnv")

MIN_EPI_CELLS = 200
MIN_REF_CELLS = 50
MAX_REF_CELLS = 500            # cap per sample, deterministic RNG below
REF_SUBSAMPLE_SEED = 42
DEFAULT_N_SAMPLES = 0          # 0 = all eligible


def build_epitype_labels(obs_df):
    """Assign NMF percentile epitype labels to epithelial cells."""
    nmf_usage = pd.read_csv(NMF_CSV, index_col=0)
    f2 = nmf_usage["Factor_2"]

    # Align
    common = obs_df.index.intersection(f2.index)
    f2 = f2.loc[common]

    # Ciliated first
    is_ciliated = obs_df.loc[common, "celltype_level2"] == "Ciliated epithelial cell"
    non_cil_f2 = f2[~is_ciliated]

    p50 = np.percentile(non_cil_f2, 50)
    p75 = np.percentile(non_cil_f2, 75)
    print(f"  NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

    labels = pd.Series("SecA", index=common)
    labels[f2 >= p50] = "Transitioning"
    labels[f2 >= p75] = "SecB"
    labels[is_ciliated] = "Ciliated"

    return labels


def collect_reference_barcodes():
    """Load per-sample ref barcodes, then cap each sample at MAX_REF_CELLS.

    Stratified cap: when a compartment contributes >1 cell, cap proportionally
    across fibroblast / T-NK / endothelial so the diploid baseline keeps all
    three compartments represented.
    """
    # Load raw by compartment
    ref_by_compartment = {name: {} for name in H5AD_REF}
    for name, path in H5AD_REF.items():
        print(f"  Loading {name} (backed)...")
        adata_ref = ad.read_h5ad(path, backed="r")
        for sid, grp in adata_ref.obs[["sample_id"]].groupby("sample_id"):
            ref_by_compartment[name].setdefault(sid, [])
            ref_by_compartment[name][sid].extend(grp.index.tolist())
        del adata_ref

    # Union of sample ids
    all_sids = set()
    for d in ref_by_compartment.values():
        all_sids.update(d.keys())

    rng = np.random.default_rng(REF_SUBSAMPLE_SEED)
    ref_by_sample = {}
    total_pre = 0
    total_post = 0

    for sid in all_sids:
        per_comp = {n: ref_by_compartment[n].get(sid, []) for n in H5AD_REF}
        n_total = sum(len(v) for v in per_comp.values())
        total_pre += n_total

        if n_total <= MAX_REF_CELLS:
            picked = []
            for v in per_comp.values():
                picked.extend(v)
        else:
            picked = []
            remaining = MAX_REF_CELLS
            comps = sorted(per_comp.items(), key=lambda x: len(x[1]))
            for i, (_n, bcs) in enumerate(comps):
                comps_left = len(comps) - i
                share = remaining // comps_left
                if len(bcs) <= share:
                    picked.extend(bcs)
                    remaining -= len(bcs)
                else:
                    idx = rng.choice(len(bcs), size=share, replace=False)
                    picked.extend([bcs[j] for j in idx])
                    remaining -= share

        ref_by_sample[sid] = picked
        total_post += len(picked)

    print(f"  Reference cells: {total_pre:,} raw -> "
          f"{total_post:,} capped at {MAX_REF_CELLS}/sample "
          f"across {len(ref_by_sample)} samples")
    return ref_by_sample


def write_sample_data(sample_id, epi_barcodes, ref_barcodes, epitype_labels,
                      epi_counts, ref_counts, gene_names, out_base):
    """Write per-sample CopyKAT input files."""
    sample_dir = os.path.join(out_base, "per_sample", sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    # Combine epithelial + reference counts
    all_counts = sparse.vstack([epi_counts, ref_counts], format="csc")
    all_barcodes = list(epi_barcodes) + list(ref_barcodes)

    # CopyKAT expects genes × cells
    counts_t = all_counts.T.tocsc()

    # Write sparse matrix (genes × cells)
    mtx_path = os.path.join(sample_dir, "counts.mtx.gz")
    with gzip.open(mtx_path, "wb") as f:
        sio.mmwrite(f, counts_t)

    # Write gene names
    genes_path = os.path.join(sample_dir, "genes.txt")
    with open(genes_path, "w") as f:
        f.write("\n".join(gene_names))

    # Write barcodes with metadata
    bc_df = pd.DataFrame({
        "barcode": all_barcodes,
        "epitype": [epitype_labels.get(b, "reference") for b in all_barcodes],
        "is_reference": [b in set(ref_barcodes) for b in all_barcodes],
    })
    bc_df.to_csv(os.path.join(sample_dir, "barcodes.csv"), index=False)

    # Write reference barcodes
    with open(os.path.join(sample_dir, "ref_barcodes.txt"), "w") as f:
        f.write("\n".join(ref_barcodes))

    return len(epi_barcodes), len(ref_barcodes)


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-sample CopyKAT inputs")
    parser.add_argument("--n-samples", type=int, default=DEFAULT_N_SAMPLES,
                        help="Number of top samples to extract (0 = all eligible)")
    args = parser.parse_args()

    os.makedirs(os.path.join(OUT_DIR, "per_sample"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "tables"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "figs"), exist_ok=True)

    limit_msg = "all eligible" if args.n_samples == 0 else f"top {args.n_samples}"
    print("=" * 65)
    print("  19_cnv — Extract per-sample CopyKAT inputs")
    print(f"  Target: {limit_msg} samples   Ref cap: {MAX_REF_CELLS}/sample")
    print("=" * 65)

    # ── 1. Load epithelial metadata ───────────────────────────
    print("\n[1] Loading epithelial h5ad...")
    t0 = time.time()
    adata_epi = ad.read_h5ad(H5AD_EPI, backed="r")
    obs_epi = adata_epi.obs.copy()
    gene_names = adata_epi.var_names.tolist()
    print(f"  {len(obs_epi):,} cells × {len(gene_names):,} genes  "
          f"({time.time() - t0:.0f}s)")

    # ── 2. Build epitype labels ───────────────────────────────
    print("\n[2] Building NMF percentile epitype labels...")
    epitype_labels = build_epitype_labels(obs_epi)

    for lab in ["SecA", "Transitioning", "SecB", "Ciliated"]:
        n = (epitype_labels == lab).sum()
        print(f"  {lab}: {n:,}")

    # ── 3. Collect reference barcodes ─────────────────────────
    print("\n[3] Collecting reference cell barcodes...")
    ref_by_sample = collect_reference_barcodes()

    # ── 4. Identify eligible samples ──────────────────────────
    print("\n[4] Identifying eligible samples...")

    epi_per_sample = obs_epi.groupby("sample_id").size()
    epi_per_sample = epi_per_sample.sort_values(ascending=False)

    manifest_rows = []
    eligible = []

    for sid in epi_per_sample.index:
        n_epi = int(epi_per_sample[sid])
        n_ref = len(ref_by_sample.get(sid, []))
        status = "eligible" if (n_epi >= MIN_EPI_CELLS and
                                n_ref >= MIN_REF_CELLS) else "skipped"

        # Count epitypes
        sample_mask = obs_epi["sample_id"] == sid
        sample_barcodes = obs_epi.index[sample_mask]
        sample_labels = epitype_labels.reindex(sample_barcodes).dropna()
        n_seca = int((sample_labels == "SecA").sum())
        n_trans = int((sample_labels == "Transitioning").sum())
        n_secb = int((sample_labels == "SecB").sum())
        n_cil = int((sample_labels == "Ciliated").sum())

        manifest_rows.append({
            "sample_id": sid,
            "study": obs_epi.loc[sample_barcodes[0], "study"] if len(sample_barcodes) > 0 else "",
            "patient_id": obs_epi.loc[sample_barcodes[0], "patient_id"] if len(sample_barcodes) > 0 else "",
            "n_epithelial": n_epi,
            "n_reference": n_ref,
            "n_secA": n_seca,
            "n_transitioning": n_trans,
            "n_secB": n_secb,
            "n_ciliated": n_cil,
            "status": status,
        })

        if status == "eligible":
            eligible.append(sid)

    manifest = pd.DataFrame(manifest_rows)
    print(f"  Total samples: {len(manifest)}")
    print(f"  Eligible (>={MIN_EPI_CELLS} epi + >={MIN_REF_CELLS} ref): {len(eligible)}")

    if args.n_samples == 0 or args.n_samples >= len(eligible):
        selected = list(eligible)
    else:
        selected = eligible[:args.n_samples]
    manifest["selected"] = manifest["sample_id"].isin(selected)
    print(f"  Selected for extraction: {len(selected)}")

    manifest_path = os.path.join(OUT_DIR, "tables", "sample_manifest.csv")
    manifest.to_csv(manifest_path, index=False)
    print(f"  Saved  {manifest_path}")

    # ── 5. Load epithelial counts into memory ─────────────────
    print("\n[5] Loading epithelial raw counts into memory...")
    t0 = time.time()

    # Read counts layer — this is ~4GB for 575k × 14k sparse int32
    counts_epi = adata_epi.layers["counts"][:]
    if not sparse.issparse(counts_epi):
        counts_epi = sparse.csr_matrix(counts_epi)
    else:
        counts_epi = counts_epi.tocsr()

    epi_barcodes_all = obs_epi.index.tolist()
    epi_barcode_to_idx = {b: i for i, b in enumerate(epi_barcodes_all)}
    print(f"  Loaded {counts_epi.shape[0]:,} × {counts_epi.shape[1]:,} "
          f"({time.time() - t0:.0f}s)")

    # Close backed epithelial h5ad
    del adata_epi

    # ── 6. Load reference counts per compartment ──────────────
    print("\n[6] Loading reference counts...")
    ref_counts_by_barcode = {}

    for name, path in H5AD_REF.items():
        print(f"  Loading {name} counts...")
        t0 = time.time()
        adata_ref = ad.read_h5ad(path, backed="r")

        # Only load counts for barcodes we need (in selected samples)
        needed_barcodes = set()
        for sid in selected:
            for b in ref_by_sample.get(sid, []):
                needed_barcodes.add(b)

        ref_obs = adata_ref.obs
        ref_barcodes_all = ref_obs.index.tolist()
        needed_idx = [i for i, b in enumerate(ref_barcodes_all)
                      if b in needed_barcodes]

        if needed_idx:
            # Read in chunks to manage memory
            ref_counts_raw = adata_ref.layers["counts"][needed_idx, :]
            if not sparse.issparse(ref_counts_raw):
                ref_counts_raw = sparse.csr_matrix(ref_counts_raw)
            else:
                ref_counts_raw = ref_counts_raw.tocsr()

            for local_i, global_i in enumerate(needed_idx):
                bc = ref_barcodes_all[global_i]
                ref_counts_by_barcode[bc] = ref_counts_raw[local_i]

        del adata_ref
        print(f"    {len(needed_idx):,} reference cells loaded ({time.time() - t0:.0f}s)")

    print(f"  Total reference cells in memory: {len(ref_counts_by_barcode):,}")

    # ── 7. Extract per-sample data ────────────────────────────
    print(f"\n[7] Extracting data for {len(selected)} samples...")
    t_extract = time.time()

    for i, sid in enumerate(selected):
        # Epithelial cells for this sample
        sample_mask = obs_epi["sample_id"] == sid
        sample_epi_barcodes = obs_epi.index[sample_mask].tolist()
        epi_idx = [epi_barcode_to_idx[b] for b in sample_epi_barcodes]
        epi_counts_sample = counts_epi[epi_idx, :]

        # Reference cells for this sample
        ref_bcs = ref_by_sample.get(sid, [])
        ref_bcs_in_mem = [b for b in ref_bcs if b in ref_counts_by_barcode]

        if ref_bcs_in_mem:
            ref_rows = sparse.vstack(
                [ref_counts_by_barcode[b] for b in ref_bcs_in_mem],
                format="csr"
            )
        else:
            ref_rows = sparse.csr_matrix((0, counts_epi.shape[1]))

        # Epitype labels dict
        ep_labels = epitype_labels.reindex(sample_epi_barcodes).to_dict()

        n_epi, n_ref = write_sample_data(
            sid, sample_epi_barcodes, ref_bcs_in_mem, ep_labels,
            epi_counts_sample, ref_rows, gene_names, OUT_DIR)

        if (i + 1) % 10 == 0 or i == 0:
            print(f"  [{i+1}/{len(selected)}] {sid}: {n_epi} epi + {n_ref} ref")

    elapsed = time.time() - t_extract
    print(f"\n  Extraction complete in {elapsed / 60:.1f} min")

    print(f"\n{'=' * 65}")
    print(f"  DONE — {len(selected)} samples extracted to {OUT_DIR}/per_sample/")
    print(f"  Next: bash scripts/19_cnv_run.sh 2 4 10    # pilot of 10")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
