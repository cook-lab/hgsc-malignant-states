#!/usr/bin/env python3
"""
CNV — extract per-sample count matrices for CopyKAT
===================================================
HGSC malignant-states atlas backend.

Extracts per-sample raw count matrices (epithelial + reference cells) for
CopyKAT CNV inference. Reference cells = fibroblast + T/NK + endothelial from
the same sample, capped at MAX_REF_CELLS/sample (stratified across the three
compartments) so the diploid baseline stays balanced and runtime bounded.
Epithelial cells are labeled with the NMF percentile epitypes
(SecA / Intermediate / SecB / Ciliated).

INPUTS (<data_root>/2026_final_atlas/celltype_h5ad/):
  - hgsc_atlas_final_epithelial.h5ad
  - hgsc_atlas_final_fibroblast.h5ad
  - hgsc_atlas_final_t_nk_cell.h5ad
  - hgsc_atlas_final_endothelial.h5ad
  - output_root/03_epithelial_nmf/11d_nmf_usage.csv  (Factor_2)

OUTPUTS (output_root/05_cnv/):
  - per_sample/<sample_id>/{counts.mtx.gz,genes.txt,barcodes.csv,ref_barcodes.txt}
  - tables/sample_manifest.csv

MANUSCRIPT PANELS: upstream of Fig 1J, SF4C, SF7 (CopyKAT CNV chain).

RUNTIME TIER: heavy (loads full epithelial counts into memory + ref slices).

SEEDING: reference subsample uses np.random.default_rng(SEED) for determinism.

Usage:
    python 01_cnv_extract.py                 # all eligible
    python 01_cnv_extract.py --n-samples 10  # pilot
"""

import argparse
import gzip
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

# ── Paths ────────────────────────────────────────────────────
_CT = lambda f: path("data_root", "2026_final_atlas", "celltype_h5ad", f)  # noqa: E731
H5AD_EPI = _CT("hgsc_atlas_final_epithelial.h5ad")
H5AD_REF = {
    "fibroblast":  _CT("hgsc_atlas_final_fibroblast.h5ad"),
    "t_nk_cell":   _CT("hgsc_atlas_final_t_nk_cell.h5ad"),
    "endothelial": _CT("hgsc_atlas_final_endothelial.h5ad"),
}
NMF_CSV = path("output_root", "03_epithelial_nmf", "11d_nmf_usage.csv")
OUT_DIR = path("output_root", "05_cnv")

MIN_EPI_CELLS = 200
MIN_REF_CELLS = 50
MAX_REF_CELLS = 500            # cap per sample
DEFAULT_N_SAMPLES = 0          # 0 = all eligible


def build_epitype_labels(obs_df):
    """Assign NMF percentile epitype labels to epithelial cells."""
    nmf_usage = pd.read_csv(NMF_CSV, index_col=0)
    f2 = nmf_usage["Factor_2"]
    common = obs_df.index.intersection(f2.index)
    f2 = f2.loc[common]

    is_ciliated = obs_df.loc[common, "celltype_level2"] == "Ciliated epithelial cell"
    non_cil_f2 = f2[~is_ciliated]
    p50 = np.percentile(non_cil_f2, 50)
    p75 = np.percentile(non_cil_f2, 75)
    print(f"  NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

    labels = pd.Series("SecA", index=common)
    labels[f2 >= p50] = "Intermediate"
    labels[f2 >= p75] = "SecB"
    labels[is_ciliated] = "Ciliated"
    return labels


def collect_reference_barcodes():
    """Per-sample ref barcodes, capped at MAX_REF_CELLS (stratified by compartment)."""
    ref_by_compartment = {name: {} for name in H5AD_REF}
    for name, p in H5AD_REF.items():
        print(f"  Loading {name} (backed)...")
        adata_ref = ad.read_h5ad(p, backed="r")
        for sid, grp in adata_ref.obs[["sample_id"]].groupby("sample_id"):
            ref_by_compartment[name].setdefault(sid, [])
            ref_by_compartment[name][sid].extend(grp.index.tolist())
        del adata_ref

    all_sids = set()
    for d in ref_by_compartment.values():
        all_sids.update(d.keys())

    rng = np.random.default_rng(SEED)
    ref_by_sample = {}
    total_pre = total_post = 0

    for sid in all_sids:
        per_comp = {n: ref_by_compartment[n].get(sid, []) for n in H5AD_REF}
        n_total = sum(len(v) for v in per_comp.values())
        total_pre += n_total

        if n_total <= MAX_REF_CELLS:
            picked = [b for v in per_comp.values() for b in v]
        else:
            picked = []
            remaining = MAX_REF_CELLS
            comps = sorted(per_comp.items(), key=lambda x: len(x[1]))
            for i, (_n, bcs) in enumerate(comps):
                share = remaining // (len(comps) - i)
                if len(bcs) <= share:
                    picked.extend(bcs)
                    remaining -= len(bcs)
                else:
                    idx = rng.choice(len(bcs), size=share, replace=False)
                    picked.extend([bcs[j] for j in idx])
                    remaining -= share
        ref_by_sample[sid] = picked
        total_post += len(picked)

    print(f"  Reference cells: {total_pre:,} raw -> {total_post:,} capped at "
          f"{MAX_REF_CELLS}/sample across {len(ref_by_sample)} samples")
    return ref_by_sample


def write_sample_data(sample_id, epi_barcodes, ref_barcodes, epitype_labels,
                      epi_counts, ref_counts, gene_names, out_base):
    """Write per-sample CopyKAT input files."""
    sample_dir = os.path.join(out_base, "per_sample", sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    all_counts = sparse.vstack([epi_counts, ref_counts], format="csc")
    all_barcodes = list(epi_barcodes) + list(ref_barcodes)
    counts_t = all_counts.T.tocsc()  # CopyKAT expects genes x cells

    with gzip.open(os.path.join(sample_dir, "counts.mtx.gz"), "wb") as f:
        sio.mmwrite(f, counts_t)
    with open(os.path.join(sample_dir, "genes.txt"), "w") as f:
        f.write("\n".join(gene_names))

    bc_df = pd.DataFrame({
        "barcode": all_barcodes,
        "epitype": [epitype_labels.get(b, "reference") for b in all_barcodes],
        "is_reference": [b in set(ref_barcodes) for b in all_barcodes],
    })
    bc_df.to_csv(os.path.join(sample_dir, "barcodes.csv"), index=False)
    with open(os.path.join(sample_dir, "ref_barcodes.txt"), "w") as f:
        f.write("\n".join(ref_barcodes))
    return len(epi_barcodes), len(ref_barcodes)


def main():
    parser = argparse.ArgumentParser(description="Extract per-sample CopyKAT inputs")
    parser.add_argument("--n-samples", type=int, default=DEFAULT_N_SAMPLES,
                        help="Number of top samples to extract (0 = all eligible)")
    args = parser.parse_args()

    os.makedirs(os.path.join(OUT_DIR, "per_sample"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "tables"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "figs"), exist_ok=True)

    limit_msg = "all eligible" if args.n_samples == 0 else f"top {args.n_samples}"
    print("=" * 65)
    print("  CNV — Extract per-sample CopyKAT inputs")
    print(f"  Target: {limit_msg} samples   Ref cap: {MAX_REF_CELLS}/sample")
    print("=" * 65)

    print("\n[1] Loading epithelial h5ad...")
    t0 = time.time()
    adata_epi = ad.read_h5ad(H5AD_EPI, backed="r")
    obs_epi = adata_epi.obs.copy()
    gene_names = adata_epi.var_names.tolist()
    print(f"  {len(obs_epi):,} cells x {len(gene_names):,} genes  ({time.time()-t0:.0f}s)")

    print("\n[2] Building NMF percentile epitype labels...")
    epitype_labels = build_epitype_labels(obs_epi)
    for lab in ["SecA", "Intermediate", "SecB", "Ciliated"]:
        print(f"  {lab}: {(epitype_labels == lab).sum():,}")

    print("\n[3] Collecting reference cell barcodes...")
    ref_by_sample = collect_reference_barcodes()

    print("\n[4] Identifying eligible samples...")
    epi_per_sample = obs_epi.groupby("sample_id").size().sort_values(ascending=False)
    manifest_rows = []
    eligible = []
    for sid in epi_per_sample.index:
        n_epi = int(epi_per_sample[sid])
        n_ref = len(ref_by_sample.get(sid, []))
        status = "eligible" if (n_epi >= MIN_EPI_CELLS and n_ref >= MIN_REF_CELLS) else "skipped"

        sample_barcodes = obs_epi.index[obs_epi["sample_id"] == sid]
        sample_labels = epitype_labels.reindex(sample_barcodes).dropna()
        manifest_rows.append({
            "sample_id": sid,
            "study": obs_epi.loc[sample_barcodes[0], "study"] if len(sample_barcodes) else "",
            "patient_id": obs_epi.loc[sample_barcodes[0], "patient_id"] if len(sample_barcodes) else "",
            "n_epithelial": n_epi, "n_reference": n_ref,
            "n_secA": int((sample_labels == "SecA").sum()),
            "n_intermediate": int((sample_labels == "Intermediate").sum()),
            "n_secB": int((sample_labels == "SecB").sum()),
            "n_ciliated": int((sample_labels == "Ciliated").sum()),
            "status": status,
        })
        if status == "eligible":
            eligible.append(sid)

    manifest = pd.DataFrame(manifest_rows)
    print(f"  Total samples: {len(manifest)}")
    print(f"  Eligible (>={MIN_EPI_CELLS} epi + >={MIN_REF_CELLS} ref): {len(eligible)}")

    selected = list(eligible) if (args.n_samples == 0 or args.n_samples >= len(eligible)) \
        else eligible[:args.n_samples]
    manifest["selected"] = manifest["sample_id"].isin(selected)
    print(f"  Selected for extraction: {len(selected)}")

    manifest_path = os.path.join(OUT_DIR, "tables", "sample_manifest.csv")
    manifest.to_csv(manifest_path, index=False)
    print(f"  Saved  {manifest_path}")

    print("\n[5] Loading epithelial raw counts into memory...")
    t0 = time.time()
    counts_epi = adata_epi.layers["counts"][:]
    counts_epi = counts_epi.tocsr() if sparse.issparse(counts_epi) else sparse.csr_matrix(counts_epi)
    epi_barcodes_all = obs_epi.index.tolist()
    epi_barcode_to_idx = {b: i for i, b in enumerate(epi_barcodes_all)}
    print(f"  Loaded {counts_epi.shape[0]:,} x {counts_epi.shape[1]:,} ({time.time()-t0:.0f}s)")
    del adata_epi

    print("\n[6] Loading reference counts...")
    ref_counts_by_barcode = {}
    for name, p in H5AD_REF.items():
        print(f"  Loading {name} counts...")
        t0 = time.time()
        adata_ref = ad.read_h5ad(p, backed="r")
        needed_barcodes = {b for sid in selected for b in ref_by_sample.get(sid, [])}
        ref_barcodes_all = adata_ref.obs.index.tolist()
        needed_idx = [i for i, b in enumerate(ref_barcodes_all) if b in needed_barcodes]
        if needed_idx:
            ref_counts_raw = adata_ref.layers["counts"][needed_idx, :]
            ref_counts_raw = ref_counts_raw.tocsr() if sparse.issparse(ref_counts_raw) \
                else sparse.csr_matrix(ref_counts_raw)
            for local_i, global_i in enumerate(needed_idx):
                ref_counts_by_barcode[ref_barcodes_all[global_i]] = ref_counts_raw[local_i]
        del adata_ref
        print(f"    {len(needed_idx):,} reference cells loaded ({time.time()-t0:.0f}s)")
    print(f"  Total reference cells in memory: {len(ref_counts_by_barcode):,}")

    print(f"\n[7] Extracting data for {len(selected)} samples...")
    t_extract = time.time()
    for i, sid in enumerate(selected):
        sample_epi_barcodes = obs_epi.index[obs_epi["sample_id"] == sid].tolist()
        epi_idx = [epi_barcode_to_idx[b] for b in sample_epi_barcodes]
        epi_counts_sample = counts_epi[epi_idx, :]

        ref_bcs_in_mem = [b for b in ref_by_sample.get(sid, []) if b in ref_counts_by_barcode]
        if ref_bcs_in_mem:
            ref_rows = sparse.vstack([ref_counts_by_barcode[b] for b in ref_bcs_in_mem],
                                     format="csr")
        else:
            ref_rows = sparse.csr_matrix((0, counts_epi.shape[1]))

        ep_labels = epitype_labels.reindex(sample_epi_barcodes).to_dict()
        n_epi, n_ref = write_sample_data(sid, sample_epi_barcodes, ref_bcs_in_mem,
                                         ep_labels, epi_counts_sample, ref_rows,
                                         gene_names, OUT_DIR)
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  [{i+1}/{len(selected)}] {sid}: {n_epi} epi + {n_ref} ref")

    print(f"\n  Extraction complete in {(time.time()-t_extract)/60:.1f} min")
    print(f"\n{'=' * 65}")
    print(f"  DONE — {len(selected)} samples extracted to {OUT_DIR}/per_sample/")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
