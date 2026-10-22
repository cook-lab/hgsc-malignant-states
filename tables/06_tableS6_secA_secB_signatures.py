#!/usr/bin/env python3
"""Supplementary Table 6 — SecA / SecB gene signatures used for polarization scoring, their NMF and pseudobulk support,
and the gene lists each downstream analysis actually used.

Sources (read-only):
  shared/signatures.yml                                   7 + 7 scoring signatures (source of truth)
  ATLAS/output/11d_epithelial_nmf/11d_nmf_loadings.csv    Factor_3 (SecA) / Factor_2 (SecB) loadings
  ATLAS/output/11e_nmf_characterization/11e_gene_classification.csv   177 SecA- / 177 SecB-specific genes (TCGA lists)
  RERUN/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv         Fig 2g pseudobulk DESeq2 (the copy Fig 6 reads;
                                                          byte-identical to the Fig 2g output in output_root/figure2/)
  ATLAS/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad (counts CSR indices) + ATLAS/output/fig_secretory_polarization/
      data/meta.parquet (schema_nmf, the Fig 2g epitypes)  -> % of cells expressing, all 575,366 cells, no subsampling
  Xenium panel: 202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv (477 genes)
  RERUN/06d_annotation_noBCAM/shared_gene_list_xenium_atlas.txt      457-gene Xenium x atlas UCell ranking space
  ATLAS/output/18_ucell_atlas/shared_gene_list.txt        atlas x organoid UCell ranking space (18b)
  ATLAS/data/cibersort_data_prev/tcga_ecotyper.txt        TCGA-OV expression used by 22d (gene presence only)
  ICON7 expr_gene_symbol.rds (config icon7_expr)          ICON7 array (gene presence only, via Rscript)
  Downstream scripts whose literal lists are parsed and asserted equal to signatures.yml (see Signature_usage sheet).
Outputs: SuppTable6_secA_secB_scoring_signatures.xlsx (+ one CSV per sheet) via supp_xlsx.write_table; the
  Scoring_signatures and Extended_lists CSVs are read by figures/supplementary/16_figS10_kpca_mouse_vs_human_programs.py.
Manuscript item(s): Table S6; Methods 457 / 14,271 genes and the 177 + 177 TCGA-OV lists.
Run:  python 06_tableS6_secA_secB_signatures.py
"""
import os
import re
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supp_xlsx import write_table  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, ".."))
ATLAS = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas")
RERUN = os.environ.get("OUTPUT_ROOT", "output")
XEN = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_xenium_analysis")
ICON7 = f"{XEN}/davids side quests/ICON7"
ORG = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202603_organoids/scripts")

F_SIG = f"{REPO}/shared/signatures.yml"
F_LOAD = f"{ATLAS}/output/11d_epithelial_nmf/11d_nmf_loadings.csv"
F_CLASS = f"{ATLAS}/output/11e_nmf_characterization/11e_gene_classification.csv"
F_DE = f"{RERUN}/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv"
F_DE_REPO = f"{REPO}/output/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv"
F_H5 = f"{ATLAS}/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad"
F_META = f"{ATLAS}/output/fig_secretory_polarization/data/meta.parquet"
F_PANEL = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv")
F_X457 = f"{RERUN}/06d_annotation_noBCAM/shared_gene_list_xenium_atlas.txt"
F_ORGSPACE = f"{ATLAS}/output/18_ucell_atlas/shared_gene_list.txt"
F_TCGA = f"{ATLAS}/data/cibersort_data_prev/tcga_ecotyper.txt"
F_ICON7 = f"{ICON7}/data/processed/expr_gene_symbol.rds"

# Scripts whose literal SecA/SecB lists were used by the final panels: (analysis, panels, script, SecA pattern, SecB pattern,
# how the list is obtained, scoring method)
S_18B = f"{ATLAS}/scripts/18b_ucell_atlas_scoring.R"
S_XEN = f"{REPO}/spatial/00_setup/00_setup.R"
S_VIS = f"{REPO}/figures/figure6/01_fig6abcd_visium_hd_bins.py"
S_ORG = f"{ORG}/00_setup.R"
S_ORG5 = f"{ORG}/04v5_secB_only_comparison.py"
S_ICON = f"{ICON7}/scripts/03_score_signatures.R"
S_11D = f"{ATLAS}/scripts/11d_epithelial_nmf.py"
S_V2 = f"{ATLAS}/scripts/18b_v2_ucell_atlas_scoring.R"

EXTRA = ["WT1", "BCAM", "MMP7"]  # named in text / 11d signature but not in the 7 + 7 scoring lists


def _vec(path, var, lang):
    txt = open(path).read()
    if lang == "R":
        m = re.search(rf'(?m)^\s*{var}\s*(?:<-|=)\s*c\(([^)]*)\)', txt)
    else:
        m = re.search(rf'{var}\s*=\s*\[([^\]]*)\]', txt)
    if m is None:
        raise KeyError(f"{var} not in {path}")
    return re.findall(r'"([^"]+)"', m.group(1))


def pct_expressing(genes):
    """% of SecA / Intermediate / SecB atlas epithelial cells (schema_nmf) with >= 1 count, all cells."""
    f = h5py.File(F_H5, "r")
    idx = [x.decode() if isinstance(x, bytes) else x for x in f["obs"][f["obs"].attrs["_index"]][:]]
    lab = pd.read_parquet(F_META, columns=["schema_nmf"])["schema_nmf"].reindex(idx).values
    var = [x.decode() if isinstance(x, bytes) else x for x in f["var"][f["var"].attrs["_index"]][:]]
    gi = {v: i for i, v in enumerate(var)}
    names = [g for g in genes if g in gi]
    lut = np.full(len(var), -1, np.int64)
    lut[[gi[g] for g in names]] = np.arange(len(names))
    grp = {"SecA": "SecA", "Transitioning": "Intermediate", "SecB": "SecB"}
    code = np.array([list(grp).index(x) if x in grp else -1 for x in lab])
    cnt = np.zeros((3, len(names)), np.int64)
    ncell = np.array([(code == k).sum() for k in range(3)])
    g = f["layers/counts"]
    ip = g["indptr"][:]
    B = 20000
    for s in range(0, len(idx), B):
        e = min(s + B, len(idx))
        ind = g["indices"][ip[s]:ip[e]]
        rows = np.repeat(np.arange(e - s), np.diff(ip[s:e + 1]))
        k = lut[ind]
        keep = k >= 0
        c = code[s:e][rows[keep]]
        ok = c >= 0
        np.add.at(cnt, (c[ok], k[keep][ok]), 1)
        if (s // B) % 5 == 0:
            print(f"  pct: {e:,}/{len(idx):,}", flush=True)
    out = pd.DataFrame({"gene": names})
    for j, nm in enumerate(grp.values()):
        out[f"pct_expr_atlas_{nm}"] = (100 * cnt[j] / ncell[j]).round(2)
    return out.set_index("gene"), dict(zip(grp.values(), ncell.tolist()))


def main():
    sig = yaml.safe_load(open(F_SIG))
    SECA, SECB = list(sig["SecA"]), list(sig["SecB"])
    assert len(SECA) == 7 and len(SECB) == 7
    assert "WT1" not in SECA and "FBXO21" in SECA

    # ---- literal lists in each downstream script, asserted equal to signatures.yml ----
    org_a = [g for g in _vec(S_ORG, "secA_genes", "R") if g != "BCAM"]   # 04_v4 / 08 apply setdiff(secA, "BCAM")
    checks = {
        "18b atlas": (_vec(S_18B, "secA_genes", "R"), _vec(S_18B, "secB_genes", "R")),
        "Visium HD": (_vec(S_VIS, "SECA", "py"), _vec(S_VIS, "SECB", "py")),
        "organoid 00_setup (minus BCAM)": (org_a, _vec(S_ORG, "secB_genes", "R")),
        "organoid 04v5 SecB": (SECA, _vec(S_ORG5, "SECB_GENES", "py")),
        "ICON7": (_vec(S_ICON, "SecA", "R"), _vec(S_ICON, "SecB", "R")),
    }
    for k, (a, b) in checks.items():
        assert sorted(a) == sorted(SECA) and sorted(b) == sorted(SECB), f"{k} list differs from signatures.yml"
    assert "signatures.yml" in open(S_XEN).read()
    canon11d_a, canon11d_b = _vec(S_11D, "SECA_GENES", "py"), _vec(S_11D, "SECB_GENES", "py")
    v2_a, v2_b = _vec(S_V2, "secA_genes", "R"), _vec(S_V2, "secB_genes", "R")

    # ---- NMF, classification, DE ----
    L = pd.read_csv(F_LOAD, index_col=0)
    f3, f2 = L.loc["Factor_3"], L.loc["Factor_2"]
    r3 = f3.rank(ascending=False, method="first").astype(int)
    r2 = f2.rank(ascending=False, method="first").astype(int)
    cls = pd.read_csv(F_CLASS, index_col=0).set_index("gene")["class"]
    ext_a, ext_b = cls.index[cls == "SecA-specific"].tolist(), cls.index[cls == "SecB-specific"].tolist()
    assert len(ext_a) == 177 and len(ext_b) == 177
    de = pd.read_csv(F_DE).set_index("gene")
    de_repo = pd.read_csv(F_DE_REPO).set_index("gene")
    assert de.equals(de_repo), "the two Fig 2g DESeq2 copies differ"

    panel = pd.read_csv(F_PANEL).set_index("gene")
    x457 = set(open(F_X457).read().split())
    orgspace = set(open(F_ORGSPACE).read().split())
    tcga_genes = set(pd.read_csv(F_TCGA, sep="\t", usecols=[0]).iloc[:, 0])
    icon7_genes = set(subprocess.run(["Rscript", "-e", f'writeLines(rownames(readRDS("{F_ICON7}")))'],
                                     capture_output=True, text=True, check=True).stdout.split())
    assert len(icon7_genes) > 19000

    allg = list(dict.fromkeys(SECA + SECB + EXTRA + ext_a + ext_b))
    pct, ncell = pct_expressing(allg)

    def row(g):
        return {
            "NMF_F3_SecA_loading": round(float(f3[g]), 5) if g in f3 else np.nan,
            "NMF_F3_rank": int(r3[g]) if g in r3 else np.nan,
            "NMF_F2_SecB_loading": round(float(f2[g]), 5) if g in f2 else np.nan,
            "NMF_F2_rank": int(r2[g]) if g in r2 else np.nan,
            "NMF_90th_pct_class": cls.get(g, "not in NMF genes"),
            "Fig2g_log2FC_SecB_vs_SecA": round(float(de.loc[g, "log2FoldChange"]), 4) if g in de.index else np.nan,
            "Fig2g_padj": float(f"{de.loc[g, 'padj']:.3g}") if g in de.index else np.nan,
            "Fig2g_baseMean": round(float(de.loc[g, "baseMean"]), 2) if g in de.index else np.nan,
            **{c: pct.loc[g, c] if g in pct.index else np.nan for c in pct.columns},
            "on_xenium_panel_477": "yes" if g in panel.index else "no",
            "xenium_probe_qc": panel.loc[g, "qc_status"] if g in panel.index else "",
        }

    # ---- Scoring_signatures ----
    rows = []
    for s, genes in (("SecA", SECA), ("SecB", SECB)):
        for g in genes:
            r = {"gene": g, "signature": s, **row(g)}
            r["in_xenium_atlas_ucell_space_457"] = "yes" if g in x457 else "no"
            r.update({
                "used_atlas_UCell_18b": "yes",
                "used_Xenium_UCell": "yes" if g in x457 else "no",
                "used_Visium_HD_zscore": "yes",
                "used_organoid_UCell": "yes" if g in orgspace else "no",
                "used_TCGA_22d": "yes" if g in (ext_a if s == "SecA" else ext_b) else "no (not factor-specific)",
                "used_ICON7_UCell": "yes" if g in icon7_genes else "no (not on array)",
            })
            rows.append(r)
    S = pd.DataFrame(rows)
    S["Fig2g_direction_as_expected"] = np.where(
        (S.signature == "SecA") & (S.Fig2g_log2FC_SecB_vs_SecA < 0) & (S.Fig2g_padj < 0.05), "yes",
        np.where((S.signature == "SecB") & (S.Fig2g_log2FC_SecB_vs_SecA > 0) & (S.Fig2g_padj < 0.05), "yes", "no"))

    # ---- genes named in text / earlier signatures but not in the scoring lists ----
    cites = {"WT1": "Results ¶49, methods ¶109, Supp Fig 6 legend ¶168, Fig 1h subtitle and SF5 label '(WT1, RCN2, PBX1)' "
                    "list WT1 as a SecA marker; it is not in any final 7-gene scoring list (FBXO21 is). On the Xenium "
                    "panel but excluded from SingleR annotation by probe QC. In the 11d NMF forced gene set and in the "
                    "177 SecA-specific (TCGA) list.",
             "BCAM": "In the 11d NMF forced SecA set and the organoid base list; dropped from all scoring ('noBCAM'). "
                     "Shared class (not factor-specific), so not in the TCGA list.",
             "MMP7": "In the 11d NMF forced SecB set; not in the 7-gene SecB scoring list. SecB-specific, so in the TCGA list. "
                     "Xenium probe-QC FAIL (excluded from annotation)."}
    X = pd.DataFrame([{"gene": g, "role": cites[g], **row(g)} for g in EXTRA])

    # ---- Extended_lists: 177 + 177 factor-specific genes (TCGA-OV Fig 7g-j) ----
    E = pd.DataFrame([{"gene": g, "list": c, **row(g), "present_in_TCGA_expression": "yes" if g in tcga_genes else "no",
                       "in_7gene_scoring_signature": "yes" if g in SECA + SECB else ""}
                      for c, gl in (("SecA-specific", ext_a), ("SecB-specific", ext_b)) for g in gl])
    E["_s"] = np.where(E["list"] == "SecA-specific", E.NMF_F3_rank, E.NMF_F2_rank)
    E = E.sort_values(["list", "_s"]).drop(columns="_s").reset_index(drop=True)
    n_tcga_a = int(((E["list"] == "SecA-specific") & (E.present_in_TCGA_expression == "yes")).sum())
    n_tcga_b = int(((E["list"] == "SecB-specific") & (E.present_in_TCGA_expression == "yes")).sum())

    # ---- Signature_usage ----
    U = pd.DataFrame([
        ("Atlas scRNA-seq UCell (SecA_UCell / SecB_UCell)", "Fig 3b/c-e atlas reference lines",
         "202601_final_atlas/scripts/18b_ucell_atlas_scoring.R -> output/18_ucell_atlas/atlas_ucell_scores.csv",
         "7 + 7 (hardcoded, identical to signatures.yml)", "UCell, maxRank 1500, raw counts, atlas x organoid gene space "
         f"({len(orgspace):,} genes)"),
        ("Atlas epitypes for Fig 1-2", "Fig 1h-m, Fig 2 (SecA / Intermediate / SecB)", "11d_nmf_usage.csv Factor_2 percentiles",
         "none (NMF Factor_2 usage, not a gene signature)", "SecA < p50, Intermediate p50-p75, SecB > p75 of Factor_2 usage"),
        ("Fig 2g highlighted markers", "Fig 2g", "figures/figure2/04_fig2g_pseudobulk_volcano.py",
         "7 + 7 (reads signatures.yml)", "pseudobulk DESeq2, 261 tumours / 522 profiles, design ~sample + epitype (paired within tumours), BH FDR"),
        ("Xenium polarization (SecA_UCell, SecB_UCell, polarization_UCell = SecB - SecA)",
         "Fig 4-7 Xenium panels, SecA/Intermediate/SecB labels, Figures S15-S18",
         "spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R (lists from spatial/00_setup/00_setup.R, which reads "
         "signatures.yml); outputs 202609_spatial_rerun/06d_annotation_noBCAM/", "7 + 7 (reads signatures.yml)",
         "UCell on the 457-gene Xenium x atlas shared space (maxRank 457); atlas secretory cells re-scored in the same "
         "space for thresholds (t_low 0.109, t_high 0.211)"),
        ("Visium HD polarization", "Fig 6a-d, Visium rows of Fig 6h", "figures/figure6/01_fig6abcd_visium_hd_bins.py",
         "7 + 7 (hardcoded, identical)", "NOT UCell: mean per-gene z of log1p CP10k, SecB minus SecA; tertiles among "
         "epithelial bins"),
        ("Patient-derived organoids / OPTO98", "Fig 3b, 3c-e", "202603_organoids/scripts/04_secretory_polarization_v4_UCell_"
         "aligned.R, 08_opto98_ucell_aligned.R, 04v5_secB_only_comparison.py (00_setup.R lists)",
         "7 + 7 (hardcoded 8 SecA incl. BCAM, BCAM removed with setdiff)", "UCell, maxRank 1500, atlas-aligned gene space"),
        ("External in vitro models", "Fig 3f-h", "202603_organoids/scripts/12_gse284507_spheroid_organoid_secAB.R",
         "7 + 7 (hardcoded)", "as scripted (also HIF and proliferation modules)"),
        ("TCGA-OV survival", "Fig 7g-j", "202601_final_atlas/scripts/22d_signature_survival.py -> output/22_tcga_deconvolution/"
         "22d_signature_scores.csv (repo port atlas/07_deconvolution_survival/05_signature_survival.py)",
         f"177 SecA-specific + 177 SecB-specific from 11e_gene_classification.csv (all present in TCGA: {n_tcga_a}/{n_tcga_b})",
         "NOT UCell: per-gene z across samples, mean per list, polarization = SecA - SecB (inverted in Fig 7); 323 samples, "
         "tertiles 108/107/108"),
        ("ICON7 (GSE140082)", "Fig 7k-n; figure_icon7_bevacizumab", "ICON7/scripts/03_score_signatures.R -> "
         "data/processed/cohort_filtered.tsv (polarization_ucell)", "7 + 7 (hardcoded, identical)",
         "UCell ScoreSignatures_UCell (default maxRank) on the gene-symbol array, z-scored, SecB_z - SecA_z; median split"),
        ("NMF gene forcing (not scoring)", "Supp Table 5", "202601_final_atlas/scripts/11d_epithelial_nmf.py",
         f"{len(canon11d_a)} SecA ({', '.join(canon11d_a)}) + {len(canon11d_b)} SecB ({', '.join(canon11d_b)})",
         "forced into the 3,000 HVGs (3,007 genes)"),
        ("Fig 1g polarization PC alignment (organoid-derived sets)", "Fig 1g",
         "variance-decomposition analysis (lib_decompose.py, outside the repo); same lists as "
         "202601_final_atlas/scripts/18b_v2_ucell_atlas_scoring.R",
         f"SecA {', '.join(v2_a)}; SecB {', '.join(v2_b)}",
         "NOT UCell: scanpy score_genes, SecB - SecA, per tumour (250 tumours)"),
        ("Mouse KPCA.B polarization", "Fig 3j", "figures/figure3/04_fig3j_kpca_polarization_7plus7.py",
         "7 + 7 mouse 1:1 orthologues (MGI): SecA Mecom, Fbxo21, Lgr5, Lpar3, Pbx1, Sox17, Rcn2; SecB Krt17, Krt19, Krt7, "
         "Lcn2, Prss22, Slpi, Tacstd2",
         "UCell (KPCA.B project implementation: maxRank 1500, ties in gene order) over all genes of the 12,429 tumour cells; "
         "SecB - SecA. Replaces the earlier 5 + 4 score (Bcam, Mecom, Pbx1, Sox17, Rcn2 vs Krt19, Krt7, Lcn2, Prss22; "
         "Spearman 0.98 with the new score)"),
    ], columns=["analysis", "figure_panels", "script / output", "gene lists", "scoring method"])

    readme = [
        "Source (read-only): scoring signatures = shared/signatures.yml (SecA: " + ", ".join(SECA) + "; SecB: " +
        ", ".join(SECB) + "). NMF loadings = 202601_final_atlas/output/11d_epithelial_nmf/11d_nmf_loadings.csv (Supp Table 5); "
        "gene classes = 11e_nmf_characterization/11e_gene_classification.csv; Fig 2g = 202609_spatial_rerun/figure2/"
        "atlas_volcano_secA_secB_pseudobulk_deseq2.csv (the copy Fig 6 reads; checked identical to "
        "hgsc-malignant-states/output/figure2/).",
        "Derivation: the 7 + 7 genes are the SecA (Factor_3) and SecB (Factor_2) markers of the 11d NMF after dropping BCAM "
        "('noBCAM'); every one is significant in the expected direction in the Fig 2g pseudobulk DESeq2 (design "
        "~sample + epitype, SecB and SecA paired within each tumour; column Fig2g_direction_as_expected; log2FC > 0 = "
        "higher in SecB). Fig2g_padj = 0 means p below the double-precision floor (< 1e-300).",
        "Sheet Scoring_signatures: the 14 scoring genes with Factor_3 / Factor_2 loading and within-factor rank (of 3,007), "
        "90th-percentile class, Fig 2g log2FC / padj / baseMean, % of atlas epithelial cells with >= 1 count in SecA / "
        f"Intermediate / SecB (schema_nmf labels as in Fig 2g; all cells: SecA {ncell['SecA']:,}, Intermediate "
        f"{ncell['Intermediate']:,}, SecB {ncell['SecB']:,}), Xenium panel membership and probe-QC status, and yes/no for "
        "each downstream analysis (Signature_usage explains each).",
        "Sheet Signature_usage: for every downstream analysis, the final script, the list it used, and the scoring method. "
        "This script parses the literal lists in the 18b, Visium, organoid and ICON7 scripts and checks they equal "
        "signatures.yml; the Xenium and Fig 2g scripts read signatures.yml directly.",
        "Sheet Extended_lists: the 177 SecA-specific and 177 SecB-specific genes (Supp Table 5 classification) that are the "
        f"TCGA-OV signatures of Fig 7g-j ({n_tcga_a}/177 and {n_tcga_b}/177 present in the TCGA expression matrix), with the "
        "same annotation columns.",
        "Sheet Other_named_markers: WT1, BCAM and MMP7, which appear in the text or earlier signatures but are not in the "
        "scoring lists.",
    ]
    write_table(6, "secA_secB_scoring_signatures",
                title="Supplementary Table 6. SecA and SecB gene signatures used for polarization scoring",
                readme=readme,
                sheets=[("Scoring_signatures", S), ("Signature_usage", U), ("Extended_lists", E),
                        ("Other_named_markers", X)])
    print(S[["gene", "signature", "NMF_F3_rank", "NMF_F2_rank", "NMF_90th_pct_class", "Fig2g_log2FC_SecB_vs_SecA",
             "Fig2g_padj", "pct_expr_atlas_SecA", "pct_expr_atlas_SecB", "Fig2g_direction_as_expected"]].to_string())


if __name__ == "__main__":
    main()
