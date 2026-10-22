#!/usr/bin/env python3
"""Supplementary Table 9 — gene sets scored on the Xenium data (37 UCell pathway modules + the custom sets the
final figures use), with panel membership and where each set is reported.

Sources (read-only):
  RERUN/9b_scoring/pathway_gene_sets_v2.csv       the 37 modules as read by the final all-cell scoring
                                                  (spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R,
                                                  UCell maxRank = 477, raw counts, all cells, 11 SFEs)
  RERUN/9b_scoring/pathway_coverage_allcells.csv  per-sample genes requested / found / probe-QC FAIL + WARN
  RERUN/9b_scoring/scoring_summary_allcells.csv   per-sample cells scored
  spatial/04_pathway_scoring/04b_pathway_gene_sets.R       inlined canonical definitions (asserted identical)
  figures/figure5/01_fig5a_epithelial_programs_build.R          Fig 5a custom epithelial sets (DDR, EPIS)
  figures/figure5/08_fig5i_macrophage_programs_build.R     Fig 5i macrophage program sets (SETS)
  figures/figure6/06_fig6g_macrophage_secB_neighbours.py Fig 6g macrophage glycolytic vs resident program (GLYCO, RESID)
  Xenium panel: 202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv
                (477 genes, in_analysis_panel = the 449 genes not excluded by probe QC), cross-checked against
                the SFE gene list RERUN/fig6_viewer/degdata/WT01_genes.txt
Outputs: SuppTable9_ucell_gene_sets.xlsx (Module_summary, Gene_sets_long; + one CSV per sheet) via supp_xlsx.write_table.
Manuscript item(s): Table S9 (37 modules, 325 gene entries, and the custom sets of Figures 5a, 5i and 6g).
Run:  python 09_tableS9_ucell_gene_sets.py
"""
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supp_xlsx import write_table  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, ".."))
RERUN = os.environ.get("OUTPUT_ROOT", "output")
F_SETS = f"{RERUN}/9b_scoring/pathway_gene_sets_v2.csv"
F_COV = f"{RERUN}/9b_scoring/pathway_coverage_allcells.csv"
F_SUMM = f"{RERUN}/9b_scoring/scoring_summary_allcells.csv"
F_CANON = f"{REPO}/spatial/04_pathway_scoring/04b_pathway_gene_sets.R"
F_5A = f"{REPO}/figures/figure5/01_fig5a_epithelial_programs_build.R"
F_5I = f"{REPO}/figures/figure5/08_fig5i_macrophage_programs_build.R"
F_6G = f"{REPO}/figures/figure6/06_fig6g_macrophage_secB_neighbours.py"
F_PANEL = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202605_epitype_manuscript/supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv")
F_SFEGENES = f"{RERUN}/fig6_viewer/degdata/WT01_genes.txt"


def _r_vec(path, name):
    """Genes of an R assignment `name <- c("A","B")` or a list entry `name = c(...)`."""
    txt = open(path).read()
    m = re.search(rf'(?m)^\s*{name}\s*(?:<-|=)\s*c\(([^)]*)\)', txt)
    if m is None:
        raise KeyError(f"{name} not found in {path}")
    return re.findall(r'"([^"]+)"', m.group(1))


def _r_list(path, name):
    txt = open(path).read()
    m = re.search(rf'{name}\s*<-\s*list\((.*?)\)\s*\n\s*\w', txt, re.S)
    body = m.group(1) + ")"
    return {k: re.findall(r'"([^"]+)"', v) for k, v in re.findall(r'(\w+)\s*=\s*c\(([^)]*)\)', body)}


def _py_list(path, name):
    m = re.search(rf'{name}\s*=\s*\[([^\]]*)\]', open(path).read())
    return re.findall(r'"([^"]+)"', m.group(1))


# Category and one-line description for the 37 modules (grouping follows the Methods)
CAT = {
    "proliferation": ("Proliferation / cell cycle", "G2/M and S-phase proliferation markers"),
    "hypoxia": ("Hypoxia / metabolism", "HIF targets and hypoxia-responsive genes"),
    "glycolysis": ("Hypoxia / metabolism", "Glycolytic enzymes, glucose and lactate transport"),
    "lactate_metabolism": ("Hypoxia / metabolism", "Lactate production and export"),
    "oxidative_stress": ("Hypoxia / metabolism", "Antioxidant / NRF2 response"),
    "rtk_ras": ("Oncogenic signalling", "Receptor tyrosine kinase / RAS-MAPK pathway"),
    "pi3k_activating": ("Oncogenic signalling", "PI3K/AKT/mTOR pathway activators"),
    "pi3k_inhibitory": ("Oncogenic signalling", "PI3K/AKT/mTOR pathway inhibitors"),
    "myc_activating": ("Oncogenic signalling", "MYC network activators"),
    "myc_inhibitory": ("Oncogenic signalling", "MYC network antagonists (MXD/MNT/MGA)"),
    "wnt_activating": ("Developmental signalling", "Wnt/beta-catenin activators"),
    "wnt_inhibitory": ("Developmental signalling", "Wnt/beta-catenin inhibitors"),
    "notch_activating": ("Developmental signalling", "Notch ligands, receptors and co-activator"),
    "notch_inhibitory": ("Developmental signalling", "Notch repressors"),
    "hippo": ("Developmental signalling", "Hippo pathway"),
    "tgfb": ("Developmental signalling", "TGF-beta / activin pathway"),
    "p53_activating": ("DNA damage / p53", "p53 pathway activators (TP53, ATM, CHEK2, CDKN2A)"),
    "p53_inhibitory": ("DNA damage / p53", "p53 inhibitors (MDM2, MDM4)"),
    "apoptosis_pro": ("Cell death / survival", "Pro-apoptotic genes"),
    "apoptosis_anti": ("Cell death / survival", "Anti-apoptotic BCL2 family"),
    "survival_antiapoptotic": ("Cell death / survival", "Anti-apoptotic and complement-protective survival genes"),
    "nfkb": ("Inflammation", "NF-kB pathway and targets"),
    "jak_stat": ("Inflammation", "Cytokine receptors and JAK/STAT"),
    "type_i_ifn": ("Interferon", "Type I interferon response"),
    "type_ii_ifn": ("Interferon", "Type II (IFN-gamma) response"),
    "chemokine": ("Immune", "Chemokines and chemokine receptors"),
    "cytotoxicity": ("Immune", "Cytotoxic lymphocyte effector genes"),
    "antigen_presentation": ("Immune", "MHC class I antigen processing and presentation"),
    "checkpoint_inhibitory": ("Immune", "Inhibitory immune checkpoints"),
    "checkpoint_costimulatory": ("Immune", "Co-stimulatory immune receptors"),
    "complement_activation": ("Immune", "Complement components"),
    "complement_regulation": ("Immune", "Complement regulators (CD55, CLU)"),
    "emt": ("Stroma / ECM", "EMT and CAF-associated matrix genes"),
    "ecm_adhesion": ("Stroma / ECM", "Cell adhesion and integrins"),
    "matrix_remodeling": ("Stroma / ECM", "Metalloproteinases and matrix-remodelling enzymes"),
    "angiogenesis": ("Vasculature", "Endothelial / angiogenesis genes"),
    "adc_targets": ("Therapeutic targets", "Antibody-drug-conjugate target antigens"),
}

# Where each set is reported (first element: written to the figure_panels column; second: not written). Letters
# follow the manuscript legends; the final_panels/ file letter is given in brackets when it differs (Fig 5 files
# are lettered a-i; the Figure 5 legend runs a-j because it splits the ROI into c/d).
USE = {
    "hypoxia": ("Fig 4f (neighbourhood mean, 50 um), Fig 4j (ROI map); Fig 5f [file 5e] x-axis 'local epithelial hypoxia'; "
                "Fig 5g [file 5f] TMA core map; Fig5supp_macrophage_programs_hypoxia x-axis; Figure S16d",
                "¶69, ¶71 (Fig 4 legend lists all 12 genes), ¶75, ¶78, ¶137, ¶182"),
    "matrix_remodeling": ("Fig 4f; Fig 4k (ROI map; final_panels only, not in the legend); Figure S16c",
                          "¶69, ¶71 (legend lists the 6 genes), ¶182"),
    "nfkb": ("Fig 5a; Fig 5d [file 5c] ROI map; Figure S17f", "¶74, ¶78, ¶184"),
    "glycolysis": ("Fig 5a; Fig 5d [file 5c] ROI map; Figure S17g", "¶74, ¶78, ¶184"),
    "proliferation": ("Fig 5a; Fig 5d [file 5c] ROI map; Figure S17a", "¶74, ¶78, ¶184"),
    "myc_activating": ("Fig 5a 'Oncogenic growth' composite; Figure S17b 'Myc'", "¶74 (MYC, MAX, MYCL), ¶184"),
    "pi3k_activating": ("Fig 5a 'Oncogenic growth' composite; Figure S17c 'PI3K'", "¶74 (PIK3CA, AKT1, MTOR), ¶184"),
    "wnt_activating": ("Fig 5a 'Oncogenic growth' composite; Figure S17d 'Wnt'", "¶74 (FZD3, CTNNB1, TCF7L2), ¶184"),
    "notch_activating": ("Fig 5a 'Oncogenic growth' composite", "¶74 (JAG1, NOTCH1, NOTCH2)"),
    # p53_activating: scored, not shown (Figure S17e plots the custom ddr_repair set)
}


def main():
    sets = pd.read_csv(F_SETS)
    assert list(sets.columns) == ["pathway", "gene"]
    mods = list(dict.fromkeys(sets.pathway))
    assert len(mods) == 37, len(mods)
    # identical to the inlined canonical definitions (spatial/04_pathway_scoring/04b_pathway_gene_sets.R)
    canon_txt = open(F_CANON).read()
    for m in mods:
        assert _r_vec(F_CANON, m) == sets.loc[sets.pathway == m, "gene"].tolist(), m
    assert canon_txt.count(" = c(") >= 37

    panel = pd.read_csv(F_PANEL)
    assert len(panel) == 477 and panel.in_analysis_panel.sum() == 449
    sfe_genes = set(open(F_SFEGENES).read().split())
    assert sfe_genes == set(panel.gene), "panel CSV and SFE gene list differ"
    on_panel, in_449 = set(panel.gene), set(panel.gene[panel.in_analysis_panel])
    qc = panel.set_index("gene")["qc_status"]

    cov = pd.read_csv(F_COV)
    summ = pd.read_csv(F_SUMM)
    assert (cov.n_found == cov.n_requested).all()

    # ---------------- custom sets used by the final figures ----------------
    custom = []  # (module, group, category, genes, method, description, figures, text)
    custom += [
        ("ddr_repair", "Custom epithelial set (Fig 5a)", "DNA damage / p53", _r_vec(F_5A, "DDR"),
         "UCell (maxRank 477, raw counts), secretory cells", "DNA damage repair machinery",
         "Fig 5a 'DNA damage repair'; Figure S17e", "¶74 (BRCA1/2, ATM, CHEK2)"),
        ("epigenetic_silencing", "Custom epithelial set (Fig 5a)", "Epigenetic", _r_vec(F_5A, "EPIS"),
         "UCell (maxRank 477, raw counts), secretory cells", "Epigenetic silencing / chromatin remodelling",
         "Fig 5a 'Epigenetic silencing'", "¶74 ('epigenetic remodelling', DNMT3A, HDAC1, KDM5A)"),
    ]
    mac = _r_list(F_5I, "SETS")
    lab = {"glycolytic": ("Glycolytic", "Macrophage glycolytic program"),
           "nfkb_activated": ("NF-kB activated", "Macrophage NF-kB activation"),
           "tissue_resident": ("Tissue-resident identity", "Tissue-resident macrophage identity"),
           "chemokine": ("Chemokine", "Macrophage CXCR3-ligand chemokines")}
    for k, g in mac.items():
        custom.append((f"macrophage_{k}", "Macrophage program (Fig 5i)", "Macrophage", g,
                       "UCell (maxRank 477, raw counts) on all labelled cells, macrophages retained",
                       lab[k][1], f"Fig 5i [file 5h] '{lab[k][0]}'; Fig5supp_macrophage_programs_hypoxia"
                       + {"glycolytic": "; Figure S18a", "nfkb_activated": "; Figure S18b"}.get(k, ""),
                       "¶75, ¶78, ¶81 (glycolytic)"))
    custom += [
        ("macrophage_glyco_score_GLYCO", "Macrophage program score (Fig 6g)", "Macrophage", _py_list(F_6G, "GLYCO"),
         "Not UCell: log1p CP10k, each gene z-scored within the sample's macrophages, set mean; score = z(GLYCO) - z(RESID)",
         "Hypoxic / glycolytic macrophage genes (positive arm)", "Fig 6g (tertiles: resident-high vs hypoxic-high)",
         "¶81, ¶86"),
        ("macrophage_glyco_score_RESID", "Macrophage program score (Fig 6g)", "Macrophage", _py_list(F_6G, "RESID"),
         "Not UCell: negative arm of the Fig 6g score (see GLYCO)", "Tissue-resident macrophage genes (negative arm)",
         "Fig 6g", "¶81, ¶86"),
    ]
    assert mac == {"glycolytic": ["LDHA", "ENO1", "PGK1", "ALDOA"], "nfkb_activated": ["ICAM1", "NFKBIA", "BIRC3"],
                   "tissue_resident": ["A2M", "C1QC"], "chemokine": ["CXCL10", "CXCL11"]}

    # ---------------- long table ----------------
    rows = []
    for m in mods:
        for i, g in enumerate(sets.loc[sets.pathway == m, "gene"], 1):
            rows.append(("37 UCell pathway modules", m, CAT[m][0], i, g))
    for m, grp, cat, genes, *_ in custom:
        for i, g in enumerate(genes, 1):
            rows.append((grp, m, cat, i, g))
    long = pd.DataFrame(rows, columns=["set_group", "module", "category", "gene_index", "gene"])
    long["on_xenium_panel_477"] = long.gene.isin(on_panel).map({True: "yes", False: "no"})
    long["not_excluded_by_probe_QC"] = long.gene.isin(in_449).map({True: "yes", False: "no"})
    long["probe_qc_status"] = long.gene.map(qc).fillna("not on panel")

    # ---------------- module summary ----------------
    srows = []
    ncells = int(summ.n_cells_scored.sum())
    for m in mods:
        g = sets.loc[sets.pathway == m, "gene"].tolist()
        c = cov[cov.pathway == m]
        fail = sorted({x for s in c.probeqc_FAIL_genes.fillna("") for x in s.split(";") if x})
        warn = sorted({x for s in c.probeqc_WARN_genes.fillna("") for x in s.split(";") if x})
        fig, txt = USE.get(m, ("", ""))
        srows.append({"set_group": "37 UCell pathway modules", "module": m, "category": CAT[m][0],
                      "description": CAT[m][1], "n_genes": len(g), "n_on_panel_477": sum(x in on_panel for x in g),
                      "n_not_excluded_by_probe_QC": sum(x in in_449 for x in g), "genes": ", ".join(g),
                      "probe_qc_FAIL_genes_kept": ", ".join(fail), "probe_qc_WARN_genes": ", ".join(warn),
                      "scoring": f"UCell AddModuleScore_UCell, maxRank 477, raw counts, all {ncells:,} cells in 11 SFEs",
                      "figure_panels": fig,
                      "reported_in_final_figures": "yes" if fig else "no - scored, not shown"})
    for m, grp, cat, genes, method, desc, fig, txt in custom:
        srows.append({"set_group": grp, "module": m, "category": cat, "description": desc, "n_genes": len(genes),
                      "n_on_panel_477": sum(x in on_panel for x in genes),
                      "n_not_excluded_by_probe_QC": sum(x in in_449 for x in genes), "genes": ", ".join(genes),
                      "probe_qc_FAIL_genes_kept": ", ".join(x for x in genes if qc.get(x) == "Excluded from annotation"),
                      "probe_qc_WARN_genes": ", ".join(x for x in genes if str(qc.get(x, "")).startswith("Pass (")),
                      "scoring": method, "figure_panels": fig,
                      "reported_in_final_figures": "no - methods only" if fig.startswith("none") else "yes"})
    S = pd.DataFrame(srows)
    n_rep = int((S.set_group.eq("37 UCell pathway modules") & S.reported_in_final_figures.eq("yes")).sum())
    n_uniq = sets.gene.nunique()

    readme = [
        "Source (read-only): the 37 modules are read from 202609_spatial_rerun/9b_scoring/pathway_gene_sets_v2.csv, the file "
        "the final all-cell scoring read (spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R; finished 2026-09-22). All "
        "other copies found (fig6_singler/inputs/, 2026_final_xenium_analysis/output/9b_scoring/) are byte-identical, the list "
        "matches the definitions inlined in spatial/04_pathway_scoring/04b_pathway_gene_sets.R exactly (checked by this script), and "
        "the 325 module-gene pairs equal the May 2026 Supplemental_Table_6_UCell_pathway_gene_sets.csv.",
        f"37 modules, {len(sets)} module-gene entries, {n_uniq} unique genes. All {n_uniq} genes are on the 477-gene Xenium panel, so no "
        "gene was lost from any module; 5 of them (BAX, CLDN3, ERBB2, MMP7, VTCN1) fail probe QC and are among the 28 genes excluded from SingleR annotation, "
        "and were kept for scoring by design (FINAL_RUN_LOG 2026-09-21).",
        f"Scoring: UCell AddModuleScore_UCell (R), maxRank = 477 (full panel), raw counts(sfe), every cell (no subsampling) of "
        f"the 11 SFEs used in the final figures (core 8 whole tissues + CHUM TMA + COEUR TMA G/H; {ncells:,} cells). "
        "Scores are stored as pathway_<module> in the SFE colData. Figures then use these per-cell scores; the Fig 4f / 5a pooled GAMs "
        "are fit on all cells with tissue-equal weights (since 2026-09-28).",
        "Custom sets used by final figures are listed as separately labelled groups: Fig 5a epithelial sets (DNA damage repair, "
        "epigenetic silencing; UCell, same settings), Fig 5i macrophage programs (glycolytic, NF-kB activated, tissue-resident "
        "identity, chemokine; UCell), the Fig 6g macrophage program score (GLYCO vs RESID; NOT UCell: z-scored log1p CP10k "
        "means). Fig 5a also plots LDHA and TNF as single genes (not sets). The T/NK exhaustion composites of the old Fig 5d "
        "are not listed: no final figure or text uses them.",
        "Sheet Gene_sets_long: set_group, module, category, gene_index (order in the definition), gene, on_xenium_panel_477, "
        "not_excluded_by_probe_QC (gene is not among the 28 probe-QC exclusions; SingleR itself used 429 genes = 457 shared with the atlas reference - 28), probe_qc_status (Pass / Warn / Fail from the Xenium panel QC table).",
        "Sheet Module_summary: n genes, n on panel, genes, probe-QC FAIL/WARN genes kept (37 modules: from the scoring coverage file "
        "pathway_coverage_allcells.csv; custom sets: FAIL = 'Excluded from annotation', WARN = 'Pass (monitored/caution)' in the panel QC table), scoring method, figure panels "
        "that use each set, and whether it is reported. Figure letters follow "
        "the manuscript legends; [file x] gives the final_panels file letter where it differs (Fig 5 files run a-i, the legend a-j).",
        f"Count check: 37 modules (verified). {n_rep} of them appear in final figures (hypoxia, matrix_remodeling, nfkb, "
        "glycolysis, proliferation, myc/pi3k/wnt/notch_activating; Figures 4, 5, S16 and S17); the other "
        f"{37 - n_rep} were scored but are not reported, p53_activating among them (Figure S17e plots the custom ddr_repair "
        "set). No final Fig 6 panel uses a 37-module score. Fig 6g's macrophage glycolytic program (9 + 4 genes) is a "
        "different set from the 4-gene Fig 5i glycolytic UCell program.",
    ]
    write_table(9, "ucell_gene_sets",
                title="Supplementary Table 9. Gene sets scored on the Xenium spatial data: 37 UCell pathway modules and "
                      "custom epithelial and macrophage sets",
                readme=readme, sheets=[("Module_summary", S), ("Gene_sets_long", long)])
    print(S[["set_group", "module", "n_genes", "n_on_panel_477", "reported_in_final_figures"]].to_string())


if __name__ == "__main__":
    main()
