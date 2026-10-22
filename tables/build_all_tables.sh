#!/usr/bin/env bash
# Rebuild every manuscript Supplementary Table (1-12) from source, in order.
# Outputs: $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable{N}_*.xlsx (+ one CSV per sheet).
# ST6 reads all 575,366 epithelial cells.
# ST8 writes its de-identification key ONLY to $OUTPUT_ROOT/private_keys/ (never to the repo or the publication folder).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}; RS=${RS:-Rscript}
run() { echo; echo "=== $*"; "$@"; }
run $PY 01_tableS1_atlas_studies.py                                   # Table 1  study characteristics
run $PY 02_tableS2-4_annotation_markers.py                            # Tables 2-4: annotation-stage Wilcoxon markers (published)
run $PY 05_tableS5_epithelial_nmf_programs.py                                # Table 5 NMF programs
run $PY 06_tableS6_secA_secB_signatures.py                           # Table 6 SecA/SecB scoring genes
run $PY 07_tableS7_xenium_gene_panel.py                                      # Table 7 Xenium panel
run $PY 08_tableS8_xenium_cohort.py                                 # Table 8 Xenium clinical cohort (de-identified)
run $PY 09_tableS9_ucell_gene_sets.py                                        # Table 9 gene sets
run $PY 10_tableS10_secb_niche_macrophage_degs.py                            # Table 10 SecB-niche macrophage DEGs
run $PY 11_tableS11_xenium_tma_cox.py                          # Table 11 Xenium TMA multivariable Cox
run $PY 12_tableS12_icon7_cox.py                                             # Table 12 ICON7 Cox
echo; echo "ALL SUPP TABLES BUILT"
