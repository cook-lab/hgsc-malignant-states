#!/usr/bin/env bash
# CellBender 0.3.2 on the KPCA-B_B2 sample of the Cell Ranger multi run (cellranger_multi_config.csv), run on
# 2025-02-05 in the raw-data folder 2025_kpca/KPCA-B. The command is copied from the first lines of
# cellbender/KPCA-B_background_removed.log; every other setting is CellBender's default.
cellbender remove-background --input sample_raw_feature_bc_matrix.h5 --output ./cellbender/KPCA-B_background_removed.h5
