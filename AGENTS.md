# AGENTS.md: working in this repository

Orientation for people and coding agents. Read [`README.md`](README.md) first. This is a code repository; the data
are deposited separately ([`data/README.md`](data/README.md)).

## Layout in one line each

- `config/`: `config.yml` holds paths, seed, cohort and polarization thresholds. `config.py` and `config.R` load it.
  `study_codes.{py,R}` translate lab identifiers to study codes.
- `shared/`: `signatures.yml` (the 7 + 7 SecA / SecB genes) and `cellassign_markers_v3.csv`.
- `atlas/`, `spatial/`: backend stages, numbered in run order. Each stage README lists its scripts with inputs and
  outputs.
- `models/`, `icon7/`, `*/as_run/`: trimmed copies of the scripts that ran in other lab project folders. The code that
  wrote the outputs the paper uses is kept as it ran; plots and analyses not in the paper were removed, and each trimmed
  copy opens with a note saying so. They keep their original project paths. Do not edit them; they record what was
  run.
- `figures/figureN/`: the scripts of one figure, named `NN_fig<number><panels>_<what>` in run order.
  `figures/supplementary/` holds Figures S1–S18 and `figures/_prep/` the shared loaders and caches.
- `tables/`: Supplementary Tables S1–S12. `build_all_tables.sh` runs them in order.
- `figures/README.md`: the map from every panel, table and reported number to its scripts.

## Conventions

- **Paths:** resolve inputs and outputs through `config/`, never with new absolute paths.
  - Python: `from config.config import CFG, obj, path, SEED`
  - R: `source("config/config.R")`, then `cfg_obj()` and `cfg_path()`
  - Inputs sit under `DATA_ROOT` (the deposit) and outputs go to `OUTPUT_ROOT` (default `./output`; use one absolute
    path for every step). Run scripts from the repository root.
  - The copies (`as_run/`, `models/`, `icon7/`) are the exception: they keep the paths of the project they ran in.
- **R locale:** R scripts need `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8`, otherwise `config.yml` is read incorrectly.
- **States:** SecA / Intermediate / SecB / Ciliated.
  - Atlas: NMF Factor 2 usage, with SecA < p50, Intermediate p50–p75 and SecB ≥ p75.
  - Xenium: UCell polarization with frozen thresholds, SecA ≤ 0.109429 and SecB ≥ 0.210501.
- **Signatures:** read SecA / SecB genes from `shared/signatures.yml`. Do not inline gene lists.
- **Spatial cell labels:** use the SingleR label tables written by
  `spatial/03_annotation_polarization/06_extract_singler.R` (`cell_label_singler`). Do not use `colData(sfe)$cell_label`.
- **TMA cohort:** the final core set (378 patients, 503 cores). It is written by `spatial/02_qc/06b_final_tma_cores.py`
  and read with `figures/_prep/final_tma_cores.py` or `spatial/00_setup/final_tma_cores.R` (`final_core_ids()`).
- **De-identification:** samples appear only as study codes (WT01–WT08, TMA-P###, FTE01–15). Never write other sample
  identifiers into code, file names or outputs.
- **Seeds:** stochastic steps use `config.yml` `seed` (42).
- **Script headers:** repository scripts open with a header on their purpose, inputs, outputs, manuscript panels and
  runtime, most as PURPOSE / INPUTS / OUTPUTS / MANUSCRIPT PANEL(S) / RUNTIME TIER. Use that format when you add one.
  The copies have their own headers, from the project they ran in.
