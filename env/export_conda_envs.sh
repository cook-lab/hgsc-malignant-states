#!/usr/bin/env bash
# ============================================================================
# Full exports of the two local conda envs that ran the repo's Python (E2, E3 in env/README.md)
# ----------------------------------------------------------------------------
# PURPOSE
#   Records every package (conda and pip) of
#     E2  the conda base env: all final-run Python and the downstream atlas steps (Python 3.12.9)
#     E3  conda env ovarian_atlas: atlas data prep, QC and Scrublet doublet scores (Python 3.12.3)
#   with `conda env export --no-builds`; the machine-specific prefix line is dropped and the env is
#   renamed so the file can be created next to an existing base env. Then rewrites the curated
#   environment.yml from E2 (env/make_environment_yml.py).
#
# INPUTS
#   the conda install that holds both envs (conda on PATH, or CONDA_EXE)
#
# OUTPUTS
#   env/E2_base_full.yml, env/E3_ovarian_atlas_full.yml, environment.yml
#
# MANUSCRIPT PANEL(S): none (environment record)
# RUNTIME TIER: fast (under a minute)
#
# USAGE (from the repo root)
#   bash env/export_conda_envs.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
CONDA="${CONDA_EXE:-conda}"
BASE="$("$CONDA" info --base)"

export_env() {   # <prefix> <output file> <env name written into the file>
  "$CONDA" env export --no-builds -p "$1" | sed -e '/^prefix:/d' -e "s/^name: .*/name: $3/" > "env/$2"
  echo "env/$2: $(grep -c '^ *- ' "env/$2") entries"
}
export_env "$BASE"                     E2_base_full.yml          hgsc-e2-base-full
export_env "$BASE/envs/ovarian_atlas"  E3_ovarian_atlas_full.yml hgsc-e3-atlas-qc-full
"$BASE/bin/python" env/make_environment_yml.py
