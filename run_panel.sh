#!/usr/bin/env bash
# ============================================================================
# run_panel.sh — uniform launcher for HGSC figure/table generator scripts.
# ----------------------------------------------------------------------------
# Used by the final reproducibility / panel-export validation. Runs ONE
# generator script in the correct environment, captures exit code + wall time,
# tees full stdout/stderr to a per-script log, and reports the output files the
# script newly wrote/updated under OUTPUT_ROOT.
#
#   Usage:  run_panel.sh <relpath-from-repo-root> [script args...]
#   e.g.    run_panel.sh figures/figure2/03_atlas_volcano_secA_secB.py
#           run_panel.sh figures/figure4/02_xenium_whole_tissue_snapshot.R all
#
# Design notes:
#   - OUTPUT_ROOT is forced ABSOLUTE so outputs land in repo/output regardless
#     of the CWD the script is launched from.
#   - R scripts use three config-resolution idioms; CWD is chosen per-pattern:
#       * relative  source(file.path("..","..",...))  -> needs CWD = script dir
#       * everything else (hardened --file= / repo-root fallback) -> CWD = repo
#   - Python scripts self-insert the repo root onto sys.path, so CWD = repo.
#   - A soft timeout (PANEL_TIMEOUT, default 3600s) guards against hangs.
# ============================================================================
set -uo pipefail

# Repo root = directory containing this script (portable; no hardcoded location).
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Interpreters: override PY/RS to point at a specific env; otherwise use what's on PATH
# (activate the project conda env, and ensure an R 4.5.x Rscript is on PATH, first).
PY="${PY:-python}"
RS="${RS:-Rscript}"
TIMEOUT_SECS="${PANEL_TIMEOUT:-3600}"

# Data location: override DATA_ROOT / ORGANOIDS_ROOT to point at a mounted drive or a
# Zenodo download (mirrors the ${VAR:-default} contract in config/config.yml).
export DATA_ROOT="${DATA_ROOT:-/Volumes/Nosepass/Epitype_Project}"
export OUTPUT_ROOT="${OUTPUT_ROOT:-$REPO/output}"
export ORGANOIDS_ROOT="${ORGANOIDS_ROOT:-$DATA_ROOT/2026_organoids}"
export PYTHONPATH="$REPO"
export MPLBACKEND="Agg"

rel="${1:-}"; shift || true
if [ -z "$rel" ]; then echo "===PANEL_RESULT==="; echo "exit=2"; echo "ERROR: no script given"; exit 2; fi
script="$REPO/$rel"
if [ ! -f "$script" ]; then echo "===PANEL_RESULT==="; echo "script=$rel"; echo "exit=127"; echo "secs=0"; echo "ERROR: script not found: $script"; exit 127; fi

sdir=$(dirname "$script")
sbase=$(basename "$script")
ext="${sbase##*.}"

logdir="$REPO/_run_logs/panels"
mkdir -p "$logdir"
safe=$(echo "$rel" | tr '/ ' '__')
log="$logdir/${safe}.log"

case "$ext" in
  py|PY|Py) interp="$PY"; runcwd="$REPO"; target="$rel" ;;
  R|r)
    if grep -qF 'source(file.path("..", ".."' "$script"; then
      interp="$RS"; runcwd="$sdir"; target="$sbase"      # relative-source idiom -> CWD = script dir
    else
      interp="$RS"; runcwd="$REPO"; target="$rel"         # hardened/repo-root idiom -> CWD = repo
    fi ;;
  *) echo "===PANEL_RESULT==="; echo "script=$rel"; echo "exit=126"; echo "secs=0"; echo "ERROR: unknown extension .$ext"; exit 126 ;;
esac

{
  echo "### run_panel.sh"
  echo "### script : $rel"
  echo "### interp : $interp"
  echo "### cwd    : $runcwd"
  echo "### args   : $*"
  echo "### started: $(date -u +%FT%TZ)"
  echo "============================================================"
} > "$log"

marker=$(mktemp)
start=$(date +%s)
( cd "$runcwd" && perl -e 'alarm shift @ARGV; exec @ARGV or die "exec failed: $!\n"' "$TIMEOUT_SECS" "$interp" "$target" "$@" ) >> "$log" 2>&1
code=$?
end=$(date +%s)
secs=$(( end - start ))

newfiles=$(find "$REPO/output" -type f -newer "$marker" ! -name '._*' 2>/dev/null | sort)
imgs=$(printf '%s\n' "$newfiles" | grep -iE '\.(png|pdf|svg)$' || true)
rm -f "$marker"

{
  echo "============================================================"
  echo "### exit: $code   secs: $secs"
  echo "### new_output_files:"
  printf '%s\n' "$newfiles" | sed '/^$/d;s/^/###   /'
} >> "$log"

echo "===PANEL_RESULT==="
echo "script=$rel"
echo "exit=$code"
echo "secs=$secs"
[ "$code" -eq 142 ] && echo "note=TIMEOUT_after_${TIMEOUT_SECS}s"
echo "log=$log"
echo "NEW_IMAGE_FILES:"
printf '%s\n' "$imgs" | sed '/^$/d'
echo "ALL_NEW_FILES:"
printf '%s\n' "$newfiles" | sed '/^$/d'
echo "STDERR_TAIL:"
tail -n 25 "$log"
exit "$code"
