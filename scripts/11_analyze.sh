#!/usr/bin/env bash
# =============================================================================
# 11_analyze.sh
# Peptide-protein analysis of all replicas -> $RUN_DIR/analysis/
# Also writes rep*/md_clean.xtc (centered, fitted, whole molecules) for movies.
# Refuses to run on unfinished replicas unless ANALYZE_PARTIAL=1.
# Extra options for analyze_complex.py can be passed in ANALYZE_ARGS.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

STAGE="11_analyze"
stage_start
log "=== Stage: $STAGE ==="
require_cmd python3

for r in $(seq 1 "$N_REPLICAS"); do
  rd="$RUN_DIR/rep${r}"
  require_file "$rd/md.xtc"
  require_file "$rd/protein.tpr"
  if [ ! -f "$rd/.production.done" ] && [ "${ANALYZE_PARTIAL:-0}" != "1" ]; then
    die "rep${r}: production not finished (ANALYZE_PARTIAL=1 analyses what exists so far)."
  fi
done

python3 "$SCRIPT_DIR/analyze_complex.py" "$RUN_DIR" ${ANALYZE_ARGS:-} || die "Analysis failed."
log "Results in $RUN_DIR/analysis/"
mark_done "$STAGE"
