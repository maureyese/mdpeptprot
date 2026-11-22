#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
# Chains stages 00 through 07 (input prep through benchmark). Stops on the
# first failure. Does NOT launch production (08) automatically — that is a
# separate, deliberate step since it commits GPU time for days.
#
# Usage:
#   export RUN_NAME=my_complex
#   export INPUT_PDB=~/velozbio/mdpeptprot/inputs/my_complex.pdb
#   export CHARMM_FF_SOURCE=~/velozbio/mdpeptprot/forcefields/charmm36-jul2022.ff.tgz
#   ./run_pipeline.sh
#
# Resume behavior: stages already marked done (see common.sh mark_done) are
# skipped unless FORCE_RERUN=1 is set, so this script is safe to re-run after
# a failure partway through.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

STAGES=(
  00_prepare_input
  01_setup_forcefield
  02_pdb2gmx
  03_solvate_ionize
  04_minimize
  05_equilibrate
  06_setup_production
  07_benchmark
)

for stage in "${STAGES[@]}"; do
  if [ "${FORCE_RERUN:-0}" != "1" ] && stage_is_done "$stage"; then
    log "Skipping $stage (already done — set FORCE_RERUN=1 to redo everything)"
    continue
  fi
  log ">>> Running $stage"
  "$SCRIPT_DIR/${stage}.sh"
done

log ""
log "Pipeline setup complete. Review:"
log "  - $RUN_DIR/logs/05e_npt_density.xvg  (density should be flat/equilibrated)"
log "  - benchmark ns/day printed above"
log "When ready, launch production with:"
log "  $SCRIPT_DIR/08_launch_production.sh"
