#!/usr/bin/env bash
# =============================================================================
# 08_launch_production.sh
# Launches (or resumes) all N_REPLICAS production runs concurrently.
#   default      : replicas run detached in the background (survive closing the terminal)
#   WAIT=1       : block until all replicas finish; fail if any is incomplete
# Safe to re-run any time: finished replicas (.production.done) and running
# ones (verified to be a real mdrun process) are skipped; stopped ones resume
# from md.cpt. FORCE_RELAUNCH=1 ignores the .production.done marker.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="08_launch_production"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_gpu

WAIT="${WAIT:-0}"
pids=()
waiting=()
for r in $(seq 1 "$N_REPLICAS"); do
  rd="$RUN_DIR/rep${r}"
  require_file "$rd/md.tpr"

  if [ -f "$rd/.production.done" ] && [ "${FORCE_RELAUNCH:-0}" != "1" ]; then
    log "rep${r}: already complete, skipping."
    continue
  fi
  if pid_file_running "$rd/.mdrun.pid"; then
    log "rep${r}: already running (PID $(cat "$rd/.mdrun.pid")), skipping."
    if [ "$WAIT" = "1" ]; then waiting+=("$rd"); fi
    continue
  fi

  offset=$(( (r - 1) * CORES_PER_REPLICA ))
  if [ -f "$rd/md.cpt" ]; then
    log "rep${r}: resuming from checkpoint (cores from ${offset})..."
  else
    log "rep${r}: starting fresh (cores from ${offset})..."
  fi

  if [ "$WAIT" = "1" ]; then
    bash "$SCRIPT_DIR/run_replica.sh" "$rd" "$offset" &
    pids+=($!)
  else
    nohup setsid bash "$SCRIPT_DIR/run_replica.sh" "$rd" "$offset" > /dev/null 2>&1 &
  fi
done

if [ "$WAIT" = "1" ]; then
  rc=0
  for p in ${pids[@]+"${pids[@]}"}; do wait "$p" || rc=1; done
  for rd in ${waiting[@]+"${waiting[@]}"}; do
    log "Waiting for $(basename "$rd") (already running)..."
    while pid_file_running "$rd/.mdrun.pid"; do sleep 30; done
    for _ in $(seq 1 30); do [ -f "$rd/.production.done" ] && break; sleep 1; done
  done
  incomplete=0
  for r in $(seq 1 "$N_REPLICAS"); do
    [ -f "$RUN_DIR/rep${r}/.production.done" ] || { log "rep${r}: NOT complete"; incomplete=1; }
  done
  { [ "$rc" -eq 0 ] && [ "$incomplete" -eq 0 ]; } || die "Production did not complete for all replicas (see rep*/md_run.log)."
  log "All replicas complete."
  mark_done "$STAGE"
else
  log "Launch requests sent. Waiting 45 s before reporting status..."
  sleep 45
  for r in $(seq 1 "$N_REPLICAS"); do
    echo "--- rep${r} ---"
    tail -n 3 "$RUN_DIR/rep${r}/md_run.log" 2>/dev/null || echo "(no log yet)"
  done
  log "Use 09_monitor.sh to check progress; 10_stop_production.sh to stop."
fi
