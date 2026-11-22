#!/usr/bin/env bash
# =============================================================================
# 10_stop_production.sh : gracefully stop all running replicas.
# Sends SIGINT (the same as one Ctrl-C) to each mdrun PID saved in rep*/.mdrun.pid. GROMACS then
# finishes the current step, writes a clean checkpoint and exits. NEVER use kill -9 on mdrun: it
# can corrupt the checkpoint. To continue later, run 08_launch_production.sh (or run_md.py) again:
# stopped replicas resume from md.cpt, finished ones are skipped.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

log "Stopping running replicas of $RUN_NAME ..."
any=0
for r in $(seq 1 "$N_REPLICAS"); do
  pf="$RUN_DIR/rep${r}/.mdrun.pid"
  if pid_file_running "$pf"; then
    log "rep${r}: running (PID $(cat "$pf")), sending SIGINT"
    kill -INT "$(cat "$pf")"
    any=1
  else
    log "rep${r}: not running"
  fi
done

if [ "$any" -eq 1 ]; then
  log "Waiting up to 120 s for a clean shutdown..."
  for _ in $(seq 1 24); do
    sleep 5
    left=0
    for r in $(seq 1 "$N_REPLICAS"); do pid_file_running "$RUN_DIR/rep${r}/.mdrun.pid" && left=1; done
    [ "$left" -eq 0 ] && break
  done
fi

bad=0
for r in $(seq 1 "$N_REPLICAS"); do
  if pid_file_running "$RUN_DIR/rep${r}/.mdrun.pid"; then
    log "rep${r}: STILL RUNNING after 120 s (SIGINT sent; look at rep${r}/md.log)"; bad=1
  else
    log "rep${r}: stopped. $(bash "$SCRIPT_DIR/09_monitor.sh" 2>/dev/null | grep "^rep${r}:" | head -n 1)"
  fi
done
[ "$bad" -eq 0 ] || die "Not all replicas stopped; do not delete anything or start another run yet."
log "All stopped. Resume with: RUN_NAME=$RUN_NAME bash scripts/08_launch_production.sh"
