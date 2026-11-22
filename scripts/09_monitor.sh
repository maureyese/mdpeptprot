#!/usr/bin/env bash
# =============================================================================
# 09_monitor.sh : progress, speed and ETA of each production replica, plus GPU status.
# Read-only; safe to run any time:  watch -n 60 "RUN_NAME=<name> bash scripts/09_monitor.sh"
# Everything comes from each replica's md.log (progress) and its process (elapsed time),
# so it is also correct for resumed or extended runs.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

echo "=== Production status: $(date) ==="
for r in $(seq 1 "$N_REPLICAS"); do
  rd="$RUN_DIR/rep${r}"
  if [ -f "$rd/.production.done" ]; then echo "rep${r}: DONE"; continue; fi
  if [ ! -f "$rd/md.log" ]; then echo "rep${r}: not started"; continue; fi
  el=0
  if pid_file_running "$rd/.mdrun.pid"; then
    el=$(ps -o etimes= -p "$(cat "$rd/.mdrun.pid")" 2>/dev/null | tr -d ' ')
  fi
  awk -v r="rep${r}" -v el="${el:-0}" '
    /^ +nsteps +=/ {ns=$3}
    /^ +dt +=/     {dt=$3}
    /^ +tinit +=/  {t0=$3}
    /Started mdrun/ {first=""}                       # a restart: measure speed from here
    /^ +Step +Time/ {getline; last=$2; if (first=="") first=$2}
    END {
      tot = t0 + ns*dt; pct = (tot > 0) ? 100*last/tot : 0
      if (el > 20 && last > first) {
        rate = (last-first)/el
        printf "%s: %.0f / %.0f ps (%.1f%%)  ~%.0f ns/day  ETA %.0f min\n", r, last, tot, pct, rate*86.4, (tot-last)/rate/60
      } else {
        printf "%s: %.0f / %.0f ps (%.1f%%)  %s\n", r, last, tot, pct, (el > 0 ? "(measuring speed...)" : "(not running, not finished)")
      }
    }' "$rd/md.log"
done
echo ""
echo "=== GPU ==="
nvidia-smi --query-gpu=utilization.gpu,memory.used,temperature.gpu,power.draw,power.limit,clocks.current.sm,clocks.max.sm --format=csv
