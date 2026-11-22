#!/usr/bin/env bash
# =============================================================================
# run_replica.sh <replica_dir> <pin_offset>
# Runs ONE production replica (fresh start, or resume from md.cpt), and writes
# <replica_dir>/.production.done only if the run really reached its end time.
# Called by 08_launch_production.sh.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
[ -f "$RUN_DIR/.gmxlib_env" ] && source "$RUN_DIR/.gmxlib_env"

rd="$1"; offset="$2"
cd "$rd"
resume=()
[ -f md.cpt ] && resume=(-cpi md.cpt -append)

"$GMX_BIN" mdrun -deffnm md -ntmpi 1 -ntomp "$CORES_PER_REPLICA" \
  -nb gpu -pme gpu -bonded gpu -update gpu \
  -pin on -pinoffset "$offset" -pinstride 1 -nstlist 100 \
  -gpu_id "$GPU_ID" -cpt 15 ${resume[@]+"${resume[@]}"} \
  > md_run.log 2>&1 &
pid=$!
echo "$pid" > .mdrun.pid          # 10_stop_production.sh sends SIGINT to this PID
rc=0
wait "$pid" || rc=$?
rm -f .mdrun.pid
if [ "$rc" -eq 0 ] && replica_complete "$rd"; then
  touch .production.done
fi
exit "$rc"
