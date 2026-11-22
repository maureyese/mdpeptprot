#!/usr/bin/env bash
# Shared helpers. Sourced by every numbered stage script after config.sh.
set -euo pipefail

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# Every stage needs $RUN_DIR/logs to exist before it can tee into it.
mkdir -p "${RUN_DIR:?RUN_DIR not set — source config.sh first}/logs"

require_file() {
  [ -f "$1" ] || die "Required file not found: $1 (did a previous stage fail?)"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Command not found: $1 (activate the conda env with GROMACS installed)"
}

# Stop early, with a clear message, if no usable NVIDIA GPU is visible.
require_gpu() {
  command -v nvidia-smi >/dev/null 2>&1 \
    || die "nvidia-smi not found: no NVIDIA driver visible (in Docker, start the container with --gpus all)"
  nvidia-smi -i "${GPU_ID:-0}" >/dev/null 2>&1 \
    || die "GPU ${GPU_ID:-0} is not usable (check GPU_ID and the driver)"
}

# Call at the start of a stage; mark_done then records the stage's duration
# (seconds) in $RUN_DIR/stage_timings.tsv.
stage_start() { _STAGE_T0=$(date +%s); }

# Marks a stage complete so run_pipeline.sh / run_md.py can skip finished
# stages on resume. Call at the very end of each stage script.
mark_done() {
  mkdir -p "$RUN_DIR/.stage_markers"
  touch "$RUN_DIR/.stage_markers/$1.done"
  if [ -n "${_STAGE_T0:-}" ]; then
    printf '%s\t%s\n' "$1" "$(( $(date +%s) - _STAGE_T0 ))" >> "$RUN_DIR/stage_timings.tsv"
  fi
}

stage_is_done() {
  [ -f "$RUN_DIR/.stage_markers/$1.done" ]
}

# True only if the PID saved in file $1 is alive AND is really an mdrun
# process (a bare "is this PID alive" check can be fooled after a reboot,
# when the OS reuses an old PID for an unrelated program).
pid_file_running() {
  local f="$1" pid cmd
  [ -f "$f" ] || return 1
  pid=$(cat "$f" 2>/dev/null) || return 1
  [ -n "$pid" ] || return 1
  kill -0 "$pid" 2>/dev/null || return 1
  cmd=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null) || return 1
  [[ "$cmd" == *mdrun* ]]
}

# Prints the atom count from line 2 of a .gro file.
gro_atom_count() {
  sed -n '2p' "$1" | tr -d '[:space:]'
}

# Prints "<last_logged_ps> <target_ps>" for a replica directory, read from its md.log.
replica_progress() {
  awk '
    /^ +nsteps +=/ {ns=$3}
    /^ +dt +=/     {dt=$3}
    /^ +tinit +=/  {t0=$3}
    /^ +Step +Time/ {getline; last=$2}
    END { printf "%.3f %.3f\n", last+0, t0+ns*dt }
  ' "$1/md.log" 2>/dev/null
}

# True if the replica's log shows it reached its target end time.
replica_complete() {
  local p
  p=$(replica_progress "$1") || return 1
  awk -v p="$p" 'BEGIN{split(p,a," "); exit !(a[2]>0 && a[1] >= a[2]-0.01)}'
}
