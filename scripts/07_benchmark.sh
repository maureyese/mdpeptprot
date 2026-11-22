#!/usr/bin/env bash
# =============================================================================
# 07_benchmark.sh
# Runs all N_REPLICAS production .tpr files concurrently for a short burst, on
# scratch copies, to measure REAL combined throughput on this GPU. Can compare
# several thread counts (BENCH_NTOMP_LIST="4 6 10") and prints the longest
# production run (ns per replica) that fits in BUDGET_H hours.
#
# NOTE: a short benchmark can look faster than sustained multi-hour throughput
# if the GPU is power/thermal limited. Use BENCH_STEPS large enough for a few
# minutes per job, and re-check with 09_monitor.sh after hours of production.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="07_benchmark"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_gpu
cd "$RUN_DIR"

BENCH_STEPS="${BENCH_STEPS:-60000}"
BENCH_NTOMP_LIST="${BENCH_NTOMP_LIST:-$CORES_PER_REPLICA}"
BUDGET_H="${BUDGET_H:-20}"           # hours available for production

BENCH_DIR="$RUN_DIR/.bench"
rm -rf "$BENCH_DIR"; mkdir -p "$BENCH_DIR"
RESULT_FILE="$RUN_DIR/07_benchmark.tsv"
printf 'threads_per_replica\treplicas\tcombined_ns_per_day\tper_replica_ns_per_day\n' > "$RESULT_FILE"

best=0; best_ntomp=""
for ntomp in $BENCH_NTOMP_LIST; do
  if [ $(( ntomp * N_REPLICAS )) -gt "$TOTAL_CORES" ]; then
    log "Skipping ${ntomp} threads: ${N_REPLICAS} x ${ntomp} exceeds ${TOTAL_CORES} cores"
    continue
  fi
  pids=()
  for r in $(seq 1 "$N_REPLICAS"); do
    require_file "$RUN_DIR/rep${r}/md.tpr"
    bd="$BENCH_DIR/t${ntomp}/rep${r}"
    mkdir -p "$bd"
    cp "$RUN_DIR/rep${r}/md.tpr" "$bd/bench.tpr"
    offset=$(( (r - 1) * ntomp ))
    (
      cd "$bd"
      "$GMX_BIN" mdrun -s bench.tpr -deffnm bench -nsteps "$BENCH_STEPS" -resethway -noconfout \
        -ntmpi 1 -ntomp "$ntomp" \
        -nb gpu -pme gpu -bonded gpu -update gpu \
        -pin on -pinoffset "$offset" -pinstride 1 -nstlist 100 -gpu_id "$GPU_ID" \
        > run.out 2>&1
    ) &
    pids+=($!)
  done
  log "Benchmark: ${N_REPLICAS} concurrent jobs, ${ntomp} threads each, ${BENCH_STEPS} steps..."
  wait "${pids[@]}" || die "A benchmark job failed; see $BENCH_DIR/t${ntomp}/rep*/run.out"

  total=0
  for r in $(seq 1 "$N_REPLICAS"); do
    perf=$(awk '/^Performance:/{v=$2} END{print v}' "$BENCH_DIR/t${ntomp}/rep${r}/bench.log")
    [[ "$perf" =~ ^[0-9.]+$ ]] || die "Could not read ns/day for rep${r} (see $BENCH_DIR/t${ntomp}/rep${r}/run.out)"
    log "  rep${r}: ${perf} ns/day"
    total=$(awk -v t="$total" -v p="$perf" 'BEGIN{print t+p}')
  done
  per=$(awk -v t="$total" -v n="$N_REPLICAS" 'BEGIN{printf "%.1f", t/n}')
  log "  combined: ${total} ns/day (${per} per replica)"
  printf '%s\t%s\t%s\t%s\n' "$ntomp" "$N_REPLICAS" "$total" "$per" >> "$RESULT_FILE"
  if awk -v t="$total" -v b="$best" 'BEGIN{exit !(t>b)}'; then best=$total; best_ntomp=$ntomp; fi
done

[ -n "$best_ntomp" ] || die "No benchmark configuration ran (check BENCH_NTOMP_LIST against ${TOTAL_CORES} cores)"
log "Best: ${best} ns/day combined with ${best_ntomp} threads per replica."

# All replicas run concurrently, each at (combined / N) ns/day.
max_ns=$(awk -v t="$best" -v h="$BUDGET_H" -v n="$N_REPLICAS" 'BEGIN{printf "%.0f", t*(h/24)/n}')
eta_h=$(awk -v ns="$PRODUCTION_NS" -v t="$best" -v n="$N_REPLICAS" 'BEGIN{printf "%.1f", ns*n/t*24}')
log "Current setting (${PRODUCTION_NS} ns x ${N_REPLICAS} replicas) would take about ${eta_h} h."
log "With a ${BUDGET_H} h production budget you could run up to ~${max_ns} ns per replica."
log "Results table: $RESULT_FILE"

mark_done "$STAGE"
