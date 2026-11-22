#!/usr/bin/env bash
# =============================================================================
# config.sh — single source of truth for the MD pipeline.
# Every stage script sources this file. Values here are DEFAULTS: any of them
# can be overridden from the environment (run_md.py sets them from a job YAML).
# =============================================================================

# ---- Project layout ---------------------------------------------------------
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="${PROJECT_ROOT:-$(dirname "$_here")}"
export RUNS_ROOT="${RUNS_ROOT:-${PROJECT_ROOT}/runs}"
export RUN_NAME="${RUN_NAME:?Set RUN_NAME, e.g. export RUN_NAME=1WCY_BD}"
export RUN_DIR="${RUNS_ROOT}/${RUN_NAME}"
export INPUT_PDB="${INPUT_PDB:-${PROJECT_ROOT}/inputs/${RUN_NAME}.pdb}"

# Force field: extracted charmm36*.ff folder, OR the .tgz (extracted into
# $RUN_DIR/ff automatically).
export CHARMM_FF_SOURCE="${CHARMM_FF_SOURCE:-${PROJECT_ROOT}/forcefields/charmm36-jul2022.ff.tgz}"
export FF_NAME="${FF_NAME:-charmm36-jul2022}"   # must match the folder name inside the tgz

# GROMACS binary. Leave as "gmx" if it is on PATH (conda env active).
export GMX_BIN="${GMX_BIN:-gmx}"

# ---- Physiological conditions ----------------------------------------------
# Defaults describe blood plasma (extracellular targets such as DPP-IV).
export PH="${PH:-7.4}"
export TEMPERATURE_K="${TEMPERATURE_K:-310}"
export ION_CONC="${ION_CONC:-0.15}"     # mol/L
export PNAME="${PNAME:-NA}"             # cation: NA (extracellular) or K (intracellular)
export NNAME="${NNAME:-CL}"

# ---- Hardware ----------------------------------------------------------------
export GPU_ID="${GPU_ID:-0}"
export TOTAL_CORES="${TOTAL_CORES:-$(nproc)}"
export N_REPLICAS="${N_REPLICAS:-3}"
_cpr=$(( TOTAL_CORES / N_REPLICAS ))
[ "$_cpr" -le "${MAX_THREADS_PER_REPLICA:-4}" ] || _cpr="${MAX_THREADS_PER_REPLICA:-4}"   # more than 4 threads per replica gave no gain (GPU-bound)
[ "$_cpr" -ge 1 ] || _cpr=1
# Cores per replica for production (replicas share the GPU). Tune with 07_benchmark.sh.
export CORES_PER_REPLICA="${CORES_PER_REPLICA:-$_cpr}"
# Cores for single-job steps (minimization / NVT / NPT).
export CORES_SOLO="${CORES_SOLO:-$TOTAL_CORES}"

# ---- Box / solvation ---------------------------------------------------------
export BOX_TYPE="${BOX_TYPE:-dodecahedron}"
export BOX_DISTANCE="${BOX_DISTANCE:-1.0}"   # nm, min distance solute-to-box-edge
export WATER_MODEL="${WATER_MODEL:-tip3p}"

# ---- Timestep / HMR ----------------------------------------------------------
# 4 fs needs hydrogen mass repartitioning AND verified stability for your
# system (test a short NVT run first). Default is the SAFE choice: 2 fs.
export USE_HMR="${USE_HMR:-0}"          # 1 = pdb2gmx -heavyh + dt=0.004
if [ "$USE_HMR" = "1" ]; then
  export DT="0.004"
  export HEAVYH_FLAG="-heavyh"
else
  export DT="0.002"
  export HEAVYH_FLAG=""
fi

# ---- Equilibration lengths ---------------------------------------------------
export MINIM_NSTEPS="${MINIM_NSTEPS:-50000}"
export NVT_TIME_PS="${NVT_TIME_PS:-100}"
export NPT_TIME_PS="${NPT_TIME_PS:-500}"

# ---- Production (screening defaults; extend only the best peptides later) ----
export PRODUCTION_NS="${PRODUCTION_NS:-30}"                   # ns per replica (placeholder, set from benchmark)
export PROD_OUTPUT_INTERVAL_PS="${PROD_OUTPUT_INTERVAL_PS:-10}"   # frame / energy / log interval

# ---- Derived step counts (do not edit) --------------------------------------
export NVT_NSTEPS=$(awk -v ps="$NVT_TIME_PS" -v dt="$DT" 'BEGIN{printf "%d", ps/dt}')
export NPT_NSTEPS=$(awk -v ps="$NPT_TIME_PS" -v dt="$DT" 'BEGIN{printf "%d", ps/dt}')
export PROD_NSTEPS=$(awk -v ns="$PRODUCTION_NS" -v dt="$DT" 'BEGIN{printf "%d", (ns*1000)/dt}')
export OUT_INTERVAL_STEPS=$(awk -v ps="$PROD_OUTPUT_INTERVAL_PS" -v dt="$DT" 'BEGIN{printf "%d", ps/dt}')
