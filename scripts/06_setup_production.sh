#!/usr/bin/env bash
# =============================================================================
# 06_setup_production.sh
# Builds N_REPLICAS production .tpr files from npt.gro. Each replica gets a
# distinct velocity seed (gen-vel=yes, continuation=no) so replicas diverge
# from the start. Only the Protein group (protein + peptide) is written to the
# trajectory; protein.tpr is the matching topology for analysis.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="06_setup_production"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_file "$RUN_DIR/npt.gro"
cd "$RUN_DIR"

for r in $(seq 1 "$N_REPLICAS"); do
  rd="$RUN_DIR/rep${r}"
  mkdir -p "$rd"
  seed=$(( r * 111 ))

  cat > "$rd/md.mdp" << MDP
integrator           = md
dt                   = ${DT}
nsteps               = ${PROD_NSTEPS}
continuation         = no
gen-vel              = yes
gen-temp             = ${TEMPERATURE_K}
gen-seed             = ${seed}
constraints          = h-bonds
constraint-algorithm = lincs
lincs-order          = 4
lincs-iter           = 1
cutoff-scheme        = Verlet
nstlist              = 100
coulombtype          = PME
rcoulomb             = 1.2
vdwtype              = cutoff
vdw-modifier         = force-switch
rvdw-switch          = 1.0
rvdw                 = 1.2
tcoupl               = V-rescale
tc-grps              = Protein Non-Protein
tau-t                = 0.1 0.1
ref-t                = ${TEMPERATURE_K} ${TEMPERATURE_K}
pcoupl               = C-rescale
pcoupltype           = isotropic
tau-p                = 2.0
ref-p                = 1.0
compressibility      = 4.5e-5
nstxout-compressed   = ${OUT_INTERVAL_STEPS}
compressed-x-grps    = Protein
nstenergy            = ${OUT_INTERVAL_STEPS}
nstlog               = ${OUT_INTERVAL_STEPS}
nstcalcenergy        = 100
MDP

  log "grompp for rep${r} (seed=${seed}, ${PRODUCTION_NS} ns, ${TEMPERATURE_K} K)..."
  "$GMX_BIN" grompp -f "$rd/md.mdp" -c npt.gro -p topol.top -o "$rd/md.tpr" \
    -maxwarn "${GROMPP_MAXWARN:-0}" 2>&1 | tee "logs/06_grompp_rep${r}.log"
  require_file "$rd/md.tpr"

  # Topology matching the protein-only trajectory (needed by the analysis step).
  echo Protein | "$GMX_BIN" convert-tpr -s "$rd/md.tpr" -o "$rd/protein.tpr" \
    > "logs/06_convert_protein_rep${r}.log" 2>&1
  require_file "$rd/protein.tpr"
done

log "All ${N_REPLICAS} replica .tpr files created."
mark_done "$STAGE"
