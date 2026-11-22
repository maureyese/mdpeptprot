#!/usr/bin/env bash
# =============================================================================
# 04_minimize.sh — energy minimization (steepest descent).
# Fails if the structure is still badly strained after minimization.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="04_minimize"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_file "$RUN_DIR/ionized.gro"
cd "$RUN_DIR"
rm -f minim.* '#'*minim*

MINIM_FMAX_FAIL="${MINIM_FMAX_FAIL:-3000}"   # kJ/mol/nm; above this after minimization = stop

cat > minim.mdp << MDP
integrator           = steep
nsteps               = ${MINIM_NSTEPS}
emtol                = 1000.0
emstep               = 0.01
cutoff-scheme        = Verlet
nstlist              = 20
coulombtype          = PME
rcoulomb             = 1.2
vdwtype              = cutoff
vdw-modifier         = force-switch
rvdw-switch          = 1.0
rvdw                 = 1.2
constraints          = h-bonds
constraint-algorithm = lincs
MDP

# grompp fails loudly on any topology/coordinate atom-count mismatch.
"$GMX_BIN" grompp -f minim.mdp -c ionized.gro -p topol.top -o minim.tpr -maxwarn 1 \
  2>&1 | tee logs/04a_grompp.log

"$GMX_BIN" mdrun -v -deffnm minim -ntmpi 1 -ntomp "$CORES_SOLO" -nb gpu -gpu_id "$GPU_ID" \
  2>&1 | tee logs/04b_mdrun.log

require_file "$RUN_DIR/minim.gro"

fmax=$(awk '/Maximum force/{v=$4} END{print v}' minim.log)
[[ "$fmax" =~ ^[0-9.eE+-]+$ ]] || die "Could not read the final maximum force from minim.log"

if grep -q "Steepest Descents converged" minim.log; then
  log "Minimization converged (Fmax ${fmax} kJ/mol/nm)."
else
  if awk -v f="$fmax" -v t="$MINIM_FMAX_FAIL" 'BEGIN{exit !(f>t)}'; then
    die "Minimization did not converge and Fmax is still ${fmax} (> ${MINIM_FMAX_FAIL}). The structure probably has severe clashes; see logs/04b_mdrun.log."
  fi
  log "WARNING: minimization did not reach Fmax < 1000 in ${MINIM_NSTEPS} steps (Fmax ${fmax}); continuing because it is below ${MINIM_FMAX_FAIL}."
fi

log "Minimization complete: $(gro_atom_count minim.gro) atoms"
mark_done "$STAGE"
