#!/usr/bin/env bash
# =============================================================================
# 03_solvate_ionize.sh
# Boxes the complex, solvates it, neutralizes and adds salt (ion names and
# concentration come from config.sh). Safe to re-run: it always restarts from
# the pristine topology saved by stage 02, because solvate/genion edit
# topol.top in place.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="03_solvate_ionize"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_file "$RUN_DIR/processed.gro"
require_file "$RUN_DIR/topol.pristine.top"
cd "$RUN_DIR"

# Always start from the pristine topology (solvate and genion modify topol.top).
cp -f topol.pristine.top topol.top
rm -f boxed.gro solvated.gro ionized.gro ions.tpr ions.mdp mdout.mdp '#'*

log "Boxing (${BOX_TYPE}, ${BOX_DISTANCE} nm)..."
"$GMX_BIN" editconf -f processed.gro -o boxed.gro -c -d "$BOX_DISTANCE" -bt "$BOX_TYPE" \
  2>&1 | tee logs/03a_editconf.log

log "Solvating..."
"$GMX_BIN" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top \
  2>&1 | tee logs/03b_solvate.log

log "Adding ions (neutralize + ${ION_CONC} M ${PNAME}/${NNAME})..."
cat > ions.mdp << 'MDP'
integrator = steep
nsteps     = 5000
emtol      = 1000.0
MDP
"$GMX_BIN" grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1 \
  2>&1 | tee logs/03c_grompp_ions.log
echo SOL | "$GMX_BIN" genion -s ions.tpr -o ionized.gro -p topol.top \
  -pname "$PNAME" -nname "$NNAME" -neutral -conc "$ION_CONC" \
  2>&1 | tee logs/03d_genion.log

require_file "$RUN_DIR/ionized.gro"

# ---- Verify: neutral system, ion counts, approximate salt concentration -----
n_atoms=$(gro_atom_count ionized.gro)
q_sys=$(grep "Total charge in system" logs/02_pdb2gmx.log | tail -1 | awk '{print $5}')
count_mol() { awk -v n="$1" '/\[ *molecules *\]/{f=1;next} f && $1==n {s+=$2} END{print s+0}' topol.top; }
n_pos=$(count_mol "$PNAME"); n_neg=$(count_mol "$NNAME"); n_wat=$(count_mol SOL)
net=$(awk -v q="$q_sys" -v p="$n_pos" -v m="$n_neg" 'BEGIN{printf "%.3f", q+p-m}')
vol=$(tail -1 ionized.gro | awk '{printf "%.1f", $1*$2*$3}')
conc=$(awk -v n="$n_neg" -v v="$vol" 'BEGIN{printf "%.3f", n/(0.6022*v)}')

log "System: ${n_atoms} atoms | water ${n_wat} | ${PNAME} ${n_pos} | ${NNAME} ${n_neg} | box ${vol} nm^3 | net charge ${net} e | ~${conc} M ${NNAME}"
if awk -v x="$net" 'BEGIN{exit !(x>0.001 || x<-0.001)}'; then
  die "System is not neutral (net charge ${net} e). Check logs/03d_genion.log."
fi

tail -6 topol.top | tee -a logs/03d_genion.log
mark_done "$STAGE"
