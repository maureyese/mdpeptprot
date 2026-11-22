#!/usr/bin/env bash
# =============================================================================
# 05_equilibrate.sh — NVT (heat to TEMPERATURE_K), then NPT (settle density).
# Heavy atoms are position-restrained (-DPOSRES) while the solvent relaxes.
# After NPT, a gate checks mean temperature, density and density drift, and
# stops the job if the system is not equilibrated (EQUIL_GATE=strict, default).
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"

STAGE="05_equilibrate"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_file "$RUN_DIR/minim.gro"
cd "$RUN_DIR"
rm -f nvt.* npt.* '#'*nvt* '#'*npt* logs/05e_npt_props.xvg

EQUIL_GATE="${EQUIL_GATE:-strict}"                          # strict | warn
EQUIL_TEMP_TOL="${EQUIL_TEMP_TOL:-3}"                       # K
EQUIL_DENSITY_DRIFT_PCT="${EQUIL_DENSITY_DRIFT_PCT:-0.3}"   # %
XTC_STEPS=$(( OUT_INTERVAL_STEPS * 10 ))                    # coarse trajectory during equilibration

MDRUN_FLAGS="-ntmpi 1 -ntomp $CORES_SOLO -nb gpu -pme gpu -bonded gpu -update gpu -pin on -gpu_id $GPU_ID"

# ---- NVT ---------------------------------------------------------------------
cat > nvt.mdp << MDP
define               = -DPOSRES
integrator           = md
dt                   = ${DT}
nsteps               = ${NVT_NSTEPS}
nstxout-compressed   = ${XTC_STEPS}
nstenergy            = ${OUT_INTERVAL_STEPS}
nstlog               = ${OUT_INTERVAL_STEPS}
continuation         = no
gen-vel              = yes
gen-temp             = ${TEMPERATURE_K}
gen-seed             = -1
constraints          = h-bonds
constraint-algorithm = lincs
lincs-order          = 4
lincs-iter           = 1
cutoff-scheme        = Verlet
nstlist              = 20
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
pcoupl               = no
nstcalcenergy        = 100
MDP

log "Running NVT (${NVT_TIME_PS} ps, dt=${DT}, ${TEMPERATURE_K} K)..."
"$GMX_BIN" grompp -f nvt.mdp -c minim.gro -r minim.gro -p topol.top -o nvt.tpr -maxwarn "${GROMPP_MAXWARN:-0}" \
  2>&1 | tee logs/05a_grompp_nvt.log
"$GMX_BIN" mdrun -deffnm nvt $MDRUN_FLAGS 2>&1 | tee logs/05b_mdrun_nvt.log
require_file nvt.gro
log "NVT done: $(gro_atom_count nvt.gro) atoms"

# ---- NPT ---------------------------------------------------------------------
cat > npt.mdp << MDP
define               = -DPOSRES
integrator           = md
dt                   = ${DT}
nsteps               = ${NPT_NSTEPS}
nstxout-compressed   = ${XTC_STEPS}
nstenergy            = ${OUT_INTERVAL_STEPS}
nstlog               = ${OUT_INTERVAL_STEPS}
continuation         = yes
gen-vel              = no
constraints          = h-bonds
constraint-algorithm = lincs
lincs-order          = 4
lincs-iter           = 1
cutoff-scheme        = Verlet
nstlist              = 20
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
refcoord-scaling     = com
nstcalcenergy        = 100
MDP

log "Running NPT (${NPT_TIME_PS} ps, dt=${DT}, ${TEMPERATURE_K} K)..."
"$GMX_BIN" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn "${GROMPP_MAXWARN:-0}" \
  2>&1 | tee logs/05c_grompp_npt.log
"$GMX_BIN" mdrun -deffnm npt $MDRUN_FLAGS 2>&1 | tee logs/05d_mdrun_npt.log
require_file npt.gro
log "NPT done: $(gro_atom_count npt.gro) atoms"

# ---- Equilibration gate -----------------------------------------------------
echo -e "Temperature\nDensity\n0\n" | "$GMX_BIN" energy -f npt.edr -o logs/05e_npt_props.xvg \
  > logs/05e_energy.log 2>&1 || die "gmx energy failed (see logs/05e_energy.log)"
require_file logs/05e_npt_props.xvg

# Means over the last 40% of NPT (T, density) and the previous 40% (density).
read -r n_pts T_LAST D_LAST D_PREV < <(awk '
  !/^[#@]/ {n++; t[n]=$2; d[n]=$3}
  END {
    if (n < 10) { print n, 0, 0, 0; exit }
    a=int(0.2*n)+1; b=int(0.6*n)
    for (i=b+1; i<=n; i++) { st+=t[i]; sl+=d[i]; nl++ }
    for (i=a;   i<=b; i++) { sp+=d[i]; np_++ }
    printf "%d %.3f %.3f %.3f\n", n, st/nl, sl/nl, sp/np_
  }' logs/05e_npt_props.xvg)
[ "$n_pts" -ge 10 ] || die "Too few NPT energy points (${n_pts}) to judge equilibration; lengthen NPT_TIME_PS."

drift=$(awk -v a="$D_LAST" -v b="$D_PREV" 'BEGIN{d=(a-b)/b*100; if(d<0)d=-d; printf "%.3f", d}')
t_off=$(awk -v a="$T_LAST" -v b="$TEMPERATURE_K" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.2f", d}')

printf '{"mean_temperature_K": %s, "target_temperature_K": %s, "mean_density_kg_m3": %s, "density_drift_pct": %s, "npt_ps": %s}\n' \
  "$T_LAST" "$TEMPERATURE_K" "$D_LAST" "$drift" "$NPT_TIME_PS" > 05_equilibration_report.json

problems=""
awk -v x="$t_off" -v t="$EQUIL_TEMP_TOL" 'BEGIN{exit !(x>t)}' \
  && problems+="mean temperature ${T_LAST} K is ${t_off} K from target; "
awk -v x="$drift" -v t="$EQUIL_DENSITY_DRIFT_PCT" 'BEGIN{exit !(x>t)}' \
  && problems+="density still drifting (${drift}% between the last two blocks); "
awk -v x="$D_LAST" 'BEGIN{exit !(x<950 || x>1100)}' \
  && problems+="density ${D_LAST} kg/m3 outside the plausible 950-1100 range; "

log "NPT summary: mean T ${T_LAST} K (target ${TEMPERATURE_K}), density ${D_LAST} kg/m3, drift ${drift}%"
if [ -n "$problems" ]; then
  if [ "$EQUIL_GATE" = "strict" ]; then
    die "Equilibration gate failed: ${problems}Try a longer NPT (NPT_TIME_PS) or set EQUIL_GATE=warn to continue anyway."
  fi
  log "WARNING: ${problems}"
fi

mark_done "$STAGE"
