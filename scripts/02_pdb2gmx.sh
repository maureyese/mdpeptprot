#!/usr/bin/env bash
# =============================================================================
# 02_pdb2gmx.sh
# 1) Decides protonation states at the target pH (PROPKA via PDB2PQR) and
#    writes them into the residue names CHARMM36 understands
#    (HSD/HSE/HSP, ASPP, GLUP, LSN). Coordinates and atom names are untouched.
# 2) Builds the GROMACS topology with pdb2gmx, non-interactively, using the
#    standard charged termini (NH3+ / COO-).
# PROTONATION=propka (default) | default (GROMACS defaults, NOT pH-aware)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

STAGE="02_pdb2gmx"
stage_start
log "=== Stage: $STAGE ==="
require_cmd "$GMX_BIN"
require_file "$RUN_DIR/input_topology.pdb"
require_file "$RUN_DIR/.gmxlib_env"
source "$RUN_DIR/.gmxlib_env"
cd "$RUN_DIR"

PROTONATION="${PROTONATION:-propka}"
FF_PATH="$GMXLIB/${FF_NAME}.ff"
require_file "$FF_PATH/aminoacids.rtp"

# Clean stale outputs so topol.top can never go out of sync with coordinates.
rm -f processed.gro topol.top topol.pristine.top posre*.itp topol_*.itp posre_*.itp \
      input_protonated.pdb 02_pdb2pqr.pqr 02_protonation.tsv '#'*

# ---- Step 1: protonation states at the target pH ---------------------------
if [ "$PROTONATION" = "propka" ]; then
  PDB2PQR_BIN=""
  for c in pdb2pqr pdb2pqr30; do
    if command -v "$c" >/dev/null 2>&1; then PDB2PQR_BIN="$c"; break; fi
  done
  [ -n "$PDB2PQR_BIN" ] || die "pdb2pqr not found. Install it with: pip install pdb2pqr   (or set PROTONATION=default, which is NOT pH-aware)"

  log "Assigning protonation states at pH ${PH} (PROPKA via ${PDB2PQR_BIN})..."
  "$PDB2PQR_BIN" --ff=AMBER --with-ph="$PH" --titration-state-method=propka \
      --keep-chain --drop-water input_topology.pdb 02_pdb2pqr.pqr \
      > logs/02a_pdb2pqr.log 2>&1 \
    || { tail -20 logs/02a_pdb2pqr.log >&2; die "pdb2pqr failed (full log: logs/02a_pdb2pqr.log)"; }
  require_file "$RUN_DIR/02_pdb2pqr.pqr"

  python3 "$SCRIPT_DIR/protonate_from_pqr.py" input_topology.pdb 02_pdb2pqr.pqr input_protonated.pdb "$FF_PATH/aminoacids.rtp" 02_protonation.tsv "$PH"
elif [ "$PROTONATION" = "default" ]; then
  log "WARNING: PROTONATION=default uses GROMACS defaults and does NOT follow the pH setting."
  cp input_topology.pdb input_protonated.pdb
else
  die "PROTONATION must be 'propka' or 'default' (got '$PROTONATION')"
fi
require_file "$RUN_DIR/input_protonated.pdb"

# ---- Step 2: topology (non-interactive; standard charged termini) ----------
"$GMX_BIN" pdb2gmx \
  -f "$RUN_DIR/input_protonated.pdb" \
  -o "$RUN_DIR/processed.gro" \
  -water "$WATER_MODEL" \
  -ff "$FF_NAME" \
  -ignh \
  $HEAVYH_FLAG \
  < /dev/null 2>&1 | tee "$RUN_DIR/logs/02_pdb2gmx.log"

require_file "$RUN_DIR/processed.gro"
require_file "$RUN_DIR/topol.top"
# Pristine copy: stage 03 edits topol.top in place and restores from this on every run.
# Make force-field includes relative (resolved through GMXLIB) so the run is portable
sed -i -E "s|\"[^\"]*/(${FF_NAME}\.ff/)|\"\1|" "$RUN_DIR"/*.top "$RUN_DIR"/*.itp 2>/dev/null || true
if grep -q '#include "/' "$RUN_DIR/topol.top"; then die "Absolute include path left in topol.top"; fi
cp -f "$RUN_DIR/topol.top" "$RUN_DIR/topol.pristine.top"

if [ "$USE_HMR" = "1" ]; then
  mass=$(awk '/\[ atoms \]/{f=1;next} /^\[/{f=0} f && $5=="HN" {print $8; exit}' topol_Protein_chain_A.itp 2>/dev/null || true)
  if [ -z "$mass" ] || awk -v m="$mass" 'BEGIN{exit !(m<2)}'; then
    die "HMR requested (USE_HMR=1) but could not confirm repartitioned HN mass (got '${mass:-none}', expected ~4.03)."
  fi
  log "HMR verified: HN mass = $mass"
fi

log "Termini chosen by pdb2gmx:"
grep -iE "terminus" logs/02_pdb2gmx.log | grep -viE "identified|select" || true
grep -E "Total charge" logs/02_pdb2gmx.log || true

log "pdb2gmx complete: $(gro_atom_count "$RUN_DIR/processed.gro") atoms"
mark_done "$STAGE"
