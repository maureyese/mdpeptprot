#!/usr/bin/env bash
# =============================================================================
# 01_setup_forcefield.sh
# Makes the CHARMM36 force field visible to GROMACS (via GMXLIB). The archive
# is extracted ONCE into a shared folder (FF_CACHE_DIR), not once per run.
# Also checks that the chosen water model and ions exist in this force field.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

STAGE="01_setup_forcefield"
stage_start
log "=== Stage: $STAGE ==="

FF_CACHE_DIR="${FF_CACHE_DIR:-${PROJECT_ROOT}/forcefields/_extracted}"

if [ -d "$CHARMM_FF_SOURCE" ]; then
  # Already-extracted *.ff folder: use its parent as the force field library.
  [ "$(basename "$CHARMM_FF_SOURCE")" = "${FF_NAME}.ff" ] \
    || die "Folder name must be ${FF_NAME}.ff (got $(basename "$CHARMM_FF_SOURCE")); fix FF_NAME or the folder name."
  FF_LIB_DIR="$(cd "$(dirname "$CHARMM_FF_SOURCE")" && pwd)"
elif [ -f "$CHARMM_FF_SOURCE" ]; then
  FF_LIB_DIR="$FF_CACHE_DIR"
  if [ ! -d "$FF_LIB_DIR/${FF_NAME}.ff" ]; then
    log "Extracting force field once into $FF_LIB_DIR ..."
    mkdir -p "$FF_LIB_DIR"
    tmp=$(mktemp -d "$FF_LIB_DIR/.extract.XXXXXX")
    case "$CHARMM_FF_SOURCE" in
      *.tgz|*.tar.gz) tar -xzf "$CHARMM_FF_SOURCE" -C "$tmp" ;;
      *.zip)          unzip -oq "$CHARMM_FF_SOURCE" -d "$tmp" ;;
      *) rm -rf "$tmp"; die "Unrecognized force field archive format: $CHARMM_FF_SOURCE" ;;
    esac
    found=$(find "$tmp" -maxdepth 2 -type d -name "${FF_NAME}.ff" -print -quit)
    [ -n "$found" ] || { rm -rf "$tmp"; die "Archive does not contain a ${FF_NAME}.ff folder: $CHARMM_FF_SOURCE"; }
    mv "$found" "$FF_LIB_DIR/" 2>/dev/null || true   # another job may have extracted it first
    rm -rf "$tmp"
  else
    log "Force field already extracted at $FF_LIB_DIR/${FF_NAME}.ff (reusing)"
  fi
else
  die "Force field source not found: $CHARMM_FF_SOURCE (set CHARMM_FF_SOURCE to the .ff folder or its .tgz)"
fi

FF_PATH="$FF_LIB_DIR/${FF_NAME}.ff"
require_file "$FF_PATH/forcefield.itp"
log "Force field ready at: $FF_PATH"

# Verify the water model and the ions requested in config exist in this force field.
[ -f "$FF_PATH/${WATER_MODEL}.itp" ] || die "Water model '${WATER_MODEL}' not found in $FF_PATH"
ion_names=$(awk '/\[ *moleculetype *\]/{f=1;next} f && !/^ *;/ && NF {print $1; f=0}' "$FF_PATH/ions.itp")
for ion in "$PNAME" "$NNAME"; do
  grep -qx "$ion" <<< "$ion_names" || die "Ion '$ion' is not defined in $FF_PATH/ions.itp"
done
log "Water model '${WATER_MODEL}' and ions '${PNAME}'/'${NNAME}' confirmed."

# GROMACS looks for *.ff folders inside GMXLIB. Later stages source this file.
export GMXLIB="$FF_LIB_DIR"
echo "export GMXLIB=\"$FF_LIB_DIR\"" > "$RUN_DIR/.gmxlib_env"
log "GMXLIB set to $FF_LIB_DIR (saved to $RUN_DIR/.gmxlib_env for later stages)"

mark_done "$STAGE"
