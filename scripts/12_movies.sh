#!/usr/bin/env bash
# =============================================================================
# 12_movies.sh : two videos of one replica (run AFTER 11_analyze.sh)
#   runs/<run>/movies/rep<N>_structure.mp4         3D binding only (title, H-bond labels, contacts table)
#   runs/<run>/movies/rep<N>_structure_plots.mp4   same 3D (left) + RMSD / distance / contacts / H-bonds (right)
# Options (environment): REP=1  MOVIE_FRAMES=300  MOVIE_FPS=10  MOVIE_ENV=mdmovie
#   MOVIE_CPUS=16-31 (restrict to these cores, e.g. while a simulation runs)
#   KEEP_FRAMES=1 (keep the PNG frames)   REUSE_FRAMES=1 (skip PyMOL, redo overlays/plots/videos from kept frames)
# PyMOL and ffmpeg live in a separate conda env (default: mdmovie), not in the simulation image.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

REP="${REP:-1}"; MOVIE_FRAMES="${MOVIE_FRAMES:-300}"; MOVIE_FPS="${MOVIE_FPS:-10}"; MOVIE_ENV="${MOVIE_ENV:-mdmovie}"
REUSE_FRAMES="${REUSE_FRAMES:-0}"
log "=== Stage: 12_movies (rep${REP}) ==="
rd="$RUN_DIR/rep${REP}"
an="$RUN_DIR/analysis"
require_file "$rd/md_clean.xtc"
require_file "$rd/complex_start.pdb"
require_file "$an/rep${REP}_timeseries.tsv"
require_file "$an/rep${REP}_hbond_pairs.json"
require_file "$an/rep${REP}_contact_pairs.json"
require_file "$RUN_DIR/00_input_report.json"

if [ -n "${MOVIE_BIN:-}" ]; then
  BIN="$MOVIE_BIN"
else
  if [ -n "${CONDA_EXE:-}" ]; then base="$(dirname "$(dirname "$CONDA_EXE")")"; else base="$(conda info --base 2>/dev/null)"; fi
  BIN="$base/envs/${MOVIE_ENV}/bin"
fi
[ -x "$BIN/pymol" ]  || die "PyMOL not found in conda env '$MOVIE_ENV' ($BIN/pymol)."
[ -x "$BIN/ffmpeg" ] || die "ffmpeg not found in conda env '$MOVIE_ENV' ($BIN/ffmpeg)."

TS=()
[ -n "${MOVIE_CPUS:-}" ] && TS=(taskset -c "$MOVIE_CPUS")

out="$RUN_DIR/movies/rep${REP}"
if [ "$REUSE_FRAMES" = "1" ]; then
  require_file "$out/states.txt"
  require_file "$out/labels.json"
  log "Re-using the rendered 3D frames in $out"
else
  rm -rf "$out"; mkdir -p "$out"
  log "Rendering 3D frames with PyMOL (ray tracing, ~2-5 s per frame)..."
  export MV_PDB="$rd/complex_start.pdb" MV_XTC="$rd/md_clean.xtc" MV_OUT="$out" MV_TARGET="$MOVIE_FRAMES"
  export MV_HB="$an/rep${REP}_hbond_pairs.json"
  ${TS[@]+"${TS[@]}"} "$BIN/pymol" -cq "$SCRIPT_DIR/pymol_render.py" > "$out/pymol.log" 2>&1 \
    || { tail -n 15 "$out/pymol.log"; die "PyMOL rendering failed (log: $out/pymol.log)"; }
  tail -n 1 "$out/pymol.log"
fi

log "Drawing plot panels, labels and contacts table..."
${TS[@]+"${TS[@]}"} python3 "$SCRIPT_DIR/plot_frames.py" "$RUN_DIR" "$REP" "$out" || die "plot_frames failed"

log "Encoding videos..."
v1="$RUN_DIR/movies/rep${REP}_structure.mp4"
v2="$RUN_DIR/movies/rep${REP}_structure_plots.mp4"
${TS[@]+"${TS[@]}"} "$BIN/ffmpeg" -y -loglevel error -framerate "$MOVIE_FPS" -i "$out/lbl_%04d.png" \
  -c:v libx264 -pix_fmt yuv420p -crf 18 "$v1" || die "ffmpeg failed (structure video)"
${TS[@]+"${TS[@]}"} "$BIN/ffmpeg" -y -loglevel error -framerate "$MOVIE_FPS" -i "$out/lbl_%04d.png" \
  -framerate "$MOVIE_FPS" -i "$out/p_%04d.png" -filter_complex "[0:v][1:v]hstack=inputs=2" \
  -c:v libx264 -pix_fmt yuv420p -crf 18 "$v2" || die "ffmpeg failed (structure+plots video)"

[ "${KEEP_FRAMES:-0}" = "1" ] || [ "$REUSE_FRAMES" = "1" ] || rm -rf "$out"
log "Videos written:"; ls -lh "$v1" "$v2"
