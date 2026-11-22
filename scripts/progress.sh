#!/usr/bin/env bash
# progress.sh : progress and ETA of each production replica. Usage: RUN_NAME=<name> bash scripts/progress.sh
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
run="${RUN_NAME:?set RUN_NAME}"
for d in runs/"$run"/rep*/; do
  r=$(basename "$d"); pid=$(cat "$d/.mdrun.pid" 2>/dev/null); el=0
  [ -n "$pid" ] && el=$(ps -o etimes= -p "$pid" 2>/dev/null | tr -d ' ')
  [ -f "$d/.production.done" ] && { echo "$r: DONE"; continue; }
  awk -v r="$r" -v el="${el:-0}" '
    /^ +nsteps +=/ {ns=$3}
    /^ +dt +=/     {dt=$3}
    /^ +tinit +=/  {t0=$3}
    /^ +Step +Time/ {getline; t=$2}
    END {
      tot = t0 + ns*dt; pct = (tot > 0) ? 100*t/tot : 0
      if (el > 0 && t > 0) {
        rate = t/el
        printf "%s: %.0f / %.0f ps (%.1f%%)  ~%.0f ns/day  ETA %.0f min\n", r, t, tot, pct, rate*86.4, (tot-t)/rate/60
      } else printf "%s: %.0f / %.0f ps (%.1f%%)  (not running)\n", r, t, tot, pct
    }' "$d/md.log" 2>/dev/null
done
