#!/bin/bash
# Run step-7 ISOMIP+ Ocean0 cases at the same time, one run per core.
# Each case uses the common ISOMIP+ settings (experiments/isoplus_common) and its own
# data.shelfice (and data, if the case folder has one) on top; the 2D and 3D input
# files written by notebook 07 (part A) are the base input folder.
#
# Usage:  bash run_isoplus_cases.sh <case> [<case> ...]
#   case   experiment folder name, e.g. isoplus_g050 (at most 4, one per core)

set -e                                                   # stop at the first error

SCRIPTS=$(cd "$(dirname "$0")" && pwd)
E=$(cd "$SCRIPTS/.." && pwd)/experiments
RUNS=$HOME/mitgcm_runs/isomip
BUILD=$RUNS/build_isoplus_fast                           # build with the ISOMIP+ copy of shelfice_thermodynamics.F
BASE=$RUNS/inputs_isoplus                                # geometry, WARM profiles and restoring mask

JOBS=""
for CASE in "$@"; do
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/$CASE" "$E/isoplus_common $E/$CASE" > "$RUNS/$CASE.log" 2>&1 &
    JOBS="$JOBS $CASE:$!"                                # remember each background job
done

STATUS=0
for job in $JOBS; do
    if wait "${job#*:}"; then
        echo "== $(date '+%H:%M') ${job%%:*} finished normally"
    else
        echo "== $(date '+%H:%M') ${job%%:*} FAILED (see $RUNS/${job%%:*}.log)"
        STATUS=1
    fi
done
exit $STATUS
