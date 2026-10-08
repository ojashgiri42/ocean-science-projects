#!/bin/bash
# Run several restoring-zone cases at the same time, one model run per core.
# Each case uses: the baseline input files, the common rbcs data.pkg, the case's own
# data and data.rbcs, and the 3D input fields written by notebook 03 (part A).
#
# Usage:  bash run_rbcs_cases.sh <field_dir> <case> [<case> ...]
#   field_dir  folder with rbcs_mask.bin, salt_34p4.bin and theta_<case>.bin
#   case       experiment folder name, e.g. rbcs_base (at most 4, one per performance core)

set -e                                                   # stop at the first error

SCRIPTS=$(cd "$(dirname "$0")" && pwd)                   # this scripts/ folder
PROJECT=$(cd "$SCRIPTS/.." && pwd)                       # the project folder
E=$PROJECT/experiments
FIELDS=$1                                                # folder with the 3D input fields
shift                                                    # the remaining arguments are the case names

BUILD=$HOME/mitgcm_runs/isomip/build_rbcs_fast           # build with pkg/rbcs and the NaN-safe bounds check
BASE=$HOME/mitgcm_runs/isomip/source_copy/input          # our copy of the original isomip inputs
RUNS=$HOME/mitgcm_runs/isomip                            # all runs live here, outside the repo

JOBS=""
for CASE in "$@"; do
    CHANGES="$E/baseline/input $E/rbcs_common $E/$CASE/input $FIELDS/rbcs_mask.bin $FIELDS/salt_34p4.bin $FIELDS/theta_$CASE.bin"
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/$CASE" "$CHANGES" > "$RUNS/$CASE.log" 2>&1 &
    JOBS="$JOBS $CASE:$!"                                # remember each background job
done

STATUS=0
for job in $JOBS; do
    if wait "${job#*:}"; then                            # wait for this job and check that it succeeded
        echo "== $(date '+%H:%M') ${job%%:*} finished normally"
    else
        echo "== $(date '+%H:%M') ${job%%:*} FAILED (see $RUNS/${job%%:*}.log)"
        STATUS=1
    fi
done
exit $STATUS
