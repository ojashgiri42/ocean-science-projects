#!/bin/bash
# Run the final step-3 cases (open-ocean restoring zone) at the same time, one run per core.
# Each case uses: the baseline input files, the common rbcs data.pkg, the case's own data,
# data.shelfice and data.rbcs, the open-ocean geometry files, and the 3D restoring fields
# written by notebook 03 (part A).
#
# Usage:  bash run_open_cases.sh <field_dir> <case> [<case> ...]
#   field_dir  folder with rbcs_mask_open.bin, salt_34p4.bin and theta_rbcs_<x>.bin
#   case       name after open_ in experiments/, e.g. base, p20, g05 (at most 4, one per core)

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
    CHANGES="$E/baseline/input $E/rbcs_common $E/open_$CASE/input $E/open_geometry"
    TFILE=$(grep "relaxTFile" "$E/open_$CASE/input/data.rbcs" | cut -d"'" -f2)   # target-temperature file named in data.rbcs
    CHANGES="$CHANGES $FIELDS/rbcs_mask_open.bin $FIELDS/salt_34p4.bin $FIELDS/$TFILE"
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/open_$CASE" "$CHANGES" > "$RUNS/open_$CASE.log" 2>&1 &
    JOBS="$JOBS open_$CASE:$!"                           # remember each background job
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
