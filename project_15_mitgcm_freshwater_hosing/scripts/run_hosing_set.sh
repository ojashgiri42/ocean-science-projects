#!/bin/bash
# Run one set of hosing experiments: a control run and three hosing strengths.
# Each hosing run has two parts: 200 years with hosing ("on"), then 300 years
# without it ("off"), restarted from the "on" run's pickup at step 792000.
# The control and the three hosing chains run at the same time: 4 model runs,
# one per performance core of the M1.
#
# Usage:  bash run_hosing_set.sh <set> <forcing_dir> <start_run_dir> ["extra"] [skip_control]
#   set            name prefix of the experiment folders, e.g. "restore" or "mixed"
#   forcing_dir    folder holding emp_hose010.bin, emp_hose030.bin, emp_hose050.bin
#   start_run_dir  run folder with the pickup the set starts from (step 720000)
#   extra          optional files copied into every run of the set (step 6 uses this)
#   skip_control   write "skip_control" to leave out the control run (step 6 reuses
#                  mixed_control from step 5)

set -e                                                   # stop at the first error

SCRIPTS=$(cd "$(dirname "$0")" && pwd)                   # this scripts/ folder
PROJECT=$(cd "$SCRIPTS/.." && pwd)                       # the project folder
SET=$1
FORCING=$(cd "$2" && pwd)
START=$3
EXTRA=$4
SKIP_CONTROL=$5

BUILD=$HOME/mitgcm_runs/hosing/build_spinup_fast         # same executable as the spin-up
BASE=$HOME/models/MITgcm/verification/tutorial_global_oce_latlon/input   # original tutorial inputs
RUNS=$HOME/mitgcm_runs/hosing                            # all runs live here, outside the repo
COMMON=$PROJECT/experiments/spinup/input                 # data.pkg and data.diagnostics from the spin-up
END_OF_HOSING=0000792000                                 # step at the end of model year 2200

run_chain () {
    CODE=$1                                              # 010, 030 or 050
    ON=${SET}_hose${CODE}_on
    OFF=${SET}_hose${CODE}_off
    # Part 1: hosing on. Later entries in the change list win, so the run's own data file
    # replaces the spin-up data file, and the hosing EmPmR file is copied in last.
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/$ON" \
        "$COMMON $EXTRA $PROJECT/experiments/$ON/input $FORCING/emp_hose$CODE.bin" "$START" > "$RUNS/$ON.log" 2>&1
    ls "$RUNS/$ON/pickup.$END_OF_HOSING.data" > /dev/null   # the "off" part needs this pickup; stop if it is missing
    # Part 2: hosing off, restarted from the pickup written at the end of part 1
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/$OFF" \
        "$COMMON $EXTRA $PROJECT/experiments/$OFF/input" "$RUNS/$ON" > "$RUNS/$OFF.log" 2>&1
}

echo "== $(date '+%H:%M') starting set '$SET'"
JOBS=""                                                  # list of "name:process-id" to wait for
if [ "$SKIP_CONTROL" != "skip_control" ]; then
    bash "$SCRIPTS/run_experiment.sh" "$BUILD" "$BASE" "$RUNS/${SET}_control" \
        "$COMMON $EXTRA $PROJECT/experiments/${SET}_control/input" "$START" > "$RUNS/${SET}_control.log" 2>&1 &
    JOBS="control:$!"                                    # remember each background job, so we can wait for it
fi
run_chain 010 &
JOBS="$JOBS hose010:$!"
run_chain 030 &
JOBS="$JOBS hose030:$!"
run_chain 050 &
JOBS="$JOBS hose050:$!"

STATUS=0
for job in $JOBS; do
    if wait "${job#*:}"; then                            # wait for this job and check that it succeeded
        echo "== $(date '+%H:%M') ${job%%:*} finished normally"
    else
        echo "== $(date '+%H:%M') ${job%%:*} FAILED (see the .log files in $RUNS)"
        STATUS=1
    fi
done
exit $STATUS
