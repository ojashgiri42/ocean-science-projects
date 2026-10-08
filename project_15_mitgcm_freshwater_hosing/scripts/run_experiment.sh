#!/bin/bash
# Set up a run directory and run MITgcm in it.
#
# Usage:  bash run_experiment.sh <build_dir> <base_input_dir> <run_dir> ["changes"] [pickup_dir]
#   build_dir       folder containing mitgcmuv
#   base_input_dir  the tutorial's original input files
#   run_dir         where to run; must be outside the git repo
#   changes         optional list (in quotes, separated by spaces) of folders or
#                   single files with our changes. They are copied in the order
#                   given, so a later entry wins over an earlier one.
#   pickup_dir      optional run folder to restart from. The pickup files for
#                   the run's starting step (nIter0 in its data file) are linked in.

set -e                                   # stop at the first error

BUILD_DIR=$1
BASE_INPUT=$2
RUN_DIR=$3
CHANGES=$4
PICKUP_DIR=$5

if [ -z "$BUILD_DIR" ] || [ -z "$BASE_INPUT" ] || [ -z "$RUN_DIR" ]; then
    echo "usage: bash run_experiment.sh <build_dir> <base_input_dir> <run_dir> [\"changes\"] [pickup_dir]"
    exit 1
fi

if [ -e "$RUN_DIR/output.txt" ]; then    # never overwrite a finished run by accident
    echo "error: $RUN_DIR already has output.txt; remove it first if you mean to rerun"
    exit 1
fi

mkdir -p "$RUN_DIR"
ln -sf "$(cd "$BASE_INPUT" && pwd)"/* "$RUN_DIR"/       # link the original inputs (no copies of big files)

for item in $CHANGES; do                                 # apply our changes in order
    if [ -d "$item" ]; then
        files="$item"/*                                  # a folder: copy every file in it
    else
        files="$item"                                    # a single file
    fi
    for f in $files; do
        rm -f "$RUN_DIR/$(basename "$f")"                # remove a link first, so the original is never overwritten
        cp "$f" "$RUN_DIR"/                              # our file replaces the original
    done
done

if [ -n "$PICKUP_DIR" ]; then
    NITER0=$(grep -E "^ *nIter0" "$RUN_DIR/data" | sed 's/.*=//' | tr -dc '0-9')   # starting step: the number after "=" in the data file
    SUFFIX=$(printf "%010d" "$NITER0")                               # pickup files use a 10-digit step number
    ls "$PICKUP_DIR"/pickup*."$SUFFIX".* > /dev/null                 # stop here if the pickup is missing
    ln -sf "$PICKUP_DIR"/pickup*."$SUFFIX".* "$RUN_DIR"/             # link only the pickup for this start
    echo "== starting from pickup $SUFFIX in $PICKUP_DIR"
fi

cp "$BUILD_DIR/mitgcmuv" "$RUN_DIR"/                     # copy the executable, so a later rebuild cannot change this run

cd "$RUN_DIR"
echo "== running in $RUN_DIR"
/usr/bin/time -p ./mitgcmuv > output.txt 2> time.txt     # run; model messages go to output.txt, timing to time.txt
tail -3 time.txt                                         # show the wall-clock time
grep -c "%MON time_secondsf" output.txt | sed 's/^/monitor outputs: /'   # how many monitor blocks were written
grep "Execution ended Normally" output.txt               # fails (and stops a chain) if the model did not end normally
if grep -q "NaN" output.txt; then                      # "ended normally" is not enough: the 0.5 Sv mixed run
    echo "error: NaN found in output.txt"                 # printed it although its fields had become NaN
    exit 1
fi
