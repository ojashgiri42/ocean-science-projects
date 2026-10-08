#!/bin/bash
# Continue a finished step-3 run from its last pickup (restart) file.
# Same input files as run_open_cases.sh; the case's data file sets nIter0 to the pickup step.
#
# Usage:  bash run_extension.sh <field_dir> <case>
#   field_dir  folder with rbcs_mask_open.bin, salt_34p4.bin and theta_rbcs_<x>.bin
#   case       name after open_ in experiments/, e.g. base_ext (continues open_base)

set -e                                                   # stop at the first error

SCRIPTS=$(cd "$(dirname "$0")" && pwd)
E=$(cd "$SCRIPTS/.." && pwd)/experiments
FIELDS=$1
CASE=$2
RUNS=$HOME/mitgcm_runs/isomip
PARENT=$RUNS/open_${CASE%_ext}                           # the run being continued: open_base for base_ext

CHANGES="$E/baseline/input $E/rbcs_common $E/open_$CASE/input $E/open_geometry"
TFILE=$(grep "relaxTFile" "$E/open_$CASE/input/data.rbcs" | cut -d"'" -f2)
CHANGES="$CHANGES $FIELDS/rbcs_mask_open.bin $FIELDS/salt_34p4.bin $FIELDS/$TFILE"
bash "$SCRIPTS/run_experiment.sh" "$RUNS/build_rbcs_fast" "$RUNS/source_copy/input" "$RUNS/open_$CASE" "$CHANGES" "$PARENT"
