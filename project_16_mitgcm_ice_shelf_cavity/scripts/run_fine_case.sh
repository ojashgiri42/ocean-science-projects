#!/bin/bash
# Run the step-6 resolution test (fine grid, half the grid spacing) for one case.
#
# Usage:  bash run_fine_case.sh <field_dir> <case>
#   field_dir  folder with the fine-grid files written by notebook 06 (part A)
#   case       experiment folder name, e.g. fine_p10

set -e                                                   # stop at the first error

SCRIPTS=$(cd "$(dirname "$0")" && pwd)                   # this scripts/ folder
E=$(cd "$SCRIPTS/.." && pwd)/experiments
F=$1                                                     # fine-grid input fields
CASE=$2
RUNS=$HOME/mitgcm_runs/isomip                            # all runs live here, outside the repo
TFILE=$(grep "relaxTFile" "$E/$CASE/input/data.rbcs" | cut -d"'" -f2)   # target-temperature file named in data.rbcs

# baseline settings, the rbcs data.pkg, the case's own files, then the fine-grid fields
CHANGES="$E/baseline/input $E/rbcs_common $E/$CASE/input"
CHANGES="$CHANGES $F/bathy_fine.bin $F/icetopo_fine.bin $F/phi0surf_fine.bin $F/rbcs_mask_fine.bin $F/salt_fine.bin $F/$TFILE"
bash "$SCRIPTS/run_experiment.sh" "$RUNS/build_fine_fast" "$RUNS/source_copy/input" "$RUNS/$CASE" "$CHANGES"
