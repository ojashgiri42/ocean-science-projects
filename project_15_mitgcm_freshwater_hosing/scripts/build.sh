#!/bin/bash
# Build the MITgcm executable (mitgcmuv) for the tutorial_global_oce_latlon setup.
#
# Usage:  bash build.sh <code_dir> <build_dir> [ieee|fast] [changes_dir]
#   code_dir     folder with the compile-time files (SIZE.h, packages.conf, ...)
#   build_dir    where to build; must be outside the git repo
#   ieee|fast    "ieee" = strict floating-point flags (best for comparing with the
#                reference output), "fast" = optimized flags (faster for long runs)
#   changes_dir  optional folder with only the compile-time files we changed;
#                they are copied over the files from code_dir

set -e                                   # stop at the first error

MITGCM_ROOT=$HOME/models/MITgcm          # MITgcm source code (cloned outside the repo)
OPTFILE=$MITGCM_ROOT/tools/build_options/darwin_arm64_gfortran   # compiler settings for an Apple-silicon Mac with gfortran

# The tutorial compiles pkg/mnc (netCDF output), so the link step needs -lnetcdf.
# Homebrew's nf-config only gives the netCDF-Fortran library folder, not the
# folder of the netCDF C library, so the first build failed with
# "ld: library 'netcdf' not found". gfortran also searches LIBRARY_PATH,
# so we add Homebrew's lib folder there. No MITgcm file is changed.
export LIBRARY_PATH=/opt/homebrew/lib:$LIBRARY_PATH

CODE_DIR=$1                              # compile-time code folder
BUILD_DIR=$2                             # build folder
MODE=${3:-fast}                          # default to the optimized build
CHANGES=$4                               # optional folder of changed code files

if [ -z "$CODE_DIR" ] || [ -z "$BUILD_DIR" ]; then
    echo "usage: bash build.sh <code_dir> <build_dir> [ieee|fast] [changes_dir]"
    exit 1
fi

CODE_DIR=$(cd "$CODE_DIR" && pwd)        # make the path absolute, because we cd below
mkdir -p "$BUILD_DIR"                    # create the build folder if needed

if [ -n "$CHANGES" ]; then
    # Combine the code folders ourselves, so it is clear which file wins:
    # start from the original code, then copy our changed files on top.
    rm -rf "$BUILD_DIR/code_combined"
    mkdir "$BUILD_DIR/code_combined"
    cp "$CODE_DIR"/* "$BUILD_DIR/code_combined"/
    cp "$CHANGES"/* "$BUILD_DIR/code_combined"/
    CODE_DIR="$BUILD_DIR/code_combined"
fi

cd "$BUILD_DIR"                          # genmake2 writes its Makefile in the current folder

if [ "$MODE" = "ieee" ]; then
    IEEE_FLAG="-ieee"                    # strict IEEE arithmetic, no aggressive optimization
else
    IEEE_FLAG=""                         # use the optimized flags from the build-options file
fi

echo "== genmake2 (code: $CODE_DIR, mode: $MODE)"
$MITGCM_ROOT/tools/genmake2 -rootdir=$MITGCM_ROOT -mods=$CODE_DIR -optfile=$OPTFILE $IEEE_FLAG > genmake2.log 2>&1   # write the Makefile
echo "== make depend"
make depend > make_depend.log 2>&1       # work out which source files depend on which
echo "== make"
make > make.log 2>&1                     # compile and link mitgcmuv

ls -l mitgcmuv                           # show that the executable exists
