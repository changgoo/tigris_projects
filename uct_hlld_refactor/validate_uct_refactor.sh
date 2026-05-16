#!/bin/bash
# Validate UCT-HLLD refactor:
#   1. Timing: hlld/lhlld show ratio~0%, hlle shows non-zero ratio (fallback active)
#   2. Numerical equivalence: hlld/lhlld inline UCT matches CalculateHLLWaveSpeed output

set -e

WORKTREE=$(cd "$(dirname "$0")" && pwd)
INFILE=inputs/mhd/athinput.linear_wave3d
OUTDIR=$WORKTREE/validate_output
mkdir -p $OUTDIR

source /etc/profile.d/modules.sh
module purge
module load anaconda3/2023.3 fftw/gcc/3.3.10 intel-mpi/gcc/2021.13 hdf5/gcc/intel-mpi/1.14.4

cd $WORKTREE

# Run args: 32x16x16 grid, fast wave, uct_hlld, serial
ARGS=(
    mesh/nx1=32 mesh/nx2=16 mesh/nx3=16
    meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=8
    time/nlim=100 time/ncycle_out=9999
    hydro/ct_method=uct_hlld
    problem/wave_flag=0
    problem/compute_error=true
    output1/dt=-1 output2/dt=-1
)

build_quiet() {  # usage: build_quiet <flux> [extra_cflag]
    local flux=$1; local cflag=${2:-}
    python configure.py -b --prob=linear_wave --coord=cartesian \
        --flux=$flux ${cflag:+--cflag=$cflag} 2>&1 | grep "Riemann solver:"
    make clean -s && make -j8 2>&1 | grep -c "\.cpp" | xargs -I{} echo "  compiled {} files"
}

echo ""
echo "================================================================="
echo "  PHASE 1: Timing check (PROFILE_UCT_WAVESPEED)"
echo "  Expected: hlld/lhlld ratio~0%, hlle ratio>0%"
echo "================================================================="

for flux in hlld lhlld hlle; do
    echo ""
    echo "----- flux = $flux -----"
    build_quiet $flux -DPROFILE_UCT_WAVESPEED
    rm -f linearwave-errors.dat
    ./bin/athena -i $INFILE "${ARGS[@]}" 2>&1 | grep "UCT profile"
done

echo ""
echo "================================================================="
echo "  PHASE 2: Numerical equivalence (hlld, lhlld)"
echo "  Compare new inline UCT vs legacy CalculateHLLWaveSpeed fallback"
echo "================================================================="

for flux in hlld lhlld; do
    echo ""
    echo "----- flux = $flux -----"

    # NEW path: solver_handles_uct=true (inline UCT)
    build_quiet $flux
    rm -f linearwave-errors.dat
    ./bin/athena -i $INFILE "${ARGS[@]}" 2>&1 | grep -E "^cycle|terminate"
    cp linearwave-errors.dat $OUTDIR/errors_${flux}_new.dat

    # LEGACY path: force solver_handles_uct=false via patch
    sed -i 's/solver_handles_uct = uct_on &&/solver_handles_uct = false; \/\/ TEMP force fallback -- orig: solver_handles_uct = uct_on \&\&/' \
        src/field/field.cpp
    make clean -s && make -j8 2>&1 | grep -c "\.cpp" | xargs -I{} echo "  compiled {} files (legacy patch)"
    rm -f linearwave-errors.dat
    ./bin/athena -i $INFILE "${ARGS[@]}" 2>&1 | grep -E "^cycle|terminate"
    cp linearwave-errors.dat $OUTDIR/errors_${flux}_legacy.dat

    # Restore field.cpp
    git checkout src/field/field.cpp

    # Compare
    if diff -q $OUTDIR/errors_${flux}_new.dat $OUTDIR/errors_${flux}_legacy.dat > /dev/null 2>&1; then
        echo "  PASS: bit-for-bit identical errors.dat"
    else
        echo "  DIFF found — new vs legacy:"
        diff $OUTDIR/errors_${flux}_new.dat $OUTDIR/errors_${flux}_legacy.dat || true
    fi
done

# Restore build to clean state (no profiling, hlld)
build_quiet hlld

echo ""
echo "================================================================="
echo "  Done. Output files in $OUTDIR/"
echo "================================================================="
