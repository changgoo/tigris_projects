# Small-box Comparison Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create four scripts in `tst/tigress_classic/` that build both TIGRIS branches, run 6 SLURM jobs (2 branches × 3 r_return values), and compare timing and cross-branch equivalence.

**Architecture:** A build script creates a git worktree for `mass-return-refactor` and compiles both branches independently. A parameterized SLURM script runs each job variant. A Python comparison script parses Athena++ stdout timing and `.hst` history files to produce timing and equivalence tables.

**Tech Stack:** bash, Python 3 + numpy, SLURM (tiger cluster), Athena++/TIGRIS, IntelMPI

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `../tigris/tst/tigress_classic/build_both.sh` | Create | Build both branches via git worktree, copy exes here |
| `../tigris/tst/tigress_classic/run_mhd_smallbox.slurm` | Create | Parameterized SLURM job (branch + r_return as args) |
| `../tigris/tst/tigress_classic/submit_tests.sh` | Create | Submit all 6 jobs |
| `../tigris/tst/tigress_classic/compare.py` | Create | Parse timing + .hst files, print tables |
| `../tigris/tst/tigress_classic/test_compare.py` | Create | pytest unit tests for compare.py parsing/comparison functions |

All paths below use `TESTDIR` to mean `../tigris/tst/tigress_classic/` (i.e., `/path/to/tigris/tst/tigress_classic/`).

---

## Task 1: Create directory and `build_both.sh`

**Files:**
- Create: `TESTDIR/build_both.sh`

- [ ] **Step 1: Create the tst/tigress_classic directory**

```bash
mkdir -p /Users/changgoo/Sources/tigris/tst/tigress_classic
```

- [ ] **Step 2: Write `build_both.sh`**

Create `TESTDIR/build_both.sh` with the following content:

```bash
#!/bin/bash
# Build tigris-master and mass-return-refactor; place exes in this directory.
# Usage: ./build_both.sh [srcbase=auto] [build_option=0]
#   build_option: 0=full clean build, 2=incremental (no make clean)
set -e

TESTDIR="$(cd "$(dirname "$0")" && pwd)"
SRCBASE=${1:-}
BUILD_OPTION=${2:-0}

# Auto-detect SRCBASE from hostname
if [ -z "$SRCBASE" ]; then
    case "$(hostname -s)" in
        tiger*|stellar*) SRCBASE="$HOME" ;;
        *)               SRCBASE="$HOME/Sources" ;;
    esac
fi

MAINDIR="$SRCBASE/tigris"
WTDIR="$MAINDIR/.worktrees/mass-return-refactor"

echo "SRCBASE=$SRCBASE  MAINDIR=$MAINDIR  TESTDIR=$TESTDIR"

# Verify MAINDIR is on tigris-master
current_branch=$(git -C "$MAINDIR" rev-parse --abbrev-ref HEAD)
if [ "$current_branch" != "tigris-master" ]; then
    echo "ERROR: $MAINDIR is on '$current_branch', expected 'tigris-master'"
    exit 1
fi

# Load modules on cluster; fall back to system g++ on mac
if command -v module &>/dev/null; then
    module purge
    module load anaconda3/2023.3 intel-oneapi/2024.2 \
        openmpi/oneapi-2024.2/4.1.6 \
        hdf5/oneapi-2024.2/openmpi-4.1.6/1.14.4 \
        fftw/oneapi-2024.2/3.3.10
    CC="icpx"
else
    CC="g++"
fi

# Create worktree if it does not already exist
if [ ! -d "$WTDIR" ]; then
    echo "Creating worktree for mass-return-refactor at $WTDIR..."
    git -C "$MAINDIR" worktree add "$WTDIR" mass-return-refactor
fi

build_one() {
    local srcdir=$1
    local exedst=$2
    cd "$srcdir"
    if [ "$BUILD_OPTION" != "2" ]; then
        ./configure.py --prob=tigress_classic --nghost=4 -fft -fb --grav=fft \
            -mpi -hdf5 -b --flux=hlld --cxx="$CC"
        make clean
    fi
    make all -j4
    cp bin/athena "$exedst"
    echo "Copied: $exedst"
}

echo "=== Building tigris-master ==="
build_one "$MAINDIR" "$TESTDIR/tigris_master_mhd.exe"

echo "=== Building mass-return-refactor ==="
build_one "$WTDIR" "$TESTDIR/tigris_mass_return_refactor_mhd.exe"

echo "=== Done ==="
ls -lh "$TESTDIR"/*.exe
```

- [ ] **Step 3: Make executable and syntax-check**

```bash
chmod +x /Users/changgoo/Sources/tigris/tst/tigress_classic/build_both.sh
bash -n /Users/changgoo/Sources/tigris/tst/tigress_classic/build_both.sh
```

Expected: no output (no syntax errors).

- [ ] **Step 4: Verify SRCBASE detection on mac**

```bash
bash -c 'source /Users/changgoo/Sources/tigris/tst/tigress_classic/build_both.sh --dry 2>&1 | head -3' || true
# Just check that hostname detection resolves without error on mac:
SRCBASE=""; case "$(hostname -s)" in tiger*|stellar*) SRCBASE="$HOME" ;; *) SRCBASE="$HOME/Sources" ;; esac; echo "SRCBASE=$SRCBASE"
```

Expected: `SRCBASE=/Users/changgoo/Sources`

- [ ] **Step 5: Commit**

```bash
git -C /Users/changgoo/Sources/tigris add tst/tigress_classic/build_both.sh
git -C /Users/changgoo/Sources/tigris commit -m "Add build_both.sh for smallbox branch comparison test"
```

---

## Task 2: `run_mhd_smallbox.slurm`

**Files:**
- Create: `TESTDIR/run_mhd_smallbox.slurm`

- [ ] **Step 1: Write `run_mhd_smallbox.slurm`**

Create `TESTDIR/run_mhd_smallbox.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=smallbox-test
#SBATCH --account=eost
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mail-type=fail
#SBATCH --mail-user=changgoo@princeton.edu
#SBATCH --output=smallbox-%j.out
#SBATCH --error=smallbox-%j.err

usage="Usage: sbatch $0 <branch> <r_return>"
BRANCH=${1:-}
R_RETURN=${2:-}

[[ -z "$BRANCH" ]]   && echo "$usage" && exit 1
[[ -z "$R_RETURN" ]] && echo "$usage" && exit 1

case "$BRANCH" in
    tigris-master)        EXE_NAME="tigris_master_mhd.exe" ;;
    mass-return-refactor) EXE_NAME="tigris_mass_return_refactor_mhd.exe" ;;
    *) echo "Unknown branch: $BRANCH"; exit 1 ;;
esac

module purge
module load anaconda3/2023.3 intel-oneapi/2024.2 \
    openmpi/oneapi-2024.2/4.1.6 \
    hdf5/oneapi-2024.2/openmpi-4.1.6/1.14.4 \
    fftw/oneapi-2024.2/3.3.10

prob=tigress_classic
PID=TIGRESS
SRCDIR=$HOME/tigris
SCRIPTDIR=$SLURM_SUBMIT_DIR        # directory from which sbatch was called
TBLDIR="$SRCDIR/inputs/tables"
COOL_TBL="tigress_coolftn.txt"
POPSYNTH_TBL="Z014_GenevaV00.txt"
INPUT="athinput.$prob"
# athinput lives one level above the tiger/ subdirectory in tigris_scripts
ATHINPUT_SRC="$HOME/tigris_scripts/$prob/$INPUT"

RUNDIR="/scratch/gpfs/EOST/$USER/tigress_classic/smallbox-test/${BRANCH}/mhd-r${R_RETURN}"

echo "BRANCH=$BRANCH  R_RETURN=$R_RETURN"
echo "EXE=$SCRIPTDIR/$EXE_NAME"
echo "RUNDIR=$RUNDIR"

params="job/problem_id=$PID time/tlim=10 \
    mesh/nx3=64 mesh/x3min=-512 mesh/x3max=512 \
    cooling/coolftn_file=$COOL_TBL \
    feedback/pop_synth_file=$POPSYNTH_TBL \
    orbital_advection/Omega0=0.0"

extra_params="perturbation/rseed=1 \
    particle1/fgas=0.7 particle1/r_return=$R_RETURN \
    gravity/solve_grav_hyperbolic_dt=false \
    mesh/mhd_outflow_bc=diode problem/beta0=1 \
    output2/dt=5 output3/dt=5 \
    hydro/neighbor_flooring=false cooling/ceiling=true \
    feedback/vmax=1.e4 hydro/dfloor=1.e-8 hydro/pfloor=1.e-8 \
    feedback/tdec_rt=1 particle1/type_ia_sn=false feedback/fbinary=0.0 \
    time/integrator=rk2 time/xorder=2 \
    mesh/ix1_bc=periodic mesh/ox1_bc=periodic"

if [ -d "$RUNDIR" ]; then
    echo "WARNING: $RUNDIR exists, cleaning..."
    rm -rf "$RUNDIR"/*
else
    mkdir -p "$RUNDIR"
fi

cd "$RUNDIR"
cp "$ATHINPUT_SRC" "$INPUT"
cp "$TBLDIR/$COOL_TBL" .
cp "$TBLDIR/$POPSYNTH_TBL" .
cp "$SCRIPTDIR/$EXE_NAME" .

set -o pipefail
srun "./$EXE_NAME" -i "$INPUT" -t 00:55:00 $params $extra_params \
    1> out.txt 2> err.txt
EXITCODE=$?
set +o pipefail

echo "EXITCODE=$EXITCODE"
exit $EXITCODE
```

- [ ] **Step 2: Syntax-check**

```bash
bash -n /Users/changgoo/Sources/tigris/tst/tigress_classic/run_mhd_smallbox.slurm
```

Expected: no output.

- [ ] **Step 3: Commit**

```bash
git -C /Users/changgoo/Sources/tigris add tst/tigress_classic/run_mhd_smallbox.slurm
git -C /Users/changgoo/Sources/tigris commit -m "Add run_mhd_smallbox.slurm for smallbox branch comparison test"
```

---

## Task 3: `submit_tests.sh`

**Files:**
- Create: `TESTDIR/submit_tests.sh`

- [ ] **Step 1: Write `submit_tests.sh`**

Create `TESTDIR/submit_tests.sh`:

```bash
#!/bin/bash
# Submit 6 SLURM jobs: 2 branches x 3 r_return values.
# Must be run from tst/tigress_classic/ so SLURM_SUBMIT_DIR resolves correctly.
# Usage: ./submit_tests.sh [--dry-run]
set -e

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

BRANCHES=(tigris-master mass-return-refactor)
R_RETURNS=(128 256 512)

for branch in "${BRANCHES[@]}"; do
    for r_return in "${R_RETURNS[@]}"; do
        if [ "$DRY_RUN" -eq 1 ]; then
            echo "DRY-RUN: sbatch run_mhd_smallbox.slurm $branch $r_return"
        else
            jobid=$(sbatch run_mhd_smallbox.slurm "$branch" "$r_return" | awk '{print $NF}')
            echo "Submitted: branch=$branch r_return=$r_return jobid=$jobid"
        fi
    done
done
```

- [ ] **Step 2: Make executable and syntax-check**

```bash
chmod +x /Users/changgoo/Sources/tigris/tst/tigress_classic/submit_tests.sh
bash -n /Users/changgoo/Sources/tigris/tst/tigress_classic/submit_tests.sh
```

Expected: no output.

- [ ] **Step 3: Dry-run test**

```bash
cd /Users/changgoo/Sources/tigris/tst/tigress_classic
bash submit_tests.sh --dry-run
```

Expected (6 lines):
```
DRY-RUN: sbatch run_mhd_smallbox.slurm tigris-master 128
DRY-RUN: sbatch run_mhd_smallbox.slurm tigris-master 256
DRY-RUN: sbatch run_mhd_smallbox.slurm tigris-master 512
DRY-RUN: sbatch run_mhd_smallbox.slurm mass-return-refactor 128
DRY-RUN: sbatch run_mhd_smallbox.slurm mass-return-refactor 256
DRY-RUN: sbatch run_mhd_smallbox.slurm mass-return-refactor 512
```

- [ ] **Step 4: Commit**

```bash
git -C /Users/changgoo/Sources/tigris add tst/tigress_classic/submit_tests.sh
git -C /Users/changgoo/Sources/tigris commit -m "Add submit_tests.sh for smallbox branch comparison test"
```

---

## Task 4: `compare.py` — parsing functions with unit tests

**Files:**
- Create: `TESTDIR/test_compare.py`
- Create: `TESTDIR/compare.py` (parsing functions only)

- [ ] **Step 1: Write `test_compare.py` with failing tests for parse functions**

Create `TESTDIR/test_compare.py`:

```python
import re
import textwrap
import numpy as np
import pytest
from compare import parse_timing, parse_hst

SAMPLE_OUT = textwrap.dedent("""\
    cycle=100 time=5.00000e-01 dt=5.00000e-03 zone-cycles/cpu_s=2.00000e+06
    cycle=200 time=1.00000e+00 dt=5.00000e-03 zone-cycles/cpu_s=4.00000e+06
""")

SAMPLE_HST = textwrap.dedent("""\
    # Athena++ history data
    # [1]=time      [2]=dt        [3]=mass      [4]=1-mom    
     0.000000e+00  5.000000e+00  1.000000e+08  1.000000e+05
     5.000000e+00  5.000000e+00  1.100000e+08  2.000000e+05
     1.000000e+01  5.000000e+00  1.200000e+08  3.000000e+05
""")


def test_parse_timing_mean(tmp_path):
    f = tmp_path / "out.txt"
    f.write_text(SAMPLE_OUT)
    result = parse_timing(str(f))
    assert result == pytest.approx(3.0e6)  # mean of 2e6 and 4e6


def test_parse_timing_missing_file(tmp_path):
    result = parse_timing(str(tmp_path / "missing.txt"))
    assert result is None


def test_parse_timing_no_cycles(tmp_path):
    f = tmp_path / "out.txt"
    f.write_text("no timing data here\n")
    assert parse_timing(str(f)) is None


def test_parse_hst_columns(tmp_path):
    f = tmp_path / "TIGRESS.hst"
    f.write_text(SAMPLE_HST)
    result = parse_hst(str(f))
    assert set(result.keys()) == {"time", "dt", "mass", "1-mom"}


def test_parse_hst_values(tmp_path):
    f = tmp_path / "TIGRESS.hst"
    f.write_text(SAMPLE_HST)
    result = parse_hst(str(f))
    np.testing.assert_allclose(result["time"], [0.0, 5.0, 10.0])
    np.testing.assert_allclose(result["mass"], [1e8, 1.1e8, 1.2e8])
    np.testing.assert_allclose(result["1-mom"], [1e5, 2e5, 3e5])


def test_parse_hst_missing_file(tmp_path):
    result = parse_hst(str(tmp_path / "missing.hst"))
    assert result is None
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd /Users/changgoo/Sources/tigris/tst/tigress_classic
python -m pytest test_compare.py -v 2>&1 | head -20
```

Expected: `ModuleNotFoundError: No module named 'compare'`

- [ ] **Step 3: Write `compare.py` with parsing functions**

Create `TESTDIR/compare.py`:

```python
#!/usr/bin/env python3
"""Compare timing and cross-branch equivalence for smallbox test runs."""
import os
import re
import sys
import numpy as np

BRANCHES = ["tigris-master", "mass-return-refactor"]
R_RETURNS = [128, 256, 512]
EQUIV_THRESHOLD = 1e-10


def parse_timing(out_file):
    """Return mean zone-cycles/cpu_s from Athena++ stdout, or None."""
    if not os.path.exists(out_file):
        return None
    pattern = re.compile(r'zone-cycles/cpu_s=(\S+)')
    values = []
    with open(out_file) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                values.append(float(m.group(1)))
    return float(np.mean(values)) if values else None


def parse_hst(hst_file):
    """Return dict of column_name -> np.array from Athena++ .hst file, or None."""
    if not os.path.exists(hst_file):
        return None
    col_pattern = re.compile(r'\[(\d+)\]=(\S+)')
    headers = {}
    with open(hst_file) as f:
        for line in f:
            if line.startswith('#') and '[1]=' in line:
                for m in col_pattern.finditer(line):
                    headers[int(m.group(1)) - 1] = m.group(2)  # 0-based
    data = np.loadtxt(hst_file, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {headers[i]: data[:, i] for i in sorted(headers)}


def get_rundir(scratchbase, branch, r_return):
    return os.path.join(
        scratchbase, "tigress_classic", "smallbox-test",
        branch, f"mhd-r{r_return}"
    )


def main():
    user = os.environ.get("USER", "unknown")
    default_scratch = f"/scratch/gpfs/EOST/{user}"
    scratchbase = sys.argv[1] if len(sys.argv) > 1 else default_scratch

    results = {}
    for branch in BRANCHES:
        for r_return in R_RETURNS:
            rundir = get_rundir(scratchbase, branch, r_return)
            if not os.path.isdir(rundir):
                print(f"WARNING: missing {rundir}")
                continue
            hst_candidates = [
                os.path.join(rundir, f)
                for f in os.listdir(rundir) if f.endswith('.hst')
            ]
            results[(branch, r_return)] = {
                "zcps": parse_timing(os.path.join(rundir, "out.txt")),
                "hst":  parse_hst(hst_candidates[0]) if hst_candidates else None,
            }

    # Table 1: Timing
    print("\n=== Timing (mean zone-cycles/cpu_s) ===")
    print(f"{'branch':<30} {'r_return':>8} {'mean_zcps':>12} {'speedup':>10}")
    print("-" * 65)
    for r_return in R_RETURNS:
        master_zcps = results.get(("tigris-master", r_return), {}).get("zcps")
        for branch in BRANCHES:
            key = (branch, r_return)
            if key not in results:
                continue
            zcps = results[key]["zcps"]
            zcps_str = f"{zcps:.3e}" if zcps else "N/A"
            speedup = (
                f"{zcps / master_zcps:.2f}x"
                if (zcps and master_zcps and branch != "tigris-master")
                else ("1.00x" if branch == "tigris-master" else "N/A")
            )
            print(f"{branch:<30} {r_return:>8} {zcps_str:>12} {speedup:>10}")

    # Table 2: Cross-branch equivalence
    print("\n=== Cross-branch equivalence (mass-return-refactor vs tigris-master) ===")
    print(f"{'r_return':>8} {'column':<20} {'max_abs_diff':>14} {'max_rel_diff':>14} {'status':>6}")
    print("-" * 68)
    for r_return in R_RETURNS:
        master  = results.get(("tigris-master", r_return), {}).get("hst")
        refactor = results.get(("mass-return-refactor", r_return), {}).get("hst")
        if master is None or refactor is None:
            print(f"{r_return:>8}  {'(missing data)'}")
            continue
        for col in master:
            if col not in refactor:
                continue
            a, b = master[col], refactor[col]
            n = min(len(a), len(b))
            abs_diff = float(np.max(np.abs(a[:n] - b[:n])))
            scale = float(np.max(np.abs(a[:n]))) + 1e-300
            rel_diff = abs_diff / scale
            status = "OK" if rel_diff <= EQUIV_THRESHOLD else "FAIL"
            print(f"{r_return:>8} {col:<20} {abs_diff:>14.3e} {rel_diff:>14.3e} {status:>6}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd /Users/changgoo/Sources/tigris/tst/tigress_classic
python -m pytest test_compare.py -v
```

Expected:
```
test_compare.py::test_parse_timing_mean          PASSED
test_compare.py::test_parse_timing_missing_file  PASSED
test_compare.py::test_parse_timing_no_cycles     PASSED
test_compare.py::test_parse_hst_columns          PASSED
test_compare.py::test_parse_hst_values           PASSED
test_compare.py::test_parse_hst_missing_file     PASSED
6 passed
```

- [ ] **Step 5: Commit**

```bash
git -C /Users/changgoo/Sources/tigris add \
    tst/tigress_classic/compare.py \
    tst/tigress_classic/test_compare.py
git -C /Users/changgoo/Sources/tigris commit -m "Add compare.py and unit tests for smallbox comparison"
```

---

## Task 5: `compare.py` — equivalence comparison tests

**Files:**
- Modify: `TESTDIR/test_compare.py` (add `compare_hst` tests)
- The `compare_hst` logic is already inside `main()` in compare.py — extract it to a testable function

- [ ] **Step 1: Add `compare_hst` tests to `test_compare.py`**

Append to `TESTDIR/test_compare.py`:

```python
from compare import compare_hst

HST_A = {
    "time": np.array([0.0, 5.0, 10.0]),
    "mass": np.array([1e8, 1.1e8, 1.2e8]),
    "1-mom": np.array([1e5, 2e5, 3e5]),
}
HST_B_IDENTICAL = {k: v.copy() for k, v in HST_A.items()}
HST_B_DIFFERENT = {
    "time": np.array([0.0, 5.0, 10.0]),
    "mass": np.array([1e8, 1.1e8, 1.3e8]),   # last value differs
    "1-mom": np.array([1e5, 2e5, 3e5]),
}


def test_compare_hst_identical():
    rows = compare_hst(HST_A, HST_B_IDENTICAL)
    for row in rows:
        assert row["status"] == "OK", f"Expected OK for {row['col']}"


def test_compare_hst_different():
    rows = compare_hst(HST_A, HST_B_DIFFERENT)
    statuses = {r["col"]: r["status"] for r in rows}
    assert statuses["mass"] == "FAIL"
    assert statuses["1-mom"] == "OK"
    assert statuses["time"] == "OK"


def test_compare_hst_length_mismatch():
    # Shorter refactor array — compare only overlapping prefix
    short = {"mass": np.array([1e8, 1.1e8])}
    rows = compare_hst({"mass": np.array([1e8, 1.1e8, 1.2e8])}, short)
    assert rows[0]["status"] == "OK"
```

- [ ] **Step 2: Run tests — expect failure on compare_hst tests**

```bash
cd /Users/changgoo/Sources/tigris/tst/tigress_classic
python -m pytest test_compare.py::test_compare_hst_identical -v 2>&1 | tail -5
```

Expected: `ImportError: cannot import name 'compare_hst'`

- [ ] **Step 3: Extract `compare_hst` function in `compare.py`**

Add this function above `main()` in `compare.py` (replace the inline logic in `main()` with a call to it):

```python
def compare_hst(master, refactor):
    """Compare two hst dicts column by column.

    Returns list of dicts with keys: col, max_abs_diff, max_rel_diff, status.
    """
    rows = []
    for col in master:
        if col not in refactor:
            continue
        a, b = master[col], refactor[col]
        n = min(len(a), len(b))
        abs_diff = float(np.max(np.abs(a[:n] - b[:n])))
        scale = float(np.max(np.abs(a[:n]))) + 1e-300
        rel_diff = abs_diff / scale
        rows.append({
            "col": col,
            "max_abs_diff": abs_diff,
            "max_rel_diff": rel_diff,
            "status": "OK" if rel_diff <= EQUIV_THRESHOLD else "FAIL",
        })
    return rows
```

Replace the inline comparison loop in `main()` Table 2 section with:

```python
        for row in compare_hst(master, refactor):
            print(f"{r_return:>8} {row['col']:<20} "
                  f"{row['max_abs_diff']:>14.3e} {row['max_rel_diff']:>14.3e} "
                  f"{row['status']:>6}")
```

- [ ] **Step 4: Run all tests — expect all pass**

```bash
cd /Users/changgoo/Sources/tigris/tst/tigress_classic
python -m pytest test_compare.py -v
```

Expected:
```
test_compare.py::test_parse_timing_mean          PASSED
test_compare.py::test_parse_timing_missing_file  PASSED
test_compare.py::test_parse_timing_no_cycles     PASSED
test_compare.py::test_parse_hst_columns          PASSED
test_compare.py::test_parse_hst_values           PASSED
test_compare.py::test_parse_hst_missing_file     PASSED
test_compare.py::test_compare_hst_identical      PASSED
test_compare.py::test_compare_hst_different      PASSED
test_compare.py::test_compare_hst_length_mismatch PASSED
9 passed
```

- [ ] **Step 5: Commit**

```bash
git -C /Users/changgoo/Sources/tigris add \
    tst/tigress_classic/compare.py \
    tst/tigress_classic/test_compare.py
git -C /Users/changgoo/Sources/tigris commit -m "Add compare_hst function and equivalence tests"
```

---

## Post-implementation: Running on tiger

After all jobs complete, collect results:

```bash
# From tst/tigress_classic/ on tiger
python compare.py
# or with explicit scratch path:
python compare.py /scratch/gpfs/EOST/$USER
```

Expected output shape:
```
=== Timing (mean zone-cycles/cpu_s) ===
branch                          r_return    mean_zcps    speedup
-----------------------------------------------------------------
tigris-master                        128    X.XXXe+06      1.00x
mass-return-refactor                 128    X.XXXe+06      X.XXx
...

=== Cross-branch equivalence (mass-return-refactor vs tigris-master) ===
r_return column               max_abs_diff   max_rel_diff status
--------------------------------------------------------------------
     128 time                  0.000e+00      0.000e+00     OK
     128 mass                  X.XXXe+XX      X.XXXe-XX     OK
...
```

All status entries should be `OK`. Any `FAIL` indicates the refactor changed physics output and needs investigation.
