# Small-box comparison test: mass-return-refactor vs tigris-master

## Purpose

Validate that the P2P routing refactor in `mass-return-refactor` produces
physically identical results to `tigris-master` across three `r_return` values,
and measure the timing difference. The three values (128, 256, 512 pc) exercise
all three routing code paths on a ±512 pc z-domain:

| r_return | path in mass-return-refactor |
|----------|------------------------------|
| 128 pc   | geometric routing (`RouteToRanks`) |
| 256 pc   | geometric routing |
| 512 pc   | full-domain broadcast (`IsFullDomain`, since 512 ≥ 0.5 × 1024) |

## Test location

`$TIGRIS/tst/tigress_classic/`  where `$TIGRIS = $SRCBASE/tigris`

## Files

| File | Role |
|------|------|
| `build_both.sh` | Build both branches via git worktree, copy exes here |
| `run_mhd_smallbox.slurm` | Parameterized SLURM job (branch + r_return as args) |
| `submit_tests.sh` | Submit all 6 jobs |
| `compare.py` | Parse timing + `.hst` files, print comparison tables |

## Domain and physics

| Parameter | Value | Notes |
|-----------|-------|-------|
| `mesh/nx3` | 64 | 16 pc/cell, same as lowres |
| `mesh/x3min/max` | ±512 pc | |
| `mesh/nx1`, `mesh/nx2` | from `athinput.tigress_classic` | unchanged |
| `time/tlim` | 10 | |
| `physics` | mhd | |
| `fgas` | 0.7 | |
| `beta` | 1 | |
| `MHDBC` | diode | |
| `output2/dt`, `output3/dt` | 5 | 2 snapshots |
| SLURM ranks | 8 (1 node) | scaled from lowres 64 by nx3 ratio |
| `--time` | 1:00:00 | |

## `build_both.sh`

```
Usage: ./build_both.sh [srcbase=auto] [build_option=0]

SRCBASE auto-detect:
  tiger / stellar → $HOME
  local mac       → $HOME/Sources

Flow:
  1. Set SRCBASE (hostname-based default or explicit arg)
  2. MAINDIR=$SRCBASE/tigris  (assumed on tigris-master)
  3. WTDIR=$MAINDIR/.worktrees/mass-return-refactor
     git -C $MAINDIR worktree add $WTDIR mass-return-refactor  (skip if exists)
  4. Load tiger modules (intel-oneapi, openmpi, hdf5, fftw)
  5. Build tigris-master in $MAINDIR:
       configure.py --prob=tigress_classic --nghost=4 -fft -fb --grav=fft
                    -mpi -hdf5 -b --flux=hlld --cxx=icpx
       make all -j4
       cp bin/athena → tst/tigress_classic/tigris_master_mhd.exe
  6. Build mass-return-refactor in $WTDIR: same flags
       cp bin/athena → tst/tigress_classic/tigris_mass_return_refactor_mhd.exe
```

Inlines configure/make logic rather than calling `build_tigress.sh` to avoid
its `$HOME/$SRC` path assumption.

## `run_mhd_smallbox.slurm`

```
Usage: sbatch run_mhd_smallbox.slurm <branch> <r_return>

RUNDIR: /scratch/gpfs/EOST/$USER/tigress_classic/smallbox-test/{branch}/mhd-r{r_return}/

EXE map:
  tigris-master            → tigris_master_mhd.exe
  mass-return-refactor     → tigris_mass_return_refactor_mhd.exe

Athena++ overrides:
  mesh/nx3=64 mesh/x3min=-512 mesh/x3max=512
  time/tlim=10
  particle1/r_return=<r_return>
  output2/dt=5 output3/dt=5
```

No auto-resubmit; tlim=10 completes in a single job.

## `submit_tests.sh`

```bash
for branch in tigris-master mass-return-refactor; do
  for r_return in 128 256 512; do
    sbatch run_mhd_smallbox.slurm $branch $r_return
  done
done
```

## `compare.py`

```
Usage: python compare.py [scratchbase=/scratch/gpfs/EOST/$USER]

Table 1 — Timing:
  branch                    r_return  mean_zcps   speedup_vs_master
  tigris-master             128       ...         1.00×
  mass-return-refactor      128       ...         X.XX×
  ...
  (zone-cycles/cpu_s parsed from out.txt)

Table 2 — Cross-branch equivalence (for each r_return):
  r_return  hst_column    max_abs_diff  max_rel_diff  status
  128       total_mass    ...           ...           OK / FAIL
  128       total_mom_x   ...           ...           OK / FAIL
  ...
  (full .hst time-series compared column by column; flag if max_rel_diff > 1e-10)
```

Absolute conservation is not expected (source terms, cooling, feedback).
The check is that both branches produce the **same** `.hst` evolution at each
`r_return`. Missing run directories are flagged as warnings.
