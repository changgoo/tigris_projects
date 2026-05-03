# FFT gravity optimization handoff

Date: 2026-05-02

## Current session update: rank-aggregated ShearingRemapper exchange

Branch:

```text
shearing-remapper
```

Latest committed baseline before this edit:

```text
57f7d7c75 Add rank-aware tags and local copies to ShearingRemapper
```

Current uncommitted implementation work:

```text
src/gravity/shearing_remap.cpp
src/gravity/shearing_remap.hpp
```

`ShearingRemapper::RollUnrollAll(dt)` now runs as a coordinated all-local-block
pipeline:

```text
FillGhostZonesAll()
FractionalShiftAll(dt)
IntegerShiftAll(dt)
```

The new ghost and integer-shift paths build deterministic `RemapSegment` lists,
group remote payloads by peer rank, sort send/receive segment lists by the same
tuple, and exchange one nonblocking aggregate payload per peer rank for each phase.
Same-rank block transfers are copied directly. Remote `MeshBlock` objects are never
dereferenced; remote segment indices are inferred from the uniform non-AMR block
layout, while local `GetBlockBuffer(gid)` is used only for local copy, pack, and
unpack operations.

The old `RollUnrollBlock()` implementation is still present in the file as an unused
reference/fallback while this phase is validated.

Build and style status after the edit:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd blockfft
```

completed and copied:

```text
/home/changgoo/tigris_scripts/tigress_classic/tiger/tigris_mhd-blockfft-nofb.exe
```

Style check:

```bash
cd /home/changgoo/tigris/tst/style
./check_athena_cpp_style.sh ../../src/gravity/shearing_remap.cpp ../../src/gravity/shearing_remap.hpp
```

completed with no errors after fixing two line-length violations.

Next validation step:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd blockfft
```

After it finishes, check `err.txt` first. Then compare `phi` between
`mhd-4pc-b1-shear-nofb-blockfft` and `mhd-4pc-b1-shear-nofb-blockfft-refactor`
with `vis/python/compare_phi.py`, and compare timing with
`vis/python/fft_gravity_timing.py`.

`compare_phi.py` needs the `pyathena` conda environment. The alias in
`~/.bashrc` is:

```bash
alias pyathena='module purge; module load anaconda3/2024.6 openmpi/gcc/4.1.6 fftw/gcc/3.3.10; conda activate pyathena'
```

Run the comparison as:

```bash
source ~/.bashrc
pyathena
python vis/python/compare_phi.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft-refactor
```

Validation job submitted after the aggregate exchange edit:

```text
JobID: 2642763
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:02:52
Executable: tigris_mhd-blockfft-nofb.exe
Output directory:
/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft-refactor
```

`err.txt` tail contained only restart MeshBlock write messages. `out.txt` ended at:

```text
Terminating on time limit
time=1.0000000000000000e+00 cycle=174
```

Phi validation:

```bash
bash -lc 'source ~/.bashrc; module purge; module load anaconda3/2024.6 openmpi/gcc/4.1.6 fftw/gcc/3.3.10; conda activate pyathena; python vis/python/compare_phi.py ... 0'
bash -lc 'source ~/.bashrc; module purge; module load anaconda3/2024.6 openmpi/gcc/4.1.6 fftw/gcc/3.3.10; conda activate pyathena; python vis/python/compare_phi.py ... 1'
```

Both `out2.00000` and `out2.00001` are bit-exact:

```text
max|phi_new - phi_ref| = 0.000000e+00
Relative L-inf = 0.000000e+00
Relative L2 = 0.000000e+00
RESULT: BIT-EXACT
```

Timing comparison against the committed BlockFFT baseline:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft-refactor
```

All-row profile mean:

| Metric | Refactor/baseline | Delta |
| --- | ---: | ---: |
| TotalMax | 1.0221 | +0.002856966 |
| ShearSourceMax | 0.7777 | -0.0043691971 |
| RetrieveMax | 1.1439 | +0.0016119027 |
| SelfGravity loop | 1.0664 | +1.8496176 |

There are 13 gas-solve rows with `TotalMax > 0.2`, mostly in FFT forward/backward
and a few retrieve/task-list outliers. Excluding those outlier rows from the new run
only, the remapper effect looks favorable:

| Metric | Refactor mean | Baseline mean | Ratio |
| --- | ---: | ---: | ---: |
| TotalMax | 0.12357399 | 0.12943773 | 0.9547 |
| ShearSourceMax | 0.015367032 | 0.019654601 | 0.7819 |
| ForwardMax | 0.051739989 | 0.052813018 | 0.9797 |
| BackwardMax | 0.053622288 | 0.053451847 | 1.0032 |
| RetrieveMax | 0.011287788 | 0.011201457 | 1.0077 |
| TaskListMax | 0.0061235791 | 0.0059948299 | 1.0215 |

Interpretation: correctness passed. The aggregated remapper substantially reduces
`ShearSourceMax`, but this single run contains large non-remap timing outliers, so
do not judge final performance from the all-row means alone. A repeat validation run
or a robust/outlier-resistant summary is needed before committing the performance
claim.

User accepted the outlier behavior as known and acceptable. The aggregate remapper
change was committed:

```text
5d46841e2 Aggregate ShearingRemapper MPI exchanges
```

## Current session update: Phase 3 FFTGravity source path

Implemented Phase 3 source-side integration, leaving the result path on the existing
FFT-block shearing remap for now.

Changed files:

```text
src/gravity/fft_gravity.cpp
src/gravity/fft_gravity.hpp
```

Implementation details:

- `FFTGravityDriver` now owns a `ShearingRemapper *premapper_` when
  `pmy_mesh_->shear_periodic` is true.
- Added `FFTGravity::LoadShearedSource(const AthenaArray<Real> &src_kij,
  LogicalLocation loc, RegionSize bsize)`.
- In the shearing solve branch:
  - density or particle-augmented density is copied into the remapper buffer in
    `(k,i,j)` layout
  - `premapper_->RollUnrollAll(-dt)` performs the source remap
  - each local MeshBlock remapped buffer is loaded into AthenaFFT input with
    `LoadShearedSource()`
  - old `ApplyShearingSource(-1.0)` is no longer used in the shearing source path
- `ApplyShearingResult(1.0)` and `RetrieveAppliedShearingResult(...)` are unchanged
  and remain the result path for Phase 3.

Build and style status:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd fft
```

completed and copied:

```text
/home/changgoo/tigris_scripts/tigress_classic/tiger/tigris_mhd-fft-nofb.exe
```

Style check:

```bash
cd /home/changgoo/tigris/tst/style
./check_athena_cpp_style.sh ../../src/gravity/fft_gravity.cpp ../../src/gravity/fft_gravity.hpp ../../src/gravity/shearing_remap.cpp ../../src/gravity/shearing_remap.hpp
```

completed with no errors.

Next validation step:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd fft
```

After the job finishes:

1. Check Slurm status and `err.txt`.
2. Compare timing against both `mhd-4pc-b1-shear-nofb-fft` and
   `mhd-4pc-b1-shear-nofb-blockfft`.
3. Phase 3 gate is `FFTGravity ShearSourceMax <= 24 ms` for 64^3.

Phase 3 validation job:

```text
JobID: 2642774
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:02:58
Output: /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-refactor
```

Timing gate passed:

```text
FFT refactor ShearSourceMax = 0.018784267 s
BlockFFT baseline ShearSourceMax = 0.019654601 s
Phase 3 gate = <= 0.024 s
```

Compared with the old FFT reference, the Phase 3 source-only mixed path was:

```text
TotalMax ratio = 0.9365
ShearSourceMax ratio = 0.6222
SelfGravity loop ratio = 0.9428
```

However, final `phi` at `out2.00001` was outside the strict tolerance:

```text
FFT reference vs FFT refactor:
max diff = 7.629395e-04
relative L-inf = 2.302312e-06
relative L2 = 3.074121e-07

BlockFFT reference/refactor vs FFT refactor:
max diff = 9.994507e-04
relative L-inf = 3.016026e-06
relative L2 = 8.971054e-07
```

This is likely because Phase 3 intentionally mixed the new MeshBlock source remap
with the old FFT-block result remap. Therefore Phase 3 was not committed separately;
Phase 4 was implemented immediately.

## Current session update: Phase 4 FFTGravity result path

Implemented Phase 4 result-side integration on top of Phase 3.

Additional changes:

- Added `FFTGravity::RetrieveShearedResult(AthenaArray<Real> &dst_kij,
  LogicalLocation loc, RegionSize bsize)`.
- The shearing result path now:
  - retrieves raw inverse FFT output into remapper buffers in `(k,i,j)` layout,
    applying `norm_factor_` at retrieval time
  - calls `premapper_->RollUnrollAll(+dt)`
  - copies remapped buffers into `phi` or `phi_gasonly`
- Old `ApplyShearingResult(1.0)` and `RetrieveAppliedShearingResult(...)` are no
  longer used in the shearing solve branch.

Build and style after Phase 4:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd fft

cd /home/changgoo/tigris/tst/style
./check_athena_cpp_style.sh ../../src/gravity/fft_gravity.cpp ../../src/gravity/fft_gravity.hpp ../../src/gravity/shearing_remap.cpp ../../src/gravity/shearing_remap.hpp
```

Both passed.

Next validation step:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd fft
```

After completion, compare timing and `phi` again. Phase 4 gate is
`RetrieveMax <= 13 ms` and full-loop `SelfGravity` within 5% of BlockFFT, unless
the failure is purely a `RetrieveMax`/profile outlier requiring further profiling.

Phase 4 validation job:

```text
JobID: 2642791
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:02:54
Output: /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-refactor
```

`err.txt` tail contained only restart MeshBlock write messages. `out.txt` ended at:

```text
Terminating on time limit
time=1.0000000000000000e+00 cycle=179
```

Timing results versus old FFT reference and BlockFFT target:

```text
FFT refactor TotalMax       = 0.13185439 s
FFT refactor ShearSourceMax = 0.020443974 s
FFT refactor RetrieveMax    = 0.012947751 s
FFT refactor SelfGravity    = 28.803994 s

Old FFT TotalMax            = 0.14386894 s
Old FFT ShearSourceMax      = 0.030190929 s
Old FFT RetrieveMax         = 0.015436346 s
Old FFT SelfGravity         = 31.6595 s

BlockFFT TotalMax           = 0.12943773 s
BlockFFT ShearSourceMax     = 0.019654601 s
BlockFFT RetrieveMax        = 0.011201457 s
BlockFFT SelfGravity        = 27.849465 s
```

Ratios:

```text
FFT refactor / old FFT:
TotalMax       = 0.9165
ShearSourceMax = 0.6772
RetrieveMax    = 0.8388
SelfGravity    = 0.9098

FFT refactor / BlockFFT:
TotalMax       = 1.0187
ShearSourceMax = 1.0402
RetrieveMax    = 1.1559
SelfGravity    = 1.0343
```

Phase 4 timing gates:

```text
RetrieveMax gate <= 13 ms: passed (12.947751 ms)
SelfGravity within 5% of BlockFFT: passed (3.43% slower)
```

Phi comparisons at `out2.00001`:

```text
Old FFT reference vs FFT refactor:
max diff = 1.907349e-03
relative L-inf = 5.755779e-06
relative L2 = 6.305215e-07
RESULT: FAIL by old-FFT strict gate

BlockFFT reference vs FFT refactor:
max diff = 1.831055e-04
relative L-inf = 5.525545e-07
relative L2 = 8.334752e-07
RESULT: WARN by script text, but relative L-inf is below 1e-6

BlockFFT refactor vs FFT refactor:
max diff = 1.831055e-04
relative L-inf = 5.525545e-07
relative L2 = 8.334752e-07
RESULT: WARN by script text, but relative L-inf is below 1e-6
```

Interpretation:

- The refactored FFT path now follows the BlockFFT-style remapper for both source
  and result.
- It no longer matches the old FFT global-row remap at strict tolerance, but it does
  match BlockFFT/BlockFFT-refactor below `1e-6` relative L-inf.
- This is consistent with intentionally migrating FFTGravity to the BlockFFT-style
  shearing remap. Treat the old-FFT comparison as an algorithm-change diagnostic,
  not a rollback trigger, unless the review insists on bitwise compatibility with
  the previous FFT global-row remap.

## Goal

Optimize `FFTGravity`/AthenaFFT self-gravity so the shearing-box `fft` gravity path
approaches the speed of `BlockFFTGravity` for the TIGRESS classic MHD benchmark.
The remaining comparison target is same meshblock size, same total block count:

```text
64^3 meshblock, 1 MeshBlock per rank, Nblocks=128
fft       vs blockfft
```

Earlier work compared `FFTGravity` with `64^3` meshblocks against `FFTGravity` with
`32^3` meshblocks, where the AthenaFFT path supports multiple MeshBlocks per rank
(`8 MB/rank`). The current implementation has made those two `FFTGravity` cases
close enough. Do not spend more optimization time on the `fft` vs `fft-8mb` gap
unless new data says it regressed.

Benchmark location:

```text
/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-*
```

The relevant benchmark variants found so far are:

```text
mhd-4pc-b1-shear-nofb-fft       # 64^3 MeshBlock, 1 MB/rank, AthenaFFT path
mhd-4pc-b1-shear-nofb-fft-8mb   # 32^3 MeshBlock, 8 MB/rank, AthenaFFT path
mhd-4pc-b1-shear-nofb-blockfft  # 64^3 MeshBlock, 1 MB/rank, BlockFFT path
```

## Build and submission commands

Build commands are run from:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
```

FFTGravity/AthenaFFT build:

```bash
./build_tigress_nofb.sh tiger mhd fft
```

BlockFFTGravity build:

```bash
./build_tigress_nofb.sh tiger mhd blockfft
```

Job submissions are run from:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
```

64^3 FFTGravity job:

```bash
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd fft
```

64^3 BlockFFTGravity job:

```bash
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd blockfft
```

32^3, 8 MeshBlocks/rank FFTGravity job:

```bash
sbatch tigress_classic_mhd_shear_nofb_8mb.slurm -i mhd fft
```

Monitor submitted jobs with:

```bash
squeue -u changgoo
```

## Current uncommitted files

Tracked code changes:

```text
src/gravity/fft_gravity.cpp
src/gravity/fft_gravity.hpp
src/gravity/block_fft_gravity.cpp
src/gravity/block_fft_gravity.hpp
```

New/untracked files observed:

```text
vis/python/fft_gravity_timing.py
CLAUDE.md
TASKLIST.md
fofc_conservation_fix.md
particle_accretion_conservation.md
src/outputs/zprof_outputs.md
tst/test_notebooks/radiative_snr_particle.ipynb
tst/test_notebooks/radiative_snr_solvers.ipynb
```

Only the gravity source/header diffs and `vis/python/fft_gravity_timing.py` appear
directly related to this optimization. The notebooks and other notes may belong to
separate work.

## Profiling instrumentation

Both gravity paths now have optional timing output controlled by new input flags:

```text
<gravity>
profile_fft_gravity = true
profile_block_fft_gravity = true
```

`FFTGravityDriver` writes:

```text
<problem_id>.fft_gravity_time.txt
```

`BlockFFTGravity` writes:

```text
<problem_id>.block_fft_gravity_time.txt
```

Each row is a comma-separated key-value record with local timings reduced over MPI:

```text
Total, Particle, Load, ShearSource, Forward, Kernel, Backward, Retrieve, TaskList
```

For each phase, the file records both `*Sum` and `*Max`; the most useful quantity
for wall-clock comparison is the rank-max mean.

The timing output also records:

```text
ncycle, time, stage, mode, gas_only, Nblocks, shear, grav_bc
```

`mode` currently exists only in the `FFTGravity` output.

## Timing analysis script

New script:

```bash
python vis/python/fft_gravity_timing.py RUN_DIR [RUN_DIR ...]
```

It looks for either:

```text
TIGRESS.fft_gravity_time.txt
TIGRESS.block_fft_gravity_time.txt
```

and optionally:

```text
TIGRESS.loop_time.txt
```

It prints per-run phase means/min/max, phase fractions of `TotalMax`, worst solve,
and pairwise ratios relative to the first run argument.

Primary benchmark comparison command:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft
```

## Benchmark results observed

Baseline: `mhd-4pc-b1-shear-nofb-fft`

The actionable comparison is `fft` vs `blockfft`, both with `Nblocks=128`.
`fft-8mb` is retained only as historical context for the multiple-MeshBlocks-per-rank
work.

Profile means, rank-max:

| Run | Rows | Nblocks | TotalMax | ShearSourceMax | ForwardMax | BackwardMax | RetrieveMax | TaskListMax |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| fft | 363 | 128 | 0.14339926 | 0.030093307 | 0.053869524 | 0.053066694 | 0.014874995 | 0.0069285554 |
| fft-8mb | 363 | 1024 | 0.14523799 | 0.031377972 | 0.056785874 | 0.051971194 | 0.01447353 | 0.0092503147 |
| blockfft | 353 | 128 | 0.12943773 | 0.019654601 | 0.052813018 | 0.053451847 | 0.011201457 | 0.0059948299 |

Loop timing means:

| Run | All | TimeIntegratorTaskList | SelfGravity | OpSplit | NewDt | Particle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fft | 122.90894 | 73.317776 | 32.211606 | 12.967965 | 0.78872541 | 0.92305729 |
| fft-8mb | 163.56059 | 83.794024 | 32.822365 | 39.5652 | 0.93008835 | 2.8402241 |
| blockfft | 116.40759 | 71.876653 | 27.849465 | 12.570259 | 0.78164353 | 0.79167347 |

Ratios versus `fft`:

| Run | TotalMax | SelfGravity loop | All loop |
| --- | ---: | ---: | ---: |
| fft-8mb / fft | 1.0128 | 1.0190 | 1.3307 |
| blockfft / fft | 0.9026 | 0.8646 | 0.9471 |

Interpretation:

- `blockfft` is about 9.7% faster than `fft` in profile `TotalMax`.
- `blockfft` is about 13.5% faster in loop `SelfGravity`.
- The biggest `blockfft` wins are `ShearSourceMax` and `RetrieveMax`.
- Forward/backward FFT costs are already similar between `fft` and `blockfft`.
- `fft-8mb` is not the remaining optimization target; it was used to validate the
  multiple-MeshBlocks-per-rank AthenaFFT work.
- `blockfft` has 353 profile rows versus 363 for the two `fft` runs, so keep the
  sample-count caveat in mind.

Same-64^3 `blockfft / fft` profile deltas:

| Phase | fft mean | blockfft mean | blockfft/fft | Delta |
| --- | ---: | ---: | ---: | ---: |
| TotalMax | 0.14339926 | 0.12943773 | 0.9026 | -0.013961525 |
| ShearSourceMax | 0.030093307 | 0.019654601 | 0.6531 | -0.010438706 |
| RetrieveMax | 0.014874995 | 0.011201457 | 0.7530 | -0.0036735381 |
| TaskListMax | 0.0069285554 | 0.0059948299 | 0.8652 | -0.00093372551 |
| ForwardMax | 0.053869524 | 0.052813018 | 0.9804 | -0.0010565064 |
| KernelMax | 0.0035135744 | 0.003513006 | 0.9998 | -0.000000568 |
| BackwardMax | 0.053066694 | 0.053451847 | 1.0073 | +0.000385153 |
| LoadMax | 0.0031398725 | 0.0049114471 | 1.5642 | +0.001771575 |
| ParticleMax | 0.0025600811 | 0.0032902972 | 1.2852 | +0.000730216 |

Optimization priority from these numbers:

1. `FFTGravity::ApplyShearingSource()` is the main gap. It is about 10.4 ms slower
   than the corresponding shearing source/remap work in `BlockFFTGravity`.
2. `FFTGravity::ApplyShearingResult()` plus `RetrieveAppliedShearingResult()` is the
   second gap. The measured `RetrieveMax` difference is about 3.7 ms.
3. Do not focus on AthenaFFT forward/backward/kernel unless new profiling contradicts
   this run; those phases are already within roughly 1 ms, and kernel is identical.
4. `LoadMax` and `ParticleMax` are actually slower in `blockfft`, so they are not the
   reason `fft` loses overall.

## Partial post-instrumentation result

While the rebuilt `64^3` `fft` benchmark was still running, the timing file had
303 rows. This is not the final benchmark result, but it already shows the current
shape of the remaining gap.

Command:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft
```

Partial `fft` vs existing `blockfft`:

| Phase | fft partial mean | blockfft mean | blockfft/fft | Delta |
| --- | ---: | ---: | ---: | ---: |
| TotalMax | 0.13709555 | 0.12943773 | 0.9441 | -0.007657817 |
| ShearSourceMax | 0.02404242 | 0.019654601 | 0.8175 | -0.004387819 |
| RetrieveMax | 0.014659152 | 0.011201457 | 0.7641 | -0.003457696 |
| ForwardMax | 0.052368123 | 0.052813018 | 1.0085 | +0.000444895 |
| BackwardMax | 0.05389835 | 0.053451847 | 0.9917 | -0.000446502 |

Relative to the previous `fft` log, `ShearSourceMax` improved from about 30.1 ms
to about 24.0 ms. The remaining total profile gap is now about 7.7 ms rather than
about 14.0 ms.

FFT shearing subphase rank-max means from the partial `fft` log:

| Subphase | Mean |
| --- | ---: |
| ShearSourceCopyMax | 0.002247699 |
| ShearSourceGhostMax | 0.011795335 |
| ShearSourceOffsetMax | 0.000003658 |
| ShearSourceRollMax | 0.002373386 |
| ShearSourceIndexMax | 0.000017613 |
| ShearSourceExchangeMax | 0.003876193 |
| ShearSourceScatterMax | 0.012801286 |
| ShearResultCopyMax | 0.002285110 |
| ShearResultGhostMax | 0.008962249 |
| ShearResultOffsetMax | 0.000002140 |
| ShearResultRollMax | 0.001873586 |
| ShearResultIndexMax | 0.000016364 |
| ShearResultExchangeMax | 0.004024694 |
| ShearRetrieveScatterMax | 0.004446027 |

Do not sum these subphase `Max` values as if they were one rank's timeline; each
subphase is reduced independently by rank max. Still, the priority is clear:

1. `FillShearingGhosts()` dominates both source and result remaps.
2. Source scatter back into AthenaFFT input is also large.
3. Result exchange and retrieve scatter are secondary.
4. Offset and row-index construction are negligible after the current changes.

## Implemented ghost-fill optimization

A follow-up optimization targeted `FillShearingGhosts()` directly:

- local/single-shear-rank ghost fills now use contiguous `std::memcpy()` per `k` slab
  instead of nested `k/i/j` scalar loops
- MPI ghost exchange now packs the lower and upper active slabs contiguously
- both y-direction receives and sends are posted together with nonblocking MPI
  (`MPI_Irecv`/`MPI_Isend` + `MPI_Waitall`) instead of two sequential
  `MPI_Sendrecv` calls
- reusable send/receive buffers are kept on `FFTGravity`

This should reduce the measured `ShearSourceGhostMax` and `ShearResultGhostMax`
hotspots without changing the remap arithmetic.

Build verification:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd fft
```

The build completed and copied:

```text
/home/changgoo/tigris_scripts/tigress_classic/tiger/tigris_mhd-fft-nofb.exe
```

The linker emitted the same LLVM gold plugin vectorization warning seen before, but
the executable was produced successfully.

## Post ghost-fill optimization benchmark

Submitted after the nonblocking/slab-copy ghost-fill change:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd fft
```

Job `2641127` ran the simulation executable successfully, but the overall Slurm job
was marked failed because the follow-up Python step failed:

```text
2641127.0  tigris_mhd-fft-nofb.exe  COMPLETED  0:0
2641127.1  python                   FAILED     1:0
```

The timing log had 363 rows and was usable. Result versus the same existing
`blockfft` baseline:

| Phase | fft mean | blockfft mean | blockfft/fft | Delta |
| --- | ---: | ---: | ---: | ---: |
| TotalMax | 0.14386894 | 0.12943773 | 0.8997 | -0.014431211 |
| ShearSourceMax | 0.030190929 | 0.019654601 | 0.6510 | -0.010536328 |
| RetrieveMax | 0.015436346 | 0.011201457 | 0.7257 | -0.004234889 |
| SelfGravity loop | 31.6595 | 27.849465 | 0.8797 | -3.8100353 |

Key FFT shearing subphase means:

| Subphase | Mean |
| --- | ---: |
| ShearSourceGhostMax | 0.020356079 |
| ShearSourceScatterMax | 0.012483911 |
| ShearSourceExchangeMax | 0.007859430 |
| ShearResultGhostMax | 0.009997791 |
| ShearResultExchangeMax | 0.004451013 |
| ShearRetrieveScatterMax | 0.004263766 |

Conclusion: the nonblocking/slab-copy ghost-fill change did not improve the full
benchmark. The ghost-fill rank-max time stayed around 20 ms for source remap and
around 10 ms for result remap. The large rank-max cost is therefore likely dominated
by communication imbalance/wait time or the exchange pattern itself, not the local
pack/unpack scalar loops.

## Current code direction in `FFTGravity`

The main optimization attempt is for shearing-periodic gravity. Based on the same-64^3
comparison, optimize these functions first:

```text
FFTGravity::ApplyShearingSource()
FFTGravity::ApplyShearingResult()
FFTGravity::RetrieveAppliedShearingResult()
FFTGravity::FillShearingGhosts()
FFTGravity::ExchangeShearingRows()
```

The matching `blockfft` reference path is:

```text
BlockFFTGravity::RollUnroll()
```

`RollUnroll()` handles fractional and integer shearing shifts in one routine using
MeshBlock-local buffers plus direct neighbor/target exchanges. `FFTGravity` currently
does a full local FFT-block copy into `shear_ghost_`, a ghost fill, a fractional roll
into `shear_roll_`, row-index construction, a y-row exchange, then a scatter. That
extra staging is the likely source of the remaining gap.

The current patch already addresses the multi-MeshBlock-per-rank retrieval overhead:
it moved result unroll from per-MeshBlock retrieval to one whole-local-FFT-block pass.

Before the patch, `RetrieveShearingResult(...)` did the shearing-coordinate unroll
and row exchange per destination MeshBlock. That meant repeated setup work while
retrieving `phi`/`phi_gasonly`.

The patch splits this into:

```cpp
void FFTGravity::ApplyShearingResult(Real dt);
void FFTGravity::RetrieveAppliedShearingResult(AthenaArray<Real> &dst, bool nu,
                                               int ngh, LogicalLocation loc,
                                               RegionSize bsize);
```

Current flow in `FFTGravityDriver::Solve()` for shearing runs:

```text
Load source
ApplyShearingSource(-1.0)
ExecuteForward()
ApplyKernel(mode)
ExecuteBackward()
ApplyShearingResult(1.0)
RetrieveAppliedShearingResult(...) for each local MeshBlock
GravityBoundaryTaskList
```

`ApplyShearingResult(1.0)` now:

- copies the full local FFT result into `shear_ghost_`
- fills shearing ghosts once
- computes fractional roll into `shear_roll_`
- builds the needed global source-row index for the whole local FFT block
- calls `ExchangeShearingRows()` once

`RetrieveAppliedShearingResult(...)` now only scatters the already-unrolled rows
into one MeshBlock destination and applies `norm_factor_` at assignment time.

## Supporting micro-optimizations already in the patch

`FFTGravity` now reuses vectors with `resize(...)` instead of `assign(..., 0.0)` in
the shearing buffers where old contents are overwritten:

```text
shear_column_
shear_ghost_
shear_roll_
```

The patch caches per-x-column shearing quantities:

```cpp
std::vector<int> shear_joffset_;
std::vector<Real> shear_eps_;
```

This avoids recomputing `ceil(s)` and `eps` inside multiple `j/k` loops.

The patch replaces per-call `std::vector<bool> in_index(Nx[1], false)` row marking
with a reusable integer marker:

```cpp
std::vector<int> shear_row_mark_;
int shear_row_mark_id_;
```

This avoids repeated allocation/clear while building `shear_row_index_`.

## Important correctness notes and risks

The shearing-result refactor changes when normalization is applied:

- old `RetrieveShearingResult(...)` copied `out_ * norm_factor_` into `shear_ghost_`
- new `ApplyShearingResult(...)` copies raw `out_` into `shear_ghost_`
- new `RetrieveAppliedShearingResult(...)` multiplies by `norm_factor_` when writing
  to `dst`

This should be equivalent for linear interpolation, but it needs a correctness check.

`ApplyShearingResult(...)` computes `shear_joffset_` and `shear_eps_` for local FFT
block `i = 0..nx1-1`. `RetrieveAppliedShearingResult(...)` indexes
`shear_joffset_[mi]`, where `mi` is the meshblock-local offset inside this FFT block.
This relies on `mi` being in the same local FFT-block coordinate system used by
`ApplyShearingResult(...)`; that matched the old `RetrieveShearingResult(...)`
offset calculation.

`shear_row_mark_id_` increments monotonically. If this code runs for extremely many
calls in one process, integer overflow would eventually break row marking. This is
unlikely in normal runs, but a robust version should reset and clear if it approaches
`std::numeric_limits<int>::max()`.

The profiling code adds file append I/O from rank 0 every gravity solve when enabled.
Keep it disabled for production timing unless intentionally benchmarking.

## Next steps

1. Build with the current patch and run a short shearing-box correctness check against
   the pre-refactor `fft` path if possible.
2. Compare potentials/accelerations or conservation-sensitive diagnostics, not only
   loop timing, because the normalization and retrieval order changed.
3. Re-run the benchmark after the `ApplyShearingResult` split to see whether
   `RetrieveMax` and `ShearSourceMax` moved toward `blockfft`.
4. If `RetrieveMax` remains high, profile inside `ExchangeShearingRows()` to separate
   row-index construction, packing, `MPI_Alltoallv`, and scatter.
5. If `ShearSourceMax` remains high, compare `FFTGravity::ApplyShearingSource()` with
   `BlockFFTGravity::RollUnroll()`; current measurements suggest this is the largest
   algorithmic gap.
6. Add subphase timers for `FFTGravity` shearing remap:
   `copy -> FillShearingGhosts -> fractional roll -> row index/sort ->
   ExchangeShearingRows -> final scatter`. This should identify whether the gap is
   memory traffic, MPI exchange, or row-index bookkeeping.
7. Consider adding a guarded overflow reset for `shear_row_mark_id_` before committing.

## Useful commands

Check uncommitted files:

```bash
git status --short
git diff --stat
```

Analyze current benchmark logs:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft
```

## Current refactor status, 2026-05-02 00:59 EDT

Active branch:

```text
shearing-remapper
```

Recent commits:

```text
c27896a07 Extract BlockFFTGravity shearing remap into ShearingRemapper
2a48419dc Carry over gravity profiling baseline
```

Phase 1 status:

- `src/gravity/shearing_remap.hpp/cpp` exists.
- `BlockFFTGravity::RollUnroll(dat, dt)` remains as a compatibility wrapper.
- The hot shearing path in `BlockFFTGravity::Solve()` now uses
  `premapper_->Buffer(pmy_block_)` directly and calls `premapper_->RollUnrollAll(dt)`
  directly, avoiding wrapper copy overhead.
- Lazy remapper buffer initialization is required. Do not move buffer allocation back
  into the constructor; doing so segfaulted because `Mesh::my_blocks` was not ready
  while MeshBlocks were being constructed.
- AMR is guarded in `ShearingRemapper` constructor.

Phase 1 validation:

- Job `2641347` initially validated the constructor crash fix but showed copy-overhead
  regression.
- After direct hot-path buffer use, job `2641382` completed successfully.
- Timing for `mhd-4pc-b1-shear-nofb-blockfft-refactor` versus original
  `mhd-4pc-b1-shear-nofb-blockfft`:

```text
TotalMax       ratio 0.9528
ShearSourceMax ratio 0.6785
RetrieveMax    ratio 0.9922
SelfGravity    ratio 1.0153
All loop        ratio 1.0076
```

This is acceptable for the Phase 1 timing gate.

Phase 2 in-progress local changes:

- Modified but not committed:

```text
src/gravity/shearing_remap.cpp
src/gravity/shearing_remap.hpp
```

- The current local Phase 2 substep is incremental rank/tag plumbing, not the final
  rank-aggregated payload exchange.
- Added `remap_phys_id_ = pm->ReserveTagPhysIDs(1)`.
- Added `CreateRemapTag(gid, bufid)` using
  `BoundaryBase::CreateBvalsMPITag(local_id, bufid, remap_phys_id_)`.
- Added `GetRank(gid)` and `GetLocalId(gid)` helpers using `ranklist/nslist`.
- Replaced `gid`-as-rank MPI destinations with rank-aware destinations.
- Added explicit same-rank copies for ghost and integer-shift receive paths.

Validation currently running:

```text
job id: 2641797
command: sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd blockfft
state at last check: RUNNING, 2 nodes, tiger-i01c4n[2-3]
```

After job `2641797` finishes:

1. Check `err.txt` in
   `/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft-refactor`.
2. Compare timing:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft-refactor
```

3. If the 1 MB/rank run completes cleanly and timing is not regressed, commit:

```bash
git add src/gravity/shearing_remap.cpp src/gravity/shearing_remap.hpp
git commit -m "Add rank-aware tags and local copies to ShearingRemapper"
```

4. Then continue Phase 2 with the real rank-aggregated metadata+payload exchange.
   The current per-segment MPI shape is only an intermediate validation step.

## 2026-05-02 Phase 5 Cleanup In Progress

Current branch: `shearing-remapper`.

Last committed implementation checkpoint:

```text
8d08581b3 Use ShearingRemapper in FFTGravity shearing paths
```

The user accepted the current timing, including the known outlier behavior, and asked
to move on. The next step is the cleanup commit planned after the FFTGravity
integration.

Local tracked changes currently in progress:

```text
src/gravity/fft_gravity.cpp
src/gravity/fft_gravity.hpp
```

Cleanup performed so far:

- Removed the old FFTGravity global-row shearing remap API and implementation:
  `ApplyShearingSource`, `ApplyShearingResult`,
  `RetrieveAppliedShearingResult`, `FillShearingGhosts`,
  `ExchangeShearingRows`, `BuildShearingRowPositions`, and
  `EnsureShearingIndexCache`.
- Removed the old FFTGravity-owned shearing row/cache buffers and
  `MPI_COMM_SHEAR`; the shearing data motion now belongs to
  `ShearingRemapper`.
- Removed the old FFT shearing subphase profile columns, which were zero after
  the new remapper path became the only shearing path.
- Kept `SetShearQuantities()` and `qomt_`; they are still live because the
  sheared gravity kernel uses `qomt_`.

Validation still needed for this cleanup:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd fft
cd /home/changgoo/tigris/tst/style
./check_athena_cpp_style.sh ../../src/gravity/fft_gravity.cpp ../../src/gravity/fft_gravity.hpp
```

If build and style pass, commit:

```bash
git add src/gravity/fft_gravity.cpp src/gravity/fft_gravity.hpp
git commit -m "Remove old FFTGravity global-row shearing remap"
```

Update after build/style:

- FFT build completed and copied
  `/home/changgoo/tigris_scripts/tigress_classic/tiger/tigris_mhd-fft-nofb.exe`.
- Style check passed:

```bash
cd /home/changgoo/tigris/tst/style
./check_athena_cpp_style.sh ../../src/gravity/fft_gravity.cpp ../../src/gravity/fft_gravity.hpp
```

32^3 / 8 MB-rank validation:

- We had a pre-refactor 32^3 baseline directory:
  `/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb`.
- No refactor 32^3 directory existed before this check.
- Submitted the missing validation job:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb_8mb.slurm -i mhd fft
```

- Job id: `2642824`.
- At submit check it was running on 2 nodes:
  `tiger-i01c12n2,tiger-i02c6n1`.
- Expected output directory from the script:
  `/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb-refactor`.

After job `2642824` finishes:

1. Check `err.txt` and `out.txt` in the `fft-8mb-refactor` directory.
2. Compare timing against `fft-8mb`:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb-refactor
```

3. Run phi comparison if both outputs exist:

```bash
source ~/.bashrc
module purge
module load anaconda3/2024.6 openmpi/gcc/4.1.6 fftw/gcc/3.3.10
conda activate pyathena
python vis/python/compare_phi.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-8mb-refactor \
  1
```

Final 32^3 validation results:

- Slurm batch `2642824` ended `FAILED 1:0` only because the post-run plotting
  step failed immediately (`2642824.1 python FAILED 1:0`).
- The simulation step completed successfully:
  `2642824.0 tigris_mh+ COMPLETED 0:0`, elapsed `00:03:49`, MaxRSS
  `686304K`.
- `out.txt` terminated normally:

```text
Terminating on time limit
time=1.0000000000000000e+00 cycle=179
```

- `err.txt` tail only showed restart MeshBlock writes.
- Timing versus pre-refactor 32^3 FFT baseline:

```text
profile Max mean:
  TotalMax       ratio 0.9221
  ShearSourceMax ratio 0.6755
  RetrieveMax    ratio 1.0622

loop timing mean:
  All            ratio 0.9814
  SelfGravity    ratio 0.9397
```

- Phi versus pre-refactor 32^3 FFT baseline shows the expected algorithm-change
  mismatch already seen in the 64^3 comparison:

```text
max|phi_new - phi_ref| = 1.907349e-03
Relative L-inf         = 5.755779e-06
Relative L2            = 6.305215e-07
```

- The more relevant meshblock-size consistency check passed bit-exactly:
  64^3 refactor (`fft-refactor`) versus 32^3 refactor (`fft-8mb-refactor`),
  `out2.00001`, `phi`:

```text
max|phi_new - phi_ref| = 0.000000e+00
Relative L-inf         = 0.000000e+00
Relative L2            = 0.000000e+00
RESULT: BIT-EXACT
```

- Timing 32^3 refactor versus 64^3 refactor:

```text
profile Max mean:
  TotalMax       ratio 1.0459
  ShearSourceMax ratio 1.1506
  RetrieveMax    ratio 1.1628

loop timing mean:
  SelfGravity    ratio 1.0704
  All            ratio 1.3360
```

Interpretation: the refactored algorithm works for multiple MeshBlocks per rank.
The 32^3/8 MB-rank gravity solve is close to the 64^3/1 MB-rank case and produces
bit-exact `phi`; remaining `All` loop overhead is outside the gravity solve.

Cleanup commit created:

```text
9e8a02c52 Remove old FFTGravity global-row shearing remap
```
