# CLAUDE.md: TIGRIS/Athena++ Project Guide

## Code style
run tst/style/check_athena_cpp_style.sh
should follow C++11 std

## Project Overview

**TIGRIS** is a custom fork of **Athena++**, a finite-volume MHD code with AMR. It adds astrophysical physics for star-forming region simulations (TIGRESS).

- **Primary repo**: PrincetonUniversity/tigris
- **Upstream**: PrincetonUniversity/athena (public)
- **Main branch**: `tigris-master` (protected — all changes via PR)
- **Language**: C++11, MPI/OpenMP, HDF5, FFTW

**Ignore**: `*_gr.cpp`, `*_sr.cpp`, `*_rel*` — general/special relativity code not used in TIGRIS.

## Build

```bash
# Configure (generates Makefile and src/defs.hpp)
./configure.py --prob=tigress_classic --nghost=4 --grav=blockfft \
  -b -fft -mpi -hdf5 --cxx=g++ \
  --lib_path=<hdf5_lib> --include=<hdf5_inc>

make clean && make all -j
# Output: bin/athena
```

Key `configure.py` flags: `--prob`, `--eos`, `--flux`, `--nghost`, `--nscalars`, `-b` (MHD), `-mpi`, `-hdf5`, `-fft`, `--grav`.

## Source Layout (non-GR)

```
src/
  hydro/          Hydro + Riemann solvers (HLLE, HLLD, LLF, Roe)
  field/ct.cpp    Constrained transport (divergence-free B)
  eos/general/    EOS: ideal.cpp, hydrogen.cpp, general_mhd.cpp
  reconstruct/    PPM (ppm.cpp), donor-cell (dc.cpp)
  scalars/        Passive scalar transport + diffusion
  bvals/          Boundary conditions (cc/, fc/, orbital/)
  fft/            Poisson solver (self-gravity)
  pgen/           Problem generators (initial conditions)
  outputs/        HDF5/ASCII I/O
  task_list/      Task scheduling framework
  particles/      ← CURRENT FOCUS
  feedback/       Stellar feedback injection
```

## Execution Flow

```
main.cpp — per cycle:
  for stage = 1..N:
    ptlist->DoTaskListOneStage()     ← RK integration (hydro, field, scalars, particles)
      Last stage only:
        INT_PAR → SEND/RECV_PAR → CREATE_PAR → SEND/RECV_GPAR → INTERACT

  popstlist->DoTaskListOneStage()    ← Operator-split (when ops_task=true)
    CREATE_PAR → SET_FBR → SEND/RECV_GPAR → INTERACT → REMOVE_PAR
    (INTERACT calls ppar->InteractWithMesh())

  Particles::ProcessNewParticles()   ← Assign unique IDs to new particles (after all stages)
```

**INT_PAR** uses primitives from the PREVIOUS RK stage (intentional — ensures deterministic results across MPI configurations).

`InteractWithMesh()` sequence:
1. `Merge()` — merge particles with overlapping AccrCVs
2. `Accrete()` → calls `TestAccrete(i)` per particle
3. `MechanicalFeedback::DoFeedback()` → calls `TestFeedbackOverlap()` then injects

## Particle Module (`src/particles/`)

### Key Files

| File | Purpose |
|------|---------|
| `particles.hpp` | Base `Particles` class + `ComplexParticles` class definitions; all flags |
| `particles.cpp` | Base implementation; `IsControlVolumeOverlap()` at ~line 412 |
| `complex_particles.cpp` | Star/feedback particles; `TestAccrete()` ~553, `TestFeedbackOverlap()` ~803 |
| `accretion.hpp/cpp` | `Accretion` class — AccrCV logic, `Accrete()` |
| `mass_return.hpp/cpp` | Stellar mass return |
| `particle_mesh.hpp/cpp` | Particle↔mesh interpolation |
| `particle_gravity.hpp/cpp` | Particle gravity |
| `particles_bvals.cpp` | Ghost particle communication |
| `particle_buffer.hpp/cpp` | MPI communication buffers |
| `dust_particles.cpp` | Dust particles |
| `tracer_particles.cpp` | Tracer particles |

### Class Hierarchy

```
Particles (base)
  ├── ComplexParticles  ← star particles; accretion + feedback
  ├── DustParticles
  └── TracerParticles
```

**ComplexParticles** key members:
- `pacc` — `Accretion*`, manages AccrCV
- `pmf` — `MechanicalFeedback*`, manages FeedbackCV
- `pmret` — `MassReturn*`
- `const int accretion` — accretion mode (0=off)
- `const bool feedback` — feedback enabled

### Particle Flags (bitmask, `particles.hpp`)

```
GROWING       = 1<<0   accretion enabled (sink particle)
FEEDBACK      = 1<<1   feedback enabled (star particle)
FEEDBACK_NOW  = 1<<2   inject feedback THIS timestep
FEEDBACK_NEXT = 1<<3   feedback postponed → inject NEXT timestep
PASSIVE       = 1<<4   no feedback even if FEEDBACK set
SN_IA         = 1<<5   Type Ia supernova
RUNAWAY       = 1<<6   runaway particle
GHOST         = 1<<7   ghost particle (on neighbor block)
```

Utilities: `SetFlag(i, FLAG)`, `ClearFlag(i, FLAG)`, `TestFlag(i, FLAG)` on particle index `i`.

## Overlap Check Logic (Current Focus)

### Control Volumes

| CV | Radius | Owner |
|----|--------|-------|
| **AccrCV** | `pacc->rctrl = 1` cell (3×3×3 cube) | Accreting particles |
| **FeedbackCV** | `pmf->pischeme->GetNOverlap()` cells (larger) | Feedback particles |

### Detection: `IsControlVolumeOverlap()` — `particles.cpp:~412`

```cpp
bool Particles::IsControlVolumeOverlap(x1,y1,z1,r1, x2,y2,z2,r2, &flag_strong)
```

- Converts positions to integer cell indices
- Compares cubic regions `[pos-r, pos+r]^3` along each axis
- Returns `true` if overlap exists; `flag_strong=true` if one CV is entirely inside the other

### Three Overlap Checks

| Check | Function | CV₁ | CV₂ | Trigger | Action |
|-------|----------|-----|-----|---------|--------|
| Accretion cancel | `TestAccrete(i)` ~553 | AccrCV of `i` | FeedbackCV of `j` | `j` has `FEEDBACK_NOW` | Cancel accretion; delete if `i` is NEW |
| Feedback delay | `TestFeedbackOverlap()` ~803 | FeedbackCV of `i` | FeedbackCV of `j` | Both have `FEEDBACK_NOW` | Postpone lower-SNR to `FEEDBACK_NEXT` |
| Merger | `TestMerger()` | AccrCV of `i` | AccrCV of `j` | Both `GROWING` | Merge particles |

### Accretion Cancellation — `TestAccrete(i)`

**Conditions (all required):**
1. Particle `i` has `GROWING` flag
2. Some neighbor `j` has `FEEDBACK_NOW` (injecting this step)
3. AccrCV of `i` overlaps FeedbackCV of `j`

**Result:**
- NEW particle (just born, `pid == NEW`) → deleted (`pid = DEL`)
- Established particle → accretion skipped this timestep

**Physical intent:** Gas inside an active feedback region should receive that energy, not be accreted away.

### Feedback Delay — `TestFeedbackOverlap()`

**Conditions (all required):**
1. Both particles `i` and `j` have `FEEDBACK_NOW`
2. Their FeedbackCVs overlap

**Result:**
- Particle with lower `snr` (supernova rate) → `FEEDBACK_NOW` cleared, `FEEDBACK_NEXT` set
- That particle injects feedback next timestep instead

**Why feedback regions cannot safely overlap (numerical constraints):**
1. **Order-dependence**: Injection is additive, but when two FeedbackCVs overlap the result depends on the ORDER of particle processing — making outcomes non-reproducible
2. **Meshblock boundary sync**: When overlap spans a meshblock boundary, ghost zones are NOT synchronized at injection time. Each block applies feedback independently using stale ghost cell data from its neighbor → inconsistent combined state across the domain

### Known Limitation (Issue filed)

At high star formation rates, many large accreting particles can cluster near active feedback regions. The current logic may:
- Repeatedly cancel accretion across many particles simultaneously
- Cascade-delay feedback (FEEDBACK_NEXT chains)
→ Suppressing star formation beyond what's physically warranted ("pessimistic" regime).

**GitHub issue**: https://github.com/PrincetonUniversity/tigris/issues/275

## Git Workflow

- **Main branch**: `tigris-master` (protected, requires 1 PR review)
- **PR target**: Always `tigris-master`, not `master`
- **Commit style**: Descriptive; reference `#N` for issue numbers
- **Interactive rebase** before merging to keep history linear

## Useful Links

- **Doxygen docs**: https://changgoo.github.io/athena
- **Wiki**: https://github.com/PrincetonUniversity/tigris/wiki
- **Athena++ paper**: https://ui.adsabs.harvard.edu/abs/2020ApJS..249....4S/abstract
- **TIGRIS Slack**: tigrepp.slack.com

---

**Last updated**: 2026-04-24
