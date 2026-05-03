# TIGRIS — Source Code Structure

TIGRIS is a private fork of Athena++ (grid-based GRMHD + AMR) extended with ISM physics.
Primary use case: **CRMHD simulations** with the `tigress_classic` problem generator.

- Base code: `$HOME/tigris` (branch: `tigris-master`)
- Scripts repo: `$HOME/tigris_scripts/tigress_classic/`
- Doxygen docs: https://changgoo.github.io/athena

## Source Layout (CRMHD-relevant)

```
src/
  athena.hpp / athena_arrays.hpp   # Core types, AthenaArray
  defs.hpp(.in)                    # Compile-time defines (generated)
  main.cpp                         # Entry point
  mesh/                            # AMR mesh, MeshBlock, time integration
  hydro/                           # MHD solver, reconstruction, Riemann
  field/                           # Magnetic field evolution (CT)
  eos/                             # Equation of state (adiabatic γ=5/3)
  cr/, cr_mg/                      # Cosmic ray transport — multigroup (streaming, diffusion)
  fft/                             # Distributed FFT infrastructure
    athena_fft.hpp/cpp             # FFTBlock — multi-meshblock FFT (cuboid, fftmpi backend)
    fft_driver.cpp                 # FFTDriver — cuboid MeshBlock grouping, plan creation
    block_fft.hpp/cpp              # BlockFFT — single-meshblock FFT (1 MB/rank, MPI only)
    fftmpi/                        # fftMPI C++ library (active backend)
    plimpton/                      # Plimpton C library (legacy, being retired)
  gravity/                         # Self-gravity (block-FFT Poisson solvers)
    fft_gravity.hpp/cpp            # FFTGravity : FFTBlock — multi-MB/rank, periodic/open/disk/shear BCs
    block_fft_gravity.hpp/cpp      # BlockFFTGravity : BlockFFT — 1 MB/rank, all BCs (legacy)
  particles/                       # Stellar/complex particles
    particles.hpp/cpp              # Base particle class
    complex_particles.cpp          # Star particles (accretion, mass return)
    accretion.cpp/hpp
    mass_return.cpp/hpp
    particle_gravity.cpp/hpp       # Particle-mesh gravity
    particle_mesh.cpp/hpp
  feedback/                        # SN/stellar feedback injection
    feedback.hpp/cpp
    classic_scheme.cpp             # Classic (resolved) feedback scheme
    subcell_scheme.cpp             # Sub-cell feedback
    injection.cpp/hpp
    pop_synth.cpp/hpp              # Population synthesis (SB99)
    pop_synth_table.cpp
  microphysics/                    # Cooling/heating, units
    cooling_function.cpp
    cooling_solver.cpp
    cooling.hpp
    units.cpp/hpp
  outputs/                         # HDF5, history, zprof, parbin outputs
  task_list/                       # Task-based operator splitting
  coordinates/                     # Cartesian shearing-box
  orbital_advection/               # Shear (Omega0, qshear)
  bvals/                           # Boundary conditions (shear-periodic, user)
  pgen/
    tigress_classic.cpp            # PRIMARY problem generator
```

## Input File Structure

Template: `$HOME/tigris_scripts/tigress_classic/athinput.tigress_classic`

Key blocks and their roles:
| Block | Purpose |
|-------|---------|
| `<mesh>` | Shearing-box domain; X1/X2 shear-periodic, X3 user BC |
| `<meshblock>` | 32³ per MPI block |
| `<time>` | rk2 integrator, cr_integrator=rk1 |
| `<hydro>` | γ=5/3, FOFC, ct_method=uct_hlld |
| `<orbital_advection>` | Shear: qshear=1, Omega0=0.028 |
| `<cooling>` | Op-split cooling, tigress cooling function |
| `<cr>` | CR physics flags (streaming, losses, self-consistent) |
| `<particle1>` | Complex (stellar) particles: feedback, accretion, mass return |
| `<particle2>` | Runaway particles (type=none by default) |
| `<feedback>` | Classic scheme, SB99 pop synthesis, fe_CR=0.1 |
| `<problem>` | ISM parameters: Sigma_gas, Sigma_star, ext_grav, turbulence |
| `<gravity>` | grav_bc=disk |

Output types: `hst`, `hdf5` (prim + uov), `parbin`, `rst`, `phase_hst`, `zprof`

## Execution Flow

Per-cycle loop in `main.cpp`:

```
for stage in [1, 2]:                         ← RK2 stages
    TimeIntegratorTaskList(stage)            ← MHD, CR, scalars, particle motion
    BlockFFT/FFTGravity::Solve(stage)        ← Poisson gravity (between stages)

OperatorSplitTaskList()                      ← once per cycle (ops_task=true):
    CREATE_PAR → SET_FBR → SEND/RECV_GPAR
    → INTERACT (Merge→Accrete→Feedback)
    → OPS_INT_COOLING → REMOVE_PAR
    → boundary comms → CONS2PRIM → USERWORK

Particles::ProcessNewParticles()             ← assign global IDs to new particles
```

`INT_PAR` (particle trajectory integration) uses primitives from the **previous** RK
stage — intentional, so results are deterministic across MPI configurations.
`USERWORK` calls `UserWorkInLoop()` from `tigress_classic.cpp`: history, phase
diagnostics, z-profile outputs.

See `TASKLIST.md` for the full task-dependency graph.

## Key Physics

- **Domain**: Stratified shearing box, typically 1×1×6 kpc (64×64×384 cells at 16 pc/cell)
- **CR transport**: multigroup (`--cr=mg`), streaming + diffusion + losses, self-consistent Alfven speed
- **Feedback**: SN energy split — `fkin` fraction kinetic, `fe_CR=0.1` to CRs, rest thermal
- **Cooling**: Operator-split with sub-cycling (`cfl_cool_sub`), TIGRESS cooling function
- **Gravity**: Block-FFT for gas self-gravity + external disk potential (`ext_grav=force`)
- **Particles**: Star particles formed by accretion; return mass + SN energy after delay

## Particle System Architecture

### Ghost Particles & Boundary Communication

Particles near MeshBlock boundaries are copied as **ghost particles** to neighboring
blocks so that operations like accretion and feedback can access data across boundaries.

Communication flow in the operator-split task list (`ops_task_list.cpp`):
```
CREATE_PAR → SET_FBR → SEND_GPAR/RECV_GPAR
  → [SEND_GPARSH/RECV_GPARSH if shear-periodic]
  → INTERACT → ... → REMOVE_PAR
```

- `SEND_GPAR`/`RECV_GPAR`: Exchange ghost particles with neighbors (P2P via `ParticleBuffer`)
- `SEND_GPARSH`/`RECV_GPARSH`: Additional shear-periodic exchange for x1-boundaries
- `INTERACT` (`InteractWithMesh()`): Merge, Accrete, then DoFeedback
- `REMOVE_PAR`: Clean up old/dead particles

### Key Files for Particle Communication

| File | Role |
|------|------|
| `particles_bvals.cpp` | `SendParticleBuffer()`, `ReceiveParticleBuffer()`, shear-periodic ghost exchange |
| `particle_buffer.hpp` | `ParticleBuffer` class — packs/unpacks particle data for MPI |
| `complex_particles.cpp` | `InteractWithMesh()`, `Accrete()`, `ExchangeGhostAccretionDelta()`, `Merge()` |
| `accretion.cpp/hpp` | Accretion logic — gas → particle mass transfer |
| `mass_return.cpp/hpp` | `MassReturn::CollectParticlesInfo()`, mass/energy return from aged stars |

### Accretion & Ghost Conservation

When a particle near a boundary accretes gas, both the **active** copy (on the origin
MeshBlock) and the **ghost** copy (on the neighbor) may accrete from their respective
cells. `ExchangeGhostAccretionDelta()` communicates ghost accretion deltas back to the
origin so that the particle's total accreted mass is correct and conserved.

Currently uses `MPI_Allgatherv` (known limitation — see GitHub issue for P2P refactor plan).

### Shear-Periodic Boundaries

Shear-periodic boundaries (x1-direction) add complexity because the x-offset between
neighboring blocks changes with time (`qshear * Omega0 * Lx * t`). Key concerns:

- Ghost particles must have their positions shifted by the shear offset
- `noverlap_` controls how many copies of each near-boundary particle are created
  (set to `NGHOST` to ensure proximity checks work correctly)
- FOFC is **skipped on MeshBlock boundary faces** (`first_order_flux_correction.cpp`)
  to avoid conservation issues at shear-periodic interfaces

## FOFC (First-Order Flux Correction)

FOFC detects unphysical states (negative density/pressure) after Riemann solves and
falls back to first-order reconstruction on affected faces. At MeshBlock boundaries
this can create conservation mismatches because neighboring blocks may apply different
flux corrections. Current fix: skip FOFC for boundary faces entirely
(`first_order_flux_correction.cpp`).

## Coding Conventions

- C++11, BSD 3-Clause license
- Follow Athena++ style guide (snake_case, Doxygen comments; tst/style/check_athena_cpp_style.sh)
- Python code linting (flake8)
- Add regression test for new functionality (tst/regression)
- Reference issue numbers in commits: `Fixes #42`
