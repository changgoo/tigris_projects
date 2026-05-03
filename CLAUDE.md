# TIGRIS — CRMHD / tigress_classic

TIGRIS is a private fork of Athena++ (grid-based GRMHD + AMR) extended with ISM physics.
Primary use case here: **CRMHD simulations** with the `tigress_classic` problem generator.

- Main branch for PRs: **`tigris-master`**
- Scripts repo: `$HOME/tigris_scripts/tigress_classic/`
- Doxygen docs: https://changgoo.github.io/athena

## Build

```bash
# From $HOME/tigris_scripts/tigress_classic/
./build_tigress.sh tiger crmhd
# Usage: ./build_tigress.sh <machine> <physics> [build_option] [src_dir] [flux]
#   machine: tiger | stellar | anvil
#   physics: crmhd | crmhd_duale | crmhd_duals | mhd | hydro | ...
#   build_option: 0=normal (default), 1=debug (g++), 2=no-clean
```

What this does for `tiger crmhd`:
1. Loads intel-oneapi/2024.2, openmpi, hdf5, fftw modules
2. Runs `configure.py` with:
   ```
   --prob=tigress_classic --nghost=4 -fft -fb --grav=blockfft -mpi -hdf5
   -b --cr=mg --flux=hlld --cxx=icpx
   ```
3. `make all -j4`
4. Copies binary to `tiger/tigris_crmhd.exe`

Key `configure.py` flags for CRMHD:
| Flag | Meaning |
|------|---------|
| `-b` | Enable MHD (magnetic fields) |
| `--cr=mg` | Cosmic ray transport via multigroup |
| `--flux=hlld` | HLLD Riemann solver |
| `-fb` | Enable feedback module |
| `--grav=blockfft` | Block-FFT self-gravity |
| `--nghost=4` | 4 ghost cells (required) |

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
  gravity/                         # Block-FFT self-gravity
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

## Key Physics Notes

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
cells. The `ExchangeGhostAccretionDelta()` function communicates ghost accretion deltas
back to the origin so that the particle's total accreted mass is correct and conserved.

Currently uses `MPI_Allgatherv` (known limitation — see GitHub issue for P2P refactor plan).

### Shear-Periodic Boundaries

Shear-periodic boundaries (x1-direction) add complexity because the x-offset between
neighboring blocks changes with time (`qshear * Omega0 * Lx * t`). Key concerns:

- Ghost particles must have their positions shifted by the shear offset
- `noverlap_` controls how many copies of each near-boundary particle are created
  (set to `NGHOST` to ensure proximity checks work correctly)
- FOFC (first-order flux correction) is **skipped on MeshBlock boundary faces**
  (`first_order_flux_correction.cpp`) to avoid conservation issues at shear-periodic
  interfaces

## FOFC (First-Order Flux Correction)

FOFC detects unphysical states (negative density/pressure) after Riemann solves and
falls back to first-order reconstruction on affected faces. At MeshBlock boundaries
this can create conservation mismatches because neighboring blocks may apply different
flux corrections. Current fix: skip FOFC for boundary faces entirely
(`first_order_flux_correction.cpp`).

## Testing

```bash
cd tst/regression
python run_tests.py                     # all regression tests
./test_mpi.sh                           # MPI-specific tests
```

## Coding Conventions
- C++11/14, BSD 3-Clause license
- Follow Athena++ style guide (snake_case, Doxygen comments)
- Add regression test for new functionality
- Reference issue numbers in commits: `Fixes #42`
