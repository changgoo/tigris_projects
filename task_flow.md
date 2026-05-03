# Task Flow: CRMHD + Shearing-Box, No Multilevel

Configuration: `rk2` integrator, `ops_task=true`, shearing-box boundary conditions.

**Active flags** (from `athinput.tigress_classic` + `configure.py --prob=tigress_classic -b --cr=mg`):

| Flag | Value |
|------|-------|
| `SHEAR_PERIODIC` | true |
| `ORBITAL_ADVECTION` | false (OAorder=0) |
| `MAGNETIC_FIELDS_ENABLED` | true |
| `CR_ENABLED` | true (multigroup, `--cr=mg`) |
| `PARTICLES` | true (complex stellar particles) |
| `ops_task` | true |
| `multilevel` | false |
| `STS_ENABLED` | false |
| `DUAL_ENERGY_ENABLED` | false |
| `SELF_GRAVITY_ENABLED` | 3 (blockfft) |

---

## Overview

Three separate task lists are orchestrated by `main.cpp` each cycle:

```
for stage in [1, 2]:               ← RK2 has 2 main stages
    TimeIntegratorTaskList(stage)
    BlockFFT::Solve(stage)         ← gravity, called between stages from main.cpp

RayTrace()                         ← if ray_tracing enabled

OperatorSplitTaskList(nstages)     ← once per cycle, after all RK2 stages

Particles::ProcessNewParticles()   ← assign global IDs to newly created particles
```

---

## 1. TimeIntegratorTaskList — 2 stages (RK2 / SSPRK Heun)

### Startup (called before each stage)

- **Stage 1 only:** zero-clear registers `u1`, `b1.x{1,2,3}f`, `u_cr1`
- `ComputeShear(time + sbeta*dt, time + ebeta*dt)` — update shear offset for this stage
- `StartReceivingSubset(all)` — post all MPI receives for main integration

### Task graph (per stage)

Tasks with no listed dependency run immediately; all others wait for their deps.

#### Flux calculation

```
DIFFUSE_HYD (dep: NONE)  ──┐
DIFFUSE_FLD (dep: NONE)  ──┴──► CALC_HYDFLX
```

- `DIFFUSE_HYD` — compute hydro viscous / conductive diffusion fluxes
- `DIFFUSE_FLD` — compute magnetic field diffusion fluxes
- `CALC_HYDFLX` — MHD Riemann solver (HLLD); applies FOFC if enabled

#### Hydro flux correction (shear-periodic)

```
CALC_HYDFLX ──► SEND_HYDFLX
CALC_HYDFLX ──► RECV_HYDFLX
RECV_HYDFLX ──► SEND_HYDFLXSH
(SEND_HYDFLX | RECV_HYDFLX) ──► RECV_HYDFLXSH ──► INT_HYD
```

- `SEND/RECV_HYDFLX` — exchange hydro fluxes at block boundaries for flux correction
- `SEND/RECV_HYDFLXSH` — shearing-box remapping of fluxes at X1 boundaries
- `INT_HYD` — update hydro conserved vars `u` using corrected fluxes

#### CR flux correction (shear-periodic)

```
CALC_HYDFLX ──► CALC_CRFLX
CALC_CRFLX ──► SEND_CRFLX
CALC_CRFLX ──► RECV_CRFLX
RECV_CRFLX ──► SEND_CRFLXSH
(SEND_CRFLX | RECV_CRFLX) ──► RECV_CRFLXSH ──► INT_CR
```

- `CALC_CRFLX` — compute CR energy/flux advection (uses velocity from `CALC_HYDFLX`)
- `SEND/RECV_CRFLXSH` — shear remapping of CR fluxes
- `INT_CR` — update CR conserved vars `u_cr`

#### Magnetic field (CT, shear-periodic)

```
CALC_HYDFLX ──► CALC_FLDFLX ──► SEND_FLDFLX ──► RECV_FLDFLX
RECV_FLDFLX ──► SEND_EMFSH
RECV_FLDFLX ──► RECV_EMFSH
RECV_EMFSH  ──► INT_FLD ──► SEND_FLD ──► RECV_FLD ──► SETB_FLD ──► SEND_FLDSH ──► RECV_FLDSH
```

- `CALC_FLDFLX` — compute EMF from Riemann solver for constrained transport
- `SEND/RECV_FLDFLX` — communicate EMF at block boundaries (required even without AMR)
- `SEND/RECV_EMFSH` — shear remapping of EMF
- `INT_FLD` — integrate `B` via constrained transport
- `SETB_FLD` — set field ghost zones; `SEND/RECV_FLDSH` — shear ghost exchange

#### Passive scalars (if NSCALARS > 0, shear-periodic)

```
(CALC_HYDFLX | DIFFUSE_SCLR) ──► CALC_SCLRFLX
CALC_SCLRFLX ──► SEND_SCLRFLX ──► RECV_SCLRFLX
               ──► RECV_SCLRFLX ──► SEND_SCLRFLXSH
(SEND_SCLRFLX | RECV_SCLRFLX) ──► RECV_SCLRFLXSH ──► INT_SCLR
```

#### Hydro source terms

```
(INT_HYD | INT_SCLR) ──► SRC_TERM
```

- `SRC_TERM` → `AddSourceTerms()`: adds gravity acceleration (from BlockFFT potential),
  Coriolis/tidal forces (shearing box), and any user-defined source terms

#### CR source terms

```
(INT_CR | SRC_TERM) ──► SRCTERM_CR
```

- `SRCTERM_CR` → `AddSourceTermsCR()`: CR streaming, diffusion, hadronic/Coulomb losses,
  gas–CR energy exchange. Uses `cr_integrator=rk1` sub-integrator internally.

#### CR boundary communication (shear-periodic)

```
SRCTERM_CR ──► SEND_CR ──► RECV_CR (dep: NONE) ──► SETB_CR ──► SEND_CRSH ──► RECV_CRSH
```

At this point `src_term = SRC_TERM | SRCTERM_CR`.

#### Hydro/scalar boundary communication (no orbital advection, shear-periodic)

```
src_term ──► SEND_HYD ──► RECV_HYD (dep: NONE) ──► SETB_HYD ──► SEND_HYDSH ──► RECV_HYDSH
src_term ──► SEND_SCLR ──► RECV_SCLR (dep: NONE) ──► SETB_SCLR ──► SEND_SCLRSH ──► RECV_SCLRSH
```

#### Boundary sync point

```
set_boundary = RECV_HYDSH | RECV_FLDSH | RECV_SCLRSH | SEND_CR | SETB_CR
```

This is the combined "all ghost zones are up to date" barrier used as dep for `CONS2PRIM`.

#### Particle integration (independent of hydro)

```
INT_PAR (dep: NONE) ──► SEND_PAR ──► RECV_PAR (dep: NONE) ──► SEND_PARSH ──► RECV_PARSH
```

- `INT_PAR` — integrate particle trajectories using primitive vars from the *previous*
  timestep. This is intentional: `CONS2PRIM` has not yet run, so primitives are
  well-defined and deterministic.
- `SEND/RECV_PARSH` — exchange particles that have crossed the shear boundary

#### Particle-mesh density (ops_task=true: no feedback here)

```
(recvpar | set_boundary) ──► SEND_PM ──► RECV_PM ──► SETB_PM ──► SEND_PMSH ──► RECV_PMSH
```

- `SEND_PM` — deposit particle mass onto mesh, communicate PM density
- When `ops_task=false`, `CREATE_PAR` / `INTERACT` (feedback) would appear here instead.
  With `ops_task=true` they are deferred to `OperatorSplitTaskList`.

#### Primitive variable recovery

```
set_boundary ──► CONS2PRIM
```

- Converts updated conserved vars to primitives; applies density/pressure floors
  and neighbor flooring if enabled.

#### Physical boundaries and CR opacity

```
(CONS2PRIM | SETB_CR | SEND_CR | recvpar | RECV_PMSH) ──► PHY_BVAL ──► CR_OPACITY
                                                                    ──► CLEAR_ALLBND
```

- `PHY_BVAL` — apply physical boundary conditions (user-defined Z boundaries, etc.)
- `CR_OPACITY` — update CR scattering coefficients / opacity after primitives are set
- `CLEAR_ALLBND` — release MPI boundary buffers for this stage
- **`USERWORK` is NOT called here** when `ops_task=true`; it runs in `OperatorSplitTaskList`

---

## 1b. Gravity — between RK2 stages (main.cpp)

After each stage, `main.cpp` calls:

```cpp
if (ptlist->CheckNextMainStage(stage))
    pmesh->my_blocks(0)->pfft->Solve(stage)   // SELF_GRAVITY_ENABLED == 3 (blockfft)
```

- Solves the Poisson equation for gas self-gravity using a block-FFT algorithm
- Deposits the resulting gravitational potential/acceleration, which is then picked
  up by `SRC_TERM → AddSourceTerms()` in the next stage

---

## 2. OperatorSplitTaskList — 1 stage, once per cycle

Runs after all RK2 stages complete. All particle-intensive and stiff operator-split
physics lives here.

### Task graph

```
CREATE_PAR (dep: NONE)               ← form new star particles from gas that meets
    │                                   accretion criteria (last stage only)
    └──► SET_FBR                     ← compute feedback region properties (e.g.,
              │                         feedback radius rfb_res, resolved/unresolved flag)
              └──► SEND_GPAR         ← send ghost (halo) particles to neighbor blocks
RECV_GPAR (dep: NONE)
(SEND_GPAR after SET_FBR, RECV_GPAR) ──► SEND_GPARSH ──► RECV_GPARSH   (shear exchange)

recvgpar ──► INTERACT                ← gas–particle interaction:
                  │                     SN energy/momentum injection (Classic scheme),
                  │                     mass return to gas, CR energy injection (fe_CR),
                  │                     accretion onto existing stars
                  │
                  └──► OPS_INT_COOLING   ← operator-split cooling/heating:
                            │               TIGRESS cooling function, subcycled with
                            │               cfl_cool_sub; includes PE heating
                            │
                            ├──► REMOVE_PAR              ← remove particles flagged for deletion
                            │
                            ├──► SEND_HYD                ← communicate updated hydro after
                            │    RECV_HYD (dep: NONE)       cooling+feedback energy deposition
                            │    SETB_HYD
                            │    SEND_HYDSH ──► RECV_HYDSH
                            │
                            └──► SEND_SCLR               ← communicate updated scalars
                                 RECV_SCLR (dep: NONE)
                                 SETB_SCLR
                                 SEND_SCLRSH ──► RECV_SCLRSH

set_boundary = RECV_HYDSH | RECV_SCLRSH
(no CR comms here; B field unchanged by cooling/feedback)

(set_boundary | REMOVE_PAR) ──► CONS2PRIM ──► PHY_BVAL ──► USERWORK ──► CLEAR_ALLBND
```

- `USERWORK` — runs `UserWorkInLoop()` from `tigress_classic.cpp`: history outputs,
  phase diagnostics, z-profile outputs, etc.

### After OperatorSplitTaskList (main.cpp)

```
Particles::ProcessNewParticles()    ← assign unique global IDs to particles created
                                       by CREATE_PAR this cycle
```

---

## Summary of Operator Splitting

| Physics | Where | Notes |
|---------|-------|-------|
| MHD fluxes + CT | TimeIntegrator, each stage | HLLD + constrained transport |
| CR advection | TimeIntegrator, each stage | Advected with gas velocity |
| CR source terms | TimeIntegrator, after `SRC_TERM` | Streaming, losses, coupling |
| Gravity (Poisson) | `main.cpp`, between stages | BlockFFT; result used in `SRC_TERM` |
| Particle motion | TimeIntegrator, independent | Uses primitives from previous step |
| PM density | TimeIntegrator, after particle recv | No feedback in integrator stages |
| Star formation | OperatorSplit | `CREATE_PAR` runs once per cycle |
| SN feedback | OperatorSplit | After star formation, before cooling |
| Cooling/heating | OperatorSplit | After feedback; subcycled internally |
| Ray tracing | `main.cpp`, after integrator | Before OperatorSplit |
| UserWork / outputs | OperatorSplit | At end, after cooling |
