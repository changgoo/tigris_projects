# Code Context: UCT Wave Speed Duplication

## What UCT-HLL/HLLD is

Constrained transport (CT) updates face-centered magnetic fields using corner EMFs.
UCT-HLL and UCT-HLLD (Mignone & Del-Zanna 2021, MDZ'21) weight the corner EMF using
four per-face coefficients (`aL`, `aR`, `dL`, `dR`) that encode how fast information
propagates from each side of the interface.

- **UCT-HLL** (MDZ'21 Eq. 32): uses only the outer HLL wave speeds λ_L, λ_R.
- **UCT-HLLD** (MDZ'21 Eq. 44–45): additionally uses HLLD intermediate states
  λ*, ρ*_L, ρ*_R, λ*_L, λ*_R for sharper weighting near Alfvén waves.

Runtime flags: `pfield->uct_on` (either UCT mode active), `pfield->uct_hlld` (HLLD mode).
Declared in `src/field/field.hpp:66`.

---

## Where the duplication lives

### `src/hydro/calculate_fluxes.cpp` — `Hydro::CalculateHLLWaveSpeed` (line 441)

Called in `CalculateFluxes` for each direction immediately after `RiemannSolver` when
`uct_on` (lines 96–99 x1, 208–211 x2, 319–322 x3). Computes from scratch:

| Quantity | Lines | Cost |
|----------|-------|------|
| `FastMagnetosonicSpeed` × 2 | ~505–506 | 2 EOS calls per cell |
| Roe-averaged state | ~483–500 | ~20 flops per cell |
| λ_L, λ_R (HLL outer speeds) | ~532–533 | 2 min/max per cell |
| λ*, ρ*_L/R, λ*_L/R (if `uct_hlld`) | ~573–581 | ~15 flops per cell |
| chit_L/R, vL/R/vstar (if `uct_hlld`) | ~583–598 | ~20 flops per cell |

Writes: `aLx(k,j,i)`, `aRx(k,j,i)`, `dLx(k,j,i)`, `dRx(k,j,i)`.

### `src/hydro/rsolvers/mhd/hlld.cpp` — `Hydro::RiemannSolver` (line 36)

Already computes in the same `(k,j)` pass (as `private` SIMD-register values):

| `spd` index | Quantity | Line |
|-------------|----------|------|
| `spd[0]` | λ_L = `min(wl[IVX]-cfl, wr[IVX]-cfr)` | 118 |
| `spd[4]` | λ_R = `max(wl[IVX]+cfl, wr[IVX]+cfr)` | 119 |
| `spd[2]` | λ* (Miyoshi & Kusano Eq. 38) | 158 |
| `ulst.d` | ρ*_L (Eq. 43) | 165 |
| `urst.d` | ρ*_R (Eq. 43) | 166 |
| `spd[1]` | λ*_L = λ* − |Bx|/√ρ*_L | 173 |
| `spd[3]` | λ*_R = λ* + |Bx|/√ρ*_R | 174 |

`src/hydro/rsolvers/mhd/lhlld.cpp` has the same structure at the same relative lines.

---

## UCT coefficient consumption

**`src/field/calculate_corner_e.cpp`** — `Field::ComputeCornerE_UCT` (line 253)

Reads `aLx1/aRx1/dLx1/dRx1` (and x2, x3 variants) to compute corner EMFs for CT update.
Arrays declared in `src/field/field.hpp:63–64`. Called after `CalculateFluxes` in the
task graph.

---

## SIMD structure in `hlld.cpp`

```cpp
// Before loop: spd aliases and field pointer setup (loop-invariant)
#pragma omp simd simdlen(SIMD_WIDTH) private(wli,wri,spd,flxi,vf)
for (int i=il; i<=iu; ++i) {
    // ... flux computation ...
    spd[0..4] computed here        // live in SIMD registers

    // existing array writes (model for UCT writes):
    ey(k,j,i) = -flxi[IBY];
    ez(k,j,i) =  flxi[IBZ];
    wct(k,j,i) = GetWeightForCT(...);

    // UCT writes will go here (Phase 2)
    // if (uct_on) { aLx(k,j,i) = ...; aRx(k,j,i) = ...; ... }
}
```

`spd` is `private` in the pragma (per-lane in registers). Writing to `aLx(k,j,i)` inside
the loop follows the same scatter-write pattern as `ey`/`ez`/`wct`. The `uct_on` guard
is loop-invariant; the compiler hoists it.

---

## Why non-HLLD solvers are unaffected

`CalculateHLLWaveSpeed` computes λ_L = `min(vl−cfl, vr−cfr)` using its own
`FastMagnetosonicSpeed` calls — independent of whatever the active flux solver computed.
For `--flux=roe`, the Roe solver uses a different single-sided speed estimate
(`a = 0.5*max(|vl|+cfl, |vr|+cfr)`, `roe_mhd.cpp:188–190`) and does not produce λ_L/λ_R
in the HLL sense. So `CalculateHLLWaveSpeed` is not duplicating work for roe; it must run.

---

## Phase 1 profiling result

Run: `--flux=hlld --ct_method=uct_hlld`, 64×32×32 grid, 16³ meshblocks, 500 cycles.

```
[UCT profile] RiemannSolver: 29.50 s  CalculateHLLWaveSpeed: 18.24 s  ratio: 61.8%
```

`CalculateHLLWaveSpeed` consumes **62%** of `RiemannSolver` wall time — far above the 5% threshold. Phase 2 is well justified.

---

## Profiling instrumentation sketch (Phase 1)

```cpp
// calculate_fluxes.cpp — per-direction block (x1 shown):
#ifdef PROFILE_UCT_WAVESPEED
static double t_rsolver = 0.0, t_hll_ws = 0.0;
double t0 = omp_get_wtime();
#endif
RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);
#ifdef PROFILE_UCT_WAVESPEED
#pragma omp atomic
t_rsolver += omp_get_wtime() - t0;
if (pmb->pfield->uct_on) {
  t0 = omp_get_wtime();
#endif
  CalculateHLLWaveSpeed(k, j, is, ie+1, IVX, b1, wl_, wr_,
                        pmb->pfield->aLx1, pmb->pfield->aRx1,
                        pmb->pfield->dLx1, pmb->pfield->dRx1);
#ifdef PROFILE_UCT_WAVESPEED
  #pragma omp atomic
  t_hll_ws += omp_get_wtime() - t0;
}
#endif
```

Print `t_rsolver` and `t_hll_ws` in the `Hydro` destructor (or `Mesh::Finalize`).
Strip all `#ifdef PROFILE_UCT_WAVESPEED` blocks after Phase 2.
