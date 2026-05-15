# UCT-HLLD Wave Speed Deduplication — Design

## Problem

When `ct_method=uct_hlld` or `ct_method=uct_hll` is active, `Hydro::CalculateHLLWaveSpeed`
(`src/hydro/calculate_fluxes.cpp:441`) computes HLL outer wave speeds (λ_L, λ_R) and —
for `uct_hlld` — HLLD intermediate states (λ*, ρ*_L, ρ*_R, λ*_L, λ*_R). When
`--flux=hlld` or `--flux=lhlld`, the Riemann solver already computes all of these as
`spd[0..4]`. The duplication includes two `FastMagnetosonicSpeed` evaluations and a full
Roe-averaged-state computation per cell per direction.

## Goal

1. **Quantify** the overhead with instrumented timing (Phase 1).
2. **Eliminate** the duplication for `--flux=hlld` and `--flux=lhlld` by inlining UCT
   coefficient computation into the solvers' existing SIMD loops (Phase 2).
3. **Preserve** `CalculateHLLWaveSpeed` as the fallback for all other flux solvers.

## Flux / CT-method support matrix

| `--flux`             | `uct_hll`                    | `uct_hlld`                   |
|----------------------|------------------------------|------------------------------|
| `hlld`               | inline in solver             | inline in solver             |
| `lhlld`              | inline in solver             | inline in solver             |
| `roe`, `hlle`, `llf` | `CalculateHLLWaveSpeed` kept | `CalculateHLLWaveSpeed` kept |

`CalculateHLLWaveSpeed` is not deleted — it remains active for all non-HLLD solver paths.

---

## Phase 1 — Profiling

Add `#ifdef PROFILE_UCT_WAVESPEED` guards in `src/hydro/calculate_fluxes.cpp` around
the `RiemannSolver` and `CalculateHLLWaveSpeed` calls for all three directions. Use
`omp_get_wtime()` with `#pragma omp atomic` accumulation into file-scope `static double`
accumulators. Print totals when the `Hydro` destructor fires (or in `Mesh::Finalize()`).

The block is defined by a single `#define PROFILE_UCT_WAVESPEED` at the top of the file
(or passed as `-DPROFILE_UCT_WAVESPEED` at compile time). No runtime cost when undefined.
After Phase 2 validation, remove the `#ifdef` blocks entirely.

**Decision threshold:** if `CalculateHLLWaveSpeed` accounts for ≥ 5% of `RiemannSolver`
wall time in a representative 3D MHD turbulence run with `--flux=hlld --ct_method=uct_hlld`,
Phase 2 is justified.

---

## Phase 2 — Inline UCT coefficients into HLLD solvers

### Interface choice

**Option A (chosen):** compute UCT coefficients inside the solvers' existing
`#pragma omp simd` loop, writing `aLx`, `aRx`, `dLx`, `dRx` directly to the Field arrays.
`spd[0..4]` are already live in SIMD registers at the write point — no extra memory round-trip.

**Option B (rejected):** store `spd[0..4]` in scratch arrays; compute coefficients in a
second pass. Adds 5 stores + 5 loads with no benefit over Option A.

The `uct_on` and `uct_hlld` booleans are loop-invariant; the compiler hoists the branch,
preserving vectorization width. This follows the existing pattern of `ey`, `ez`, and `wct`
writes inside the same SIMD loop.

### Changes required

**`src/hydro/rsolvers/mhd/hlld.cpp` and `lhlld.cpp`**

Before the SIMD loop: set `AthenaArray<Real>&` aliases for the direction-appropriate
Field arrays based on `ivx` (loop-invariant — no SIMD impact):

```cpp
Field *pfield = pmy_block->pfield;
AthenaArray<Real> &aLx = (ivx==IVX) ? pfield->aLx1 :
                         (ivx==IVY) ? pfield->aLx2 : pfield->aLx3;
// ... similarly for aRx, dLx, dRx
```

Inside the SIMD loop, after `spd[0..4]` are computed, add an `if (uct_on)` block that
computes MDZ '21 Eq. 32 (`uct_hll` path) or Eq. 44–45 (`uct_hlld` path) and writes the
four coefficient arrays. The arithmetic is identical to what `CalculateHLLWaveSpeed`
already does in the `uct_hlld` branch (`calculate_fluxes.cpp:557–598`).

**`src/hydro/calculate_fluxes.cpp`**

Replace the unconditional `if (pmb->pfield->uct_on) CalculateHLLWaveSpeed(...)` with a
guard that skips the call when the active solver is hlld or lhlld (solver-type flag TBD
during implementation — see notes below). For all other solvers, behaviour is unchanged.

Remove `#ifdef PROFILE_UCT_WAVESPEED` blocks after Phase 2 is validated.

### Solver-type flag (implementation note)

A lightweight mechanism is needed to tell `CalculateFluxes` whether the active solver
already wrote the UCT coefficients. Options to resolve during implementation:

- Add `bool pfield->solver_handles_uct` set during `Field::Init` based on `rsolver_method`.
- Or check `phydro->rsolver_method` (if that member exists) directly in `CalculateFluxes`.

The writing-plans step will pin down which is cleaner given the existing init flow.

---

## Files changed

| File | Change |
|------|--------|
| `src/hydro/calculate_fluxes.cpp` | Phase 1 profiling guards; Phase 2 skip for hlld/lhlld |
| `src/hydro/rsolvers/mhd/hlld.cpp` | Phase 2 inline UCT coefficients in SIMD loop |
| `src/hydro/rsolvers/mhd/lhlld.cpp` | Phase 2 inline UCT coefficients in SIMD loop |
| `src/hydro/hydro.hpp` or `src/field/field.hpp` | Phase 2 solver-type flag (if needed) |

## References

- Mignone & Del-Zanna (2021) — MDZ '21 — Eq. 32 (UCT-HLL), Eq. 44–45 (UCT-HLLD)
- Miyoshi & Kusano (2005) — Eq. 38, 43, 51 (HLLD intermediate states)
