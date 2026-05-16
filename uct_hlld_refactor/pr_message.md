# Inline UCT wave-speed coefficients into HLLD/LHLLD SIMD loops

## Summary

- **Problem:** When `ct_method=uct_hll` or `ct_method=uct_hlld` is active with `--flux=hlld` or `--flux=lhlld`, `Hydro::CalculateHLLWaveSpeed` recomputes HLL outer wave speeds (λ_L, λ_R) and HLLD intermediate states (λ\*, ρ\*_L/R, λ\*_L/R) that the Riemann solver already computed as `spd[0..4]` in the same SIMD pass. Profiling showed this redundant call consumed **62%** of `RiemannSolver` wall time on a 64×32×32 test (29.5 s vs 18.2 s).
- **Fix:** Write the UCT coefficients (`aLx`, `aRx`, `dLx`, `dRx`) directly inside the existing `#pragma omp simd` loop in `hlld.cpp` and `lhlld.cpp`, using `spd[0..4]` already live in SIMD registers (MDZ '21 Eq. 32 / Eq. 44–45). A new `Field::solver_handles_uct` compile-time flag (set via `RIEMANN_SOLVER` macro) skips the redundant `CalculateHLLWaveSpeed` call in `CalculateFluxes`. The function is preserved unchanged for all other solvers (`hlle`, `roe`, etc.).
- **Measured speedup:** `CalculateHLLWaveSpeed` reduced from 18.2 s → 0.0 s for `--flux=hlld`; ~29% reduction in combined wave-speed work per cycle.

## Validation

- `mhd_linwave_ucthll` regression test: **PASS** (hlld + uct_hlld, 3D linear wave convergence with SMR)
- Timing checks (serial, 32×16×16, 100 cycles):
  - `hlld`: CalculateHLLWaveSpeed ratio **0.0%** ✓
  - `lhlld`: CalculateHLLWaveSpeed ratio **0.0%** ✓
  - `hlle`: CalculateHLLWaveSpeed ratio **102%** ✓ (fallback still active)
- Numerical equivalence vs forced-fallback path:
  - `hlld`: **bit-for-bit identical** L1 errors ✓
  - `lhlld`: L1 errors identical; 1-ULP difference in one max-norm value (expected from SIMD re-ordering) ✓

## Files changed

| File | Change |
|------|--------|
| `src/hydro/rsolvers/mhd/hlld.cpp` | Inline UCT coefficient writes in SIMD loop |
| `src/hydro/rsolvers/mhd/lhlld.cpp` | Inline UCT coefficient writes in SIMD loop |
| `src/field/field.hpp` | Add `solver_handles_uct` member |
| `src/field/field.cpp` | Set `solver_handles_uct` from `RIEMANN_SOLVER` macro |
| `src/hydro/calculate_fluxes.cpp` | Skip `CalculateHLLWaveSpeed` when `solver_handles_uct` |
| `src/hydro/hydro.hpp` / `hydro.cpp` | Add `~Hydro()` destructor declaration |

References: Mignone & Del-Zanna (2021) Eq. 32, 44–45; Miyoshi & Kusano (2005) Eq. 38, 43.
