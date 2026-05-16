# UCT Wave Speed Deduplication — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate redundant wave-speed computation in `CalculateHLLWaveSpeed` when `--flux=hlld` or `--flux=lhlld` by inlining UCT coefficient writes into those solvers' existing SIMD loops.

**Architecture:** Phase 1 adds compile-time profiling guards (`#ifdef PROFILE_UCT_WAVESPEED`) to quantify overhead. Phase 2 inlines UCT coefficients into the HLLD/LHLLD solver SIMD loops using `spd[0..4]`/`ulst.d`/`urst.d` already live in registers, and adds a `solver_handles_uct` flag so `CalculateFluxes` skips the redundant `CalculateHLLWaveSpeed` call. `CalculateHLLWaveSpeed` is preserved unchanged for all other solvers.

**Tech Stack:** C++11, Athena++/TIGRIS; OpenMP SIMD (`#pragma omp simd`); `omp_get_wtime()` for profiling.

**Worktree:** All changes go in the `uct-hlld-refactor` worktree of `/home/changgoo/tigris`.

---

## File map

| File | Change |
|------|--------|
| `src/hydro/calculate_fluxes.cpp` | Phase 1 profiling guards; Phase 2 conditional `CalculateHLLWaveSpeed` skip |
| `src/hydro/hydro.cpp` | Phase 1 timing print in destructor |
| `src/hydro/rsolvers/mhd/hlld.cpp` | Phase 2 inline UCT writes in SIMD loop |
| `src/hydro/rsolvers/mhd/lhlld.cpp` | Phase 2 inline UCT writes in SIMD loop |
| `src/field/field.hpp` | Add `bool solver_handles_uct` (line 67) |
| `src/field/field.cpp` | Set `solver_handles_uct` in constructor |

---

## Task 1 — Phase 1: Add profiling guards to `calculate_fluxes.cpp`

**Files:**
- Modify: `src/hydro/calculate_fluxes.cpp` (three RiemannSolver call sites + top of function)
- Modify: `src/hydro/hydro.cpp` (destructor)

- [ ] **Step 1: Add `#include <omp.h>` guard near top of `calculate_fluxes.cpp`**

After existing includes, add:
```cpp
#ifdef PROFILE_UCT_WAVESPEED
#include <omp.h>
#endif
```

- [ ] **Step 2: Add static accumulators at the top of `Hydro::CalculateFluxes`**

Right after the opening brace of `Hydro::CalculateFluxes(...)`, add:
```cpp
#ifdef PROFILE_UCT_WAVESPEED
  static double t_rsolver = 0.0, t_hll_ws = 0.0;
  double t0;
#endif
```

- [ ] **Step 3: Wrap the x1 `RiemannSolver` + `CalculateHLLWaveSpeed` calls (lines 94–99)**

Replace:
```cpp
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);

      if (pmb->pfield->uct_on)
        CalculateHLLWaveSpeed(k, j, is, ie+1, IVX, b1, wl_, wr_,
                              pmb->pfield->aLx1, pmb->pfield->aRx1,
                              pmb->pfield->dLx1, pmb->pfield->dRx1);
```
With:
```cpp
#ifdef PROFILE_UCT_WAVESPEED
      t0 = omp_get_wtime();
#endif
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);
#ifdef PROFILE_UCT_WAVESPEED
      #pragma omp atomic
      t_rsolver += omp_get_wtime() - t0;
#endif
      if (pmb->pfield->uct_on) {
#ifdef PROFILE_UCT_WAVESPEED
        t0 = omp_get_wtime();
#endif
        CalculateHLLWaveSpeed(k, j, is, ie+1, IVX, b1, wl_, wr_,
                              pmb->pfield->aLx1, pmb->pfield->aRx1,
                              pmb->pfield->dLx1, pmb->pfield->dRx1);
#ifdef PROFILE_UCT_WAVESPEED
        #pragma omp atomic
        t_hll_ws += omp_get_wtime() - t0;
#endif
      }
```

- [ ] **Step 4: Apply identical wrapping to x2 (`IVY`, `b2`, `aLx2`) and x3 (`IVZ`, `b3`, `aLx3`)**

Same pattern as Step 3. For x2 the call is inside `for (int j=js; j<=je+1; ...)`. For x3 it is inside `for (int k=ks; k<=ke+1; ...)`. Update the array names and direction indices accordingly.

- [ ] **Step 5: Print totals in `Hydro::~Hydro()` in `src/hydro/hydro.cpp`**

At the end of the destructor body add:
```cpp
#ifdef PROFILE_UCT_WAVESPEED
  if (t_rsolver > 0.0) {
    std::printf("[UCT profile] RiemannSolver: %.4f s  CalculateHLLWaveSpeed: %.4f s"
                "  ratio: %.1f%%\n",
                t_rsolver, t_hll_ws, 100.0*t_hll_ws/t_rsolver);
  }
#endif
```
`t_rsolver` and `t_hll_ws` are the file-scope statics from `calculate_fluxes.cpp`; add `extern` declarations at the top of `hydro.cpp` inside `#ifdef PROFILE_UCT_WAVESPEED`:
```cpp
#ifdef PROFILE_UCT_WAVESPEED
extern double t_rsolver, t_hll_ws; // defined in calculate_fluxes.cpp
#endif
```

- [ ] **Step 6: Commit**
```bash
git add src/hydro/calculate_fluxes.cpp src/hydro/hydro.cpp
git commit -m "Add PROFILE_UCT_WAVESPEED timing instrumentation"
```

---

## Task 2 — Phase 1: Build and measure

Use the existing UCT-HLL linear wave regression test as the representative run. The test file is at `tst/regression/scripts/tests/mhd/mhd_linwave_ucthll.py` — it configures with `flux=hlld`, builds, and runs 3D linear wave convergence with `hydro/ct_method=uct_hlld`.

- [ ] **Step 1: Configure and build with profiling macro**

From the TIGRIS source root (or worktree root):
```bash
python configure.py -b --prob=linear_wave --coord=cartesian --flux=hlld \
    --cxxflags="-DPROFILE_UCT_WAVESPEED"
make -j8
```
Expected: clean build, no new warnings.

- [ ] **Step 2: Run the UCT linear wave test**
```bash
cd tst/regression
python run_tests.py scripts/tests/mhd/mhd_linwave_ucthll.py
```
Expected: test passes AND stdout includes a line like:
```
[UCT profile] RiemannSolver: 12.3 s  CalculateHLLWaveSpeed: 0.9 s  ratio: 7.3%
```

- [ ] **Step 3: Decision gate**

- If ratio **< 5%**: record the numbers in `uct_hlld_refactor/context.md` and stop — Phase 2 is not cost-justified.
- If ratio **≥ 5%**: proceed to Task 3.

---

## Task 3 — Add `solver_handles_uct` flag to `Field`

**Files:**
- Modify: `src/field/field.hpp:67`
- Modify: `src/field/field.cpp` (Field constructor)

- [ ] **Step 1: Add member to `field.hpp`**

Change line 67 from:
```cpp
  bool uct_on, uct_hlld;
```
To:
```cpp
  bool uct_on, uct_hlld, solver_handles_uct;
```

- [ ] **Step 2: Set flag in the `Field` constructor in `field.cpp`**

In `Field::Field(...)`, immediately after the lines that set `uct_on` and `uct_hlld`, add:
```cpp
  std::string rsolver = pin->GetOrAddString("hydro", "rsolver", "hlld");
  solver_handles_uct = uct_on && (rsolver == "hlld" || rsolver == "lhlld");
```
Verify the pin key by searching: `grep -n "rsolver\|GetString" src/hydro/hydro.cpp` and use the exact key found there.

- [ ] **Step 3: Commit**
```bash
git add src/field/field.hpp src/field/field.cpp
git commit -m "Add solver_handles_uct flag to Field"
```

---

## Task 4 — Inline UCT coefficients in `hlld.cpp`

`spd[0]`=λ_L, `spd[4]`=λ_R, `spd[2]`=λ*, `ulst.d`=ρ*_L, `urst.d`=ρ*_R, `spd[1]`=λ*_L, `spd[3]`=λ*_R are all live after line 171.

**Files:**
- Modify: `src/hydro/rsolvers/mhd/hlld.cpp` (before SIMD loop ~line 53, and after line 171 inside loop)

- [ ] **Step 1: Add `field.hpp` include**

Near the top of `hlld.cpp`, after the existing Athena++ headers:
```cpp
#include "../../field/field.hpp"
```

- [ ] **Step 2: Add direction-aware array aliases before `#pragma omp simd` (~line 55)**

After `Real dt = pmy_block->pmy_mesh->dt;` and before `#pragma omp simd`, insert:
```cpp
#if MAGNETIC_FIELDS_ENABLED
  Field *pfield = pmy_block->pfield;
  const bool uct_on   = pfield->uct_on;
  const bool uct_hlld_flag = pfield->uct_hlld;
  AthenaArray<Real> &aLx = (ivx==IVX) ? pfield->aLx1 :
                            (ivx==IVY) ? pfield->aLx2 : pfield->aLx3;
  AthenaArray<Real> &aRx = (ivx==IVX) ? pfield->aRx1 :
                            (ivx==IVY) ? pfield->aRx2 : pfield->aRx3;
  AthenaArray<Real> &dLx = (ivx==IVX) ? pfield->dLx1 :
                            (ivx==IVY) ? pfield->dLx2 : pfield->dLx3;
  AthenaArray<Real> &dRx = (ivx==IVX) ? pfield->dRx1 :
                            (ivx==IVY) ? pfield->dRx2 : pfield->dRx3;
#endif
```

- [ ] **Step 3: Insert UCT coefficient block after `spd[3]` (line 171) inside the SIMD loop**

```cpp
    // UCT coefficients — MDZ '21 Eq. 32 (uct_hll) overwritten by Eq. 44 (uct_hlld)
    if (uct_on) {
      Real lambda_L = spd[0];
      Real lambda_R = spd[4];
      constexpr Real eps = 1e-9;
      Real ap = lambda_R > 0.0 ? lambda_R : 0.0;
      Real am = lambda_L < 0.0 ? lambda_L : 0.0;
      if ((ap == 0.0) && (am == 0.0)) { am = lambda_R; ap = lambda_L; }
      Real asum = ap + std::abs(am);
      if (ap/asum < eps)                   ap = 0.0;
      else if (std::abs(am)/asum < eps)    am = 0.0;
      asum = ap + std::abs(am);
      // MDZ '21 Eq. 32
      aLx(k,j,i) = ap               / asum;
      aRx(k,j,i) = std::abs(am)     / asum;
      dLx(k,j,i) = ap*std::abs(am)  / asum;
      dRx(k,j,i) = ap*std::abs(am)  / asum;

      if (uct_hlld_flag) {
        Real lstar   = spd[2];
        Real lstar_L = spd[1];
        Real lstar_R = spd[3];
        Real vL    = (lstar_L + lambda_L) / (std::abs(lstar_L) + std::abs(lambda_L));
        Real vR    = (lstar_R + lambda_R) / (std::abs(lstar_R) + std::abs(lambda_R));
        Real vstar = (lstar_R + lstar_L)  / (std::abs(lstar_R) + std::abs(lstar_L));
        Real chit_L = (wli[IVX]-lstar)*(lambda_L-lstar)/(lstar_L+lambda_L-2.0*lstar);
        Real chit_R = (wri[IVX]-lstar)*(lambda_R-lstar)/(lstar_R+lambda_R-2.0*lstar);
        if (std::abs(lstar_R - lstar_L) < eps*std::abs(lambda_R - lambda_L)) vstar = 0.0;
        // MDZ '21 Eq. 44
        aLx(k,j,i) = (1.0 + vstar)/2.0;
        aRx(k,j,i) = (1.0 - vstar)/2.0;
        dLx(k,j,i) = 0.5*(vL-vstar)*chit_L + 0.5*(std::abs(lstar_L) - vstar*lstar_L);
        dRx(k,j,i) = 0.5*(vR-vstar)*chit_R + 0.5*(std::abs(lstar_R) - vstar*lstar_R);
      }
    }
```

- [ ] **Step 4: Commit**
```bash
git add src/hydro/rsolvers/mhd/hlld.cpp
git commit -m "Inline UCT coefficient writes into hlld SIMD loop"
```

---

## Task 5 — Inline UCT coefficients in `lhlld.cpp`

`lhlld.cpp` has the identical `spd` layout at the same relative lines (spd[0]/spd[4] ~116–117, spd[2] ~155, `ulst.d`/`urst.d` ~162–163, spd[1]/spd[3] ~170–171). The LHLLD shock-detector `th` affects only `spd[2]`; `spd[0]`, `spd[4]`, `spd[1]`, `spd[3]` use the same expressions as HLLD.

**Files:**
- Modify: `src/hydro/rsolvers/mhd/lhlld.cpp` (same two insertion points as Task 4)

- [ ] **Step 1: Add `field.hpp` include** (same as Task 4 Step 1)

- [ ] **Step 2: Add direction-aware array aliases before `#pragma omp simd` (~line 53)**

After `CalculateVelocityDifferences(k, j, il, iu, ivx, dvn, dvt);` and before `#pragma omp simd`, insert the identical block from Task 4 Step 2.

- [ ] **Step 3: Insert UCT coefficient block after `spd[3]` (line 171)**

Insert the identical block from Task 4 Step 3.

- [ ] **Step 4: Commit**
```bash
git add src/hydro/rsolvers/mhd/lhlld.cpp
git commit -m "Inline UCT coefficient writes into lhlld SIMD loop"
```

---

## Task 6 — Skip `CalculateHLLWaveSpeed` for hlld/lhlld in `calculate_fluxes.cpp`

**Files:**
- Modify: `src/hydro/calculate_fluxes.cpp` (three `uct_on` guard sites)

- [ ] **Step 1: Update x1 guard (line 96)**

Change:
```cpp
      if (pmb->pfield->uct_on)
```
To:
```cpp
      if (pmb->pfield->uct_on && !pmb->pfield->solver_handles_uct)
```

- [ ] **Step 2: Apply the same one-line change to the x2 and x3 guards**

Same edit — `!pmb->pfield->solver_handles_uct` added to each `if (pmb->pfield->uct_on)` that guards a `CalculateHLLWaveSpeed` call.

- [ ] **Step 3: Commit**
```bash
git add src/hydro/calculate_fluxes.cpp
git commit -m "Skip CalculateHLLWaveSpeed for hlld/lhlld when UCT is inline"
```

---

## Task 7 — Build and regression test

- [ ] **Step 1: Clean build (without profiling flag)**
```bash
make clean && make -j8
```
Expected: compiles with no new warnings or errors.

- [ ] **Step 2: Run MHD regression tests**
```bash
python tst/regression/run_tests.py --config=<config> mhd_linwave mhd_blast mhd_rotor
```
Expected: all PASS.

- [ ] **Step 3: Run a UCT-specific test**
```bash
python tst/regression/run_tests.py --config=<config> field_loop
```
Expected: PASS. If a dedicated `ct_method=uct_hlld` test exists, run it too.

- [ ] **Step 4: Verify numerical equivalence**

Run the same problem twice — once with the old code path (`solver_handles_uct=false` forced) and once with the new — and diff the outputs. Differences should be exactly zero (same floating-point operations).

- [ ] **Step 5: Commit if all pass**
```bash
git commit -m "UCT wave speed deduplication: validated against regression suite"
```

---

## Task 8 — Remove profiling guards

- [ ] **Step 1: Remove all `#ifdef PROFILE_UCT_WAVESPEED` blocks from `calculate_fluxes.cpp`**

Remove the three timing wrappers (Steps 3–4 of Task 1), the static accumulator block (Step 2), and the `#include <omp.h>` guard (Step 1).

- [ ] **Step 2: Remove timing print and `extern` declarations from `hydro.cpp`**

Remove the blocks added in Task 1 Step 5.

- [ ] **Step 3: Build**
```bash
make -j8
```
Expected: clean build.

- [ ] **Step 4: Final commit**
```bash
git add src/hydro/calculate_fluxes.cpp src/hydro/hydro.cpp
git commit -m "Remove PROFILE_UCT_WAVESPEED instrumentation"
```
