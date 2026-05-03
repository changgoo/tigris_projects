# FOFC Conservation Violation: Diagnosis and Fix

## Problem

First Order Flux Correction (FOFC) causes mass conservation violation.
This occurs even in pure hydro.

## Diagnosis

### What was checked and found correct

1. **`ApplyFOFC_Hydro`** (`first_order_flux_correction.cpp:445`): Piecewise
   constant reconstruction is correct. Cell-centered primitives are directly
   used as left/right states. Velocity rotation for x2 (IVX->IVY->IVZ->IVX)
   and x3 directions follows standard Athena++ convention. Flux index mapping
   back to global momentum components (IM1/IM2/IM3) is correct.

2. **`SingleStateLLF_Hydro`** (`rsolvers/hydro/llf_single_state.cpp:34`):
   Standard Rusanov/LLF formula is correctly implemented:
   - Wave speed: `a = 0.5 * max(|v_L|+c_L, |v_R|+c_R)`
   - Physical fluxes `fl`, `fr` and conserved state difference `du` correct
   - Final flux: `flx = 0.5*(fl + fr) - a*du`

3. **Within a single MeshBlock**: Conservation is guaranteed because adjacent
   cells share the same flux array entry. When two flagged cells both write to
   a shared face, they compute the same LLF flux (same primitive states from
   `w`), so the last write is harmless.

### Root cause: flux mismatch at MeshBlock / periodic boundaries

The conservation violation occurs because **FOFC can replace the flux at a
boundary face on one side but not the other**, creating different flux values
for the same physical face.

**Mechanism:**

1. Each MeshBlock independently predicts the next-step state (`utest_`) using
   `AddFluxDivergence` (line 71), then flags bad cells via
   `ConservedToPrimitiveTest` (line 115).

2. The FOFC flag depends on the full stencil of fluxes around each cell. Two
   blocks can independently disagree on whether their respective boundary
   cells are bad. Block A may flag cell `ie` while Block B does not flag cell
   `is` (or vice versa).

3. If Block A flags `ie`, it replaces `x1flux(*,k,j,ie+1)` with LLF. Block B
   keeps the original high-order `x1flux(*,k,j,is)`. These represent the
   **same physical face** but now have different values.

4. `AddFluxDivergence` (line 57 of `add_flux_divergence.cpp`) uses each
   block's own flux array:
   ```
   dflx = area(i+1)*x1flux(n,k,j,i+1) - area(i)*x1flux(n,k,j,i)
   ```
   Mass leaving Block A != mass entering Block B -> **conservation violated**.

**Same issue for single-block periodic BCs:** Face `is` and face `ie+1` are
the same physical face but stored as separate array entries. If FOFC replaces
one but not the other, the flux sum no longer telescopes:
```
sum(F[i+1] - F[i]) = F[ie+1] - F[is] != 0
```

**No flux correction exists for this case.** From the task graph
(`time_integrator.cpp:953-964`), for same-level blocks without shear-periodic,
the chain is simply `CALC_HYDFLX -> INT_HYD` with no flux correction step.
`SEND_HYDFLX`/`RECV_HYDFLX` only exist for multilevel (AMR) or
`SHEAR_PERIODIC`, and even then they reconcile AMR fine/coarse mismatches,
not FOFC mismatches.

## Fix options considered

### Option 1: Always flag same-level boundary cells (like `fofc_shear`)

Similar to the existing `fofc_shear` mechanism (lines 118-129) which tags
`is`/`ie` cells at shear boundaries, extend this to all same-level face
neighbors. Both blocks would then apply FOFC and compute identical LLF fluxes.

- Pro: ~20 lines, no communication, guaranteed conservation
- Con: More FOFC applications than necessary at all block boundaries

### Option 2: Skip boundary face replacement (current fix)

When FOFC flags a boundary cell, replace only the interior face fluxes and
leave the boundary face (face `is` or `ie+1`) untouched. Both blocks keep the
original high-order flux at the shared face -> conservation preserved.

- Pro: No extra FOFC, no communication, minimal code change
- Con: If the boundary face flux is the problematic one, FOFC loses
  effectiveness for that face (interior faces still get corrected)

### Option 3: Add boundary communication for FOFC flags

Exchange FOFC flags at block boundaries after computing them but before
applying FOFC. If either side flags a boundary cell, both apply FOFC.

- Pro: Most correct -- FOFC applied only when needed, conservation guaranteed
- Con: Requires new MPI communication, task list reordering, significant
  complexity in the Athena++ boundary infrastructure

## Current fix: Option 2 (skip boundary face replacement)

### Changes in `first_order_flux_correction.cpp`

**`ApplyFOFC_Hydro`** (pure hydro, loop range `is` to `ie`):

Added boolean flags to skip flux replacement at MeshBlock boundary faces:
```cpp
bool at_ix1 = (i == pmb->is);   // left x1 boundary face at is
bool at_ox1 = (i == pmb->ie);   // right x1 boundary face at ie+1
bool at_ix2 = (j == pmb->js);   // (inside f2 block)
bool at_ox2 = (j == pmb->je);
bool at_ix3 = (k == pmb->ks);   // (inside f3 block)
bool at_ox3 = (k == pmb->ke);
```

Each "replace fluxes" block is wrapped:
```cpp
// replace fluxes at i (skip at MeshBlock boundary face)
if (!at_ix1) {
    x1flux(IDN,k,j,i) = flx[IDN];
    ...
}
```

The LLF flux computation is still performed (just not written), so there is
no change to the code flow.

**`ApplyFOFC_MHD`** (MHD builds, loop range `is-1` to `ie+1`):

Uses face-position-based checks because the extended loop means multiple cells
can touch the same boundary face:
```cpp
// Face is is touched by cell is (left face) and cell is-1 (right face)
// Face ie+1 is touched by cell ie (right face) and cell ie+1 (left face)
bool skip_x1_at_i   = (i == pmb->is) || (i == pmb->ie + 1);
bool skip_x1_at_ip1 = (i == pmb->is - 1) || (i == pmb->ie);
```

Same pattern for x2 (`skip_x2_at_j/jp1`) and x3 (`skip_x3_at_k/kp1`).
Each skip guard also covers the EMF (`e3x1_`, `e2x1_`, etc.), CT weight
(`w_x1f`, etc.), and UCT wave speed (`aLx1`, `aRx1`, etc.) at that face.

### Why this preserves conservation

- Both blocks independently compute the same high-order flux at their shared
  boundary face (from identical ghost zone data after exchange).
- FOFC no longer replaces this flux on either side, so both blocks use the
  same value -> the flux telescope sum cancels at the boundary.
- Interior faces of the flagged cell still get LLF replacement, so the cell
  benefits from partial FOFC correction.

### Limitation

If the problematic flux is at the boundary face itself, this fix does not
correct it. The cell gets a "hybrid" update: LLF on interior faces,
high-order on the boundary face. If testing shows this is insufficient,
Option 1 (always flag boundary cells) would be the next step.
