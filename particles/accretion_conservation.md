# Particle Accretion Mass Conservation

## Background

Sink particle accretion uses a control volume (CV) reset scheme: a 3x3x3 cell region
centered on the particle is extrapolated inward from a 5x5x5 padding shell, and the
mass difference is credited to the particle. When the CV overlaps ghost zones, a
**ghost particle** is created on the neighbor MeshBlock to ensure the same physical
cells are consistently reset on both blocks.

### How accretion works (single particle, no boundary overlap)

1. `CopyControlVolumeFromMesh(cons, ..., buf0)` — save 3x3x3 cells before reset
2. `ResetControlVolume(cons, ..., buf1)`:
   - `ConservedToPrimitive`: read 5x5x5 from `cons` (3x3x3 + 1-cell padding) into buffer
   - `ExtrapolateInward`: fill interior from outermost padding shell
   - `PrimitiveToConserved`: write extrapolated 3x3x3 back to `cons`
3. `CopyControlVolumeFromMesh(cons, ..., buf1)` — save 3x3x3 cells after reset
4. `delta_reset = sum(buf1 - buf0) * dvol` over all 27 cells
5. `particle.mass += -delta_reset(IDN)`

Conservation holds trivially: `delta_grid + delta_particle = 0`.

### How it works at periodic boundaries (ghost particle scheme)

- Active particle on Block A modifies `cons` for all 27 CV cells (active + ghost on A)
- Ghost particle on Block B independently resets the same 27 physical cells
- After boundary exchange: A's ghost cells get B's active values, B's ghost gets A's active
- For conservation, ghost particle must produce **identical** extrapolation to active particle

This works for pure periodic boundaries because the ghost zone data is a direct copy
of the neighbor's active cells — both particles see identical 5x5x5 padding data and
produce identical extrapolation results.

---

## Issue 1: Diagnostic Indexing Bug (fixed)

**File:** `src/particles/accretion.cpp`, `Accrete()` function

The diagnostic split of `delta_reset` into active vs ghost cell contributions used
wrong cell indices. Buffer index `i` maps to mesh cell `ic-rctrl+i`, but
`CheckInMeshBlock` was called with `ic+i` — off by `rctrl` in each dimension.

**Fix:** Changed `CheckInMeshBlock(i+ic, j+jc, k+kc, 0)` to
`CheckInMeshBlock(i+ic-rctrl, j+jc-rctrl, k+kc-rctrl, 0)`.

This was purely a diagnostic bug — the total `delta_reset` was always correct.

---

## Issue 2: Shear-Periodic Boundary Conservation Violation (fixed)

### Root cause

At shear-periodic (X1) boundaries, ghost zone data is **interpolated** due to the
non-integer Y-shift. This means the 5x5x5 padding data seen by the ghost particle on
the shear neighbor differs from what the active particle sees. Different padding →
different extrapolation → different `delta_reset` values for the overlapping cells.

The old scheme credited the active particle with the **full** `delta_reset` (active +
ghost cells). But after boundary exchange, the ghost cells on the active particle's
block are overwritten with the neighbor's active cell values — which reflect the
**ghost particle's** extrapolation, not the active particle's. This mismatch breaks
mass conservation.

### When it triggers

Whenever a particle's CV overlaps with a shear-periodic boundary. The violation
magnitude depends on how different the extrapolation contexts are, which varies with
the local gas structure near the shear boundary.

### Fix: Ghost delta return communication

Split `delta_reset` into active-cell and ghost-cell contributions. Each particle type
is handled differently:

- **Active particle**: credits itself only with `delta_reset_active` (mass change in
  its own block's active cells). These are the cells whose values persist after
  boundary exchange.
- **Ghost particle**: modifies the grid on the neighbor block (this always happens).
  Its `delta_reset_active` (mass change in the neighbor's active cells) is stored and
  communicated back to the active particle's block.

After the accretion loop, `ExchangeGhostAccretionDelta()` uses `MPI_Allgatherv` to
broadcast all ghost deltas, and each block applies corrections to its matching active
particles.

**Total particle mass credit** = `delta_active(own block)` + `sum(delta_active from ghost particles on neighbors)`

This equals the total grid mass change (only active cell modifications persist), so
conservation holds regardless of how different the extrapolation contexts are.

### Files modified

| File | Changes |
|------|---------|
| `src/particles/accretion.hpp` | `Accrete()` signature: added `delta_vars_active` output parameter |
| `src/particles/accretion.cpp` | `Accrete()`: computes `delta_reset_active` separately; fixed diagnostic indexing |
| `src/particles/complex_particles.cpp` | `AccreteFromSingleParticle()`: active particles use `delta_reset_active`; ghost particles store delta for return. Added `ExchangeGhostAccretionDelta()`. |
| `src/particles/particles.hpp` | Added `ExchangeGhostAccretionDelta()` declaration and ghost delta storage vectors |

### Conservation argument

After all accretion and boundary exchange:

- Block A active cells: modified by active particle. Delta = `delta_active(A)`.
- Block B active cells: modified by ghost particle. Delta = `delta_active(ghost on B)`.
- Ghost cells on both blocks: overwritten by boundary exchange (no net effect).
- Particle mass change: `delta_active(A)` + `delta_active(ghost on B)` (via exchange).
- Grid mass change: `delta_active(A)` + `delta_active(ghost on B)`.
- `grid_change + particle_change = 0`. QED.

---

## Issue 3: Cross-Block Overlapping Particles (investigated, not observed)

When two particles on different blocks have overlapping CVs near their shared boundary,
the sequential processing order could cause inconsistency: the ghost particle sees
already-modified `cons` data from the other particle's accretion. This was investigated
with overlap detection diagnostics but never triggered in test runs — it requires two
accreting particles within ~2 cells of each other straddling a block boundary, which is
extremely rare given the particle merger logic.

---

## Issue 4: `pid=NEW` Ambiguity in Ghost Delta Exchange (fixed)

### Root cause

During the cycle when new particles are created (before `ProcessNewParticles` assigns
unique IDs), multiple particles can have `pid=NEW=-1`. The original
`ExchangeGhostAccretionDelta()` matched ghost corrections to active particles by pid
only. When multiple NEW particles existed simultaneously, **every** ghost correction
was applied to **every** active particle with `pid=NEW`, instead of only to the
correct matching particle.

### When it triggers

Whenever two or more particles are created in the same cycle near boundaries that
produce ghost copies. Observed in practice at ncycle=98 with two simultaneous NEW
particles — one at a pure periodic Z boundary and one at a shear-periodic X boundary.

### Fix: Position-based matching for `pid=NEW`

The ghost delta exchange buffer was extended from `[pid, delta...]` to
`[pid, xp, yp, zp, delta...]`. For particles with unique pids (>= 0), matching
uses pid only (fast path). For `pid=NEW`, matching additionally requires position
proximity after unwrapping periodic and shear-periodic offsets:

- **Pure periodic (X2, X3)**: ghost and active have identical `(xp, yp, zp)`.
  Direct float comparison within half-cell tolerance.
- **Shear-periodic (X1)**: ghost position differs by `(±Lx, ∓deltay, 0)` where
  `deltay = fmod(qomL * time, Ly)`. Position difference is corrected for the shear
  offset before comparison.
- **Combined (shear + periodic Z)**: both corrections applied sequentially.

Tolerance: half cell width (`0.5 * dx`), which is much smaller than the minimum
separation between distinct particle creation sites.

---

## Known Limitations

- **MPI collective in per-MeshBlock task**: `ExchangeGhostAccretionDelta()` uses
  `MPI_Allgatherv` inside the INTERACT task, which runs per-MeshBlock. This requires
  all ranks to have the same number of MeshBlocks. This is always true for uniform
  mesh (no AMR), which is the standard TIGRIS configuration. Same pattern as
  `MassReturn::CollectParticlesInfo()`.
