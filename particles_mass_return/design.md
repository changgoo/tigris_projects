# Mass Return P2P — Design

**Date**: 2026-05-15
**Related issue**: [#269](https://github.com/PrincetonUniversity/tigris/issues/269)
**In-progress PR**: [#285](https://github.com/PrincetonUniversity/tigris/pull/285)

---

## Scope

Extend the mesh-level P2P mass-return path in PR #285 so that `r_return` can be
arbitrarily large instead of limited to the immediate-neighbor exchange range. Add a
regression test that validates conservation for multiple `r_return` values and MPI
decompositions.

Two problems are solved independently:

1. **General geometric routing** — send particle records only to ranks whose MeshBlocks
   intersect the return sphere, using the global block map already available on every
   rank.
2. **Test problem** — a controlled single-particle setup with a known deposition pattern
   and a mass-conservation assertion.

---

## Current State (PR #285)

PR #285 replaces `MPI_Allgatherv` with a mesh-level P2P path:

- `MassReturn::ReturnMassFromParticles(Mesh*)` is called once per cycle from
  `Mesh::UserWorkInLoop()` across all three problem generators.
- Particle records are sent to **all immediate-neighbor ranks** via
  `ExchangeRealBuffers`, then deposited, then deposited-totals are returned to source
  ranks and applied to owner particles.
- `r_return == 0` and `r_return > immediate_neighbor_range` are rejected with fatal
  errors.

The PR is functionally correct within its stated constraints. The two fatal errors become
the extension points for this design.

---

## Design 1: General Routing via Global Block Map

### Data structure

Build a flat block-bounds array once per `ReturnMassFromParticles(Mesh*)` call:

```cpp
struct BlockBounds {
  Real x1min, x1max, x2min, x2max, x3min, x3max;
  int rank;
};
std::vector<BlockBounds> block_map; // size = pm->nbtotal
```

Every rank has `pm->loclist[i]` (logical location: level, lx1, lx2, lx3) and
`pm->ranklist[i]`. Physical bounds for block `i` at level `l` are computed from
`pm->mesh_size` and `pm->nrbx1/2/3`:

```
dx1 = (mesh_size.x1max - mesh_size.x1min) / (nrbx1 * 2^(l - root_level))
x1min_i = mesh_size.x1min + lx1 * dx1
x1max_i = x1min_i + dx1
```

No MPI communication is needed to build this map.

### Sphere-box intersection

For each `ParticleData pd` in the local particle list (originals and periodic/shear-periodic
copies alike), test against every entry in `block_map`:

```
d² = sum over active dims of max(0, x_block_min - xp)² + max(0, xp - x_block_max)²
intersects if d² <= r_return²
```

This is O(1) per (particle, block) pair and O(n_particles × n_blocks) total — acceptable
because `n_blocks` is small in typical TIGRIS runs.

### Per-rank send buffers

Replace the current "same buffer to all neighbors" with per-rank packing:

```
rank_to_particles: map<int, vector<ParticleData>>

for each pd in local_particles (including copies):
  seen_ranks: set<int>
  for each block in block_map:
    if sphere(pd, r_return) intersects block.bounds:
      if block.rank not in seen_ranks:
        rank_to_particles[block.rank].push_back(pd)
        seen_ranks.insert(block.rank)

affected_ranks = sorted keys of rank_to_particles
send_buffers[n] = pack(rank_to_particles[affected_ranks[n]])
```

`ExchangeRealBuffers` already accepts per-rank buffers — only the population step
changes. The totals return exchange reuses the same `affected_ranks` list.

### Large `r_return` guard

Before the routing loop, check whether `r_return` exceeds half the domain length in any
periodic direction. If so, every block is affected regardless of particle position. In
that case, set `affected_ranks` to all ranks and send the full local particle list to
each rank — no geometric loop needed. This avoids an infinite-copy-generation problem and
is correct because at that radius every block receives deposits.

```cpp
bool full_domain = false;
if (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic)
  full_domain |= (r_return >= 0.5 * (pm->mesh_size.x1max - pm->mesh_size.x1min));
// repeat for x2, x3
```

### Periodic and shear-periodic copies

`CollectParticlesInfo()` already generates periodic position copies (`icopy=true`, same
`pid`). Routing uses these copies directly: a copy at the shifted position naturally
intersects blocks near the periodic image of the return sphere. Shear-periodic copies
carry the correct y-shift. No new copy-generation logic is needed.

The deposited-totals map is keyed by `pid`, so deposits from copies and originals
accumulate under the same owner particle. Conservation is preserved.

### Removal of the fatal error

Remove the `r_return > p2p_range` fatal error. Remove `ImmediateNeighborReturnRange`.
Replace `AddNeighborRanks` with the block-map-based affected-rank computation.

---

## Design 2: Test Problem

**This test is a separate PR branched from `tigris-master`, independent of the P2P
refactoring in PR #285.** It checks physics correctness only — no assumptions about which
communication path is active. The same test suite runs unchanged against the old
`MPI_Allgatherv` implementation and the new P2P implementation, making it suitable for
direct performance comparison between the two.

### Goals

- Validate mass conservation (particle mass loss == grid mass gain) to roundoff.
- Validate correct spatial deposition pattern for finite `r_return`.
- Exercise same-rank multi-MeshBlock boundaries and cross-rank boundaries.
- Cover three radius regimes: sub-block, multi-block, and full-domain.
- Serve as a timing baseline to compare `Allgatherv` vs. P2P performance at scale.

### Setup

Use a new pgen (or extend `particle_complex.cpp`) with:

- A single sink particle placed at the domain center with known initial mass `M`.
- Particles are flagged to return all mass in a single cycle (or over a short schedule).
- No accretion, no feedback — only mass return.
- Domain: periodic in all directions (shear-periodic variant for shearing-box coverage).
- Grid initialized to uniform density `ρ₀`.
- `r_return` and MPI decomposition are test parameters, not hard-coded in the pgen.

### Test matrix

| Run | `r_return` | MPI decomposition | What it exercises |
|-----|-----------|-------------------|-------------------|
| T1  | `< 0.5 × block_width` | 1 rank | single-block local deposition |
| T2  | `1.5 × block_width` | 2×2×1 | cross-block same-rank deposition |
| T3  | `1.5 × block_width` | 2×2×1, 2 ranks | cross-rank deposition |
| T4  | `= full domain` | 2×2×2, 4 ranks | full-domain path |
| T5  | `1.5 × block_width` | shearing box 2×2×1 | shear-periodic copy routing |

### Conservation assertion

At the end of each run, check:

```
Δm_particle = M_particle_initial - M_particle_final
Δm_grid     = M_grid_final - M_grid_initial
assert |Δm_particle - Δm_grid| / M < ε  (ε ~ machine epsilon × nsteps)
```

Read from history output or a custom diagnostic, not from in-memory state, so the check
exercises the full I/O path and is identical on both the `Allgatherv` and P2P builds.

### Spatial deposition check (T1 only)

For the single-block case, verify that no mass is deposited outside the sphere of radius
`r_return` centered on the particle position. Read cell-by-cell density from an HDF5
output, subtract the initial uniform density, and assert that cells outside the sphere
have zero change. This check is implementation-agnostic: it tests the deposit geometry,
not the communication pattern.

### Performance comparison workflow

Run the full T1–T5 matrix on both the `tigris-master` branch (Allgatherv) and the P2P
branch, capturing wall-clock time per cycle from the log output. The test pgen and
assertions are identical; only the binary differs. This gives a direct apples-to-apples
timing comparison without any instrumentation in the test itself.

---

## Grid Synchronization Note

Mass return modifies hydro and scalar state in `Mesh::UserWorkInLoop()`, which runs after
`PHY_BVAL` in the operator-split task list. Grid updates therefore propagate to neighbor
ghost zones at the start of the next cycle. Any same-cycle consumer that reads returned
mass from ghost zones would see stale values.

The first PR should explicitly state this ordering in the PR description and confirm that
no same-cycle consumer reads the returned hydro/scalar state before the next cycle's
boundary exchange. If one does, a targeted ghost refresh will be needed after
`UserWorkInLoop()`.

---

## Files to Change

| File | Change |
|------|--------|
| `src/particles/mass_return.hpp` | Add `BlockBounds` struct and `BuildBlockMap` static helper; remove `ImmediateNeighborReturnRange` declaration |
| `src/particles/mass_return.cpp` | Replace `AddNeighborRanks` + neighbor broadcast with block-map routing; remove fatal error for large radii |
No changes to `ops_task_list.cpp`, `particles_bvals.cpp`, or the pgen files beyond what
PR #285 already has.

The test problem (Design 2) lives in a separate PR branched from `tigris-master` and is
not listed here. Its files will be determined when that PR is planned.

---

## Out of Scope

- Accretion-delta P2P (separate project in `particles_accdelta_p2p/`).
- `r_return == 0` global return.
- AMR-aware block-map computation (uniform refinement assumed; assert at startup).
- Multi-species mass return.
