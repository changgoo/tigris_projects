# Mass Return P2P General Routing — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the immediate-neighbor-only routing in `MassReturn::ReturnMassFromParticles(Mesh*)` with a geometry-based block-map approach that correctly handles any `r_return`, including radii spanning multiple MeshBlocks or the full domain.

**Architecture:** Build a `BlockBounds` map from `pm->loclist` (available on every rank without MPI) to compute physical bounds of all blocks. For each local particle (including periodic copies), test its return sphere against the map to find affected remote ranks. Exchange a symmetric rank list via `MPI_Allgatherv` so `ExchangeRealBuffers` stays deadlock-free. Everything else in the exchange pipeline is unchanged.

**Tech Stack:** C++11, TIGRIS (Athena++ fork), Open MPI. All changes are in `../tigris/src/particles/`.

**Branch:** `mass-return-refactor` (PR #285 branch, already checked out at `../tigris`)

---

## File Structure

| File | Change |
|------|--------|
| `../tigris/src/particles/mass_return.cpp` | Add `BlockBounds` struct, `BuildBlockMap`, `SphereIntersectsBox`, `IsFullDomain`, `RouteToRanks`, `SymmetricRanks`; remove `AddNeighborRanks` and `ImmediateNeighborReturnRange`; replace routing in `ReturnMassFromParticles(Mesh*)` |
| `../tigris/src/particles/mass_return.hpp` | No changes required |

---

## Task 1: Add BlockBounds struct and BuildBlockMap

**Files:**
- Modify: `../tigris/src/particles/mass_return.cpp` (anonymous namespace, after the existing helpers)

The anonymous namespace in `mass_return.cpp` currently ends around line 71 with `}  // namespace`. Add before the closing brace.

- [ ] **Step 1: Add BlockBounds struct and BuildBlockMap inside the anonymous namespace**

Locate the closing `} // namespace` of the anonymous namespace (after `AddMassReturnTotal`). Insert before it:

```cpp
struct BlockBounds {
  Real x1min, x1max, x2min, x2max, x3min, x3max;
  int rank;
};

// Build physical bounds for every block from the mesh's global block list.
// Uses pm->loclist (logical locations) and pm->ranklist — both replicated on
// all ranks, so no MPI communication is needed.
// Assumes uniform grid spacing (x1rat == x2rat == x3rat == 1).
std::vector<BlockBounds> BuildBlockMap(Mesh *pm) {
  std::vector<BlockBounds> map(pm->nbtotal);
  const Real x1min = pm->mesh_size.x1min;
  const Real x2min = pm->mesh_size.x2min;
  const Real x3min = pm->mesh_size.x3min;
  const Real Lx1 = pm->mesh_size.x1max - x1min;
  const Real Lx2 = pm->mesh_size.x2max - x2min;
  const Real Lx3 = pm->mesh_size.x3max - x3min;

  for (int i = 0; i < pm->nbtotal; ++i) {
    const LogicalLocation &loc = pm->loclist[i];
    // At level l, nrbx1 * 2^l blocks span x1 uniformly.
    const std::int64_t nx1 = static_cast<std::int64_t>(pm->nrbx1) << loc.level;
    const std::int64_t nx2 = static_cast<std::int64_t>(pm->nrbx2) << loc.level;
    const std::int64_t nx3 = static_cast<std::int64_t>(pm->nrbx3) << loc.level;
    map[i].x1min = x1min + static_cast<Real>(loc.lx1) * Lx1 / nx1;
    map[i].x1max = x1min + static_cast<Real>(loc.lx1 + 1) * Lx1 / nx1;
    map[i].x2min = x2min + static_cast<Real>(loc.lx2) * Lx2 / nx2;
    map[i].x2max = x2min + static_cast<Real>(loc.lx2 + 1) * Lx2 / nx2;
    map[i].x3min = x3min + static_cast<Real>(loc.lx3) * Lx3 / nx3;
    map[i].x3max = x3min + static_cast<Real>(loc.lx3 + 1) * Lx3 / nx3;
    map[i].rank = pm->ranklist[i];
  }
  return map;
}
```

- [ ] **Step 2: Verify the build compiles cleanly**

```bash
cd ../tigris && make -j4 2>&1 | grep -E "error:|warning:" | head -20
```

Expected: no errors. Warnings from unrelated files are OK.

---

## Task 2: Add sphere-box intersection and routing helpers

**Files:**
- Modify: `../tigris/src/particles/mass_return.cpp` (anonymous namespace, after `BuildBlockMap`)

- [ ] **Step 1: Add SphereIntersectsBox, IsFullDomain, and RouteToRanks inside the anonymous namespace**

Add immediately after `BuildBlockMap`:

```cpp
// Returns true if the sphere (center x1p,x2p,x3p, radius r) overlaps box.
bool SphereIntersectsBox(Real x1p, Real x2p, Real x3p, Real r,
                         const BlockBounds &box) {
  Real d2 = 0.0;
  auto clamp_dist2 = [](Real p, Real lo, Real hi) -> Real {
    if (p < lo) return (lo - p) * (lo - p);
    if (p > hi) return (p - hi) * (p - hi);
    return 0.0;
  };
  d2 += clamp_dist2(x1p, box.x1min, box.x1max);
  d2 += clamp_dist2(x2p, box.x2min, box.x2max);
  d2 += clamp_dist2(x3p, box.x3min, box.x3max);
  return d2 <= r * r;
}

// Returns true when r_return wraps more than halfway around a periodic dimension,
// meaning every block in that dimension is affected regardless of particle position.
bool IsFullDomain(Mesh *pm, Real r_return) {
  auto is_periodic = [](BoundaryFlag f) {
    return f == BoundaryFlag::periodic || f == BoundaryFlag::shear_periodic;
  };
  const RegionSize &ms = pm->mesh_size;
  if (is_periodic(pm->mesh_bcs[BoundaryFace::inner_x1]) &&
      r_return >= 0.5 * (ms.x1max - ms.x1min)) return true;
  if (ms.nx2 > 1 && is_periodic(pm->mesh_bcs[BoundaryFace::inner_x2]) &&
      r_return >= 0.5 * (ms.x2max - ms.x2min)) return true;
  if (ms.nx3 > 1 && is_periodic(pm->mesh_bcs[BoundaryFace::inner_x3]) &&
      r_return >= 0.5 * (ms.x3max - ms.x3min)) return true;
  return false;
}

// For each remote rank whose blocks intersect any particle's return sphere,
// builds a per-rank list of ParticleData records to send.
// Excludes Globals::my_rank; uses r_return for the sphere radius.
// Each (particle, rank) pair is added at most once.
std::map<int, std::vector<ParticleData>> RouteToRanks(
    const std::vector<ParticleData> &particles, Real r_return,
    const std::vector<BlockBounds> &block_map) {
  std::map<int, std::vector<ParticleData>> result;
  for (const ParticleData &pd : particles) {
    std::set<int> seen;
    for (const BlockBounds &box : block_map) {
      if (box.rank == Globals::my_rank) continue;
      if (seen.count(box.rank)) continue;
      if (SphereIntersectsBox(pd.x1, pd.x2, pd.x3, r_return, box)) {
        result[box.rank].push_back(pd);
        seen.insert(box.rank);
      }
    }
  }
  return result;
}
```

- [ ] **Step 2: Add the required include for `<set>` if not already present**

Check the top of `mass_return.cpp` — `<set>` is already included (added by PR #285). No change needed.

- [ ] **Step 3: Verify the build compiles cleanly**

```bash
cd ../tigris && make -j4 2>&1 | grep -E "error:|warning:" | head -20
```

Expected: no errors.

---

## Task 3: Add SymmetricRanks helper (MPI only)

**Files:**
- Modify: `../tigris/src/particles/mass_return.cpp` (anonymous namespace, inside `#ifdef MPI_PARALLEL`)

`ExchangeRealBuffers` posts `MPI_Irecv` from every rank in its list before sending. If rank A sends to rank B but B does not list A, B never posts a receive for A and the `MPI_Waitall` deadlocks. `SymmetricRanks` uses a small `MPI_Allgatherv` of rank-list integers (not particle data) to make every pair symmetric.

- [ ] **Step 1: Add SymmetricRanks inside the anonymous namespace, after RouteToRanks**

```cpp
#ifdef MPI_PARALLEL
// Exchanges each rank's send-rank list via Allgatherv (integers only, not particle
// data) and returns the union of send_ranks and the ranks that listed this rank in
// their own send lists. This makes the exchange symmetric so ExchangeRealBuffers
// does not deadlock when send/receive sets differ between ranks.
std::vector<int> SymmetricRanks(const std::vector<int> &send_ranks) {
  const int my_rank = Globals::my_rank;
  const int nranks = Globals::nranks;

  int my_count = static_cast<int>(send_ranks.size());
  std::vector<int> all_counts(nranks);
  MPI_Allgather(&my_count, 1, MPI_INT,
                all_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> displs(nranks, 0);
  for (int i = 1; i < nranks; ++i)
    displs[i] = displs[i - 1] + all_counts[i - 1];
  const int total = displs[nranks - 1] + all_counts[nranks - 1];

  std::vector<int> all_lists(total);
  MPI_Allgatherv(const_cast<int*>(send_ranks.data()), my_count, MPI_INT,
                 all_lists.data(), all_counts.data(), displs.data(), MPI_INT,
                 MPI_COMM_WORLD);

  std::set<int> sym_set(send_ranks.begin(), send_ranks.end());
  for (int r = 0; r < nranks; ++r) {
    if (r == my_rank) continue;
    for (int i = displs[r]; i < displs[r] + all_counts[r]; ++i) {
      if (all_lists[i] == my_rank) {
        sym_set.insert(r);
        break;
      }
    }
  }
  sym_set.erase(my_rank);
  return std::vector<int>(sym_set.begin(), sym_set.end());
}
#endif
```

- [ ] **Step 2: Verify the build compiles cleanly**

```bash
cd ../tigris && make -j4 2>&1 | grep -E "error:|warning:" | head -20
```

Expected: no errors.

---

## Task 4: Replace routing in ReturnMassFromParticles

**Files:**
- Modify: `../tigris/src/particles/mass_return.cpp`

Remove `AddNeighborRanks` and `ImmediateNeighborReturnRange` from the anonymous namespace and replace the routing logic in `ReturnMassFromParticles(Mesh*)`.

- [ ] **Step 1: Delete AddNeighborRanks and ImmediateNeighborReturnRange from the anonymous namespace**

Remove these two functions entirely (lines ~46–69 in the current file):

```cpp
void AddNeighborRanks(Mesh *pm, std::vector<int> &neighbor_ranks) {
  std::set<int> ranks;
  ...
}

Real ImmediateNeighborReturnRange(MeshBlock *pmb) {
  ...
}
```

- [ ] **Step 2: Replace ReturnMassFromParticles(Mesh *pm) with the new implementation**

Replace the entire body of `MassReturn::ReturnMassFromParticles(Mesh *pm)` (lines ~110–209) with:

```cpp
void MassReturn::ReturnMassFromParticles(Mesh *pm) {
  std::vector<MassReturn*> handlers;
  AddMassReturnHandlers(pm, handlers);
  if (handlers.empty()) return;

  std::vector<MassReturn*> deposit_handlers;
  std::set<MeshBlock*> deposit_blocks;
  for (MassReturn *handler : handlers) {
    if (deposit_blocks.insert(handler->pmy_block).second)
      deposit_handlers.push_back(handler);
  }

  for (MassReturn *handler : handlers) {
    if (handler->r_return == 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in MassReturn::ReturnMassFromParticles" << std::endl
          << "r_return == 0 global mass return is not supported by the P2P "
          << "mass-return path. Set r_return > 0." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // Collect all local particles including periodic/shear-periodic copies.
  std::vector<ParticleData> local_particles;
  for (MassReturn *handler : handlers) {
    std::vector<ParticleData> particles = handler->CollectParticlesInfo();
    local_particles.insert(local_particles.end(), particles.begin(), particles.end());
  }

  // Use the maximum r_return across all handlers for geometric routing.
  // This is conservative: particles with smaller r_return may be sent to a few
  // extra ranks, but no rank that needs a deposit will be missed.
  Real r_max = 0.0;
  for (MassReturn *handler : handlers) r_max = std::max(r_max, handler->r_return);

  // Determine which remote ranks need to receive particle records.
  std::map<int, std::vector<ParticleData>> rank_to_particles;
  if (IsFullDomain(pm, r_max)) {
    // r_return wraps past the midpoint of a periodic dimension: every rank
    // is affected, so broadcast all records to all other ranks.
    for (int r = 0; r < Globals::nranks; ++r) {
      if (r != Globals::my_rank)
        rank_to_particles[r] = local_particles;
    }
  } else {
    std::vector<BlockBounds> block_map = BuildBlockMap(pm);
    rank_to_particles = RouteToRanks(local_particles, r_max, block_map);
  }

  std::vector<int> send_ranks;
  for (const auto &kv : rank_to_particles) send_ranks.push_back(kv.first);

#ifndef MPI_PARALLEL
  if (local_particles.empty()) return;
#endif

  std::vector<ParticleData> particles_to_deposit(local_particles);

#ifdef MPI_PARALLEL
  // Make rank lists symmetric: if A sends to B, B must also list A so that
  // ExchangeRealBuffers can post matching Irecv/Isend on both sides.
  std::vector<int> symmetric_ranks = SymmetricRanks(send_ranks);

  std::vector<std::vector<Real>> send_record_buffers(symmetric_ranks.size());
  for (std::size_t n = 0; n < symmetric_ranks.size(); ++n) {
    auto it = rank_to_particles.find(symmetric_ranks[n]);
    if (it != rank_to_particles.end())
      PackParticleData(it->second, send_record_buffers[n]);
    // else: send an empty buffer — this rank is in our recv list but has
    // nothing to send to symmetric_ranks[n].
  }

  std::vector<std::vector<Real>> recv_record_buffers;
  ExchangeRealBuffers(symmetric_ranks, send_record_buffers, recv_record_buffers,
                      kMassReturnRecordCountTag, kMassReturnRecordDataTag);

  std::vector<std::map<int, Real>> totals_to_source(symmetric_ranks.size());
  for (std::size_t n = 0; n < recv_record_buffers.size(); ++n) {
    std::vector<ParticleData> remote_particles;
    UnpackParticleData(recv_record_buffers[n], remote_particles);
    for (MassReturn *handler : deposit_handlers) {
      for (const ParticleData &pd : remote_particles) {
        Real mass_added = handler->ReturnMassFromOneParticle(pd);
        AddMassReturnTotal(totals_to_source[n], pd.pid, mass_added);
      }
    }
  }
#endif

  std::map<int, Real> local_totals;
  for (MassReturn *handler : deposit_handlers) {
    for (const ParticleData &pd : particles_to_deposit) {
      Real mass_added = handler->ReturnMassFromOneParticle(pd);
      AddMassReturnTotal(local_totals, pd.pid, mass_added);
    }
  }

#ifdef MPI_PARALLEL
  std::vector<std::vector<Real>> send_total_buffers(symmetric_ranks.size());
  for (std::size_t n = 0; n < totals_to_source.size(); ++n)
    PackMassReturnTotals(totals_to_source[n], send_total_buffers[n]);

  std::vector<std::vector<Real>> recv_total_buffers;
  ExchangeRealBuffers(symmetric_ranks, send_total_buffers, recv_total_buffers,
                      kMassReturnTotalCountTag, kMassReturnTotalDataTag);

  for (const std::vector<Real> &buffer : recv_total_buffers) {
    std::map<int, Real> remote_totals;
    UnpackMassReturnTotals(buffer, remote_totals);
    for (const auto &entry : remote_totals) local_totals[entry.first] += entry.second;
  }
#endif

  for (MassReturn *handler : handlers) {
    handler->ApplyReturnedMassToLocalParticles(local_totals);
  }
}
```

- [ ] **Step 3: Verify the build compiles cleanly**

```bash
cd ../tigris && make -j4 2>&1 | grep -E "error:|warning:" | head -20
```

Expected: no errors.

- [ ] **Step 4: Confirm AddNeighborRanks and ImmediateNeighborReturnRange are gone**

```bash
grep -n "AddNeighborRanks\|ImmediateNeighborReturnRange" ../tigris/src/particles/mass_return.cpp
```

Expected: no output.

- [ ] **Step 5: Confirm no Allgatherv on particle data (only the rank-list Allgatherv in SymmetricRanks remains)**

```bash
grep -n "Allgatherv\|Allgather" ../tigris/src/particles/mass_return.cpp
```

Expected: only the two `MPI_Allgather` / `MPI_Allgatherv` calls inside `SymmetricRanks`.

---

## Task 5: Run regression test and commit

**Files:** none — validation only, then commit

- [ ] **Step 1: Run the style checker**

```bash
cd ../tigris && bash tst/style/check_athena_cpp_style.sh src/particles/mass_return.cpp
```

Expected: no style errors.

- [ ] **Step 2: Run the existing complex_particles regression test (serial)**

```bash
cd ../tigris/tst/regression && ./run_tests.py \
  -c=--include=/opt/homebrew/opt/boost/include \
  -c=--cflag=-std=c++14 \
  particles/complex_particles
```

Expected: `Summary: 1 out of 1 test passed`

- [ ] **Step 3: Run the regression test with MPI (2 ranks) to exercise cross-rank routing**

```bash
cd ../tigris/tst/regression && ./run_tests.py \
  -c=--include=/opt/homebrew/opt/boost/include \
  -c=--cflag=-std=c++14 \
  --mpirun="mpirun -np 2" \
  particles/complex_particles
```

Expected: `Summary: 1 out of 1 test passed`

- [ ] **Step 4: Smoke-test large r_return**

The athinput is at `../tigris/inputs/particles/athinput.particle_complex`. It has
`r_return = 32` and meshblock `nx1=32` over domain `x1 ∈ [-64, 64]` → block width = 64.
Set `r_return = 70` (> one block width, triggers `IsFullDomain` since 70 > 0.5×128=64)
to exercise the full-domain broadcast path without hitting the old fatal error.

```bash
cd ../tigris
cp inputs/particles/athinput.particle_complex /tmp/athinput.particle_complex.bak
sed -i 's/r_return = 32/r_return = 70/' inputs/particles/athinput.particle_complex
grep "r_return" inputs/particles/athinput.particle_complex  # confirm change
mpirun -np 2 bin/athena -i inputs/particles/athinput.particle_complex \
  job/problem_id=pcomplex time/nlim=2 output1/dt=-1 2>&1 | \
  grep -E "FATAL|SinkParticle|mass return"
cp /tmp/athinput.particle_complex.bak inputs/particles/athinput.particle_complex
```

Expected: no `FATAL ERROR` line; one or more `[SinkParticle] pid=... total mass return=...`
lines confirming deposition occurred across ranks.

- [ ] **Step 5: Commit**

```bash
cd ../tigris && git add src/particles/mass_return.cpp
git commit -m "$(cat <<'EOF'
Extend mass return P2P routing to arbitrary r_return via block map

Replace immediate-neighbor-only routing with a geometry-based block-map
approach: build physical bounds of all blocks from pm->loclist (no MPI),
test each particle's return sphere against the map, exchange a symmetric
rank list via Allgatherv of integers, then use the existing
ExchangeRealBuffers for particle data. Removes the fatal error for
r_return > one block width and the ImmediateNeighborReturnRange guard.
Full-domain wrapping (r_return >= half domain width) falls back to
sending to all ranks.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Notes

- `r_return == 0` (global return) remains a fatal error — it is out of scope for this PR.
- All handlers are assumed to have the same `r_return` in practice (single parameter in athinput). Routing uses `r_max` across handlers as a conservative bound.
- AMR is not supported: `BuildBlockMap` assumes uniform spacing (`x1rat == 1`). If AMR is needed in the future, replace the `loc.lx1 * Lx1 / nx1` formula with `Mesh::GetBlockBounds` or equivalent.
- The `SymmetricRanks` Allgatherv exchanges only integers (rank IDs), not particle data. Its O(nranks) cost is negligible compared to the particle data exchange.
