# P2P Ghost Return Refactor — Implementation Plan (rev 2, obsolete)

*Updated after advisor review. Changes from rev 1 are marked ★.*

## Files changed

| File | Change |
|------|--------|
| `src/particles/particles.hpp` | Add `origin_rank_`, `origin_gid_` (sentinel −1); `DeltaBuffer` struct; `delta_send_/recv_`; declare `InitDeltaBuffers()` |
| `src/particles/particles_bvals.cpp` | `FlushReceiveBuffer`: record origin on ghost flush; `InitDeltaBuffers()`: allocate + tag delta buffers with `PAR_DELTA_TAG_BIT` |
| `src/particles/complex_particles.hpp` | ★ Add `ghost_accretion_origin_rank_/gid_`; declare updated `ExchangeGhostAccretionDelta()` |
| `src/particles/complex_particles.cpp` | ★ Populate origin arrays in `Accrete()`; replace `ExchangeGhostAccretionDelta()` body |
| `src/particles/mass_return.hpp` | Declare `CollectParticlesInfo()`; ★ split `ReturnMassFromParticles()` into collect + deposit |
| `src/particles/mass_return.cpp` | Replace Allgatherv for `r_return > 0`; ★ fix `r_return == 0` call structure |
| `src/particles/particles.hpp` (tags) | ★ Add `PAR_DELTA_TAG_BIT`, `PAR_MR_TAG_BIT` constants |
| `src/pgen/tigress_classic.cpp` | ★ Move `CollectParticlesInfo` for `r_return == 0` from `MeshBlock::UserWorkInLoop` to `Mesh::UserWorkInLoop` |

---

## Step 1 — Track origin on ghost particles

### particles.hpp

Add alongside `pid`, `flag` in the `Particles` class:

```cpp
AthenaArray<int> origin_rank_;  // rank of active copy; -1 for active particles
AthenaArray<int> origin_gid_;   // global block ID of active copy; -1 for active particles
```

★ Initialize to `−1` in `UpdateCapacity()` for all newly allocated slots so a stray
read on an active particle traps rather than uses garbage.

### particles_bvals.cpp — FlushReceiveBuffer (line 629)

Change signature:
```cpp
void Particles::FlushReceiveBuffer(ParticleBuffer& recv, bool ghost,
                                   int src_rank = -1, int src_gid = -1)
```

After the existing copy loop:
```cpp
if (ghost && src_rank >= 0) {
  for (int k = npartot; k < npartot + nprecv; ++k) {
    origin_rank_(k) = src_rank;
    origin_gid_(k)  = src_gid;
  }
}
```

Update two call sites:
- `ReceiveFromNeighbors` (line 517): pass `nb.snb.rank, nb.snb.gid`
- Shearing-periodic receive (line 1112): pass `snb.rank, snb.gid`
  ★ `SimpleNeighborBlock` has both `rank` and `gid` fields — this works directly.
  Add a comment at that call site to prevent future removal.

---

## Step 2 — ★ Fix index tracking for accretion deltas (was blocking bug)

The old plan read `origin_gid_(npar_ + g)` to find the origin of the g-th ghost
accretion entry. This is wrong: `g` indexes `ghost_accretion_pids_` (only accreting
ghosts), not the ghost particle storage array.

**Fix**: capture origin at push time in `Accrete()`.

### complex_particles.hpp — add parallel arrays

```cpp
// alongside ghost_accretion_pids_, ghost_accretion_xp_, etc.:
std::vector<int> ghost_accretion_origin_rank_;
std::vector<int> ghost_accretion_origin_gid_;
```

### complex_particles.cpp — AccreteFromSingleParticle() / Accrete()

Where the existing code pushes to `ghost_accretion_pids_[g]` for ghost particle at
storage index `k`, add:

```cpp
ghost_accretion_origin_rank_.push_back(origin_rank_(k));
ghost_accretion_origin_gid_.push_back(origin_gid_(k));
```

Clear both arrays alongside the others at the start of `ExchangeGhostAccretionDelta()`.

---

## Step 3 — ★ Fixed delta buffer infrastructure (tag collision fix)

### Tag constants — particles.hpp

```cpp
// High-bit flags that clear the (lid<<11 | bufid<<5 | ipar<<2) range.
// MPI_TAG_UB >= 32767 guaranteed; MPICH/OpenMPI provide >= 2^23-1.
constexpr int PAR_DELTA_TAG_BIT = 1 << 20;  // accretion delta return channel
constexpr int PAR_MR_TAG_BIT    = 1 << 21;  // mass-return neighbor channel
```

Add a startup assertion in `Particles::InitParticleBvals()`:
```cpp
#ifdef MPI_PARALLEL
int tag_ub, flag;
MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
if (!flag || tag_ub < (PAR_MR_TAG_BIT | (1<<12)))
  ATHENA_ERROR("MPI_TAG_UB too small for P2P delta tags");
#endif
```

### DeltaBuffer struct — particles.hpp (Particles protected section)

```cpp
struct DeltaBuffer {
  std::vector<Real> data;      // flat: [entry_size * nentries]
  std::vector<Real> recv_data;
  int nentries   = 0;
  int nrecv      = 0;
  int entry_size = 0;
  int tag        = -1;   // base send/recv tag for this channel (uses +0 count, +1 data)
  int src_rank   = -1;

  void Clear() { nentries = 0; data.clear(); }
  void Append(const Real* entry) {
    data.insert(data.end(), entry, entry + entry_size);
    ++nentries;
  }
};
std::vector<DeltaBuffer> delta_send_;   // indexed by nb.bufid
std::vector<DeltaBuffer> delta_recv_;
```

### InitDeltaBuffers() — particles_bvals.cpp

```cpp
void Particles::InitDeltaBuffers(int entry_size, int tag_bit) {
  const int n = pbval_->nneighbor;
  delta_send_.assign(n, DeltaBuffer{});
  delta_recv_.assign(n, DeltaBuffer{});

  for (int i = 0; i < n; ++i) {
    NeighborBlock& nb   = pbval_->neighbor[i];
    SimpleNeighborBlock& snb = nb.snb;

    delta_send_[i].entry_size = entry_size;
    delta_send_[i].src_rank   = snb.rank;
    // ★ tag_bit places tags far above the existing (lid<<11|bufid<<5|ipar<<2) range
    delta_send_[i].tag = tag_bit | (snb.lid<<11) | (nb.targetid<<5) | (ipar_<<2);

    delta_recv_[i].entry_size = entry_size;
    delta_recv_[i].src_rank   = snb.rank;
    delta_recv_[i].tag = tag_bit | (pmy_block->lid<<11) | (nb.bufid<<5) | (ipar_<<2);
  }
}
```

Call from `ComplexParticles` constructor:
```cpp
InitDeltaBuffers(4 + NHYDRO + NSCALARS, PAR_DELTA_TAG_BIT);
```

`MassReturn` will call a separate init (in the `MassReturn` constructor) with `PAR_MR_TAG_BIT`.

---

## Step 4 — Replace ExchangeGhostAccretionDelta()

```cpp
void ComplexParticles::ExchangeGhostAccretionDelta() {
  const int nvar = NHYDRO + NSCALARS;
  const int entry_size = 4 + nvar;

  // --- Pack: group deltas by origin neighbor (★ use origin arrays, not storage index) ---
  for (int i = 0; i < pbval_->nneighbor; ++i)
    delta_send_[i].Clear();

  for (int g = 0; g < (int)ghost_accretion_pids_.size(); ++g) {
    int org_gid  = ghost_accretion_origin_gid_[g];   // ★ captured at accretion time
    int org_rank = ghost_accretion_origin_rank_[g];
    int bufid = -1;
    for (int i = 0; i < pbval_->nneighbor; ++i) {
      if (pbval_->neighbor[i].snb.gid == org_gid) { bufid = i; break; }
    }
    if (bufid < 0) continue;

    Real entry[entry_size];
    entry[0] = static_cast<Real>(ghost_accretion_pids_[g]);
    entry[1] = ghost_accretion_xp_[g];
    entry[2] = ghost_accretion_yp_[g];
    entry[3] = ghost_accretion_zp_[g];
    for (int v = 0; v < nvar; ++v)
      entry[4 + v] = ghost_accretion_deltas_[g](v);
    delta_send_[bufid].Append(entry);
  }

#ifdef MPI_PARALLEL
  std::vector<int> send_counts(pbval_->nneighbor, 0);
  std::vector<int> recv_counts(pbval_->nneighbor, 0);
  std::vector<MPI_Request> reqs;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    int dst = delta_send_[i].src_rank;
    send_counts[i] = delta_send_[i].nentries;

    // ★ same-rank: direct copy (avoids deadlock, no MPI loopback needed)
    if (dst == Globals::my_rank) {
      int recv_bufid = pbval_->neighbor[i].targetid;
      delta_recv_[recv_bufid].recv_data = delta_send_[i].data;
      delta_recv_[recv_bufid].nrecv     = delta_send_[i].nentries;
      continue;
    }
    MPI_Request rq;
    MPI_Isend(&send_counts[i], 1, MPI_INT, dst,
              delta_send_[i].tag, my_comm, &rq);  reqs.push_back(rq);
    MPI_Irecv(&recv_counts[i], 1, MPI_INT, dst,
              delta_recv_[i].tag, my_comm, &rq);  reqs.push_back(rq);
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
  reqs.clear();

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    int dst = delta_send_[i].src_rank;
    if (dst == Globals::my_rank) continue;

    if (send_counts[i] > 0) {
      MPI_Request rq;
      MPI_Isend(delta_send_[i].data.data(), send_counts[i] * entry_size,
                MPI_ATHENA_REAL, dst, delta_send_[i].tag + 1, my_comm, &rq);
      reqs.push_back(rq);
    }
    if (recv_counts[i] > 0) {
      delta_recv_[i].recv_data.resize(recv_counts[i] * entry_size);
      delta_recv_[i].nrecv = recv_counts[i];
      MPI_Request rq;
      MPI_Irecv(delta_recv_[i].recv_data.data(), recv_counts[i] * entry_size,
                MPI_ATHENA_REAL, dst, delta_recv_[i].tag + 1, my_comm, &rq);
      reqs.push_back(rq);
    }
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
#else
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    delta_recv_[i].recv_data = delta_send_[i].data;
    delta_recv_[i].nrecv     = delta_send_[i].nentries;
  }
#endif

  // --- Apply: unchanged matching logic, now per-neighbor buffer ---
  RegionSize& ms = pmy_mesh_->mesh_size;
  Real Lx = ms.x1len, Ly = ms.x2len, Lz = ms.x3len;
  Real tol = 0.5 * pmy_block->pcoord->dx1f(0);
  Real shear_deltay = 0.0;
  if (pmy_mesh_->shear_periodic)
    shear_deltay = std::fmod(qomL_ * pmy_mesh_->time, Ly);

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    if (delta_recv_[i].nrecv == 0) continue;
    const Real* buf = delta_recv_[i].recv_data.data();
    for (int e = 0; e < delta_recv_[i].nrecv; ++e) {
      const Real* base = buf + e * entry_size;
      // ... (identical pid/position matching + mass application from old code)
    }
    delta_recv_[i].nrecv = 0;
  }

  // Clear origin arrays for next cycle
  ghost_accretion_origin_rank_.clear();
  ghost_accretion_origin_gid_.clear();
}
```

---

## Step 5 — ★ Fix r_return == 0 scope gap (was blocking constraint)

**Root cause confirmed**: `ReturnMassFromParticles()` is called from
`MeshBlock::UserWorkInLoop()` in `tigress_classic.cpp` (line 1061) — a per-MeshBlock
task. The Allgatherv inside is called `nblocks_per_rank` times per rank. With
non-uniform MeshBlocks this deadlocks.

**Fix**: split collect from deposit. The collection (Allgatherv or P2P) runs once per
rank; the deposit runs per MeshBlock.

### mass_return.cpp — CollectParticlesInfo() returns to a Mesh-level cache

```cpp
// New: collect once per rank from Mesh::UserWorkInLoop (global return only)
std::vector<ParticleData> MassReturn::CollectGlobalParticlesInfo();  // Allgatherv
// Existing (now P2P only for r_return > 0): called per MeshBlock
std::vector<ParticleData> MassReturn::CollectLocalParticlesInfo();   // neighbor P2P

// ReturnMassFromParticles() now takes a pre-collected list for r_return==0:
void MassReturn::ReturnMassFromParticles();                    // r_return > 0 (as before)
void MassReturn::ReturnMassFromParticles(                      // r_return == 0
    const std::vector<ParticleData>& global_particles);
```

### tigress_classic.cpp — Mesh::UserWorkInLoop() (line 1215)

```cpp
void Mesh::UserWorkInLoop() {
  // Collect mass-return info once per rank for global-return (r_return == 0) sims.
  // This runs after all MeshBlock::UserWorkInLoop() calls complete.
  for (int b = 0; b < nblocal; ++b) {
    MeshBlock* pmb = my_blocks(b);
    for (Particles* ppar : pmb->ppars) {
      if (ComplexParticles* pspar = dynamic_cast<ComplexParticles*>(ppar)) {
        if (pspar->mass_return && pspar->pmret->r_return == 0) {
          // CollectGlobalParticlesInfo is a collective — called once across all blocks
          // on all ranks by iterating only the first block (others share the result).
          if (b == 0) global_mr_info_ = pspar->pmret->CollectGlobalParticlesInfo();
          pspar->pmret->ReturnMassFromParticles(global_mr_info_);
        }
      }
    }
  }
}
```

In practice, collecting on block 0 and sharing via a Mesh-level `global_mr_info_` member
(a `std::vector<ParticleData>`) ensures the Allgatherv fires exactly once per rank per cycle.

**For r_return > 0**: no change to tigress_classic.cpp. `MeshBlock::UserWorkInLoop()`
calls `ReturnMassFromParticles()` as before; `CollectLocalParticlesInfo()` (P2P) replaces
the Allgatherv and is safe to call per-MeshBlock.

### ★ Assert r_return < block_size in MassReturn constructor

```cpp
// mass_return.cpp — MassReturn constructor
if (r_return > 0) {
  int min_block_cells = std::min({pmy_block->block_size.nx1,
                                   pmy_block->block_size.nx2,
                                   pmy_block->block_size.nx3});
  if (r_return >= min_block_cells - 1)
    ATHENA_ERROR("r_return exceeds MeshBlock size; P2P mass return is invalid");
}
```

---

## Step 6 — Replace CollectParticlesInfo() for r_return > 0

Uses the same `DeltaBuffer` infrastructure initialized with `PAR_MR_TAG_BIT`. Entry
format matches `PackParticleData` output (flat Real array per particle).

```cpp
std::vector<ParticleData> MassReturn::CollectLocalParticlesInfo() {
  // (existing local collection logic — unchanged)
  std::vector<ParticleData> particles_info = CollectOwnParticles();

#ifdef MPI_PARALLEL
  std::vector<Real> send_buf;
  PackParticleData(particles_info, send_buf);
  int local_count = static_cast<int>(send_buf.size());
  int n = pmy_par->pbval_->nneighbor;

  // Exchange counts
  std::vector<int> recv_counts(n, 0);
  std::vector<MPI_Request> reqs;
  for (int i = 0; i < n; ++i) {
    int dst = mr_delta_send_[i].src_rank;
    if (dst == Globals::my_rank) {
      // same-rank: copy directly
      int recv_bufid = pmy_par->pbval_->neighbor[i].targetid;
      mr_delta_recv_[recv_bufid].recv_data = send_buf;
      recv_counts[recv_bufid] = local_count;
      continue;
    }
    MPI_Request rq;
    MPI_Isend(&local_count, 1, MPI_INT, dst,
              mr_delta_send_[i].tag, pmy_par->my_comm, &rq);  reqs.push_back(rq);
    MPI_Irecv(&recv_counts[i], 1, MPI_INT, dst,
              mr_delta_recv_[i].tag, pmy_par->my_comm, &rq);  reqs.push_back(rq);
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
  reqs.clear();

  // Exchange data
  for (int i = 0; i < n; ++i) {
    int dst = mr_delta_send_[i].src_rank;
    if (dst == Globals::my_rank) continue;
    if (local_count > 0) {
      MPI_Request rq;
      MPI_Isend(send_buf.data(), local_count, MPI_ATHENA_REAL,
                dst, mr_delta_send_[i].tag + 1, pmy_par->my_comm, &rq);
      reqs.push_back(rq);
    }
    if (recv_counts[i] > 0) {
      mr_delta_recv_[i].recv_data.resize(recv_counts[i]);
      MPI_Request rq;
      MPI_Irecv(mr_delta_recv_[i].recv_data.data(), recv_counts[i], MPI_ATHENA_REAL,
                dst, mr_delta_recv_[i].tag + 1, pmy_par->my_comm, &rq);
      reqs.push_back(rq);
    }
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

  // Unpack
  for (int i = 0; i < n; ++i) {
    if (recv_counts[i] > 0)
      UnpackParticleData(mr_delta_recv_[i].recv_data, particles_info);
  }
#endif

  return particles_info;
}
```

---

## Implementation sequence

1. Add `ghost_accretion_origin_rank_/gid_`; populate in `Accrete()` (fixes index bug).
2. Add `PAR_DELTA_TAG_BIT`/`PAR_MR_TAG_BIT` constants + `MPI_TAG_UB` assertion (fixes tag bug).
3. Add `origin_rank_/origin_gid_` with sentinel −1; update `FlushReceiveBuffer`.
4. Implement `InitDeltaBuffers()`; call from constructors.
5. Implement `ExchangeGhostAccretionDelta()` P2P body with explicit same-rank path.
6. Add `r_return < block_size` assertion to `MassReturn` constructor.
7. Implement `CollectLocalParticlesInfo()` P2P for `r_return > 0`.
8. Split `CollectGlobalParticlesInfo()` + move to `Mesh::UserWorkInLoop()` for `r_return == 0`.
9. Run regression and conservation checks.

---

## Verification

1. Build MPI: `configure.py --prob=tigress_classic -b -mpi --cr=mg && make all -j4`
2. Build serial: `configure.py --prob=tigress_classic -b && make all -j4`
3. Regression: `cd tst/regression && python run_tests.py scripts/tests/par/`
4. Short CRMHD run at 1, 2, 4, 8 ranks; diff `.hst` files — mass must be conserved
5. Grep check: `grep -n Allgatherv src/particles/complex_particles.cpp src/particles/mass_return.cpp`
   — must find Allgatherv ONLY in `CollectGlobalParticlesInfo()` (r_return == 0)
6. ★ Conservation unit test: single sink particle near a 3-rank corner; verify
   total accreted mass matches grid mass removed across all ranks

## Known limitation (carried forward)

`MPI_Allreduce(total_mass_return)` in `ReturnMassFromParticles()` remains. It is a
genuine global sum and is not per-MeshBlock (runs once per returning particle per cycle).
Replace in a follow-up with per-neighbor accumulation if needed.
