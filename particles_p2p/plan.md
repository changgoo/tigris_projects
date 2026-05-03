# P2P Ghost Return Refactor — Implementation Plan

## Files changed

| File | Change |
|------|--------|
| `src/particles/particles.hpp` | Add `origin_rank_`, `origin_gid_`; `DeltaBuffer` struct; `delta_send_/recv_` vectors; declare `InitDeltaBuffers()` |
| `src/particles/particles_bvals.cpp` | `FlushReceiveBuffer`: record origin on ghost flush; `InitDeltaBuffers()`: allocate + tag delta buffers |
| `src/particles/complex_particles.hpp` | Declare `SendGhostAccretionDelta()`, `RecvApplyGhostAccretionDelta()` |
| `src/particles/complex_particles.cpp` | Replace `ExchangeGhostAccretionDelta()` body with P2P version |
| `src/particles/mass_return.hpp` | Declare `CollectFromNeighbors()` |
| `src/particles/mass_return.cpp` | Replace Allgatherv in `CollectParticlesInfo()` for `r_return > 0`; keep for `r_return == 0` |

---

## Step 1 — Track origin on ghost particles

### particles.hpp

Add alongside `pid`, `flag` in the `Particles` class:

```cpp
AthenaArray<int> origin_rank_;  // rank that owns the active copy (ghost particles only)
AthenaArray<int> origin_gid_;   // global block ID of the active copy (ghost particles only)
```

Resize in `UpdateCapacity()` alongside the other arrays.

### particles_bvals.cpp — FlushReceiveBuffer (line 629)

Change signature from:
```cpp
void Particles::FlushReceiveBuffer(ParticleBuffer& recv, bool ghost)
```
to:
```cpp
void Particles::FlushReceiveBuffer(ParticleBuffer& recv, bool ghost,
                                   int src_rank = -1, int src_gid = -1)
```

After the existing copy loop, add:
```cpp
if (ghost && src_rank >= 0) {
  for (int k = npartot; k < npartot + nprecv; ++k) {
    origin_rank_(k) = src_rank;
    origin_gid_(k)  = src_gid;
  }
}
```

Update the two call sites:
- `ReceiveFromNeighbors` (line 517): pass `nb.snb.rank, nb.snb.gid`
- Shearing-periodic receive (line 1112): pass `snb.rank, snb.gid`

---

## Step 2 — Delta buffer infrastructure

### particles.hpp — add to Particles protected section

```cpp
struct DeltaBuffer {
  std::vector<Real> data;           // flat: [entry_size * nentries]
  std::vector<Real> recv_data;      // pre-allocated receive buffer
  int nentries  = 0;
  int nrecv     = 0;
  int entry_size = 0;
  int tag       = -1;               // base tag for this channel
  int src_rank  = -1;               // remote rank
  MPI_Request req_s = MPI_REQUEST_NULL;  // count send
  MPI_Request req_d = MPI_REQUEST_NULL;  // data send
  MPI_Request req_rc = MPI_REQUEST_NULL; // count recv
  MPI_Request req_rd = MPI_REQUEST_NULL; // data recv
  int recv_count = 0;

  void Clear() { nentries = 0; data.clear(); }
  void Append(const Real* entry) {
    data.insert(data.end(), entry, entry + entry_size);
    ++nentries;
  }
};
std::vector<DeltaBuffer> delta_send_;  // [pbval_->nneighbor], indexed by nb.bufid
std::vector<DeltaBuffer> delta_recv_;
```

### particles_bvals.cpp — InitDeltaBuffers() (new function, called from InitParticleBvals)

```cpp
void Particles::InitDeltaBuffers(int entry_size) {
  const int n = pbval_->nneighbor;
  delta_send_.resize(n);
  delta_recv_.resize(n);

  // Tag offset: shift ipar into an unused range to avoid collisions with
  // existing ghost tags (which use ipar << 2 with bits 0-1 = 0)
  const int delta_ipar = ipar_ + 16;  // 16 > max particle types; adjust as needed

  for (int i = 0; i < n; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    SimpleNeighborBlock& snb = nb.snb;

    delta_send_[i].entry_size = entry_size;
    delta_send_[i].src_rank   = snb.rank;
    // Tag: origin receives with (my_lid, my_bufid), ghost sends to match
    delta_send_[i].tag = (snb.lid<<11) | (nb.targetid<<5) | (delta_ipar<<2);

    delta_recv_[i].entry_size = entry_size;
    delta_recv_[i].src_rank   = snb.rank;
    delta_recv_[i].tag = (pmy_block->lid<<11) | (nb.bufid<<5) | (delta_ipar<<2);
  }
}
```

Call `InitDeltaBuffers(4 + NHYDRO + NSCALARS)` from the `ComplexParticles` constructor
(or lazily on first use). For `MassReturn`, call with `PackParticleData` entry size.

---

## Step 3 — Replace ExchangeGhostAccretionDelta()

### complex_particles.cpp

Replace the existing body with:

```cpp
void ComplexParticles::ExchangeGhostAccretionDelta() {
  const int nvar = NHYDRO + NSCALARS;
  const int entry_size = 4 + nvar;

  // --- Pack: group ghost deltas by origin neighbor ---
  for (int i = 0; i < pbval_->nneighbor; ++i)
    delta_send_[i].Clear();

  for (int g = 0; g < (int)ghost_accretion_pids_.size(); ++g) {
    int org_gid = origin_gid_(npar_ + g);  // ghost particles start at npar_
    int bufid = -1;
    for (int i = 0; i < pbval_->nneighbor; ++i) {
      if (pbval_->neighbor[i].snb.gid == org_gid) { bufid = i; break; }
    }
    if (bufid < 0) continue;  // shouldn't happen

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
  // --- Send counts + data to origin ranks ---
  std::vector<int> send_counts(pbval_->nneighbor), recv_counts(pbval_->nneighbor);
  std::vector<MPI_Request> reqs;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    send_counts[i] = delta_send_[i].nentries;
    int dst = delta_send_[i].src_rank;
    if (dst == Globals::my_rank) continue;
    MPI_Request rq;
    MPI_Isend(&send_counts[i], 1, MPI_INT, dst,
              delta_send_[i].tag, my_comm, &rq);
    reqs.push_back(rq);
    MPI_Irecv(&recv_counts[i], 1, MPI_INT, dst,
              delta_recv_[i].tag, my_comm, &rq);
    reqs.push_back(rq);
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
  reqs.clear();

  // --- Send/receive data payloads ---
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    int dst = delta_send_[i].src_rank;
    if (dst == Globals::my_rank) continue;

    if (send_counts[i] > 0) {
      MPI_Request rq;
      MPI_Isend(delta_send_[i].data.data(),
                send_counts[i] * entry_size, MPI_ATHENA_REAL,
                dst, delta_send_[i].tag + 1, my_comm, &rq);
      reqs.push_back(rq);
    }
    if (recv_counts[i] > 0) {
      delta_recv_[i].recv_data.resize(recv_counts[i] * entry_size);
      MPI_Request rq;
      MPI_Irecv(delta_recv_[i].recv_data.data(),
                recv_counts[i] * entry_size, MPI_ATHENA_REAL,
                dst, delta_recv_[i].tag + 1, my_comm, &rq);
      reqs.push_back(rq);
      delta_recv_[i].nrecv = recv_counts[i];
    }
  }
  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
#else
  // Serial: copy local deltas to recv side for same-rank case
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    delta_recv_[i].recv_data = delta_send_[i].data;
    delta_recv_[i].nrecv     = delta_send_[i].nentries;
  }
#endif

  // --- Apply: same matching logic as before, per neighbor ---
  // (Mesh geometry for position matching — unchanged from old code)
  RegionSize& ms = pmy_mesh_->mesh_size;
  Real Lx = ms.x1len, Ly = ms.x2len, Lz = ms.x3len;
  Real tol = 0.5 * pmy_block->pcoord->dx1f(0);
  Real shear_deltay = 0.0;
  if (pmy_mesh_->shear_periodic)
    shear_deltay = std::fmod(qomL_ * pmy_mesh_->time, Ly);

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    int total_entries = delta_recv_[i].nrecv;
    if (total_entries == 0) continue;
    const Real* buf = delta_recv_[i].recv_data.data();

    for (int e = 0; e < total_entries; ++e) {
      const Real* base = buf + e * entry_size;
      // ... (identical matching + application loop from old ExchangeGhostAccretionDelta)
    }
    delta_recv_[i].nrecv = 0;
  }
}
```

The same-rank case (ghost and active on the same MPI rank, different MeshBlocks) is
handled by the MPI path when `dst == Globals::my_rank` — skip the MPI sends and copy
directly from `delta_send_[i].data` to `delta_recv_[i].recv_data` in the serial branch.
Or handle explicitly with a same-rank check in the MPI branch.

---

## Step 4 — Replace CollectParticlesInfo() for r_return > 0

### mass_return.cpp

```cpp
std::vector<ParticleData> MassReturn::CollectParticlesInfo() {
  std::vector<ParticleData> particles_info;
  // ... (existing local collection logic unchanged) ...

  int num_particles_local = static_cast<int>(particles_info.size());
  // ... (existing diagnostic logging unchanged) ...

#ifdef MPI_PARALLEL
  if (r_return > 0) {
    // P2P: send only to neighbors; only they can have cells within r_return
    std::vector<Real> send_buf;
    PackParticleData(particles_info, send_buf);
    int local_count = static_cast<int>(send_buf.size());
    int n = pbval_->nneighbor;
    std::vector<int> recv_counts(n);
    std::vector<MPI_Request> reqs;

    for (int i = 0; i < n; ++i) {
      int dst = pbval_->neighbor[i].snb.rank;
      if (dst == Globals::my_rank) continue;
      // Use mass_return delta buffers (separate tag namespace from accretion delta)
      MPI_Request rq;
      MPI_Isend(&local_count, 1, MPI_INT, dst, mr_tag_send_[i], my_comm, &rq);
      reqs.push_back(rq);
      MPI_Irecv(&recv_counts[i], 1, MPI_INT, dst, mr_tag_recv_[i], my_comm, &rq);
      reqs.push_back(rq);
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.clear();

    for (int i = 0; i < n; ++i) {
      int dst = pbval_->neighbor[i].snb.rank;
      if (dst == Globals::my_rank) continue;
      if (local_count > 0) {
        MPI_Request rq;
        MPI_Isend(send_buf.data(), local_count, MPI_ATHENA_REAL,
                  dst, mr_tag_send_[i] + 1, my_comm, &rq);
        reqs.push_back(rq);
      }
      if (recv_counts[i] > 0) {
        std::vector<Real> rbuf(recv_counts[i]);
        MPI_Request rq;
        MPI_Irecv(rbuf.data(), recv_counts[i], MPI_ATHENA_REAL,
                  dst, mr_tag_recv_[i] + 1, my_comm, &rq);
        // store rbuf for later unpack...
        reqs.push_back(rq);
      }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    // Unpack received neighbor particle data
    for (int i = 0; i < n; ++i) {
      if (recv_counts[i] > 0)
        UnpackParticleData(recv_bufs[i], particles_info);
    }
  } else {
    // r_return == 0: global return — all ranks genuinely need the full list
    // Keep MPI_Allgatherv (particles are deposited everywhere proportionally)
    std::vector<Real> send_buffer;
    PackParticleData(particles_info, send_buffer);
    int local_count = static_cast<int>(send_buffer.size());
    std::vector<int> counts(Globals::nranks);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> displs(Globals::nranks, 0);
    std::partial_sum(counts.begin(), counts.end() - 1, displs.begin() + 1);
    int total_count = std::accumulate(counts.begin(), counts.end(), 0);
    if (total_count == 0) return particles_info;
    std::vector<Real> recv_buffer(total_count);
    MPI_Allgatherv(send_buffer.data(), local_count, MPI_ATHENA_REAL,
                   recv_buffer.data(), counts.data(), displs.data(),
                   MPI_ATHENA_REAL, MPI_COMM_WORLD);
    UnpackParticleData(recv_buffer, particles_info);
  }
#endif

  return particles_info;
}
```

Tags `mr_tag_send_[i]` / `mr_tag_recv_[i]` use a second distinct ipar offset
(e.g., `ipar + 32`) to avoid collision with the accretion delta tags.

---

## Verification

1. Build: `configure.py --prob=tigress_classic -b -mpi --cr=mg && make all -j4`
2. Build without MPI: `configure.py --prob=tigress_classic -b && make all -j4`
3. Regression: `cd tst/regression && python run_tests.py scripts/tests/par/`
4. Short CRMHD run at 1, 2, 4, 8 ranks; compare `.hst` files for mass conservation
5. Confirm collectives removed: `grep -n Allgatherv src/particles/complex_particles.cpp src/particles/mass_return.cpp`
   — should find Allgatherv only in `r_return == 0` branch of mass_return.cpp

## Known limitation carried forward

The `MPI_Allreduce` for `total_mass_return` in `ReturnMassFromParticles()` remains.
This is a genuine global sum (all blocks contribute deposited mass back to the owner).
It can be replaced in a follow-up by: owner sends `mret` to neighbors → neighbors send
back `mass_deposited` → owner accumulates locally. But that changes the algorithm more
significantly and is deferred.
