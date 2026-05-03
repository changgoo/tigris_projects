# P2P Ghost Return Refactor — Context

**Issue**: [#269](https://github.com/PrincetonUniversity/tigris/issues/269)
**Goal**: Replace `MPI_Allgatherv` in `ExchangeGhostAccretionDelta()` and
`MassReturn::CollectParticlesInfo()` with neighbor-only point-to-point communication.

---

## Problems with the current pattern

Both functions follow the same three-step structure:
1. Pack local data into a flat buffer
2. `MPI_Allgather` counts + `MPI_Allgatherv` data to all ranks
3. Scan the full received array to find particles owned by this rank

**O(nranks) scaling**: Every rank sends to and receives from every other rank.
In practice only neighboring ranks share ghost particles, so almost all of this
traffic is wasted.

**Uniform MeshBlocks/rank requirement**: `MPI_Allgatherv` is a collective — all ranks
must call it the same number of times. Both calls happen inside the per-MeshBlock
INTERACT task (`OperatorSplitTaskList`). If ranks ever own different numbers of
MeshBlocks (AMR, load balancing), the collective deadlocks.

---

## Current code locations

### ExchangeGhostAccretionDelta — complex_particles.cpp:676

Called at the end of the accretion phase inside `InteractWithMesh()`. When a ghost
particle on Block B accretes gas, it stores the delta (mass/momentum change) locally:

```
ghost_accretion_pids_      — pid of the particle (or NEW=-1)
ghost_accretion_xp_/yp_/zp_ — position for pid=NEW disambiguation
ghost_accretion_deltas_    — AthenaArray<Real> delta per NHYDRO+NSCALARS vars
```

The function then packs, Allgatherv, and scans to apply corrections to active particles
on the owning rank. For `pid=NEW`, position matching (with shear-periodic unwrapping)
disambiguates multiple simultaneously new particles.

Buffer entry layout: `[pid, xp, yp, zp, delta[0..NHYDRO+NSCALARS-1]]`

```
entry_size = 4 + NHYDRO + NSCALARS
```

### CollectParticlesInfo — mass_return.cpp:87

Called at the start of `ReturnMassFromParticles()`. Each rank collects info about
its own active (non-ghost) particles that are due to return mass. Packs via
`PackParticleData`, Allgatherv, Unpack. The returned list is used by every rank to
deposit mass within the return radius.

**Two return modes controlled by `r_return`:**
- `r_return > 0`: geometrically local — only blocks within `r_return` of the particle
  position are affected. P2P is straightforward here.
- `r_return == 0` (global): mass is spread across ALL cells proportionally. Every rank
  genuinely needs the particle list. Allgatherv is appropriate here and will be kept.

After `ReturnMassFromParticles()` deposits mass, a second collective:
```cpp
MPI_Allreduce(total_mass_return, npar_global, MPI_SUM, MPI_COMM_WORLD)
```
sums up each particle's total deposited mass across ranks so the owner can subtract it
from the particle. This Allreduce is also valid/necessary and is NOT being replaced.

---

## Existing P2P infrastructure that can be reused

### ParticleBuffer (particle_buffer.hpp)
Packs/unpacks full particle data (intprop + realprop + auxprop) for MPI.
Used for `SEND_PAR/RECV_PAR` (active migration) and `SEND_GPAR/RECV_GPAR` (ghost copy).

### SendParticleBuffer / ReceiveParticleBuffer (particles_bvals.cpp:424–490)
Use `MPI_Send` (count) + `MPI_Isend` (ibuf + rbuf). Tags are computed as:
```
tag = (dest_lid<<11) | (dest_bufid<<5) | (ipar<<2)
```
Sub-tags +1 and +2 carry the int and real payloads. So each channel uses tags T, T+1, T+2.

### Neighbor list (pbval_->neighbor[i], i < pbval_->nneighbor)
Each `NeighborBlock` has:
- `nb.snb.rank`, `nb.snb.gid` — remote rank and global block ID
- `nb.bufid` — local buffer index for this neighbor (indexes send_[], recv_[], etc.)
- `nb.targetid` — the bufid that the remote uses for me

### FlushReceiveBuffer (particles_bvals.cpp:629)
Copies particle data from a received buffer into the local particle arrays.
Currently does NOT record which neighbor the ghost came from — this is the key gap.

---

## Missing pieces

1. **Origin tracking on ghost particles**: Ghost particles don't know which rank/block
   they came from. Need `origin_rank_` and `origin_gid_` arrays (populated at ghost
   flush time) so `ExchangeGhostAccretionDelta` knows where to send each delta.

2. **Delta buffers**: Lightweight per-neighbor flat buffers for the scalar return data
   (much simpler than full `ParticleBuffer` — no int fields needed, just Real entries).

3. **Distinct MPI tags**: The new delta-return messages need tags that don't collide with
   existing ghost-particle tags. Using `ipar + MAX_PARTICLE_TYPES` as the ipar field in
   the existing tag formula creates a collision-free namespace.

4. **`FlushReceiveBuffer` signature change**: Must be extended to receive the
   `NeighborBlock&` (or rank/gid) so origin info can be recorded.

---

## What does NOT change

- `ghost_accretion_pids_/xp_/yp_/zp_/deltas_` storage — kept as-is
- The pid / position matching logic for applying corrections — identical loop, just
  iterates per-neighbor buffer instead of one global array
- `CopyPeriodicPositions`, `ReturnMassFromOneParticle`, `ReturnMassFromOneParticleGlobal`
- The `MPI_Allreduce` for `total_mass_return` in `ReturnMassFromParticles()`
- The `MPI_Allgatherv` for `r_return == 0` in `CollectParticlesInfo()`
- Task list structure — both replacements stay within INTERACT (no new tasks)
