# P2P Ghost Return Refactor — Context

**Issue**: [#269](https://github.com/PrincetonUniversity/tigris/issues/269)
**Goal**: Replace `MPI_Allgatherv` in `ExchangeGhostAccretionDelta()` and
`MassReturn::CollectParticlesInfo()` with communication that scales to multiple
MeshBlocks per rank and arbitrary boundary combinations.

---

## Problems with the current pattern

Both functions follow the same three-step structure:
1. Pack local data into a flat buffer
2. `MPI_Allgather` counts + `MPI_Allgatherv` data to all ranks
3. Scan the full received array to find particles owned by this rank

**O(nranks) scaling**: Every rank sends to and receives from every other rank.
In practice only neighboring ranks share ghost particles, so almost all of this
traffic is wasted.

**Uniform MeshBlocks/rank requirement**: `MPI_Allgatherv` is a collective -- all ranks
must call it the same number of times. Both calls happen inside the per-MeshBlock
INTERACT task (`OperatorSplitTaskList`). If ranks ever own different numbers of
MeshBlocks (AMR, load balancing), the collective deadlocks.

Non-uniform MeshBlocks/rank is a useful future capability, but it is not the main target
for the first refactor because FFT gravity configurations are unlikely to use it. The
first priority is multi-MeshBlock/rank with normal uniform ownership, same-rank routing,
and finite-radius mass return across affected MeshBlocks.

**Task-flow problem**: Mass return is currently routed through `USERWORK` as a
temporary hack so it can read updated ghost zones. That is not an acceptable final
location. `USERWORK` runs after cooling, boundary updates, primitive recovery, physical
boundaries, and diagnostics. The refactor must give mass return an explicit
operator-split task-list phase before cooling.

For the prioritized `r_return > 0` path, do not add a new full hydro/scalar ghost
refresh before mass return. Instead, route returning-particle records to every affected
MeshBlock and deposit only into active cells on each receiver. This preserves the
pre-cooling physics order and lets the existing post-cooling boundary communication
synchronize updated active zones.

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
  position are affected. P2P is valid only if routing uses the actual MeshBlock
  neighbor/boundary geometry. `r_return` is a physical radius, so this can involve one
  MeshBlock, several adjacent MeshBlocks, or more than one neighbor layer if the radius
  is large enough.
- `r_return == 0` (global): mass is spread across ALL cells proportionally. This mode
  is unlikely to be useful and is not the priority for this refactor. If supporting it
  adds significant complexity, explicitly reject it at runtime and revisit in a later
  PR.

After `ReturnMassFromParticles()` deposits mass, a second collective:
```cpp
MPI_Allreduce(total_mass_return, npar_global, MPI_SUM, MPI_COMM_WORLD)
```
sums up each particle's total deposited mass across ranks so the owner can subtract it
from the particle. This reduction is globally meaningful, but the old per-MeshBlock
call pattern is not valid for multi-MeshBlock/rank. It must be moved into the same
rank-cooperative mass-return commit phase, replaced with owner-directed P2P totals, or
left unsupported for `r_return == 0` in the first implementation.

---

## Required task-flow target

Mass return must become a coherent operator-split phase, not a `USERWORK` side effect.
The target ordering is:

```
recvgpar
  -> INTERACT_PRE_MR        // merge, accretion, feedback decisions/deposits as needed
  -> MASS_RETURN_COLLECT    // once per rank: local inventory + P2P/global exchange
  -> MASS_RETURN_DEPOSIT    // per MeshBlock: deposit active zones on this block
  -> MASS_RETURN_COMMIT     // once per rank: return deposited totals to owners
  -> OPS_INT_COOLING
```

The exact task names can follow Athena++ conventions, but the ownership split is fixed:
collect/exchange/commit run once per rank; deposit runs per MeshBlock.

`Particles::ProcessNewParticles(pmesh, ipar)` is already a mesh-level operation: it
counts `pid == NEW` particles on every local MeshBlock, all-reduces counts indexed by
global block ID, and assigns IDs in deterministic gid order. This is compatible with
multiple MeshBlocks per rank and must remain outside per-MeshBlock task execution.
The new mass-return task runs before this ID assignment, so mass-return collection must
ignore or reject `pid < 0` particles. Accretion-delta return is the only phase that may
communicate `pid == NEW`, and it must keep position/shear matching until IDs are assigned.

---

## Boundary-condition requirements

The implementation must support arbitrary active boundary combinations, especially:

- shear-periodic in x/y,
- disk, outflow, or open boundaries in z,
- same-rank neighboring MeshBlocks,
- multiple MeshBlocks per rank. Non-uniform MeshBlocks/rank should be supported if it
  falls out naturally, but can be deferred if it adds substantial complexity.

Routing must use `pbval_->neighbor`, global block IDs, ranks, buffer IDs, and existing
boundary/shear transforms. Do not infer communication partners from coordinate wrapping
alone. `pid=NEW` matching must remain position-aware and shear-aware. For `r_return > 0`,
the affected-block set is defined by physical overlap between the return sphere and
MeshBlock domains, not by a hard-coded one-neighbor stencil unless that stencil is
guarded by an explicit physical-radius assertion.

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

2. **Return/exchange buffers**: Lightweight buffers for accretion deltas and mass-return
   particle records. Same-rank delivery needs a rank-local mailbox or rank-level exchange
   manager, not a copy into the sender object's receive vector.

3. **Distinct MPI tags**: The new delta-return messages need tags that don't collide with
   existing ghost-particle tags. High-bit tag namespaces are acceptable only with a
   startup assertion that the complete base tag plus channel offset fits under
   `MPI_TAG_UB`.

4. **`FlushReceiveBuffer` signature change**: Must be extended to receive the
   `NeighborBlock&` (or rank/gid) so origin info can be recorded.

5. **Task-list hook**: Add explicit mass-return task phases before cooling. Do not
   implement the final refactor by calling mass return from `USERWORK`.

---

## What does NOT change

- `ghost_accretion_pids_/xp_/yp_/zp_/deltas_` storage — kept as-is
- The pid / position matching logic for applying corrections — identical loop, just
  iterates per-neighbor buffer instead of one global array
- `CopyPeriodicPositions`, `ReturnMassFromOneParticle`, `ReturnMassFromOneParticleGlobal`
  numerical kernels, except for call-site/argument changes needed by the new task flow
- `r_return == 0` physics semantics if it is implemented later; the first P2P refactor
  may reject this mode
- The conservation requirement that owners subtract exactly the mass/energy deposited

---

## Rules for AI agents

1. Do not add or keep mass-return physics in `USERWORK`.
2. Do not call collectives or neighbor exchanges from independent per-MeshBlock loops.
3. Treat same-rank MeshBlock-to-MeshBlock delivery as real communication through an
   explicit local mailbox or rank-level exchange object.
4. Preserve boundary-condition generality; do not assume all boundaries are periodic.
5. Do not change `Particles::ProcessNewParticles` into a per-MeshBlock task; it is the
   mesh-level barrier that makes `pid == NEW` unique after operator-split physics.
6. For `r_return > 0`, do not depend on fresh mass-return ghost zones. Deposit active
   zones on all affected MeshBlocks and rely on the existing post-cooling boundary sync.
