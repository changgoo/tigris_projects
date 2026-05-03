# P2P Ghost Return Refactor — Final Implementation Plan

This plan supersedes `plan_rev3_obsolete.md`. It keeps the rev-3 architecture and
restores the concrete rev-2 implementation details needed for coding.

## Design Rules

1. Do not put mass-return physics in `USERWORK`.
2. Do not call MPI collectives or neighbor exchanges from independent per-MeshBlock loops.
3. All new exchange/deposit/apply phases in `OperatorSplitTaskList` are per-MeshBlock
   tasks. Do not introduce mesh-level or rank-level tasks inside the task list.
4. Same-rank MeshBlocks communicate through the same per-MeshBlock channel abstraction
   as cross-rank neighbors. If MPI loopback is not robust, use a small local mailbox
   hidden inside that channel, not a rank-level task.
5. Keep `Particles::ProcessNewParticles` mesh-level and post-operator-split.
6. Exclude `pid < 0` from mass-return collect. `pid == NEW` is allowed only in the
   accretion-delta path with position- and shear-aware matching.
7. Route through `pbval_->neighbor`, gids, ranks, `targetid`, and existing boundary/shear
   transforms. Do not infer partners from coordinate wrapping alone.
8. Prioritize `r_return > 0`. If `r_return == 0` adds significant complexity, reject it
   at runtime and defer global return to a follow-up PR.
9. Support multiple MeshBlocks per rank. Non-uniform MeshBlocks/rank is desirable, but
   not a first-PR blocker when combined with FFT gravity constraints; defer it if it
   materially complicates the implementation.

## Files Changed

| File | Change |
|------|--------|
| `src/task_list/ops_task_list.hpp` | Add task IDs and declarations for per-MeshBlock P2P exchange phases |
| `src/task_list/ops_task_list.cpp` | Register per-MeshBlock exchange/deposit/commit tasks before cooling |
| `src/particles/particles.hpp` | Add origin arrays, tag constants, and per-MeshBlock channel declarations |
| `src/particles/particles_bvals.cpp` | Preserve origin fields and record source rank/gid on ghost flush |
| `src/particles/complex_particles.hpp/cpp` | Split accretion from delta exchange/apply; capture ghost origins at push time |
| `src/particles/mass_return.hpp/cpp` | Split collect, deposit, and commit; remove per-block collectives |
| `src/pgen/tigress_classic.cpp` | Remove temporary mass-return `USERWORK` hook |

## PR Split

### PR 1: Accretion Delta P2P Infrastructure

This PR should be small enough to review independently and should not move mass return
out of `USERWORK` yet.

- Add per-MeshBlock channel buffers, channel tag namespaces, and `MPI_TAG_UB` checks.
- Add and preserve `origin_rank_` / `origin_gid_` through particle migration, ghost
  receive, deletion swaps, and capacity changes.
- Extend `FlushReceiveBuffer` so ghost particles know source rank and global block ID.
- Split `ComplexParticles::InteractWithMesh()` into `INTERACT_PRE_MR`,
  `SEND/RECV_ACCDELTA`, `ACCRETION_DELTA_APPLY`, and `FEEDBACK_INJECT`.
- Replace `ExchangeGhostAccretionDelta()` Allgatherv with owner-routed per-MeshBlock
  P2P delivery.

PR 1 must leave the existing mass-return behavior unchanged except where shared
infrastructure is added. The acceptance test is that accretion delta return no longer
uses per-MeshBlock collectives and still works with multiple MeshBlocks per rank,
including same-rank neighbors and shear-periodic ghost particles.

### PR 2: Mass Return Task Flow

This PR consumes the PR 1 infrastructure and removes the temporary `USERWORK` hook.

- Add finite-radius `SEND/RECV_MR_RECORDS`, `MASS_RETURN_DEPOSIT_ACTIVE`,
  `SEND/RECV_MR_TOTALS`, and `MASS_RETURN_COMMIT_APPLY` tasks.
- Route `r_return > 0` records to all represented neighbor-stencil MeshBlocks whose
  active domains overlap the physical return sphere.
- Deposit only active zones on receiving MeshBlocks and return deposited totals to the
  owner block before cooling.
- Reject or defer `r_return == 0` if it adds significant complexity.
- Remove `src/pgen/tigress_classic.cpp` mass-return `USERWORK` call.

PR 2 must not add a pre-mass-return hydro/scalar ghost refresh. Existing post-cooling
hydro/scalar communication remains the synchronization point for active-zone deposits.

## Final Task Graph

For the prioritized `r_return > 0` path, do not add a mass-return hydro/scalar ghost
refresh. The extra refresh would be equivalent to another full hydro/scalar boundary
communication, including shear variants, and would be expensive. Instead, route each
returning-particle record to every MeshBlock whose active domain overlaps the physical
return sphere. Each receiving MeshBlock deposits only into its active cells.

The old kernel writes through `CheckInMeshBlock(..., NGHOST)`. Refactor that behavior:
ghost-zone deposition is replaced by owner-block active-zone deposition plus a deposited
total returned to the particle owner. The existing post-cooling hydro/scalar boundary
communication then synchronizes the updated active zones before `CONS2PRIM`.

```text
CREATE_PAR -> SET_FBR -> SEND/RECV_GPAR[SH]
  -> INTERACT_PRE_MR
  -> SEND_ACCDELTA / RECV_ACCDELTA
  -> ACCRETION_DELTA_APPLY
  -> FEEDBACK_INJECT
  -> SEND_MR_RECORDS / RECV_MR_RECORDS
  -> MASS_RETURN_DEPOSIT_ACTIVE
  -> SEND_MR_TOTALS / RECV_MR_TOTALS
  -> MASS_RETURN_COMMIT_APPLY
  -> OPS_INT_COOLING
  -> existing SEND/RECV/SETB_HYD and SEND/RECV/SETB_SCLR tasks
  -> REMOVE_PAR -> CONS2PRIM -> PHY_BVAL -> USERWORK -> CLEAR_ALLBND
```

`TaskList::DoTaskListOneStage` executes tasks independently per MeshBlock. Therefore
the new communication phases must be per-MeshBlock P2P tasks, like existing particle
and hydro boundary tasks. They must not contain collectives and must not wait for "all
local MeshBlocks" to reach the same point.

Do not place mass return after the existing post-cooling communication to reuse that
sync point. That would move particle mass return after cooling and would either change
physics ordering or require another communication before later consumers. The finite-
radius path should run before cooling and avoid pre-MR ghost synchronization by using
active-zone-only deposition on all affected MeshBlocks.

If a later `r_return == 0` implementation needs ghost-zone density/phase information
for global weighting, handle it separately; this is another reason to defer global
return when it complicates the first PR.

## Task Split for `INTERACT`

| Step | Ownership | Depends On | Action |
|------|-----------|------------|--------|
| `INTERACT_PRE_MR` | per MeshBlock | `recvgpar` | `Merge()` then `AccreteLocal()`; accretion writes local/ghost cells and records ghost deltas, but does not exchange/apply them |
| `SEND_ACCDELTA` | per MeshBlock | `INTERACT_PRE_MR` | pack this block's ghost deltas into neighbor/owner buffers and send counts/payloads |
| `RECV_ACCDELTA` | per MeshBlock | matching receives posted/arrived | receive delta counts/payloads addressed to this block |
| `ACCRETION_DELTA_APPLY` | per MeshBlock | `RECV_ACCDELTA` | apply delivered deltas to active particles using existing pid/position/shear matching |
| `FEEDBACK_INJECT` | per MeshBlock | delta apply | current `pmf->DoFeedback(pmy_block, this)` unchanged |
| `SEND_MR_RECORDS` | per MeshBlock | feedback inject | send eligible `r_return > 0` records to affected neighbor-stencil blocks |
| `RECV_MR_RECORDS` | per MeshBlock | matching receives posted/arrived | receive finite-radius mass-return records addressed to this block |
| `MASS_RETURN_DEPOSIT_ACTIVE` | per MeshBlock | `RECV_MR_RECORDS` | deposit active zones only and stage deposited totals |
| `SEND_MR_TOTALS` | per MeshBlock | deposit active | send deposited totals back to owner blocks |
| `RECV_MR_TOTALS` | per MeshBlock | matching receives posted/arrived | receive deposited totals for local active particles |
| `MASS_RETURN_COMMIT_APPLY` | per MeshBlock | `RECV_MR_TOTALS` | owner blocks subtract returned totals from active particles |

Apply accretion deltas before feedback so particle masses and flags are final for any
feedback logic that inspects particle state.

## Concrete Implementation Details

### Tag Constants

```cpp
constexpr int PAR_DELTA_TAG_BIT = 1 << 20;  // accretion delta return channel
constexpr int PAR_MR_TAG_BIT    = 1 << 21;  // mass-return neighbor channel
```

In `Particles::InitParticleBvals()`:

```cpp
int tag_ub, flag;
MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
if (!flag || tag_ub < (PAR_MR_TAG_BIT | (1<<12)))
  ATHENA_ERROR("MPI_TAG_UB too small for P2P delta tags");
```

Also compute the maximum existing base tag and assert `base_tag < PAR_DELTA_TAG_BIT`
and `PAR_MR_TAG_BIT + base_tag + 1 <= tag_ub`.

### Origin Fields

```cpp
AthenaArray<int> origin_rank_;
AthenaArray<int> origin_gid_;
```

Initialize all new slots in `UpdateCapacity()`:

```cpp
origin_rank_(k) = -1;
origin_gid_(k)  = -1;
```

Preserve these fields through `ReindexOneParticleAndClear`, active migration, deletion
swaps, ghost append, and capacity changes.

Change `FlushReceiveBuffer`:

```cpp
void Particles::FlushReceiveBuffer(ParticleBuffer& recv, bool ghost,
                                   int src_rank = -1, int src_gid = -1)
```

Pass `nb.snb.rank, nb.snb.gid` in `ReceiveFromNeighbors`; pass `snb.rank, snb.gid` in
the shear-periodic receive. Add a comment there that `SimpleNeighborBlock` carries both
fields and the arguments are required for ghost-origin tracking.

### Accretion Delta Origins

Add beside `ghost_accretion_pids_`:

```cpp
std::vector<int> ghost_accretion_origin_rank_;
std::vector<int> ghost_accretion_origin_gid_;
```

When pushing a ghost accretion entry at storage index `k`:

```cpp
ghost_accretion_origin_rank_.push_back(origin_rank_(k));
ghost_accretion_origin_gid_.push_back(origin_gid_(k));
```

This fixes the index-misalignment bug: the vector index `g` counts accreting ghosts,
not ghost storage order.

### `r_return` Geometry

`r_return` is a physical radius, not a cell-count limit. The implementation must not
assert `r_return < min(meshblock.nx*)`. Local return means "finite-radius return", not
"single-MeshBlock return."

For PR 2, route each returning particle only through the existing neighbor stencil,
including face/edge/corner neighbors and shear-periodic neighbors represented in the
particle boundary infrastructure. This supports finite-radius return across the local
block and represented neighbors. General multi-hop/non-neighbor block lookup is deferred.

Add an explicit runtime assertion on the physical radius:

```cpp
if (r_return > max_supported_exchange_distance)
  ATHENA_ERROR("r_return exceeds neighbor-stencil exchange range; "
               "multi-hop routing is deferred");
```

where `max_supported_exchange_distance` is derived from physical block extents and the
actual neighbor stencil. Do not express this as a raw cell-count comparison.

## Per-MeshBlock Exchange Design

Use per-MeshBlock flat buffers modeled after existing particle boundary buffers, not a
rank-level exchange manager. Each channel has send/receive buffers indexed by
`pbval_->neighbor[i].bufid` and uses tags derived from destination lid, destination
buffer id, particle type, and a high-bit channel namespace.

Required operations per channel:

```cpp
void ClearChannelBuffers(int channel);
void PackToNeighbor(int bufid, const Real* data, int n);
void SendChannel(int channel);     // send counts then payloads to each neighbor
bool ReceiveChannel(int channel);  // receive counts/payloads; true when complete
```

Same-rank neighboring MeshBlocks are still real communication partners. In flat MPI,
the first implementation may use the same MPI send/receive path and tags as cross-rank
neighbors. If MPI loopback proves fragile, add a small local mailbox keyed by
`(dst_gid, channel)` inside the per-MeshBlock channel implementation, but do not promote
the exchange to a rank-level task.

## Accretion Delta Return

Record layout remains:

```text
[pid, xp, yp, zp, delta[0..NHYDRO+NSCALARS-1]]
```

For each ghost-delta entry, route to `ghost_accretion_origin_rank_[g]` and
`ghost_accretion_origin_gid_[g]`. The owner block drains by gid and applies the existing
matching logic. For `pid == NEW`, keep position matching with periodic/shear unwrapping
and half-cell tolerance. Do not silently drop missing-origin records.

## Mass Return Phases

### `MASS_RETURN_COLLECT`

Each local MeshBlock contributes active, non-ghost particles due for return. Only
`pid >= 0` is eligible; `pid == NEW` or `pid == DEL` with positive return is an error or
an explicit diagnostic skip.

- `r_return > 0`: route records to every represented neighbor-stencil MeshBlock whose
  physical domain can overlap the return sphere, using block geometry,
  neighbor/boundary metadata, and shear transforms. Multi-hop routing beyond this
  stencil is deferred and guarded by `max_supported_exchange_distance`.
- `r_return == 0`: optional/deferred. If this path stays simple, gather the returning
  particle list in a task-list-safe way because every rank deposits global return. If it
  complicates the per-MeshBlock task flow, reject with `ATHENA_ERROR` and document a
  follow-up.

### `MASS_RETURN_DEPOSIT`

Each MeshBlock deposits from its delivered records. The first implementation should focus
on an active-zone-only variant of `ReturnMassFromOneParticle` for `r_return > 0`.
Receiving blocks still filter geometrically, but they must not write to ghost zones.
Reuse `ReturnMassFromOneParticleGlobal` only if global return remains supported. Each
block records deposited totals:

```text
[owner_rank, owner_gid, pid, xp, yp, zp, deposited_mass_or_vars...]
```

Position fields remain for debugging and any future `pid == NEW` guard, but mass-return
records must have `pid >= 0`.

### `MASS_RETURN_COMMIT`

For `r_return > 0`, each depositing MeshBlock sends deposited totals back to the owner
MeshBlock with the per-MeshBlock `SEND_MR_TOTALS`/`RECV_MR_TOTALS` channel. The owner
block applies all totals addressed to its active particles in `MASS_RETURN_COMMIT_APPLY`.

If `r_return == 0` is kept, use one vector `MPI_Allreduce`, not one reduction per
particle. Build a deterministic list of collected `pid >= 0` values during
`MASS_RETURN_COLLECT` or allocate by global `max_pid + 1` if memory is acceptable. Each
rank contributes the deposited total at that particle's index; owners subtract their
entry from the active particle. Negative sentinel IDs cannot appear as keys.

If this global path makes the first refactor materially harder, add an input/runtime
guard:

```cpp
if (r_return == 0)
  ATHENA_ERROR("Global mass return is temporarily unsupported by P2P mass return");
```

No `MPI_Allreduce(total_mass_return, ...)` may remain inside a per-MeshBlock method.

## Boundary Requirements

Support shear-periodic x/y, disk/outflow/open z, same-rank neighbors, cross-rank
neighbors, multiple MeshBlocks per rank, and finite-radius return regions that span
multiple MeshBlocks. Non-uniform MeshBlocks/rank should work if the per-MeshBlock P2P
design naturally supports it, but it can be deferred because FFT gravity configurations are
unlikely to use it. Outflow/open/disk boundaries do not create phantom periodic partners.
Periodic/shear image positions must use existing boundary transforms, not manual
all-direction wrapping.

## Implementation Sequence

Follow the PR split above.

For PR 1:

1. Add per-MeshBlock channel buffers and safe tag assertions.
2. Add/preserve origin fields and `FlushReceiveBuffer` source arguments.
3. Add per-MeshBlock accretion-delta task IDs/functions in `ops_task_list.hpp/cpp`.
4. Split `ComplexParticles::InteractWithMesh()` into pre-MR accretion, delta exchange,
   delta apply, and feedback injection.
5. Convert accretion-delta return to owner-routed per-MeshBlock exchange.
6. Verify no per-MeshBlock collective remains in accretion-delta return.

For PR 2:

1. Add per-MeshBlock mass-return task IDs/functions in `ops_task_list.hpp/cpp`.
2. Split mass return into collect/deposit/commit and exclude `pid < 0`; make finite-
   radius deposit active-zone-only.
3. Implement `r_return > 0` neighbor-stencil physical-overlap routing with a physical
   radius guard for deferred multi-hop routing.
4. Either implement `r_return == 0` as a simple follow-on path or add a clear runtime
   guard that defers global return.
5. Keep non-uniform MeshBlocks/rank support if it falls out naturally; otherwise document
   the uniform-ownership assumption and open a follow-up.
6. Remove the `USERWORK` mass-return hook.
7. Run verification below.

## Verification

1. Build MPI and serial `tigress_classic`.
2. Run at least one multi-MeshBlock/rank case with uniform MeshBlock ownership.
3. Treat non-uniform MeshBlocks/rank as optional verification. If easy, run a 2-rank
   non-uniform case; otherwise document it as deferred.
4. Run shear-periodic x/y with disk/outflow z boundaries.
5. Test `r_return > 0` cases where the return sphere is contained in one MeshBlock,
   crosses one MeshBlock boundary, and crosses an edge/corner within the neighbor
   stencil.
6. Check deterministic results at fixed nranks. Do not require bitwise identity across
   different nranks because P2P aggregation changes floating-point summation order.
7. Grep particle mass-return code for `Allgatherv` and `Allreduce`; no collective may
   remain inside per-MeshBlock task-list execution.
8. Create multiple NEW particles on different MeshBlocks of the same rank and verify
   `ProcessNewParticles` assigns unique IDs after operator-split physics.
9. Test particles near a shear-periodic corner and near a vertical disk/outflow boundary.
10. If `r_return == 0` is deferred, verify the runtime guard fails early with a clear
    message.
11. Verify `r_return > max_supported_exchange_distance` fails early with a clear
    multi-hop-routing-deferred message.
