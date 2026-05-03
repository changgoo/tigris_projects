# P2P Ghost Return Refactor — Final Implementation Plan

This plan supersedes `plan_rev3_obsolete.md`. It keeps the rev-3 architecture and
restores the concrete rev-2 implementation details needed for coding.

## Design Rules

1. Do not put mass-return physics in `USERWORK`.
2. Do not call MPI collectives or neighbor exchanges from independent per-MeshBlock loops.
3. Rank-level phases collect/exchange/commit once per rank; deposit/apply phases run per MeshBlock.
4. Same-rank MeshBlocks communicate through a rank-level mailbox, not by writing into
   the sender object's receive vector.
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
| `src/task_list/ops_task_list.hpp` | Add task IDs and declarations for split interaction, mass-return, and refresh tasks |
| `src/task_list/ops_task_list.cpp` | Register new dependencies and task functions before cooling |
| `src/particles/particles.hpp` | Add origin arrays, tag constants, and rank exchange/mailbox declarations |
| `src/particles/particles_bvals.cpp` | Preserve origin fields and record source rank/gid on ghost flush |
| `src/particles/complex_particles.hpp/cpp` | Split accretion from delta exchange/apply; capture ghost origins at push time |
| `src/particles/mass_return.hpp/cpp` | Split collect, deposit, and commit; remove per-block collectives |
| `src/pgen/tigress_classic.cpp` | Remove temporary mass-return `USERWORK` hook |

## Final Task Graph

The deposit kernels read/write hydro variables in ghost zones. `ReturnMassFromOneParticle`
writes through `CheckInMeshBlock(..., NGHOST)`, and the global path reads ghost-zone
density/phase when `return_to_warm` is enabled. Therefore keep a real hyd/scalar ghost
refresh before mass-return deposit.

```text
recvgpar
  -> INTERACT_PRE_MR
  -> ACCRETION_DELTA_EXCHANGE      // once per rank
  -> ACCRETION_DELTA_APPLY         // per MeshBlock
  -> FEEDBACK_INJECT               // per MeshBlock
  -> MR_SEND_HYD/MR_RECV_HYD/MR_SETB_HYD
  -> MR_SEND_HYDSH/MR_RECV_HYDSH   // if shear_periodic
  -> MR_SEND_SCLR/MR_RECV_SCLR/MR_SETB_SCLR
  -> MR_SEND_SCLRSH/MR_RECV_SCLRSH // if NSCALARS > 0 and shear_periodic
  -> MASS_RETURN_COLLECT           // once per rank
  -> MASS_RETURN_DEPOSIT           // per MeshBlock
  -> MASS_RETURN_COMMIT            // once per rank
  -> OPS_INT_COOLING
```

The refresh uses hydro and scalar boundary machinery with the same face/edge/corner and
shear coverage as the post-cooling update. If `return_to_warm` depends on phase labels
derived from primitives, refresh the fields required by `pslt->AssignPhase/CheckPhase`
before `MASS_RETURN_DEPOSIT`.

## Task Split for `INTERACT`

| Step | Ownership | Depends On | Action |
|------|-----------|------------|--------|
| `INTERACT_PRE_MR` | per MeshBlock | `recvgpar` | `Merge()` then `AccreteLocal()`; accretion writes local/ghost cells and records ghost deltas, but does not exchange/apply them |
| `ACCRETION_DELTA_EXCHANGE` | once per rank | all local `INTERACT_PRE_MR` complete | pack ghost deltas from all local blocks; route to owner rank/gid |
| `ACCRETION_DELTA_APPLY` | per MeshBlock | exchange complete | apply delivered deltas to active particles using existing pid/position/shear matching |
| `FEEDBACK_INJECT` | per MeshBlock | delta apply | current `pmf->DoFeedback(pmy_block, this)` unchanged |
| `MASS_RETURN_*` | mixed | feedback inject + MR ghost refresh | collect/deposit/commit mass return before cooling |

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

For `r_return > 0`, route each returning particle to every MeshBlock whose physical
domain can intersect the sphere of radius `r_return` around the particle, including
same-rank MeshBlocks and remote MeshBlocks. A single particle may therefore deposit on
one MeshBlock, several neighboring MeshBlocks, or more than the immediate neighbor
stencil if `r_return` is large enough.

If the first implementation only supports the existing ghost-neighbor stencil, add an
explicit runtime assertion on the physical radius:

```cpp
r_return <= max_supported_exchange_distance
```

where `max_supported_exchange_distance` is derived from physical block extents and the
actual exchange stencil. Do not express this as a raw cell-count comparison.

## Exchange Manager Design

Use Option A: a Mesh-level `RankExchangeManager` owned by `Mesh` or the particle
subsystem and reused by `ComplexParticles` and `MassReturn`. It owns phase-local remote
send buffers and a same-rank mailbox:

```cpp
using MailboxKey = std::pair<int, int>;  // (dst_gid, channel)
std::map<MailboxKey, std::vector<Real>> local_mailbox_;
```

Required interface:

```cpp
void BeginPhase(int channel, int entry_size);
void Post(int dst_rank, int dst_gid, int channel, const Real* data, int n);
void Exchange(int channel);  // count exchange then payload exchange for remote ranks
void Drain(int my_gid, int channel, std::vector<Real>& out);
void EndPhase(int channel);
```

`Post` appends directly to `local_mailbox_` when `dst_rank == Globals::my_rank`; remote
records are grouped by destination rank/gid/channel and sent once per phase. `Drain`
returns all records for `(my_gid, channel)` and clears that mailbox entry. The mailbox
is cleared at `BeginPhase`.

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

- `r_return > 0`: route records to every MeshBlock whose physical domain can overlap
  the return sphere, using block geometry, neighbor/boundary metadata, and shear
  transforms. This may include multiple MeshBlocks per particle; it is not restricted
  to purely local cells or a single immediate neighbor.
- `r_return == 0`: optional/deferred. If this path stays simple, gather the returning
  particle list once per rank because every rank deposits global return. If it complicates
  the rank-level task split, reject with `ATHENA_ERROR` and document a follow-up.

### `MASS_RETURN_DEPOSIT`

Each MeshBlock deposits from its delivered records. The first implementation should focus
on `ReturnMassFromOneParticle` for `r_return > 0`. Reuse `ReturnMassFromOneParticleGlobal`
only if global return remains supported. Receiving blocks still filter geometrically.
Each block records deposited totals:

```text
[owner_rank, owner_gid, pid, xp, yp, zp, deposited_mass_or_vars...]
```

Position fields remain for debugging and any future `pid == NEW` guard, but mass-return
records must have `pid >= 0`.

### `MASS_RETURN_COMMIT`

For `r_return > 0`, aggregate local deposited totals and return them to owner rank/gid
with the rank exchange manager.

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
multiple MeshBlocks. Non-uniform MeshBlocks/rank should work if the rank-level design
naturally supports it, but it can be deferred because FFT gravity configurations are
unlikely to use it. Outflow/open/disk boundaries do not create phantom periodic partners.
Periodic/shear image positions must use existing boundary transforms, not manual
all-direction wrapping.

## Implementation Sequence

1. Add task IDs/functions in `ops_task_list.hpp/cpp` for the split graph above.
2. Add `RankExchangeManager` and safe tag assertions.
3. Add/preserve origin fields and `FlushReceiveBuffer` source arguments.
4. Split `ComplexParticles::InteractWithMesh()` into pre-MR accretion, delta exchange,
   delta apply, and feedback injection.
5. Convert accretion-delta return to owner-routed exchange.
6. Split mass return into collect/deposit/commit and exclude `pid < 0`.
7. Implement `r_return > 0` physical-overlap routing across all affected MeshBlocks.
8. Either implement `r_return == 0` as a once-per-rank vector reduction or add a clear
   runtime guard that defers global return.
9. Keep non-uniform MeshBlocks/rank support if it falls out naturally; otherwise document
   the uniform-ownership assumption and open a follow-up.
10. Remove the `USERWORK` mass-return hook.
11. Run verification below.

## Verification

1. Build MPI and serial `tigress_classic`.
2. Run at least one multi-MeshBlock/rank case with uniform MeshBlock ownership.
3. Treat non-uniform MeshBlocks/rank as optional verification. If easy, run a 2-rank
   non-uniform case; otherwise document it as deferred.
4. Run shear-periodic x/y with disk/outflow z boundaries.
5. Test `r_return > 0` cases where the return sphere is contained in one MeshBlock,
   crosses one MeshBlock boundary, crosses an edge/corner, and spans more than one
   neighbor layer if that radius is supported.
6. Check deterministic results at fixed nranks. Do not require bitwise identity across
   different nranks because P2P aggregation changes floating-point summation order.
7. Grep particle mass-return code for `Allgatherv` and `Allreduce`; any remaining
   collective must be in a once-per-rank phase.
8. Create multiple NEW particles on different MeshBlocks of the same rank and verify
   `ProcessNewParticles` assigns unique IDs after operator-split physics.
9. Test particles near a shear-periodic corner and near a vertical disk/outflow boundary.
10. If `r_return == 0` is deferred, verify the runtime guard fails early with a clear
    message.
