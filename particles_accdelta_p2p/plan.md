# Accretion Delta P2P — Plan

## Design Target

Replace the `ExchangeGhostAccretionDelta()` collective with per-MeshBlock point-to-point
communication while keeping the rest of `INTERACT` close to its current structure.

Keep the implementation constrained by three coding principles:

- DRY: reuse existing particle boundary packing, matching, and shear-position handling.
- KISS: add the smallest task split that lets feedback wait for corrected particles.
- YAGNI: do not build mass-return infrastructure or a generic exchange framework in
  this project.

Target graph:

```text
CREATE_PAR
  -> SET_FBR
  -> SEND_GPAR / RECV_GPAR
  -> SEND_GPARSH / RECV_GPARSH
  -> INTERACT_PRE_FEEDBACK
  -> SEND_ACCDELTA / RECV_ACCDELTA
  -> ACCRETION_DELTA_APPLY
  -> FEEDBACK
  -> OPS_INT_COOLING
```

The exact task names can follow existing Athena++ style. The important separation is
that feedback runs only after all accretion deltas addressed to the local MeshBlock have
been received and applied.

---

## Files to Change

| File | Planned change |
|------|----------------|
| `src/task_list/ops_task_list.hpp` | Add task IDs for accretion-delta send, receive, apply, and feedback split if needed |
| `src/task_list/ops_task_list.cpp` | Register the minimal `INTERACT` split before cooling |
| `src/particles/particles.hpp` | Add ghost-origin arrays and accretion-delta buffer/channel declarations |
| `src/particles/particles_bvals.cpp` | Record ghost source rank/gid during receive flush |
| `src/particles/complex_particles.hpp/cpp` | Split pre-feedback interaction from feedback; replace Allgatherv with owner-routed delta exchange |

---

## Step 1 — Verify Current Call Boundaries

In the TIGRIS checkout, inspect:

```text
ComplexParticles::InteractWithMesh()
ComplexParticles::ExchangeGhostAccretionDelta()
Particles::FlushReceiveBuffer()
OperatorSplitTaskList::AddTask()
```

Record the current internal ordering of merge, accretion, delta exchange, and feedback
before editing. The implementation should preserve that ordering except for extracting
feedback after the new delta-apply task.

---

## Step 2 — Add Ghost Origin Tracking

Add per-particle integer arrays:

```cpp
origin_rank_
origin_gid_
```

For active local particles, initialize these to the local owner rank/gid or `-1` if the
codebase's conventions make that clearer. For received ghost particles, store the source
rank and source global block ID from the neighbor metadata used by particle boundary
communication.

Preserve origin fields through:

- capacity growth,
- particle swaps and deletion,
- active particle migration,
- ghost append/clear,
- shear-periodic ghost receive.

Acceptance check: after `SEND/RECV_GPAR[SH]`, every ghost particle that can accrete has
a valid origin rank/gid.

---

## Step 3 — Split `INTERACT` Minimally

Refactor the current interaction into two callable pieces:

| Piece | Scope | Contents |
|-------|-------|----------|
| `INTERACT_PRE_FEEDBACK` | per MeshBlock | current merge and accretion work; record ghost accretion deltas; do not call feedback |
| `FEEDBACK` | per MeshBlock | existing feedback call path, with particle state already corrected by delta apply |

Avoid broad cleanup. The desired change is a task boundary, not a redesign of particle
interaction.

---

## Step 4 — Replace the Delta Collective

Change ghost-delta records from global broadcast to owner-directed records.

Each record should include the existing matching data:

```text
pid
xp, yp, zp
delta[0..NHYDRO+NSCALARS-1]
origin_rank
origin_gid
```

When packing from ghost-accretion storage, use the recorded origin for the accreting
ghost. Do not infer the owner from current position alone.

Send records only to their owner block. For same-rank owner blocks, use the same channel
interface as remote delivery. A local mailbox is acceptable inside the channel, but the
task graph should still see this as `SEND_ACCDELTA` / `RECV_ACCDELTA`.

---

## Step 5 — Apply Before Feedback

`ACCRETION_DELTA_APPLY` runs on the owner MeshBlock after `RECV_ACCDELTA`.

Use the same matching semantics as the current code:

- match by `pid` when `pid >= 0`,
- keep position matching for `pid == NEW`,
- preserve shear-aware position handling.

Then make `FEEDBACK` depend on `ACCRETION_DELTA_APPLY`.

For multi-MeshBlock/rank, this dependency is per MeshBlock. It guarantees that feedback
for block `B` sees all deltas delivered to block `B`. It does not wait for unrelated
blocks.

---

## Step 6 — MPI Tags and Buffers

Use a distinct tag namespace for the accretion-delta channel. Before relying on high-bit
tag offsets, check `MPI_TAG_UB` at startup and fail early if the tag range is too small.

Keep the payload format simple and close to the current packed record. This project
should not introduce a general mass-return communication manager unless the existing
particle boundary code already has a natural place for it.

---

## Validation

Minimum checks:

- `rg "Allgatherv|Allgather" src/particles/complex_particles.cpp` shows no collective
  in the accretion-delta path.
- Existing accretion conservation tests still pass.
- A multi-MeshBlock/rank test exercises ghost accretion across same-rank and cross-rank
  block boundaries.
- A shear-periodic setup exercises `pid == NEW` matching and position unwrapping.
- Feedback output is unchanged for a baseline where no ghost accretion occurs.

---

## Out of Scope for This PR

- finite-radius mass return,
- global `r_return == 0` mass return,
- moving mass return task location,
- adding hydro or scalar boundary refreshes,
- changing feedback physics.
