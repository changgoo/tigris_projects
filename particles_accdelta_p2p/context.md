# Accretion Delta P2P — Context

**Related issue**: [#269](https://github.com/PrincetonUniversity/tigris/issues/269)

**Goal**: Replace the `MPI_Allgatherv` path in
`ComplexParticles::ExchangeGhostAccretionDelta()` with local owner-directed
communication, while keeping the existing `INTERACT` task mostly intact.

This project is intentionally narrower than the older combined `particles_p2p/` plan.
It handles only the ghost-particle accretion correction needed before feedback. Mass
return is a separate project.

---

## Current Behavior

`ComplexParticles::InteractWithMesh()` currently groups the particle-gas work into one
operator-split task:

```text
CREATE_PAR -> SET_FBR -> SEND/RECV_GPAR[SH] -> INTERACT -> OPS_INT_COOLING
```

Within `INTERACT`, accretion can be evaluated using ghost particles. If a ghost particle
accretes, the MeshBlock that owns the ghost must return the particle mass and momentum
change to the MeshBlock that owns the active particle.

The current code records ghost-accretion entries and then calls
`ExchangeGhostAccretionDelta()`, which uses:

```text
MPI_Allgather counts
MPI_Allgatherv packed delta records
scan all records and apply those belonging to this rank
```

This is expensive because every rank receives every rank's records. It also happens
inside a per-MeshBlock task, so the collective call count depends on local MeshBlock
count.

---

## Desired Scope

Keep `INTERACT` recognizable:

```text
INTERACT_PRE_FEEDBACK
  -> SEND/RECV_ACCDELTA
  -> ACCRETION_DELTA_APPLY
  -> FEEDBACK
```

`INTERACT_PRE_FEEDBACK` should run the current merge/accretion work, but stop before
feedback injection. `FEEDBACK` should remain the existing feedback implementation as
much as possible.

The accretion-delta exchange corrects only particle properties:

- particle mass,
- particle momentum,
- any existing particle-side bookkeeping already corrected by the old delta path.

It must not update hydro or scalar grid data. Hydro changes from accretion are already
performed where accretion is evaluated; the delta return is for the active particle
state on its owner block.

---

## Ordering Answer: Does `RECV` Dependency Protect Feedback?

For the multi-MeshBlock case, yes, if the task dependencies are local to the receiving
MeshBlock and the receive task is defined correctly:

```text
INTERACT_PRE_FEEDBACK -> SEND_ACCDELTA
RECV_ACCDELTA -> ACCRETION_DELTA_APPLY -> FEEDBACK
```

`FEEDBACK` on a MeshBlock may run after `ACCRETION_DELTA_APPLY` on that same MeshBlock.
`ACCRETION_DELTA_APPLY` may run only after `RECV_ACCDELTA` has completed all inbound
delta messages addressed to that MeshBlock, including same-rank MeshBlock messages and
shear-periodic messages.

This is not a rank-wide or mesh-wide barrier. It is sufficient only because feedback
should use the local block's active particle state after all deltas targeting that
block have been applied. If future feedback logic reads particle state from another
MeshBlock directly, this dependency would no longer be enough.

---

## Required Code Knowledge

Relevant upstream locations, based on existing notes:

| Area | Path / symbol |
|------|---------------|
| Operator-split graph | `src/task_list/ops_task_list.cpp`, `OperatorSplitTaskList` |
| Current interaction task | `src/particles/complex_particles.cpp`, `InteractWithMesh()` |
| Current delta collective | `src/particles/complex_particles.cpp`, `ExchangeGhostAccretionDelta()` |
| Ghost particle exchange | `src/particles/particles_bvals.cpp`, `SendParticleBuffer()`, `ReceiveParticleBuffer()`, `FlushReceiveBuffer()` |
| Particle boundary metadata | `src/particles/particles.hpp`, `ParticleBuffer`, `pbval_->neighbor` |

Before implementation, verify exact line numbers in the current TIGRIS checkout.

---

## Main Missing Piece

Ghost particles need origin metadata. When a ghost particle is received, record the
source rank and source global MeshBlock ID. Then accretion delta records can be routed
back to the owning block instead of broadcast to all ranks.

The origin metadata must survive particle array resizing, swaps, deletion, migration,
and ghost flush operations. Do not depend on particle storage index alone.

---

## Non-Goals

- Do not refactor mass return in this project.
- Do not move mass return out of or into `INTERACT` here.
- Do not add hydro/scalar communication for accretion deltas.
- Do not redesign feedback.
- Do not make `Particles::ProcessNewParticles()` a per-MeshBlock task.
- Do not require a global barrier before feedback.

---

## Agent Rules

1. Cross-check the current TIGRIS source before claiming exact behavior.
2. Keep the `INTERACT` split minimal: accretion before exchange, feedback after apply.
3. Treat same-rank MeshBlock delivery as real communication; a local mailbox is allowed
   only as an implementation detail.
4. Route by recorded origin rank/gid and boundary metadata, not by coordinate guesses.
5. Preserve the existing `pid == NEW` position/shear matching behavior until global IDs
   are assigned after operator-split physics.
