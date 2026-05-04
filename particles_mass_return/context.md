# Mass Return Communication — Context

**Related issue**: [#269](https://github.com/PrincetonUniversity/tigris/issues/269)

**Goal**: Plan mass return as a project separate from accretion-delta P2P. This project
can be done after the operator-split task-list work and can use `Mesh::UserWorkInLoop()`
as the mesh-level place to coordinate communication.

This document intentionally does not reuse the older combined `particles_p2p/` design,
which tried to solve accretion-delta return and mass return inside one task-flow
refactor.

---

## Current Behavior

Existing notes identify `MassReturn::CollectParticlesInfo()` as using:

```text
MPI_Allgather counts
MPI_Allgatherv particle records
```

For finite-radius return, each rank then scans returned particle records and deposits
mass on cells within `r_return`. The current implementation also needs a way to subtract
the amount actually deposited from the owner particle.

The old combined plan tried to move mass return into explicit operator-split tasks before
cooling. This made the project large because it coupled mass-return geometry, deposited
total return, active-zone-only deposition, ghost-zone synchronization, and accretion
delta ordering.

---

## New Project Boundary

Mass return is separate from the `INTERACT` refactor.

The accretion-delta project should only ensure that feedback sees corrected particle
mass and momentum. It should not solve mass return.

The mass-return project can then focus on:

- where mass return is called,
- how returning-particle records are exchanged,
- how deposited totals are returned to owner particles,
- how grid updates are synchronized after deposition.

---

## Proposed Scheduling Direction

Call mass return from `Mesh::UserWorkInLoop()` or an equivalent mesh-level hook after
the operator-split task list.

Why this is simpler:

- `Mesh::UserWorkInLoop()` is mesh-level, so a collective or rank-cooperative exchange
  can be called once per rank per cycle instead of once per MeshBlock.
- It avoids refactoring `INTERACT` beyond the accretion-delta/feedback split.
- It decouples mass return from the delicate ordering between accretion and feedback.

Known consequence:

- This places mass-return grid updates after operator-split cooling in the current task
  graph. If that changes physics ordering relative to the original intended model, the
  PR should state the tradeoff explicitly and validate diagnostics that are sensitive to
  returned mass, momentum, energy, and scalar fields.

---

## Communication Issues to Solve

There are two distinct communication problems:

1. **Particle record distribution**: each rank or MeshBlock that may receive returned
   mass needs the active particle records relevant to its active grid cells.
2. **Deposited-total return**: the particle owner must learn how much mass, momentum,
   energy, and scalar content was deposited so the active particle can be corrected
   consistently.

These should be designed independently from accretion-delta return. Accretion deltas
only correct particle mass and momentum before feedback; mass return modifies hydro and
scalar grid state.

---

## Return Modes

| Mode | Meaning | Planning stance |
|------|---------|-----------------|
| `r_return > 0` | finite-radius return around each particle | Primary target |
| `r_return == 0` | global return over all cells | Defer, keep collective once per rank, or reject if it complicates the first PR |

For `r_return > 0`, do not assume one neighbor layer is always enough. `r_return` is a
physical radius. Any P2P design must either route to all blocks whose active domains
overlap the return sphere or assert a documented maximum supported radius.

---

## Relevant Upstream Areas

| Area | Path / symbol |
|------|---------------|
| Mass return implementation | `src/particles/mass_return.cpp`, `MassReturn::ReturnMassFromParticles()` |
| Particle record collection | `src/particles/mass_return.cpp`, `CollectParticlesInfo()` |
| Problem generator hook | `src/pgen/tigress_classic.cpp`, `Mesh::UserWorkInLoop()` |
| Operator-split ordering reference | `src/task_list/ops_task_list.cpp`, `OperatorSplitTaskList` |
| Existing notes | `particles/accretion_conservation.md`, `task_flow.md` |

Before implementation, verify exact symbols and call sites in the current TIGRIS
checkout.

---

## Open Design Questions

- Does mass return intentionally need to occur before cooling, or is post-operator-split
  `Mesh::UserWorkInLoop()` acceptable for the target science runs?
- After mass-return deposition in `Mesh::UserWorkInLoop()`, which existing boundary
  communication is guaranteed to synchronize updated hydro/scalar state before the next
  consumer?
- Should the first PR keep `MPI_Allgatherv` but call it once per rank from
  `Mesh::UserWorkInLoop()`, or replace it with P2P immediately?
- Should `r_return == 0` be rejected to keep the finite-radius project small?

---

## Agent Rules

1. Do not mix this project with accretion-delta P2P.
2. Do not claim `Mesh::UserWorkInLoop()` preserves pre-cooling ordering; verify and state
   the actual ordering.
3. Do not infer affected blocks from periodic coordinate wrapping alone. Use MeshBlock
   geometry and boundary metadata.
4. Keep the first implementation focused on finite-radius return unless global return is
   explicitly required.
5. Any remaining collective must be mesh-level or rank-level, not inside an independent
   per-MeshBlock task.
