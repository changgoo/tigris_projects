# Mass Return Communication — Plan

## Design Target

Implement mass return as its own project after the accretion-delta/feedback split is
settled. Prefer a mesh-level call from `Mesh::UserWorkInLoop()` or an equivalent
post-operator-split hook so communication happens once per rank per cycle.

The first implementation should target `r_return > 0` finite-radius return. Treat
`r_return == 0` as a separate decision.

Keep the implementation constrained by three coding principles:

- DRY: reuse existing mass-return kernels and particle communication helpers.
- KISS: first separate collection, deposition, totals, and owner correction before
  adding new routing complexity.
- YAGNI: do not support global return, multi-hop routing, or task-list phases unless
  a concrete run requires them.

---

## Files to Change

| File | Planned change |
|------|----------------|
| `src/particles/mass_return.hpp/cpp` | Split collection, deposition, and owner correction so communication can be changed independently |
| `src/pgen/tigress_classic.cpp` | Call mass return from `Mesh::UserWorkInLoop()` or a mesh-level wrapper |
| `src/particles/particles.hpp` | Add any owner/deposited-total metadata needed by the mass-return exchange |
| `src/particles/particles_bvals.cpp` | Reuse particle boundary routing helpers if implementing P2P |
| `src/task_list/ops_task_list.cpp` | No change expected for the first mesh-level plan, except documentation checks |

---

## Step 1 — Confirm Scheduling

Verify the current task order in TIGRIS:

```text
OperatorSplitTaskList
  -> cooling
  -> hydro/scalar communication
  -> CONS2PRIM
  -> PHY_BVAL
  -> USERWORK
```

If mass return is called from `Mesh::UserWorkInLoop()`, document that it runs after
operator-split cooling in the current task graph. This is the main physics-ordering
tradeoff of the simpler plan.

Acceptance requirement: the PR description must explicitly say whether this post-cooling
ordering is intended or temporary.

---

## Step 2 — Separate Mass-Return Phases in Code

Refactor mass return into callable phases without changing communication yet:

| Phase | Purpose |
|-------|---------|
| collect active particle records | find local active particles eligible for return |
| distribute records | make records visible to blocks/ranks that may deposit |
| deposit on grid | update hydro/scalar active cells |
| collect deposited totals | compute what each owner particle must subtract |
| apply owner correction | update particle mass/momentum/energy bookkeeping |

Keep this split internal to `mass_return.cpp` at first. Do not create task-list phases
unless later evidence shows the mesh-level hook is insufficient.

---

## Step 3 — First Communication Choice

Choose one of two first-PR strategies.

### Option A: Mesh-Level Collective Cleanup

Keep `MPI_Allgatherv` temporarily, but call it once per rank from `Mesh::UserWorkInLoop()`
instead of from per-MeshBlock work.

Benefits:

- smallest code change,
- fixes collective call-count risk for multi-MeshBlock/rank,
- provides a stable baseline before P2P routing.

Costs:

- still O(nranks),
- does not solve scaling in issue #269.

### Option B: Finite-Radius P2P

Route returning-particle records only to blocks whose active domains overlap the return
sphere.

Benefits:

- addresses the scaling problem directly,
- avoids broadcasting all records to all ranks.

Costs:

- requires affected-block geometry,
- requires deposited-total return from receivers to owner blocks,
- needs careful boundary and same-rank handling.

Recommended sequence: implement Option A only if the goal is to reduce risk before the
full P2P project. Otherwise start with Option B, but keep it independent from accretion
delta work.

---

## Step 4 — Finite-Radius P2P Design

For each active particle due for return:

1. Compute the physical return sphere from particle position and `r_return`.
2. Find local or neighbor MeshBlocks whose active domains overlap that sphere.
3. Send the particle return record to those blocks.
4. On each receiver, deposit only into active cells.
5. Stage the deposited totals per owner particle.
6. Send deposited totals back to the owner block/rank.
7. Apply the particle correction on the owner.

Use existing neighbor and boundary metadata where possible. Do not hard-code periodic
wrapping as the routing rule.

If only the represented particle-boundary neighbor stencil is supported, add a runtime
assertion that `r_return` is within the supported exchange distance. The assertion should
state that larger radii need a wider block lookup or multi-hop routing.

---

## Step 5 — Grid Synchronization

Mass return modifies hydro and scalar grid state. The plan must identify the next
communication that makes those updates visible to neighboring MeshBlocks.

If mass return runs from `Mesh::UserWorkInLoop()` after `PHY_BVAL`, then the next normal
boundary synchronization may be in the following cycle. Verify whether any same-cycle
consumer reads returned grid quantities before that synchronization.

Do not assume accretion-delta ordering answers this question. Accretion deltas correct
particles before feedback; mass return changes grid variables.

---

## Step 6 — `r_return == 0`

Pick one policy for the first PR:

- reject `r_return == 0` with a clear runtime error,
- keep the old global collective but call it once per rank from the mesh-level hook,
- implement a separate global algorithm.

Do not let global return force complexity into the finite-radius P2P design.

---

## Validation

Minimum checks:

- `rg "Allgatherv|Allreduce" src/particles/mass_return.cpp` confirms collectives are
  either removed or called only from the mesh-level path.
- Multi-MeshBlock/rank run does not call collectives once per MeshBlock.
- Finite-radius return conserves total returned mass and momentum to roundoff.
- A boundary test covers same-rank neighbor blocks and cross-rank neighbor blocks.
- A shear-periodic setup confirms routing/deposition near the shearing boundary.
- Diagnostics sensitive to returned hydro/scalar quantities are compared before and
  after the scheduling change.

---

## Out of Scope

- Accretion-delta P2P.
- Refactoring feedback.
- Adding new operator-split mass-return tasks unless the mesh-level hook proves
  insufficient.
- Supporting arbitrary multi-hop `r_return` without a documented geometry/routing plan.
