# Plan Review 4 — Final Pre-Coding Review

**Verdict**: All five gaps from the prior review are closed. The de-prioritizations
(rule 8: defer `r_return==0`; rule 9: defer non-uniform MeshBlocks/rank) are well-formed
— each carries an explicit `ATHENA_ERROR` escape, so the first PR ships with a defined
scope. One new blocking issue surfaces because rev 4 is now concrete enough to expose it;
two secondary items need one-line decisions. After those three fixes the plan is ready
for Codex.

---

## Prior-gap status

| Gap | Status |
|-----|--------|
| G1 — Rev 2 concrete code (tags, origin arrays, `FlushReceiveBuffer` signature, parallel origin vectors) | ✅ Closed |
| G2 — Exchange manager design | ✅ Closed (Option A, full interface, mailbox structure) |
| G3 — `REFRESH_MR_GHOSTS` conditional | ✅ Closed — dropped in favor of active-zone-only deposit + existing post-cooling sync |
| G4 — `INTERACT` task split | ✅ Closed — full table with per-MeshBlock vs rank-level ownership |
| G5 — `r_return==0` commit reduction form | ✅ Closed — vector `MPI_Allreduce` specified; `ATHENA_ERROR` deferral path present |

Note on the `r_return < min_block_cells` assertion: review_final suggested it; rev 4
removes it and replaces it with a physical-radius model. This is a **correction**, not
a regression. The cell-count assertion would have been wrong once multi-MeshBlock deposit
is supported.

---

## Blocking issue — Rank-level tasks inside a per-MeshBlock task list

Athena++'s `OperatorSplitTaskList` runs per-MeshBlock. There is no built-in primitive
for "all local MeshBlocks complete task X before any proceeds to Y." Yet rev 4 registers
`ACCRETION_DELTA_EXCHANGE`, `MASS_RETURN_COLLECT`, and `MASS_RETURN_COMMIT` as
once-per-rank tasks in `ops_task_list.cpp` without naming the mechanism.

The existing codebase handles rank-level work by calling it **outside** the task list:
- `BlockFFT/FFTGravity::Solve()` — called from `main.cpp` between RK2 stages
- `Particles::ProcessNewParticles()` — called from `main.cpp` after `OperatorSplitTaskList`

**Recommended mechanism: option B — split the operator-split list**

Split `OperatorSplitTaskList` at the rank-cooperative boundaries and interleave
mesh-level calls from `main.cpp` (or `Mesh::OperatorSplitWork()`):

```text
// main.cpp, once per cycle, after all RK2 stages:
OperatorSplitTaskListPhase1()   // INTERACT_PRE_MR (per MeshBlock) stops here
pmesh->AccretionDeltaExchange() // rank-level: RankExchangeManager::Exchange(DELTA_CHANNEL)
OperatorSplitTaskListPhase2()   // ACCRETION_DELTA_APPLY → FEEDBACK_INJECT →
                                //   MASS_RETURN_DEPOSIT (per MeshBlock) stops here
pmesh->MassReturnExchange()     // rank-level: exchange deposited totals
OperatorSplitTaskListPhase3()   // MASS_RETURN_COMMIT apply → OPS_INT_COOLING → ...
Particles::ProcessNewParticles()
```

The task-graph diagram in the plan's "Final Task Graph" section should be redrawn as
three separate task-list passes with labeled `main.cpp` calls between them.

The files-changed table must add a `main.cpp` row.

Option A (atomic counter: each block increments a rank counter; the last block runs the
exchange) is also viable but requires a signal-back mechanism so the other blocks know
to proceed. Option B matches the two existing precedents and is simpler to reason about.

**Action**: pick B (or A with a full explanation), update the task graph, and add
`main.cpp` to the files-changed table.

---

## Secondary concern — `r_return` and the neighbor stencil

Rev 4 says "a single particle may deposit on one MeshBlock, several neighboring
MeshBlocks, or more than the immediate neighbor stencil if `r_return` is large enough"
and then offers the hedge "if the first implementation only supports the existing
ghost-neighbor stencil, add a runtime assertion." This leaves Codex a choice that should
be made now.

Non-neighbor routing requires querying block geometry globally — information that is not
in `pbval_->neighbor`. That is out of scope for the first PR.

**Recommendation**: commit to neighbor-stencil-only for PR 1 with a physical-radius
assertion:

```cpp
if (r_return > max_supported_exchange_distance)
  ATHENA_ERROR("r_return exceeds neighbor-stencil exchange range; "
               "multi-hop routing is deferred");
```

where `max_supported_exchange_distance` is computed from physical block extents and
NGHOST at initialization. State "neighbor-stencil-only routing" explicitly in the PR 1
scope section, and open a follow-up for the general case.

**Action**: remove the hedge from the plan and replace it with the scoped commitment
above.

---

## Minor item — Verify `ReindexOneParticleAndClear`

The "Origin Fields" section says to preserve origin fields through
`ReindexOneParticleAndClear`. This function name is not in any context or code-structure
document in this repo. Before Codex writes code that preserves origins through it,
verify the name exists in the upstream source (grep `$HOME/tigris/src/particles/`).
If wrong, substitute the correct compaction/deletion function name.

**Action**: grep for the function before implementation; correct the plan if needed.

---

## What is correct and must not change

- Two-phase active-zone-only deposit + post-cooling hyd/scalar sync (G3 resolution) is
  cleaner and cheaper than a ghost refresh before deposit.
- De-prioritizations (rules 8, 9) are well-formed; keep both with their `ATHENA_ERROR`
  guards.
- `INTERACT` split ordering: deltas applied before feedback is correct for any feedback
  logic that inspects particle mass or flags.
- Vector `MPI_Allreduce` form for `r_return==0` commit (if kept) with deterministic
  `pid >= 0` indexing.

---

## Instructions for finalizing the plan

Three edits to `plan.md`, in order:

1. **Pick a rank-level task mechanism** (blocking): choose option B (split task list +
   `main.cpp` calls), redraw the "Final Task Graph" as three passes with labeled
   `main.cpp` calls between them, and add `main.cpp` to the files-changed table.

2. **Specify mailbox thread safety** (secondary): add one sentence to "Exchange Manager
   Design" naming the chosen strategy (mutex on `Post` or per-block staging merged in
   `Exchange`).

3. **Commit to neighbor-stencil-only for PR 1** (secondary): replace the `r_return`
   stencil hedge in the `r_return` Geometry section with the scoped assertion above.

Also: before Codex begins coding, grep for `ReindexOneParticleAndClear` in the upstream
source and fix the name in the plan if needed.

After these four steps the plan is ready to hand to Codex with no remaining design
decisions deferred to implementation.
