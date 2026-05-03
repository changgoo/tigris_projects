# Plan Review 3 — Final Pre-Implementation Review

**Verdict**: Rev 3 is architecturally correct. All structural decisions from reviews 1 and
2 are resolved: mass return runs before cooling in an explicit operator-split task phase,
communication is rank-level, same-rank delivery uses an explicit mailbox, and boundary
conditions are first-class. However, rev 3 dropped most of rev 2's concrete code when
it invalidated rev 2's task-flow placement. The result is a plan that is sound in
design but underspecified for implementation. Five gaps must be filled before handing to
Codex; none requires a design change.

---

## Gap 1 — Rev 2's concrete fixes are not carried forward

Review 1 found two blocking bugs (tag collision, ghost-index misalignment). Rev 2 fixed
them with specific code. Review 2 invalidated rev 2's *task-flow placement*, not those
fixes. Rev 3 re-describes both bugs at high level but loses the implementation. Codex
will rediscover and re-fix them, at risk of getting them wrong again.

The following rev 2 content survives intact and must be merged into rev 3:

**Tag constants and startup assertion** (from rev 2, Step 3):
```cpp
constexpr int PAR_DELTA_TAG_BIT = 1 << 20;  // accretion delta return channel
constexpr int PAR_MR_TAG_BIT    = 1 << 21;  // mass-return neighbor channel
```
```cpp
// In Particles::InitParticleBvals()
int tag_ub, flag;
MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
if (!flag || tag_ub < (PAR_MR_TAG_BIT | (1<<12)))
  ATHENA_ERROR("MPI_TAG_UB too small for P2P delta tags");
```

**Parallel origin arrays captured at push time** (from rev 2, Step 2 — fixes the
index-misalignment bug; the key insight is that `g` indexes accreting ghosts, not ghost
storage order):
```cpp
// complex_particles.hpp — alongside ghost_accretion_pids_:
std::vector<int> ghost_accretion_origin_rank_;
std::vector<int> ghost_accretion_origin_gid_;
```
```cpp
// complex_particles.cpp — in Accrete(), when pushing ghost accretion at storage index k:
ghost_accretion_origin_rank_.push_back(origin_rank_(k));
ghost_accretion_origin_gid_.push_back(origin_gid_(k));
```

**Sentinel initialization** (from rev 2, Step 1):
```cpp
// particles.hpp — in UpdateCapacity(), for all newly allocated slots:
origin_rank_(k) = -1;
origin_gid_(k)  = -1;
```

**`FlushReceiveBuffer` signature change** (from rev 2, Step 1):
```cpp
void Particles::FlushReceiveBuffer(ParticleBuffer& recv, bool ghost,
                                   int src_rank = -1, int src_gid = -1)
```
Two call sites: `ReceiveFromNeighbors` (pass `nb.snb.rank, nb.snb.gid`) and the
shear-periodic receive (pass `snb.rank, snb.gid`). Add a comment at the shear-periodic
call site (line 1112 in `particles_bvals.cpp`) that `SimpleNeighborBlock` carries both
fields — do not remove these arguments as "redundant" in future simplification passes.

**`r_return` block-size assertion** (from rev 2, Step 5):
```cpp
// mass_return.cpp — MassReturn constructor
if (r_return > 0) {
  int min_block_cells = std::min({pmy_block->block_size.nx1,
                                   pmy_block->block_size.nx2,
                                   pmy_block->block_size.nx3});
  if (r_return >= min_block_cells - 1)
    ATHENA_ERROR("r_return exceeds MeshBlock size; P2P mass return is invalid");
}
```

**Action**: Merge these five items into rev 3 verbatim under a new section
"Concrete implementation details (from rev 2)."

---

## Gap 2 — Rank-level exchange manager is not architected

Rev 3 Step 3 says "create a small exchange helper" that supports remote P2P, same-rank
mailbox delivery, count-then-payload exchange, and separate channels. None of these are
spelled out. The same-rank case is the one that previously produced a bug (rev 2 wrote
into the sender's own `delta_recv_` via `targetid`, which does not deliver data to the
destination MeshBlock object).

**Required design decision before coding:**

The exchange manager must be a rank-level object (not per-MeshBlock) so it can route
across MeshBlocks on the same rank. There are two viable locations:

- **Option A**: A `RankExchangeManager` singleton or Mesh-level member, parameterized
  by channel. Both `ComplexParticles` and `MassReturn` post to it; it aggregates all
  local-MeshBlock contributions before sending.
- **Option B**: One `DeltaExchange` object per `Particles` instance (held by the base
  class), initialized with channel-specific tag bits, and capable of routing same-rank
  sends via a Mesh-level mailbox pointer.

In either case the interface must expose exactly two operations:

```cpp
// Post a record destined for (dst_rank, dst_gid) on channel ch.
void Post(int dst_rank, int dst_gid, int ch, const Real* data, int n);

// Drain all records delivered to (my_gid, ch) this phase.
// Returns a flat buffer; caller unpacks.
void Drain(int my_gid, int ch, std::vector<Real>& out);
```

The mailbox for same-rank delivery is a `std::map<(dst_gid, ch), std::vector<Real>>`
written by `Post` and read by `Drain` on the same rank within the same phase. It is
cleared at the start of each exchange phase.

**Action**: Add a short "Exchange manager design" subsection to rev 3 that picks one
option and shows the two-operation interface with the same-rank mailbox structure.

---

## Gap 3 — `REFRESH_MR_GHOSTS` is a conditional without a decision

Rev 3 says "hyd/scalar boundary refresh *if* mass return reads ghost zones." This hedge
is the unresolved root cause of mass return being in `USERWORK` in the first place.
Until the deposit kernels are actually read, the task graph shape is unknown.

**Required upstream read before writing the plan:**

Open `src/particles/mass_return.cpp` and read `ReturnMassFromOneParticle()` and
`ReturnMassFromOneParticleGlobal()`. Determine whether the deposit kernel reads any
gas conserved/primitive variable (for weighting or skip-logic) from ghost cells, or
whether it only *writes* to cells owned by the current MeshBlock.

- If deposit only writes local cells: **drop `REFRESH_MR_GHOSTS`**. The task graph
  becomes `INTERACT_PRE_MR → [accretion delta exchange] → MASS_RETURN_COLLECT →
  MASS_RETURN_DEPOSIT → MASS_RETURN_COMMIT → OPS_INT_COOLING`.
- If deposit reads ghost-zone data: **keep `REFRESH_MR_GHOSTS`** and specify which
  fields (hyd? scalar? both?), what neighbor coverage (face only? edge? corner?), and
  whether the shear-periodic variant (`SEND_HYDSH/RECV_HYDSH`) is needed. A "reuse
  existing machinery" instruction is not sufficient — the task-dependency graph for a
  full hyd+scalar refresh with shear is six tasks deep.

**Action**: Read the deposit kernel, resolve the conditional, and update the task graph
in rev 3 with the definitive shape.

---

## Gap 4 — `INTERACT_PRE_MR` scope is undefined

Rev 3 names `INTERACT_PRE_MR` but does not say how the existing `INTERACT` task splits.
Current `InteractWithMesh()` runs Merge → Accrete → DoFeedback. The accretion-delta
exchange is now a rank-level phase, not an internal call inside `Accrete()`. Without
an explicit split, `complex_particles.cpp:InteractWithMesh()` does not know what to
call when, and `ops_task_list.cpp` does not know what tasks to register.

**Required split (to be confirmed by reading `InteractWithMesh()`):**

```text
INTERACT_PRE_MR
  = Merge → Accrete (deposit on local cells + push entries to ghost_accretion arrays)
  [rank cooperative — accretion-delta exchange:
    pack per-MeshBlock, aggregate at rank level, send/recv, deliver to owners]
ACCRETION_DELTA_APPLY   (per MeshBlock — apply received deltas to local active particles)
DoFeedback              (per MeshBlock — unchanged)
MASS_RETURN_COLLECT     (once per rank)
...
```

The key question is whether `DoFeedback` requires deltas to be applied first. If not,
it can overlap with the exchange. If yes, `ACCRETION_DELTA_APPLY` must precede it.

**Action**: Read `InteractWithMesh()` and confirm the split. Add a "Task split for
INTERACT" table to rev 3 showing each sub-step, whether it is per-MeshBlock or
rank-level, and its predecessor dependency.

---

## Gap 5 — `r_return == 0` commit reduction shape

Rev 3 says "once-per-rank global reduction is acceptable" for the deposited-total
commit. Two implementations are possible:

- **Vector form (correct)**: one `MPI_Allreduce` of a length-N `Real` array where
  `N = npar_total` and the key is a globally consistent per-particle index. Each rank
  contributes its deposited amount for each particle; the owner reads its own entry.
- **Per-particle form (slow path in disguise)**: N separate `MPI_Allreduce` calls,
  one per returning particle. This was the old pattern.

The vector form requires a stable per-particle index (pid or a temporary global index
assigned at collect time). `pid == NEW` particles are excluded from mass-return collect,
so all collected particles have stable `pid >= 0`. The index can simply be `pid`
modulo the flat array, with a size determined by the global max pid.

**Action**: Add one paragraph to rev 3 under `MASS_RETURN_COMMIT` specifying the
vector-Allreduce form, the indexing scheme, and confirming that `pid < 0` sentinels
cannot appear as array keys.

---

## Minor items

**File names in changed-files table**: `src/task_list/` should be
`src/task_list/ops_task_list.cpp` (and likely the header where task enums are declared).
Check whether there is a companion `ops_task_list.hpp` that needs the new task IDs.

**Verification step 4 — determinism expectation**: State explicitly that determinism
across runs at *fixed* nranks is required, but bitwise identity across *different*
nranks is not expected (P2P aggregation changes floating-point summation order in
deposit). This is the same expectation as for the existing ghost-particle exchange.

**Non-uniform MeshBlock test**: The headline claim of this refactor is removing the
uniform-MeshBlocks/rank requirement. If `tigress_classic` cannot currently produce a
non-uniform layout (e.g., no AMR), the verification plan must name what it tests instead
— for example, a 2-rank run where rank 0 owns 2 MeshBlocks and rank 1 owns 1, produced
by specifying `<mesh>` and `<meshblock>` sizes that do not divide evenly — or
acknowledge this as untested in the first PR and add a follow-up issue.

---

## What is correct and must not change

- Operator-split task phases before cooling; collect/exchange/commit at rank level,
  deposit per MeshBlock.
- `pid < 0` exclusion from mass-return collect; `pid == NEW` kept only in the
  accretion-delta path with position + shear-aware matching.
- `MPI_Allreduce` retained for `r_return == 0` deposited totals, called once per rank
  in `MASS_RETURN_COMMIT`.
- `ProcessNewParticles` stays mesh-level and runs after operator-split physics.
- Routing via `pbval_->neighbor` + `targetid`; no coordinate-wrapping shortcuts.
- Record layout `[pid, xp, yp, zp, delta[0..NHYDRO+NSCALARS-1]]` unchanged.
- Existing pid/position + shear-periodic matching logic for applying deltas unchanged.

---

## Instructions for finalizing the plan

Complete these steps in order before beginning any code:

1. **Merge rev 2 concrete code** (Gap 1): Add a "Concrete implementation details"
   section to `plan.md` containing the five items listed under Gap 1 verbatim.

2. **Read two upstream files** (Gaps 3, 4): Open `mass_return.cpp`
   (`ReturnMassFromOneParticle`, `ReturnMassFromOneParticleGlobal`) and
   `complex_particles.cpp` (`InteractWithMesh`). These reads resolve Gaps 3 and 4
   and take less than an hour.

3. **Resolve Gap 3** (REFRESH_MR_GHOSTS): Based on the upstream read, either drop
   `REFRESH_MR_GHOSTS` from the task graph or specify it fully (fields, coverage,
   shear variant). Update the task graph diagram in `plan.md`.

4. **Resolve Gap 4** (INTERACT split): Add the task-split table to `plan.md` showing
   each sub-step of `INTERACT_PRE_MR`, its ownership (per-MeshBlock vs rank-level),
   and its predecessor.

5. **Add exchange manager sketch** (Gap 2): Add a short "Exchange manager design"
   subsection to `plan.md` with the chosen architecture (option A or B), the two-
   operation interface, and the same-rank mailbox structure.

6. **Add commit reduction spec** (Gap 5): Add one paragraph under `MASS_RETURN_COMMIT`
   in `plan.md` specifying the vector-Allreduce form and indexing scheme.

7. **Fix minor items** in `plan.md`: name `ops_task_list.cpp` in the files table,
   clarify the determinism expectation in verification step 4, and address the
   non-uniform MeshBlock test coverage.

After steps 1–7 are complete, `plan.md` (rev 4) is ready to hand to Codex. No further
design review is needed — all remaining work is implementation.
