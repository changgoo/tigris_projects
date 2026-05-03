# P2P Ghost Return Refactor — Implementation Plan (rev 3)

This plan supersedes `plan_rev2_obsolete.md`. The main change is task-flow ownership:
mass return must be an operator-split physics phase, not a `USERWORK` workaround, and
all communication phases must support multiple MeshBlocks per rank.

## Design rules

1. Do not put mass-return physics in `USERWORK`.
2. Do not call MPI collectives or neighbor exchanges from independent per-MeshBlock loops.
3. Collect/exchange/commit phases run once per rank; deposit phases run per MeshBlock.
4. Same-rank neighboring MeshBlocks communicate through an explicit local mailbox or
   rank-level exchange manager, never by writing into the sender object's receive vector.
5. Boundary routing must use existing neighbor lists, block IDs, ranks, buffer IDs, and
   shear/boundary transforms. Do not assume all directions are periodic.
6. `pid=NEW` matching remains position-aware and shear-aware.

## Files changed

| File | Change |
|------|--------|
| `src/task_list/` operator-split task list | Add explicit mass-return collect/deposit/commit hooks before cooling |
| `src/particles/particles.hpp` | Add ghost origin arrays; declare rank-level exchange/mailbox helpers; define tag namespaces |
| `src/particles/particles_bvals.cpp` | Record ghost origins during receive flush; initialize safe channel tags |
| `src/particles/complex_particles.hpp/cpp` | Capture ghost-accretion origins at push time; replace delta Allgatherv with owner-routed return |
| `src/particles/mass_return.hpp/cpp` | Split mass return into collect, exchange, deposit, and commit phases |
| `src/pgen/tigress_classic.cpp` | Remove temporary `USERWORK` mass-return hook after task-list hook exists |

## Target operator-split flow

Mass return needs fresh ghost-zone information but must still run before cooling. Add a
named phase inside `OperatorSplitTaskList`:

```text
recvgpar
  -> INTERACT_PRE_MR        // merge, accretion, feedback decisions/deposits as needed
  -> REFRESH_MR_GHOSTS      // hyd/scalar boundary refresh if mass return reads ghost zones
  -> MASS_RETURN_COLLECT    // once per rank: all local MeshBlocks contribute records
  -> MASS_RETURN_DEPOSIT    // per MeshBlock: deposit using exchanged particle records
  -> MASS_RETURN_COMMIT     // once per rank: return deposited totals to particle owners
  -> OPS_INT_COOLING
```

`REFRESH_MR_GHOSTS` can reuse existing hydro/scalar boundary communication machinery,
but it must be a real task dependency. This replaces the temporary reason for calling
mass return from `USERWORK`.

## Step 1 — Ghost origin tracking for accretion deltas

Add to `Particles`:

```cpp
AthenaArray<int> origin_rank_;  // active-copy rank; -1 for active/local particles
AthenaArray<int> origin_gid_;   // active-copy gid; -1 for active/local particles
```

Initialize new capacity to `-1`. Preserve these fields through particle append,
compaction, deletion swaps, migration receives, and ghost receives. In
`FlushReceiveBuffer`, when `ghost == true`, set origins from the source
`NeighborBlock` or `SimpleNeighborBlock`.

In `ComplexParticles`, add:

```cpp
std::vector<int> ghost_accretion_origin_rank_;
std::vector<int> ghost_accretion_origin_gid_;
```

When a ghost particle at storage index `k` pushes a ghost-accretion delta, push
`origin_rank_(k)` and `origin_gid_(k)` at the same time. Clear these vectors only after
the exchange has packed/applied the entries.

## Step 2 — Safe communication channels

Use channel offsets above the existing particle-buffer tag range:

```cpp
constexpr int PAR_DELTA_TAG_BIT = 1 << 20;
constexpr int PAR_MR_TAG_BIT    = 1 << 21;
```

At initialization, compute the maximum possible base tag from local lid/bufid/ipar
usage and assert:

```cpp
base_tag < PAR_DELTA_TAG_BIT
PAR_MR_TAG_BIT + base_tag + 1 <= MPI_TAG_UB
```

Do not rely on `ipar + MAX_PARTICLE_TYPES`; those bits overlap the existing `bufid`
field.

## Step 3 — Rank-level exchange manager

Create a small exchange helper for flat `Real` records. It should support:

- remote sends/receives by `(dst_rank, dst_gid, channel)`,
- same-rank delivery through a mailbox keyed by destination gid and channel,
- count exchange followed by payload exchange,
- separate channels for accretion deltas, mass-return particle records, and deposited totals.

The helper can still use `pbval_->neighbor` for local block adjacency, but ownership is
rank-level: one rank posts one coherent exchange per phase, after all local MeshBlocks
have contributed their outgoing records.

## Step 4 — Replace `ExchangeGhostAccretionDelta()`

Packing remains per `ComplexParticles` object, but sending is deferred to the
rank-level accretion-delta exchange phase.

Record layout stays:

```text
[pid, xp, yp, zp, delta[0..NHYDRO+NSCALARS-1]]
```

For each ghost-accretion entry, route to `ghost_accretion_origin_rank_[g]` and
`ghost_accretion_origin_gid_[g]`. The destination block applies entries with the
existing pid/position matching logic. For `pid=NEW`, keep the shear-periodic unwrapping
and half-cell tolerance.

Do not silently drop a delta whose origin block is not in the local neighbor list.
That is an error unless the particle was already deleted by a documented path.

## Step 5 — Refactor mass return into collect/deposit/commit

### `MASS_RETURN_COLLECT` — once per rank

Each local MeshBlock contributes active, non-ghost particles due for mass return.
Then:

- `r_return > 0`: send records only to ranks/blocks whose domains may overlap the
  return region, using boundary-aware geometry and existing neighbor metadata.
- `r_return == 0`: use a global collective once per rank, because every rank genuinely
  needs every returning particle.

For local return, assert the physical return radius cannot reach beyond the exchange
stencil. This must be based on physical block extents/cell widths, not only raw cell
counts, unless `r_return` is confirmed to be cell units.

### `MASS_RETURN_DEPOSIT` — per MeshBlock

Each MeshBlock deposits from the records delivered to it. The geometric deposit kernel
may reuse `ReturnMassFromOneParticle()` and `ReturnMassFromOneParticleGlobal()`, but
the input list is now pre-collected. Receiving blocks must still filter geometrically
before deposit.

Each block records deposited totals per owner particle:

```text
[owner_rank, owner_gid, pid, xp, yp, zp, deposited_mass_or_vars...]
```

The position fields are required for `pid=NEW` disambiguation.

### `MASS_RETURN_COMMIT` — once per rank

Aggregate deposited totals across all local MeshBlocks, then return totals to particle
owners. Use owner-directed P2P for `r_return > 0`. For `r_return == 0`, a once-per-rank
global reduction is acceptable. The owner subtracts exactly the deposited total from
the active particle.

No `MPI_Allreduce(total_mass_return, ...)` may remain inside a per-MeshBlock method.

## Step 6 — Boundary-condition handling

All routing and matching must support:

- shear-periodic x/y boundaries,
- disk/outflow/open z boundaries,
- combined shear plus vertical boundary cases,
- same-rank and cross-rank neighbors,
- non-uniform MeshBlocks per rank.

Use existing boundary infrastructure for image positions and neighbor identity. Avoid
manual "wrap every coordinate periodically" logic. For outflow/open/disk boundaries,
only route to actual MeshBlocks that can receive deposited mass; do not create phantom
periodic partners.

## Step 7 — Remove temporary hook

After the task-list mass-return phase is working, remove the temporary
`Mesh::UserWorkInLoop()` or `MeshBlock::UserWorkInLoop()` mass-return call. `USERWORK`
should return to diagnostics/history/output work only.

## Implementation sequence

1. Add the operator-split mass-return task hook and dependencies.
2. Add rank-level exchange/mailbox helper and safe tag assertions.
3. Add and preserve ghost origin fields through particle storage operations.
4. Convert accretion-delta return to owner-routed exchange.
5. Split mass return into collect, deposit, and commit phases.
6. Implement local `r_return > 0` routing with boundary-aware overlap checks.
7. Implement once-per-rank global paths for `r_return == 0`.
8. Remove the `USERWORK` mass-return hack.
9. Run multi-MeshBlock/rank and boundary-combination tests.

## Verification

1. Build MPI and serial `tigress_classic`.
2. Run with at least two MeshBlocks per rank and non-uniform MeshBlocks per rank if the
   driver can create that layout.
3. Run shear-periodic x/y with disk/outflow z boundaries.
4. Compare 1, 2, 4, and 8 rank short runs for conserved mass and deterministic particle
   mass histories where expected.
5. Grep for `Allgatherv` and `Allreduce` in particle mass-return code. Any remaining
   collective must be in a once-per-rank task phase.
6. Test a particle near a shear-periodic corner and a particle near a vertical boundary.
