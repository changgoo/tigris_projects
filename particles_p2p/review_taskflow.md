# Plan Review 2 — Task Flow and Multi-MeshBlock/Rank Compatibility

**Verdict**: Rev 2 fixes the original tag-collision and ghost-index bugs in spirit, but
it still is not safe to implement. The plan must be re-centered around the task flow:
mass return is physics, not output work, and it must run inside the operator-split task
list with explicit communication phases that are called once per rank.

## Blocking issues

### 1. `USERWORK` is the wrong hook

Calling mass return from `Mesh::UserWorkInLoop()` was a temporary workaround to see
fresh ghost-zone data. It is not a valid final task-flow location. `USERWORK` runs
after cooling, boundary exchange, `CONS2PRIM`, physical boundaries, and output hooks.
Mass/energy returned by particles must participate in the physics sequence before
cooling and before diagnostics.

**Required fix**: introduce an explicit operator-split mass-return phase. If mass return
requires fresh mesh ghost zones after accretion/feedback, add a hyd/scalar boundary
refresh immediately before the mass-return deposit phase, not an output-time call.

### 2. Per-MeshBlock collectives are still forbidden

The original implementation assumed one MeshBlock per rank. The new implementation
must support arbitrary local MeshBlock counts. No `MPI_Allgatherv`, `MPI_Allreduce`,
or neighbor P2P exchange may be called independently from each block if another rank
may call it a different number of times.

**Required fix**: split mass return into rank-cooperative phases:
1. collect returning particles across all local MeshBlocks,
2. exchange particle records once per rank,
3. deposit independently on each local MeshBlock,
4. accumulate deposited totals across local MeshBlocks,
5. return deposited totals to owners once per rank.

### 3. Same-rank neighbor delivery is not a local vector copy

Rev 2 writes same-rank sends into the sender's receive buffer using `targetid`. That
does not deliver data to the destination MeshBlock object. With multiple MeshBlocks
per rank, same-rank communication needs an explicit rank-local mailbox keyed by
destination gid, or it must be handled by the same rank-level exchange manager that
routes remote messages.

### 4. Boundary combinations must be first-class

The refactor must work for shear-periodic x/y combinations and disk/outflow/open z
boundaries. Do not assume periodic wrapping in every direction. For each particle
record, route by actual neighbor/block geometry and boundary transforms, not by a
hard-coded periodic distance. Shear-periodic matching must preserve the existing
position unwrapping used for `pid=NEW`.

### 5. `r_return == 0` still needs a global path, but not per block

Global return genuinely needs every rank to know the active returning particle list.
Keeping Allgatherv is acceptable only if it is called once per rank from the mass-return
task hook. The deposited-total reduction has the same rule.

### 6. NEW particle IDs are assigned after operator-split physics

`Particles::ProcessNewParticles(pmesh, ipar)` is not the weak point: it is already a
mesh-level routine that loops over every local MeshBlock, all-reduces counts by global
block ID, and assigns deterministic unique IDs in gid order. That pattern is compatible
with multiple MeshBlocks per rank.

The refactor hazard is ordering. The proposed mass-return task runs before
`ProcessNewParticles`, so `pid == NEW` is still ambiguous during mass-return collect.
Mass return must exclude or error on `pid < 0` particles. Accretion-delta return may
continue to use `pid == NEW`, but only with position-aware and shear-aware matching.

## Required task-flow shape

Add mass return as a named operator-split subphase before cooling:

```
recvgpar
  -> INTERACT_PRE_MR        // merge, accretion, feedback decisions/deposits as needed
  -> REFRESH_MR_GHOSTS      // hyd/scalar boundary refresh if mass return reads ghost zones
  -> MASS_RETURN_COLLECT    // once per rank: local particle inventory + P2P/global exchange
  -> MASS_RETURN_DEPOSIT    // per MeshBlock: deposit using exchanged records
  -> MASS_RETURN_COMMIT     // once per rank: return deposited totals to particle owners
  -> OPS_INT_COOLING
```

Names can change to match Athena++ task-list style, but the ownership rules cannot.
Collect/exchange/commit are rank-cooperative. Deposit is per MeshBlock.

## Agent rules for implementation

1. Do not add physics to `USERWORK`.
2. Do not call MPI collectives inside a per-MeshBlock loop.
3. Treat same-rank, different-MeshBlock delivery as real communication through a local
   mailbox or rank-level exchange object.
4. Preserve arbitrary boundary behavior. Use existing neighbor lists, boundary flags,
   and shear transforms; do not infer neighbors from coordinate wrapping alone.
5. Keep `pid=NEW` matching position-aware and shear-aware.
6. Keep `Particles::ProcessNewParticles` mesh-level and post-operator-split.
7. Verify with at least one multi-MeshBlock/rank run and one shear-periodic + vertical
   disk/outflow configuration.
