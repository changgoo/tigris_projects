# Plan Review — P2P Ghost Return Refactor

**Verdict**: Plan is structurally correct but contains two implementation-blocking bugs
(tag collision, index misalignment) and one scope gap (`r_return==0` may not remove
the uniform-MeshBlocks/rank constraint). Address all three before writing code.

---

## Blocking Issue 1 — Tag formula collision

The plan proposes `delta_ipar = ipar + 16` and reuses:
```
tag = (lid<<11) | (bufid<<5) | (delta_ipar<<2)
```

This **collides with existing tags**. `bufid` ranges 0–55, so `bufid<<5` covers the same
bit positions as the shifted `delta_ipar`. Concretely: `delta_ipar=16` gives
`delta_ipar<<2 = 64`, which equals `bufid=2, ipar=0` in the existing formula.

**Fix**: use a high bit-flag that clears the entire `(lid, bufid, ipar)` range:

```cpp
constexpr int PAR_DELTA_TAG_BIT = 1 << 20;  // above any lid<<11 for typical lid values
delta_send_[i].tag = PAR_DELTA_TAG_BIT | (snb.lid<<11) | (nb.targetid<<5) | (ipar<<2);
delta_recv_[i].tag = PAR_DELTA_TAG_BIT | (pmy_block->lid<<11) | (nb.bufid<<5) | (ipar<<2);
```

Use a separate constant (e.g., `PAR_MR_TAG_BIT = 1<<21`) for the mass-return channel
so the two new channels don't collide with each other either.

Note: MPI standard guarantees `MPI_TAG_UB >= 32767`; MPICH and OpenMPI both provide
`MPI_TAG_UB = 2^23-1` or larger, so bit 20 is safe in practice. Add a startup assertion:
```cpp
int tag_ub; int flag;
MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
ATHENA_REQUIRE(tag_ub >= (1<<21), "MPI_TAG_UB too small for P2P delta tags");
```

---

## Blocking Issue 2 — Index misalignment in pack step

The plan's pack step writes:
```cpp
int org_gid = origin_gid_(npar_ + g);   // WRONG
```
`g` indexes `ghost_accretion_pids_` (only ghost particles that accreted), but
`origin_gid_(npar_ + g)` indexes the g-th ghost particle in storage order. These are
not the same: ghost particles are stored contiguously from `npar_` onward, but only a
subset of them have accretion entries, and their storage positions are arbitrary.

**Fix**: capture origin at the moment the accretion entry is pushed in `Accrete()`.
Add parallel arrays to `ComplexParticles`:

```cpp
std::vector<int> ghost_accretion_origin_rank_;
std::vector<int> ghost_accretion_origin_gid_;
```

In `Accrete()`, when pushing to `ghost_accretion_pids_[g]` for ghost particle at index `k`:
```cpp
ghost_accretion_origin_rank_.push_back(origin_rank_(k));
ghost_accretion_origin_gid_.push_back(origin_gid_(k));
```

Then in `ExchangeGhostAccretionDelta()`:
```cpp
int org_gid = ghost_accretion_origin_gid_[g];   // correct
int org_rank = ghost_accretion_origin_rank_[g];
```

No neighbor-list lookup is needed at exchange time — both rank and gid are already known.

---

## Scope Gap — `r_return == 0` may not fix the uniform-MB/rank constraint

Issue #269's primary motivation is removing the uniform-MeshBlocks/rank requirement.
That constraint exists because `MPI_Allgatherv` is called from a per-MeshBlock task.
The plan **retains Allgatherv for `r_return == 0`** without confirming whether
`ReturnMassFromParticles()` (which calls `CollectParticlesInfo()`) is called from a
per-MeshBlock task.

**Action required**: trace the call chain. If `ReturnMassFromParticles()` is called from
`InteractWithMesh()` (as suggested by the INTERACT task graph in `task_flow.md`), the
constraint persists for any run using global mass return.

Two options:
1. **Hoist** the `r_return == 0` Allgatherv out of the per-MeshBlock task by calling it
   once per rank (requires a small ops-task-list restructure).
2. **Replace** with: each rank applies only its own particles' global return;
   use `MPI_Allreduce` over per-particle `mret` values (one global sum, but no longer
   per-MeshBlock if called once per cycle).

If standard CRMHD runs only use `r_return > 0`, document this explicitly as a known
limitation of the plan and add a comment in the code.

---

## Detail Issues

### `r_return < block_size` assumption (issue 4)
The geometric "only neighbors are affected" argument for `r_return > 0` holds only
when `r_return < min(Nx1,Nx2,Nx3) - 1` (in cells). Add a runtime assertion in the
`MassReturn` constructor or `InitMassReturn()`:
```cpp
ATHENA_REQUIRE(r_return < pmy_block->block_size.nx1 - 1,
  "r_return exceeds MeshBlock size — P2P mass return is not valid");
```

### Sentinel initialization for `origin_rank_/origin_gid_` (issue 5)
Initialize both arrays to `-1` in `UpdateCapacity()`. Active particles never get these
fields set; a stray read should trap, not silently use uninitialized memory.

### Same-rank handling needs explicit pseudocode (issue 6)
The plan says "handle explicitly with a same-rank check." Spell it out:
```cpp
if (dst == Globals::my_rank) {
  // Both ghost and active are on this rank (different MeshBlocks).
  // Read directly: the receiver's bufid for this rank is nb.targetid.
  int recv_bufid = pbval_->neighbor[i].targetid;
  delta_recv_[recv_bufid].recv_data = delta_send_[i].data;
  delta_recv_[recv_bufid].nrecv = delta_send_[i].nentries;
  continue;
}
```
Without this, the code either deadlocks (MPI_Send to self with no matching recv) or
relies on MPI_THREAD_MULTIPLE semantics that may not be guaranteed.

### `FlushReceiveBuffer` call site for shear-periodic (issue 7)
The shear-periodic call site (line 1112) provides a `SimpleNeighborBlock` (not
`NeighborBlock`), but `SimpleNeighborBlock` has both `rank` and `gid` fields, so
passing them works. Add a comment at that call site to prevent future "simplification."

---

## Strengths

The approach is architecturally right: tracking origin at ghost-flush time (not in the
wire format) avoids changing the packed buffer layout. Per-neighbor `DeltaBuffer` is
appropriate — the payload is scalars only, far lighter than `ParticleBuffer`. The
decision to keep `MPI_Allreduce(total_mass_return, ...)` is correct: that one IS a
genuinely global sum. Keeping the pid/position matching logic unchanged means the
conservation-correctness argument from `accretion_conservation.md` carries over intact.

---

## Recommended implementation sequence

1. Fix index tracking: add `ghost_accretion_origin_rank_/gid_` and populate in `Accrete()`.
2. Fix tag formula: use `PAR_DELTA_TAG_BIT` approach; add startup assertion on `MPI_TAG_UB`.
3. Verify call chain of `ReturnMassFromParticles()`; decide `r_return==0` strategy.
4. Add `r_return < block_size` assertion.
5. Implement `InitDeltaBuffers()` + `FlushReceiveBuffer` origin recording.
6. Implement `ExchangeGhostAccretionDelta()` P2P body (with same-rank path spelled out).
7. Implement `CollectParticlesInfo()` P2P path for `r_return > 0`.
8. Add conservation unit test (sink particle near 3-rank corner) before merging.
