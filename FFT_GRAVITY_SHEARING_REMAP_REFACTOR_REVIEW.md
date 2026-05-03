# FFT gravity shearing remap refactor — plan review

Date: 2026-05-01

## Diagnosis confirmation

The timing data unambiguously supports the plan's conclusion. `ShearSourceMax` (~10 ms
gap) and `RetrieveMax` (~4 ms gap) account for essentially all of the remaining
`blockfft / fft` total delta. Forward/Kernel/Backward are already within 1 ms of each
other. Local micro-optimizations (slab memcpy, nonblocking exchange) failed to move
the needle, confirming the bottleneck is architectural — the global-row exchange path
is the wrong abstraction and must be replaced.

---

## Plan strengths

- Correct identification that the remaining gap is in shearing remap, not FFT transforms.
- Phased migration with a fallback flag (`use_meshblock_shearing_remap`) is the right
  risk posture: allows benchmarking each step independently.
- Shared `ShearingRemapper` used by both `FFTGravity` and `BlockFFTGravity` avoids
  permanently duplicating the remap logic.
- `RollUnrollAll(dt)` bulk-call API solves the same-rank target buffer ownership
  problem that one-at-a-time calls cannot.
- Staged segment approach for integer shift (plan section "Implementation details")
  is the right abstraction.

---

## Issues and proposed solutions

### 1. MPI design for integer shift: start with rank-aggregated payloads, not per-segment tags

The plan recommends "one-message-per-segment with hashed tag, then optimize later."
With 8 MB/rank × 64 i-columns × up to 4 segments/i-column × 2 directions, this can
produce hundreds of small messages per remap call. Tag-hashing risks collisions and
hits MPI tag-limit constraints on some implementations.

**Solution:** Start the implementation with one `MPI_Isend`/`MPI_Irecv` per
destination rank:

1. Reserve a base physics tag from `Mesh::ReserveTagPhysIDs()`.
2. Post all receives first (avoids deadlock with blocking sends in `RollUnroll`).
3. Pack a metadata header (count of segments, `src_gid`/`dst_gid`/`src_i`/...  per
   segment) followed by the concatenated payload into a single buffer per destination
   rank.
4. Post all sends.
5. `MPI_Waitall`.
6. Unpack.

This is more code than the per-segment approach but is the design that scales to many
MeshBlocks/rank without tag pressure. Do not defer this to a later optimization.

### 2. Same-rank, different-MeshBlock is a correctness requirement, not an optimization

The current `BlockFFTGravity::RollUnroll()` only short-circuits to local copy when
`Ngrids == 0` (same-block self-copy). When `target_gid != my_gid` but
`ranklist[target_gid] == Globals::my_rank`, the existing code falls through to
`MPI_Send`/`MPI_Recv` with `destination_rank == MPI_COMM_WORLD rank == source_rank`.
This is a self-send: valid in most MPI implementations but serialized and undefined in
terms of buffering requirements.

The new helper must explicitly test `ranklist[target_gid] == Globals::my_rank` for
every segment target and copy directly from `src_block.roll_buf` to `dst_block.dat`
without MPI. Make this an acceptance criterion for Phase 2 (see Phase-gate criteria
below).

### 3. `phi` vs `phi_gasonly` handling is unspecified in the new flow

The existing code calls `RetrieveAppliedShearingResult` twice with different destination
arrays but reads the same pre-unrolled FFT output both times. If the proposed
`RetrieveShearedResult(buf, loc, bsize)` is also called twice, the helper buffer must
be preserved unchanged between the two calls (i.e., `RollUnrollAll(+dt)` must not be
called again between them).

**Required clarification in the implementation:** either (a) run the full result path
once into a dedicated `phi` buffer and once into `phi_gasonly`, with `RollUnrollAll`
called only once after both FFT retrieves, or (b) call `RetrieveShearedResult` twice
sharing the same remapper buffer (re-using the already-unrolled result). Option (a) is
cleaner and uses the existing helper structure. The plan should explicitly choose one.

### 4. Normalization placement: apply inside `RetrieveShearedResult`

The optimization context notes that the previous refactor introduced a normalization
timing risk (raw `out_` copied before norm applied). The new result-path flow creates
the same risk again: FFT result → raw MB-local buffer → helper remap → caller copies
to `phi`.

**Fix:** Apply `norm_factor_` inside `RetrieveShearedResult` at the point of reading
from the FFT output array, before writing to the MB-local buffer. This makes the
remapper work on physically normalized data and matches the original
`RetrieveShearingResult` behavior. The remapper then needs no knowledge of
normalization.

For disk BC where `norm_factor_ == 1.0`, behavior is unchanged.

### 5. Phase 2 validation must include multi-MB/rank performance gate

The plan's phase order is: Phase 1 (extract) → Phase 2 (rank-aware) → Phase 3-4
(FFTGravity). Phase 2 makes `BlockFFTGravity` rank-aware but does not validate that
this is faster than the original code for the multi-MB/rank case.

The multi-MB/rank case is where the new helper might actually perform *worse* than the
old global-row `ExchangeShearingRows` (many small per-segment messages vs. one large
`MPI_Alltoallv`). If that regression exists, Phases 3-4 will inherit it.

**Required gate at end of Phase 2:**

- Run `mhd-4pc-b1-shear-nofb-fft-8mb` (32³, 8 MB/rank) with the new helper path.
- Compare `ShearSourceMax` and `RetrieveMax` against the current `blockfft` baseline
  and the old `fft-8mb` baseline.
- Do not proceed to Phase 3 if the helper is slower than the current global-row
  exchange for the 8 MB/rank case.

### 6. Phase 3 needs a quantitative checkpoint

Phase 3 ends with "Benchmark `ShearSourceMax`" but has no defined success criterion.

**Add explicit gate:** If `ShearSourceMax` does not drop to within 20% of
`BlockFFTGravity ShearSourceMax` (i.e., ≤ 24 ms, current blockfft is 19.7 ms, current
fft is 30.2 ms) after Phase 3, stop and profile the new helper's remap subphases
before starting Phase 4. Otherwise Phase 4 is done on faith.

### 7. Alternative architecture: integer-shift baked into the FFT load scatter (deferred)

Background for clarity: the shearing remap has two sub-steps — (a) a *fractional* shift
that moves each y-column by a small sub-cell amount (done by `RemapFluxPlm/Ppm` locally),
and (b) an *integer* shift that moves whole y-columns from one MeshBlock to another
(done by the current `MPI_Send`/`Recv` loop over `i`-columns). The integer shift is
what causes the per-column MPI traffic in `RollUnroll`.

The alternative idea is: instead of sending the integer-shifted data to the target
MeshBlock over MPI, one could skip that step and instead tell the AthenaFFT load
routine to *read* from a y-offset starting position when packing data into the global
FFT input buffer. Because `LoadShearedSource` already knows which global FFT cell each
MeshBlock cell maps to, it could add the integer-shift offset to that mapping. The
result is that the integer-shifted data arrives in the right place in the FFT buffer
directly at load time, with no separate MPI round.

Why this is not straightforward: the AthenaFFT scatter maps from MeshBlock-local
indices to pencil-decomposed FFT buffer indices. The integer-shift offset is different
for each i-column and may cross MeshBlock boundaries (i.e., the source row lands on a
different rank's FFT buffer). The scatter would need to handle cross-rank destinations
for individual rows, which is essentially the same MPI problem as the explicit integer
shift but hidden inside the scatter. Modifying AthenaFFT/fftMPI for this is invasive
and risky.

**Recommendation:** Do not pursue this alternative before Phase 3. Revisit only if
Phases 3–4 do not close the gap to within the Phase-gate targets above. No action
required.

### 8. Validation tolerances need explicit numbers

"Physically appropriate tolerances" is not actionable for a correctness gate.

**Recommended validation tolerances:**

- Serial single-rank: `max(|phi_new − phi_old|) / max(|phi_old|) < 1e-12` for first
  solve after restart.
- MPI N-rank: relative L∞ < 1e-10 (reduction-order non-determinism).
- Run ≥ 100 cycles, monitor `Egrav` in history file for drift exceeding 0.1% of the
  baseline run.
- Phase 1 `BlockFFTGravity` extraction: compare bit-exact against pre-refactor
  `blockfft` output at 1 MB/rank for ≥ 10 cycles at the same process count. Any bit
  difference indicates a bug in the extraction, not acceptable FP non-determinism.

### 9. Buffer ownership constraint: assert non-AMR in constructor

`RollUnrollAll` requires the helper to own all local block buffers for the full solve
duration. Memory is acceptable (~12 MB/rank at 8 × 64³ blocks). However, this
ownership model is invalid under load balancing or AMR regridding, which reassign
MeshBlocks between ranks mid-run.

**Add to `ShearingRemapper` constructor:**

```cpp
if (pm->adaptive) {
  std::stringstream msg;
  msg << "ShearingRemapper: AMR not supported.";
  ATHENA_ERROR(msg);
}
```

Add a `// TODO(AMR): requires re-init on load balance` comment near the buffer
allocation. Without this guard, a future AMR user will get silent corruption.

### 10. Phase 1 acceptance criterion: require bit-exact diff

The plan says "confirm blockfft benchmark unchanged." Benchmark timing alone does not
catch extraction bugs that corrupt physics slowly.

**Stronger criterion for Phase 1:** After extraction, run the same 10-cycle input
at the same MPI process count, redirect both pre- and post-refactor `phi` output to
HDF5, and confirm `max(|phi_pre − phi_post|) == 0.0` exactly. If extraction moved any
loop or index, bit identity will fail immediately rather than after a long benchmark
run.

---

## Resolved decisions

**Decision A — Helper location: `src/gravity/`** (confirmed)

Files will be placed at:

```text
src/gravity/shearing_remap.hpp
src/gravity/shearing_remap.cpp
```

If the helper later proves useful to non-gravity modules, it can be moved to
`src/orbital_advection/` at that time.

**Decision B — AthenaFFT integer-shift-via-scatter alternative: deferred**

See Issue 7 above for a plain-language explanation. The alternative is not well-suited
as a first approach; it is complex and may not actually avoid the cross-rank MPI work.
Revisit only if Phases 3–4 do not hit their performance targets.

---

## Branch and commit workflow

Start the implementation on a fresh branch that includes all current uncommitted
changes from the `optimize-athenafft-shear-column` branch. This keeps the profiling
instrumentation and timing analysis work intact while the refactor starts cleanly.

```bash
# From optimize-athenafft-shear-column, with current changes already staged or unstaged
git checkout -b shearing-remapper
# Commit current tracked changes (gravity src/hpp, timing script) as a baseline
git add src/gravity/fft_gravity.cpp src/gravity/fft_gravity.hpp \
        src/gravity/block_fft_gravity.cpp src/gravity/block_fft_gravity.hpp \
        vis/python/fft_gravity_timing.py
git commit -m "Carry over profiling instrumentation and shear remap micro-opts as baseline"
```

Make a commit at the end of each Phase (and at major sub-steps within a Phase).
Suggested commit points:

| Commit point | Suggested message prefix |
|---|---|
| Phase 1: `ShearingRemapper` file created, `BlockFFTGravity` calls helper | `Extract BlockFFTGravity shearing remap into ShearingRemapper` |
| Phase 1: bit-exact validation confirmed | `Add validation result: Phase 1 bit-exact vs pre-refactor blockfft` |
| Phase 2: rank-aware communication + local copy path | `Add rank-aware communication and local copy to ShearingRemapper` |
| Phase 2: multi-MB/rank benchmark confirmed | `Phase 2 validated: no regression at 8 MB/rank` |
| Phase 3: `LoadShearedSource` + source path wired | `Use ShearingRemapper for FFTGravity source remap` |
| Phase 4: `RetrieveShearedResult` + result path wired | `Use ShearingRemapper for FFTGravity result remap` |
| Phase 5: old FFT-block remap code removed | `Remove old FFTGravity global-row shearing remap` |

Commit after each phase even if the next phase starts immediately. This gives a clear
rollback point if a benchmark or correctness check fails.

---

## Recommended phase-gate criteria

| Phase | Exit criterion |
|-------|---------------|
| 1 | Bit-exact `phi` diff vs pre-refactor `blockfft` at 1 MB/rank, ≥10 cycles |
| 2 | `BlockFFTGravity` timing unchanged at 1 MB/rank; helper not slower than old `fft-8mb` for 8 MB/rank `ShearSourceMax` |
| 3 | `FFTGravity ShearSourceMax` within 20% of `blockfft ShearSourceMax` (≤24 ms vs 19.7 ms baseline) |
| 4 | `FFTGravity RetrieveMax` within 20% of `blockfft RetrieveMax` (≤13 ms vs 11.2 ms baseline); full loop `SelfGravity` within 5% of `blockfft` |
| 5 | No regressions in `fft-8mb` case; old FFT-block remap code removed; profiling columns updated |
