# FFT gravity shearing remap refactor plan

Date: 2026-05-02

Revision note: `FFT_GRAVITY_SHEARING_REMAP_REFACTOR_REVIEW.md` has been reviewed.
The review is accepted as design guidance. The plan below incorporates its main
changes: rank-aggregated MPI payloads from the start, explicit same-rank copy
correctness requirements, clarified `phi`/`phi_gasonly` result handling, normalization
placement, AMR guard, explicit correctness tolerances, and phase-gate criteria.
The plan also incorporates the May 2 benchmark convention: current outputs under
`/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc*shear*` are the reference
datasets, and future post-refactor runs will use the `-refactor` suffix for direct
timing and field comparisons.

## Summary

The remaining `fft_gravity` performance gap relative to `block_fft_gravity` is not in
AthenaFFT forward/backward transforms. It is in the shearing-coordinate remap around
the FFT solve.

The current `FFTGravity` implementation remaps the whole local FFT block through a
custom global-y-row exchange. That works for multiple MeshBlocks per rank, but it is
the wrong abstraction and is consistently slower than the MeshBlock-local remap used
by `BlockFFTGravity`.

The fundamental fix is to introduce a reusable, rank-aware, MeshBlock-local shearing
remap helper, then use it from both gravity solvers:

```text
BlockFFTGravity -> shared shearing remapper
FFTGravity      -> shared shearing remapper before/after AthenaFFT
```

The helper must support both:

```text
1 MeshBlock/rank
multiple MeshBlocks/rank
```

The `1 MeshBlock/rank` case should not be a special permanent solution. It may be
used temporarily for validation, but the final design should be general.

## Current evidence

Benchmark target:

```text
/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-*
```

Primary comparison:

```text
mhd-4pc-b1-shear-nofb-fft       # 64^3 MeshBlock, FFTGravity/AthenaFFT
mhd-4pc-b1-shear-nofb-blockfft  # 64^3 MeshBlock, BlockFFTGravity
```

Reference/output naming convention:

| Case | Reference output | Post-refactor output | Purpose |
|---|---|---|---|
| 64^3, 1 MeshBlock/rank, FFTGravity | `mhd-4pc-b1-shear-nofb-fft` | `mhd-4pc-b1-shear-nofb-fft-refactor` | Direct timing and `phi`/`phi_gas` correctness comparison |
| 64^3, 1 MeshBlock/rank, BlockFFTGravity | `mhd-4pc-b1-shear-nofb-blockfft` | optional `mhd-4pc-b1-shear-nofb-blockfft-refactor` | Performance target and Phase 1 extraction validation |
| 32^3, 8 MeshBlocks/rank, FFTGravity | current 8 MB/rank `fft` output | corresponding `fft-refactor` output | Multi-MeshBlock/rank regression gate |

Use `-refactor` exactly for new runs. Avoid alternate spellings so the comparison
scripts can pair old and new outputs automatically.

Build and submission command reference:

```bash
cd /home/changgoo/tigris_scripts/tigress_classic
./build_tigress_nofb.sh tiger mhd fft
./build_tigress_nofb.sh tiger mhd blockfft

cd /home/changgoo/tigris_scripts/tigress_classic/tiger
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd fft       # 64^3 FFTGravity
sbatch tigress_classic_mhd_shear_nofb.slurm -i mhd blockfft  # 64^3 BlockFFTGravity
sbatch tigress_classic_mhd_shear_nofb_8mb.slurm -i mhd fft   # 32^3, 8 MeshBlocks/rank FFTGravity
```

Recent full-run timing after instrumentation:

| Phase | fft mean | blockfft mean | blockfft/fft | Delta |
| --- | ---: | ---: | ---: | ---: |
| TotalMax | 0.14386894 | 0.12943773 | 0.8997 | -0.014431211 |
| ShearSourceMax | 0.030190929 | 0.019654601 | 0.6510 | -0.010536328 |
| RetrieveMax | 0.015436346 | 0.011201457 | 0.7257 | -0.004234889 |
| ForwardMax | 0.054657691 | 0.052813018 | 0.9663 | -0.001844673 |
| KernelMax | 0.0034896872 | 0.003513006 | 1.0067 | +0.000023319 |
| BackwardMax | 0.052616581 | 0.053451847 | 1.0159 | +0.000835267 |

Important observations:

- Forward/kernel/backward FFT costs are already comparable.
- `ShearSourceMax` is the dominant gap.
- `RetrieveMax` is the secondary gap.
- Offset and row-index construction are negligible after earlier caching changes.
- Attempts to optimize `FFTGravity::FillShearingGhosts()` locally, including
  contiguous slab copies and nonblocking two-direction MPI exchange, did not improve
  the full benchmark.

Conclusion:

The problem is the current `FFTGravity` shearing remap architecture, not a small
local loop implementation detail.

## Current architecture

### BlockFFTGravity

In `src/gravity/block_fft_gravity.cpp`, the shearing path does:

```text
copy MeshBlock rho/rhosum into roll_var(k,i,j)
RollUnroll(roll_var, -dt)
copy roll_var active zone into FFT buffer
ExecuteForward()
ApplyKernel()
ExecuteBackward()
copy inverse FFT result into roll_var(k,i,j)
RollUnroll(roll_var, +dt)
copy roll_var active zone into phi
GravityBoundaryTaskList
```

`BlockFFTGravity::RollUnroll()` operates on a MeshBlock-shaped buffer and performs:

1. y ghost fill for fractional interpolation
2. fractional shearing shift using `OrbitalAdvection::RemapFluxPlm/Ppm`
3. integer y shift using MeshBlock target locations

Limitation:

`RollUnroll()` currently uses `pleaf->GetGid()` as the MPI rank in places like:

```cpp
sendto_id = pleaf->GetGid();
MPI_Send(..., sendto_id, MPI_COMM_WORLD);
```

That is only correct when:

```text
global MeshBlock id == MPI rank
```

So the current implementation is naturally tied to `1 MeshBlock/rank`.

### FFTGravity

In `src/gravity/fft_gravity.cpp`, the shearing path currently does:

```text
for each local MeshBlock:
    LoadSource(rho/rhosum into FFT block)

ApplyShearingSource(-1.0)
ExecuteForward()
ApplyKernel(mode)
ExecuteBackward()
ApplyShearingResult(+1.0)

for each local MeshBlock:
    RetrieveAppliedShearingResult(phi)

GravityBoundaryTaskList
```

The expensive operations are in:

```text
FFTGravity::ApplyShearingSource()
FFTGravity::ApplyShearingResult()
FFTGravity::RetrieveAppliedShearingResult()
FFTGravity::FillShearingGhosts()
FFTGravity::ExchangeShearingRows()
```

This path remaps the whole local FFT block and uses custom global y-row exchanges.
That is the abstraction to remove.

## Desired architecture

Introduce a reusable shearing remap component that operates on MeshBlock-local buffers
and uses rank-aware MeshBlock communication.

Proposed initial location:

```text
src/gravity/shearing_remap.hpp
src/gravity/shearing_remap.cpp
```

Rationale:

- The immediate consumers are gravity solvers.
- It can reuse gravity-specific buffer shapes without disturbing core orbital
  advection code.
- If the abstraction becomes broadly useful, it can later move under
  `src/orbital_advection/`.

Alternative location:

```text
src/orbital_advection/shearing_remap.*
```

This is conceptually attractive because the remap uses `OrbitalAdvection::RemapFlux*`,
but it risks coupling orbital-advection internals to gravity-specific buffer ownership.
Start in `src/gravity/` unless reviewers prefer otherwise.

## Proposed API

The helper must coordinate all local MeshBlocks in one call. A one-block-at-a-time
API is not sufficient because the integer shift can target a different MeshBlock on
the same rank, and the helper needs access to both source and destination buffers.

Recommended class shape:

```cpp
class ShearingRemapper {
 public:
  ShearingRemapper(Mesh *pm, ParameterInput *pin);

  void InitializeLocalBuffers();
  AthenaArray<Real> &Buffer(int local_id);
  AthenaArray<Real> &Buffer(MeshBlock *pmb);

  // Remap all local block buffers together. Each buffer is shaped as (k, i, j),
  // matching BlockFFTGravity::roll_var.
  void RollUnrollAll(Real dt);

 private:
  struct BlockBuffer {
    MeshBlock *pmb;
    AthenaArray<Real> dat;      // (k, i, j), includes y ghosts
    AthenaArray<Real> roll_buf; // same active-zone shape for simplicity
    AthenaArray<Real> pflux;
  };

  struct RemapSegment {
    int src_gid;
    int dst_gid;
    int src_i;
    int dst_i;
    int src_j_start;
    int dst_j_start;
    int count_j;
    int count_k;
    int src_rank;
    int dst_rank;
  };

  Mesh *pmy_mesh_;
  std::vector<BlockBuffer> blocks_;
  std::vector<std::vector<RemapSegment>> send_segments_by_rank_;
  std::vector<std::vector<RemapSegment>> recv_segments_by_rank_;
  std::vector<std::vector<Real>> send_payloads_by_rank_;
  std::vector<std::vector<Real>> recv_payloads_by_rank_;
};
```

Constructor guard:

```cpp
if (pm->adaptive) {
  std::stringstream msg;
  msg << "ShearingRemapper: AMR not supported.";
  ATHENA_ERROR(msg);
}
```

Add a nearby TODO:

```cpp
// TODO(AMR): reinitialize buffers and communication plans after load balancing.
```

## Generalized communication model

The central change is to replace `gid == rank` assumptions with `ranklist/nslist`.

Given a target logical location:

```cpp
MeshBlockTree *pleaf = pmy_mesh_->tree.FindMeshBlock(target_loc);
int target_gid = pleaf->GetGid();
int target_rank = pmy_mesh_->ranklist[target_gid];
int target_lid = target_gid - pmy_mesh_->nslist[target_rank];
```

For same-rank targets:

```cpp
if (target_rank == Globals::my_rank) {
  MeshBlock *target = pmy_mesh_->FindMeshBlock(target_gid);
  // copy directly into target block's remap buffer
}
```

For remote targets:

```cpp
MPI_Isend(..., target_rank, tag, MPI_COMM_WORLD, ...);
MPI_Irecv(..., source_rank, tag, MPI_COMM_WORLD, ...);
```

The hard part is not just replacing rank ids. The integer-shift remap moves columns
between MeshBlocks. For same-rank target blocks, the helper needs access to the
target buffer that will receive the shifted data.

This suggests the final helper should remap all local MeshBlocks in a coordinated
call, not only one block at a time.

## Caller flow

The caller owns physics-specific loading/storing. The helper owns the remap buffers
and communication.

```cpp
for each local MeshBlock:
    fill remapper.Buffer(pmb)

remapper.RollUnrollAll(dt)

for each local MeshBlock:
    consume remapper.Buffer(pmb)
```

This solves the same-rank target issue because the remapper owns or indexes every
local block's remap buffer during the integer shift.

## Implementation details

### Buffer shape

Use the same shape as `BlockFFTGravity`:

```text
dat(k, i, j)
```

Allocate with:

```cpp
dat.NewAthenaArray(ncells3, ncells1, ncells2);
```

This includes ghost cells in all dimensions. The remapper only requires active
`k/i` and y ghost cells.

### Fractional shift

Keep the existing method:

```cpp
yshear = -qshear * Omega0 * pmb->pcoord->x1v(i) * dt;
joffset = static_cast<int>(std::ceil(yshear/dx2));
eps = std::fmod(yshear, dx2)/dx2;
osgn = (joffset > 0) ? 1 : 0;
shift0 = (joffset > 0) ? -1 : 0;

if (xorder <= 2)
  porb->RemapFluxPlm(...);
else
  porb->RemapFluxPpm(...);

roll_buf(k,i,j) = dat(k,i,j) - (pflux(j+1) - pflux(j));
```

Do not change interpolation behavior in the first refactor.

### y ghost fill

Current `BlockFFTGravity::RollUnroll()` always fills y ghosts, despite the comment
saying one direction can skip it. Preserve behavior first.

For a coordinated multi-block helper:

1. For each local block, identify upper/lower y neighbor logical locations.
2. If neighbor is local, copy from neighbor buffer directly.
3. If neighbor is remote, exchange with neighbor rank.

Use rank-aware target lookup:

```cpp
int gid = pleaf->GetGid();
int rank = pmy_mesh_->ranklist[gid];
```

Tags need to distinguish:

```text
operation: ghost fill vs integer shift
source gid
destination gid
direction
possibly i-column for integer shift
```

Do not rely on a single tag `0` once multiple blocks per rank are supported.

Use `Mesh::ReserveTagPhysIDs()` if possible, or follow existing boundary variable tag
patterns.

### integer shift

This is the hardest part to generalize cleanly.

Current `BlockFFTGravity::RollUnroll()` loops over each active `i` and sends two
segments depending on `joffset` and `joverlap`.

For multi-block support:

- A source segment may target another block on the same rank.
- Multiple local source blocks may target the same remote rank.
- Multiple i-columns require different target blocks.

Implementation strategy:

1. Preserve the exact segment decomposition from `BlockFFTGravity::RollUnroll()`.
2. Replace immediate blocking send/receive with an exchange plan:
   - build list of send segments
   - build list of expected receive segments
   - post receives
   - pack and post sends
   - wait
   - unpack
3. For local targets, apply the copy directly without MPI.

Segment structure:

```cpp
struct RemapSegment {
  int src_gid;
  int dst_gid;
  int src_i;
  int dst_i;
  int src_j_start;
  int dst_j_start;
  int count_j;
  int count_k;
  int src_rank;
  int dst_rank;
};
```

Because the shift is along y only, `src_i == dst_i` in logical MeshBlock-local
coordinates for same-size blocks.

Pack order can match current code:

```text
for k:
  for j:
    buffer[counter++] = roll_buf(k,i,j)
```

### local copy path

When `dst_rank == Globals::my_rank`, copy directly from source block's `roll_buf`
to destination block's `dat`.

This is the behavior `BlockFFTGravity` currently misses when there are multiple
MeshBlocks per rank.

### MPI path

Do not use `target_gid` as MPI rank. Use `ranklist[target_gid]`.

Start with rank-aggregated payloads. Do not implement a per-segment MPI message path.
With 8 MeshBlocks/rank, per-segment messages can create hundreds of small messages
per remap and introduce MPI tag pressure.

Required design:

1. Reserve a base physics tag with `Mesh::ReserveTagPhysIDs()`.
2. Build send segment metadata grouped by destination rank.
3. Exchange metadata counts per rank.
4. Send metadata arrays per rank.
5. Allocate receive payloads from metadata.
6. Pack one concatenated payload per destination rank.
7. Post all receives before sends.
8. Post all sends.
9. `MPI_Waitall`.
10. Unpack by segment metadata.

This is more code than the current `BlockFFTGravity::RollUnroll()` communication,
but it is the correct design for multiple MeshBlocks/rank and avoids tag collisions.

## FFTGravity integration

### Source path

Replace this shearing source flow:

```cpp
for each local MeshBlock:
    LoadSource(rho/rhosum)

pfg->ApplyShearingSource(-1.0);
```

with:

```cpp
Real dt = ShearTimeShift(stage);

for each local MeshBlock:
    AthenaArray<Real> &buf = remapper.Buffer(pmb);
    fill buf(k,i,j) from rho/rhosum(k,j,i), including y ghosts

remapper.RollUnrollAll(-dt);

for each local MeshBlock:
    pfg->LoadShearedSource(remapper.Buffer(pmb), pmb->loc, pmb->block_size);
```

Add a new `FFTGravity` method:

```cpp
void FFTGravity::LoadShearedSource(const AthenaArray<Real> &src_kij,
                                   LogicalLocation loc, RegionSize bsize);
```

This writes directly to the FFT input buffer in physical layout:

```cpp
in_[GetIndex(mi, mj, mk, f_in_)] = {src_kij(k, i, j), 0.0};
```

where `src_kij` indices are `(k,i,j)`.

### Result path

Replace this flow:

```cpp
pfg->ApplyShearingResult(1.0);
for each local MeshBlock:
    pfg->RetrieveAppliedShearingResult(phi)
```

with:

```cpp
for each local MeshBlock:
    pfg->RetrieveShearedResult(remapper.Buffer(pmb), pmb->loc, pmb->block_size);

remapper.RollUnrollAll(+dt);

for each local MeshBlock:
    copy remapper.Buffer(pmb)(k,i,j) into phi(k,j,i) or phi_gasonly(k,j,i)
```

Add:

```cpp
void FFTGravity::RetrieveShearedResult(AthenaArray<Real> &dst_kij,
                                       LogicalLocation loc, RegionSize bsize);
```

This reads raw inverse FFT result into `(k,i,j)` buffer and applies normalization
at the same point as the current code:

```cpp
dst_kij(k,i,j) = std::real(out_[GetIndex(mi,mj,mk,b_out_)]) * norm_factor_;
```

For disk BC, current shearing path uses `norm_factor_ == 1.0`; preserve this.

`phi` versus `phi_gasonly` handling:

- The solver computes only one gravity potential for a given `Solve(..., gas_only)`
  call.
- If `gas_only == false`, retrieve into the remapper buffers, remap once, and copy
  into `phi`.
- If `gas_only == true`, retrieve into the remapper buffers, remap once, and copy
  into `phi_gasonly`.
- Do not call `RollUnrollAll(+dt)` twice for the same inverse FFT result.

Normalization rule:

- Apply `norm_factor_` inside `RetrieveShearedResult()` when reading from the FFT
  output array into the MeshBlock-local remapper buffer.
- The remapper should operate only on physically normalized data and must not know
  about FFT normalization.

### Profiling

Keep current high-level columns:

```text
ShearSourceMax
RetrieveMax
```

Add new subphase names for the new helper:

```text
MBRemapGhost
MBRemapFractional
MBRemapIntegerLocal
MBRemapIntegerMPI
MBRemapPack
MBRemapUnpack
```

or keep simpler first:

```text
MeshBlockShearSourceRemap
MeshBlockShearResultRemap
```

Current subphase columns from the old FFT-block path should remain only while the
fallback exists.

## BlockFFTGravity integration

Phase 1 commits to the final multi-block API from the start. Do not ship a
single-block intermediate API in Phase 1 and reshape it in Phase 2 — that would
partially undo the "separate extraction risk from generalization risk" purpose.

The public `BlockFFTGravity::RollUnroll(dat, dt)` method is preserved so existing call
sites compile, but its body copies into and out of the remapper-owned buffer:

```cpp
void BlockFFTGravity::RollUnroll(AthenaArray<Real> &dat, Real dt) {
  auto &buf = remapper_.Buffer(pmy_block_);
  // copy caller-owned dat into remapper buffer
  for (int k=ks; k<=ke; k++)
    for (int i=is; i<=ie; i++)
      for (int j=js-NGHOST; j<=je+NGHOST; j++)
        buf(k,i,j) = dat(k,i,j);
  remapper_.RollUnrollAll(dt);
  // copy result back into caller-owned dat
  for (int k=ks; k<=ke; k++)
    for (int i=is; i<=ie; i++)
      for (int j=js; j<=je; j++)
        dat(k,i,j) = buf(k,i,j);
}
```

The extra copies are trivial at 1 MeshBlock/rank and avoid an API rewrite in Phase 2.

After both gravity solvers migrate to the helper, consider deleting the wrapper.

## Backward compatibility and fallback

During migration, keep the old FFT-block shearing path behind a runtime or compile-time
fallback:

```text
gravity/use_meshblock_shearing_remap = true/false
```

Default recommendation during development:

```text
false
```

After correctness and benchmark validation:

```text
true
```

The fallback lets us compare:

```text
old FFT-block remap
new MeshBlock-local remap
blockfft
```

## Validation plan

### Correctness

Run at least:

1. 64^3, 1 MeshBlock/rank, shearing, disk BC:
   - old FFT-block remap vs new MeshBlock remap
   - compare `phi` and `phi_gas`/`phi_gasonly` fields directly between:

```text
/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft
/scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-refactor
```

2. 32^3, 8 MeshBlocks/rank, shearing, disk BC:
   - ensure new rank-aware remap works with same-rank target blocks
   - compare `ShearSourceMax` and `RetrieveMax` against both the current
     `blockfft` baseline and the old `fft-8mb` global-row remap baseline
3. blockfft path:
   - ensure helper preserves current blockfft results

Desired checks:

```text
max(abs(phi_refactor - phi_reference))
L2(phi_refactor - phi_reference)
max(abs(phi_gas_refactor - phi_gas_reference))
L2(phi_gas_refactor - phi_gas_reference)
history/self-gravity consistency
no conservation regressions in short run
```

Explicit tolerances:

- Phase 1 extraction of `BlockFFTGravity::RollUnroll()` must be bit-exact at
  1 MeshBlock/rank for at least 10 cycles:

```text
max(|phi_pre - phi_post|) == 0.0
```

- Serial/single-rank `FFTGravity` validation:

```text
max(|phi_new - phi_old|) / max(|phi_old|) < 1e-12
max(|phi_gas_new - phi_gas_old|) / max(|phi_gas_old|) < 1e-12
```

- MPI validation:

```text
relative L_inf < 1e-10
```

Apply the MPI relative L_inf tolerance separately to `phi` and `phi_gas`/`phi_gasonly`.

- Run at least 100 cycles and monitor `Egrav` in the history file. Drift relative to
  the baseline should stay below 0.1%.

### Performance

Use:

```bash
python vis/python/fft_gravity_timing.py \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-fft-refactor \
  /scratch/gpfs/EOST/changgoo/tigress_classic/mhd-4pc-b1-shear-nofb-blockfft
```

For the 32^3, 8 MeshBlocks/rank regression gate, compare the current 8 MB/rank
reference output against the matching `-refactor` output. The exact basename should
come from the revised Slurm scripts, but the comparison must include both directories
explicitly in the timing command.

Target:

```text
FFTGravity ShearSourceMax should approach BlockFFTGravity ShearSourceMax.
FFTGravity RetrieveMax should approach BlockFFTGravity RetrieveMax.
Forward/Kernel/Backward should remain comparable.
```

Current gap to eliminate:

```text
ShearSourceMax: ~10 ms
RetrieveMax:    ~4 ms
```

## Suggested implementation phases

### Phase 1: Extract behavior without changing communication

- Create `src/gravity/shearing_remap.hpp/cpp`.
- Move the core of `BlockFFTGravity::RollUnroll()` into helper with minimal changes.
- Keep `gid == rank` behavior for this phase.
- Convert `BlockFFTGravity` to call helper.
- Confirm blockfft output is bit-exact versus pre-refactor at 1 MeshBlock/rank for
  at least 10 cycles.

Purpose:

Separate extraction risk from generalization risk.

### Phase 2: Add rank-aware communication

- Replace `pleaf->GetGid()` as MPI rank with `ranklist[target_gid]`.
- Add local-copy path when target rank is current rank.
- Add rank-aggregated metadata and payload exchange.
- Validate with multiple MeshBlocks/rank.
- Gate: run the 32^3, 8 MeshBlocks/rank case and compare `ShearSourceMax` and
  `RetrieveMax` against the current `blockfft` baseline and the old `fft-8mb`
  global-row remap baseline. Do not proceed to Phase 3 if the helper is slower
  than the old `fft-8mb` path for either remap phase.

Purpose:

Make the helper generally correct.

### Phase 3: Use helper in FFTGravity source path

- Add `FFTGravity::LoadShearedSource()`.
- In shearing solve, fill MeshBlock-local remap buffers from rho/rhosum.
- Call helper with `-dt`.
- Load sheared buffers into FFT input.
- Keep result path old initially.
- Benchmark `ShearSourceMax`.
- For the 64^3 case, write the new run to the `-refactor` output directory so timing
  and field comparisons can be paired with the current reference output.
- Gate: stop after Phase 3 if `FFTGravity ShearSourceMax` is not within 20% of
  `BlockFFTGravity ShearSourceMax`. With the current baseline, the target is
  approximately:

```text
ShearSourceMax <= 24 ms
```

- If this gate fails, profile the new helper subphases before implementing Phase 4.

Purpose:

Attack the largest gap first.

### Phase 4: Use helper in FFTGravity result path

- Add `FFTGravity::RetrieveShearedResult()`.
- Retrieve inverse FFT result to MeshBlock-local remap buffers.
- Call helper with `+dt`.
- Copy remapped buffer to `phi`/`phi_gasonly`.
- Compare 64^3 `phi` and `phi_gas`/`phi_gasonly` against the current
  `mhd-4pc-b1-shear-nofb-fft` reference output.
- Benchmark `RetrieveMax`.
- Gate: `FFTGravity RetrieveMax` should be within 20% of `BlockFFTGravity RetrieveMax`.
  With the current baseline, target:

```text
RetrieveMax <= 13 ms
```

- Full-loop `SelfGravity` should be within 5% of `blockfft`.

Purpose:

Remove the second major gap.

### Phase 5: Remove old FFT-block remap

After correctness and performance validation:

- Remove or deprecate:

```text
FFTGravity::ApplyShearingSource()
FFTGravity::ApplyShearingResult()
FFTGravity::RetrieveAppliedShearingResult()
FFTGravity::FillShearingGhosts()
FFTGravity::ExchangeShearingRows()
```

- Remove old subphase timers specific to global-row exchange.
- Keep profiling for the new helper.

## Risks and open questions

1. **MPI tag design**
   - Use rank-aggregated metadata+payload exchange from the start.
   - Do not rely on per-segment hashed tags.
   - `Mesh::ReserveTagPhysIDs()` exists at `src/mesh/mesh.hpp:334`; use it to
     reserve the physics tag IDs needed for ghost-fill and integer-shift messages.
   - The existing tag pattern is `BoundaryValues::CreateBvalsMPITag(lid, bufid, phys)`;
     follow the same convention. Tag reservation happens in the `ShearingRemapper`
     constructor.
   - The metadata exchange round (step 3 in the 10-step MPI protocol) may be
     eliminatable because every rank can derive its receive plan locally from the
     global shear time and mesh layout. Keep the explicit metadata round in the
     first implementation for clarity; eliminate it later if it shows up in profiling.

2. **Same-rank target copies**
   - Requires coordinated all-local-block remap, not isolated one-block calls.
   - The helper should own or index all local remap buffers.

3. **AMR/multilevel**
   - Current benchmark is no multilevel.
   - The helper can initially require uniform non-AMR mesh and assert/fallback for
     multilevel.

4. **Threading**
   - Current benchmark is MPI-only.
   - If OpenMP is used later, helper-owned shared buffers need thread-safe design.

5. **Numerical equivalence**
   - Reordering communication may change bitwise results.
   - Validate with tolerances and physical diagnostics.

6. **Code ownership**
   - The helper lives in `src/gravity/` initially, but it uses orbital advection remap
     methods. Reviewers may prefer a location under `src/orbital_advection/`.

7. **Buffer layout transpose cost**
   - The remapper uses `(k,i,j)` but `rho`/`rhosum`/`phi` are `(k,j,i)`.
   - Source-fill and result-store loops are real memory-traffic work.
   - Add `MBRemapLoad` and `MBRemapStore` to the subphase profiling set so a
     regression in those transpose loops does not hide inside `ShearSourceMax` or
     `RetrieveMax`.

8. **Phase 4 `SelfGravity` 5% gate**
   - Current values: `blockfft`=27.85s, `fft`=32.21s. The 5% target is ≤29.24s.
   - If Phase 4 hits its `RetrieveMax` target but misses the loop `SelfGravity`
     gate, this indicates non-FFT gravity work (boundary task list, etc.), not a
     Phase 4 bug. In that case: profile rather than roll back Phase 4.

9. **`pflux` per-block allocation**
   - `pflux` is a 1D scratch used in the fractional-shift inner loop, written and
     read in the same `(k,i)` iteration. One shared `pflux` across all blocks in
     `RollUnrollAll` is sufficient. Do not restructure the `BlockBuffer` struct for
     this now; note it for a later cleanup pass.

## Phase-gate criteria

| Phase | Exit criterion |
|-------|---------------|
| 1 | Bit-exact `phi` diff vs pre-refactor `blockfft` at 1 MeshBlock/rank, >=10 cycles |
| 2 | `BlockFFTGravity` timing unchanged at 1 MeshBlock/rank; helper not slower than old `fft-8mb` for 8 MeshBlocks/rank `ShearSourceMax` or `RetrieveMax`; compare both against current `blockfft` baseline |
| 3 | 64^3 `-refactor` run produced; `FFTGravity ShearSourceMax` within 20% of `blockfft` (`<= 24 ms` with current baseline) |
| 4 | 64^3 `phi` and `phi_gas`/`phi_gasonly` match the current `fft` reference within tolerance; `FFTGravity RetrieveMax` within 20% of `blockfft` (`<= 13 ms` with current baseline); full-loop `SelfGravity` within 5% of `blockfft` |
| 5 | No regression in `fft-8mb`; old FFT-block remap code removed; profiling columns updated |

## Branch and commit workflow

Start from the current optimization branch with profiling instrumentation present:

```bash
git checkout -b shearing-remapper
git add src/gravity/fft_gravity.cpp src/gravity/fft_gravity.hpp \
        src/gravity/block_fft_gravity.cpp src/gravity/block_fft_gravity.hpp \
        vis/python/fft_gravity_timing.py
git commit -m "Carry over profiling instrumentation and shear remap micro-opts as baseline"
```

Commit at the end of each phase. Suggested commit points:

| Commit point | Suggested message prefix |
|---|---|
| Phase 1 extraction | `Extract BlockFFTGravity shearing remap into ShearingRemapper` |
| Phase 1 validation | `Add validation result: Phase 1 bit-exact vs pre-refactor blockfft` |
| Phase 2 communication | `Add rank-aware communication and local copy to ShearingRemapper` |
| Phase 2 validation | `Phase 2 validated: no regression at 8 MB/rank` |
| Phase 3 source path | `Use ShearingRemapper for FFTGravity source remap` |
| Phase 4 result path | `Use ShearingRemapper for FFTGravity result remap` |
| Phase 5 cleanup | `Remove old FFTGravity global-row shearing remap` |

## Recommendation

Do not keep optimizing the current `FFTGravity` global-row shearing remap. The data
does not support that path.

Proceed with the shared MeshBlock-local, rank-aware shearing remapper. Extract first,
generalize second, then migrate `FFTGravity` source and result paths in separate
steps so each change can be benchmarked and validated independently.
