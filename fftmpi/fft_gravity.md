# FFTGravity and ShearingRemapper

This document describes the current AthenaFFT gravity layer in the TIGRIS checkout at
`$HOME/tigris`: `FFTGravityDriver`, `FFTGravity`, and the `ShearingRemapper` helper used
for shearing-periodic gravity. The lower-level FFT wrapper stack is documented in
[`athena_fft_classes.md`](athena_fft_classes.md).

---

## 1. Files and Class Hierarchy

**Files**:

- `src/gravity/fft_gravity.hpp`
- `src/gravity/fft_gravity.cpp`
- `src/gravity/shearing_remap.hpp`
- `src/gravity/shearing_remap.cpp`

```text
FFTGravityDriver : public FFTDriver     (one per simulation)
  ├── FFTGravity : public FFTBlock      (one per MPI rank; cuboid FFT block)
  └── ShearingRemapper                  (only for shearing-periodic meshes)
```

`FFTGravity` inherits `FFTBlock`'s cuboid multi-MeshBlock/rank layout and overrides
`ApplyKernel`, `ExecuteForward`, and `ExecuteBackward` for gravity boundary conditions.
`FFTGravityDriver::Solve(stage, mode, gas_only)` owns the high-level flow: load source,
run FFTs, apply the kernel or Green's function, retrieve the result, then apply gravity
boundary tasks.

`BlockFFTGravity` has been restored to the `origin/tigris-master` implementation and is
kept as the legacy 1 MeshBlock/rank path. Shared gravity boundary helpers are duplicated
in `block_fft_gravity.hpp` and `fft_gravity.hpp`; they are intentionally not declared in
`gravity.hpp`.

---

## 2. Boundary Conditions

| BC | Selection | Algorithm | Extra data |
|----|-----------|-----------|------------|
| Periodic | `gbflag == periodic`, non-shearing mesh | `LoadSource` → forward FFT → `ApplyKernel(mode)` → backward FFT → `RetrieveResult` | none |
| Open | `gbflag == open` | 8-parity convolution with a cell-averaged Green's function on a 2× domain | `grf_`, `pf3dgrf_` |
| Disk | `gbflag == disk` | Direct physical-layout source load, dedicated z-pencil transform, even/odd vertical split, disk kernel, inverse transform | `in_e_`, `in_o_`, `pf3d_disk_` or serial FFTW plans |
| Shearing-periodic | `Mesh::shear_periodic` | MeshBlock-local roll/unroll through `ShearingRemapper`, then normal AthenaFFT solve with shearing kernel | `premapper_`, `qomt_` |

Open BC is rejected for shearing-periodic meshes in `FFTGravityDriver` because those
boundary conditions are incompatible.

---

## 3. FFTGravityDriver

### Construction

`FFTGravityDriver(Mesh *pm, ParameterInput *pin)`:

1. Reads `four_pi_G_`, gravity profiling options, and `grav_bc`.
2. Creates `ShearingRemapper` only when `pmy_mesh_->shear_periodic` is true.
3. Creates one `FFTGravity` as `pmy_fb`, using the cuboid FFT block location and size from
   `FFTDriver`.
4. Sets normalization:
   - periodic and shearing: `four_pi_G_/gcnt_`
   - open: ignored by `RetrieveOBCResult`, which applies `1/(8*gcnt_)`
   - disk: `1.0`, with disk normalization handled in `ExecuteBackward`
5. Calls `QuickCreatePlan()`.
6. Initializes BC-specific structures: `InitGreen(four_pi_G_)` for open BC or
   `InitDiskBC(four_pi_G_)` for disk BC.
7. Creates `GravityBoundaryTaskList`.

### ShearTimeShift

`ShearTimeShift(stage)` computes elapsed time from the nearest shearing-periodic remap time:

```cpp
time_int = mesh time + ebeta[stage-1]*dt
qomL = qshear_*Omega_0_*Lx1_
tn = n*Lx2_/qomL
return time_int - tn
```

The result is used as `dt` for roll/unroll. Source loading uses `-dt`; result retrieval
uses `+dt` to undo the coordinate transform.

---

## 4. FFTGravity Members

| Member | Meaning |
|--------|---------|
| `gbflag` | Gravity boundary condition: `periodic`, `disk`, or `open` |
| `grf_` | k-space Green's function buffer for open BC |
| `in_e_`, `in_o_` | Even/odd slow-pencil buffers for disk BC |
| `I_` | Complex unit |
| `four_pi_G_`, `Lx3_`, `qomt_` | Disk/shearing physical parameters |
| `slow_nx1_`, `slow_nx2_`, `slow_nx3_` | Local dimensions of the disk slow-pencil layout |
| `slow_ilo_`, `slow_jlo_`, `slow_klo_` | Global offsets of that slow-pencil layout |
| `pf3dgrf_` | MPI-only 2×-domain FFT3d for the open-BC Green's function |
| `pf3d_disk_` | MPI-only disk-BC transform from physical block layout to z-pencil layout |
| `permute_disk_` | Fixed to 2 for z-pencil output in `pf3d_disk_` |
| `fplan_xy_`, `bplan_xy_`, `fplan_z_`, `bplan_z_` | Serial FFTW plans for disk BC |

`SetShearQuantities(qomt)` stores `q*Omega*dt` for the shearing-periodic kernel.

---

## 5. Periodic Path

For non-shearing periodic gravity, `FFTGravityDriver::Solve()` uses the base `FFTBlock`
load/retrieve methods:

1. For every local MeshBlock, load gas density or gas+particle density with `LoadSource`.
2. `ExecuteForward()`.
3. `ApplyKernel(mode)`, where `mode=0` uses the discrete finite-difference kernel and
   `mode=1` uses the continuous Poisson kernel.
4. `ExecuteBackward()`.
5. Retrieve into `phi` or `phi_gasonly` with `RetrieveResult`.

`ApplyKernel()` stores the kernel-multiplied result in `in_`, which is the input to the
backward FFT.

---

## 6. Open Boundary Path

Open BC uses an 8-parity convolution to represent isolated boundaries with periodic FFTs.

### InitGreen

`InitGreen(four_pi_G)`:

1. Creates `pf3dgrf_` for a `2*Nx[0] × 2*Nx[1] × 2*Nx[2]` domain with the same output
   permutation as the main forward plan.
2. Allocates `grf_` large enough for the 2× brick input and all fftMPI intermediate
   pencil layouts.
3. Fills the real-space cell-averaged Green's function using `_GetIGF`, folded into the
   `[-Nx, Nx-1]` range.
4. Forward-FFTs `grf_` in place.

### Source and Result Indexing

`LoadOBCSource(src, nu, ngh, loc, bsize, px, py, pz)` is called once per MeshBlock and per
parity pass. It computes the MeshBlock's zero-based offset inside the owning `FFTBlock`:

```cpp
is = loc.lx1*bsize.nx1 - loc_.lx1*bsize_.nx1
js = loc.lx2*bsize.nx2 - loc_.lx2*bsize_.nx2
ks = loc.lx3*bsize.nx3 - loc_.lx3*bsize_.nx3
```

The local FFTBlock coordinates `mi,mj,mk` are used for `GetIndex(..., f_in_)`, while
global full-domain indices `gi,gj,gk = FFTBlock offset + mi,mj,mk` are used in the parity
phase:

```cpp
phase = PI*(gi*px/Nx[0] + gj*py/Nx[1] + gk*pz/Nx[2])
in_[idx] = rho * exp(-i*phase)
```

`MultiplyGreen(px, py, pz)` multiplies `out_` by `grf_` at the corresponding 2× parity
offset. `RetrieveOBCResult()` uses the same MeshBlock-to-FFTBlock mapping and accumulates
all eight parity passes into `phi` with normalization `1/(8*gcnt_)`.

---

## 7. Disk Boundary Path

Disk BC is decomposition-independent by avoiding assumptions about the standard
`fplan_->pf3d` output layout.

1. `LoadDiskSource()` loads directly into `orig_idx_`, the physical x-fast layout.
2. `InitDiskBC()` creates a dedicated transform:
   - MPI: `pf3d_disk_` with `permute_disk_=2`, producing z-pencil output.
   - Serial: batched FFTW plans for xy and z transforms.
3. `ExecuteForward()` dispatches to the disk path when `gbflag == disk`, remapping to the
   slow-pencil layout and splitting vertical parity into `in_e_` and `in_o_`.
4. `ApplyKernel()` calls `ComputeDiskKernel()` for each global mode.
5. `ExecuteBackward()` combines even/odd pieces, runs the inverse disk transform, and
   applies disk normalization.

The disk kernel uses global slow-pencil offsets (`slow_ilo_`, `slow_jlo_`, `slow_klo_`) so
the answer does not depend on the MPI process grid.

---

## 8. Shearing-Periodic Path

The current shearing path no longer assembles global density vectors or uses
`MPI_Allreduce`. It keeps data MeshBlock-local and uses `ShearingRemapper` to apply the
same roll/unroll concept originally implemented in `BlockFFTGravity::RollUnroll()`.

### Solve Flow

For `pmy_mesh_->shear_periodic`:

1. Compute `dt = ShearTimeShift(stage)`.
2. Store `qomt_ = qshear_*Omega_0_*dt` in `FFTGravity`.
3. Copy each MeshBlock's density into `premapper_->Buffer(pmb)` in `(k,i,j)` order. If
   particle gravity is active and `gas_only` is false, add particle density first.
4. Call `premapper_->RollUnrollAll(-dt)` to transform source density to shearing
   coordinates.
5. For each MeshBlock, call `FFTGravity::LoadShearedSource()` to load the remapped buffer
   into `in_` with `GetIndex(..., f_in_)`.
6. Run `ExecuteForward()`, `ApplyKernel(mode)`, and `ExecuteBackward()`.
7. Retrieve raw inverse FFT output into each remap buffer with `RetrieveShearedResult()`.
8. Call `premapper_->RollUnrollAll(dt)` to transform potential back to physical
   coordinates.
9. Copy the remapped buffers to `phi` or `phi_gasonly`.

`LoadShearedSource()` and `RetrieveShearedResult()` use the same cuboid offset calculation
as `FFTBlock::LoadSource()`, but their external buffer order is `(k,i,j)` to match the
shearing remap helper.

---

## 9. ShearingRemapper

`ShearingRemapper` owns one `BlockBuffer` per local MeshBlock:

| Member | Meaning |
|--------|---------|
| `pmb` | Owning MeshBlock |
| `dat` | Active field in `(k,i,j)` layout, including y ghost slots |
| `roll_buf` | Fractionally shifted field before integer y remap |
| `pflux` | Scratch flux array for PLM/PPM orbital remap |

The helper lazily allocates buffers through `InitializeLocalBuffers()` and exposes them
through `Buffer(local_id)` or `Buffer(MeshBlock*)`.

### RollUnrollAll

`RollUnrollAll(dt)` performs three phases:

1. `FillGhostZonesAll()`
2. `FractionalShiftAll(dt)`
3. `IntegerShiftAll(dt)`

The sign of `dt` selects roll vs. unroll. The source path uses `-dt`; the result path uses
`+dt`.

### FillGhostZonesAll

The fractional shift needs neighboring y cells. For every local destination MeshBlock,
`FillGhostZonesAll()` fills:

- lower-y ghost zones from the active upper-y edge of the previous y block
- upper-y ghost zones from the active lower-y edge of the next y block

Same-rank neighbors are copied directly. Remote transfers are expressed as `RemapSegment`
records and packed into one message per remote rank. `SegmentLess()` provides a
deterministic order so sender and receiver agree on the packed layout without sending
per-segment metadata.

### FractionalShiftAll

For each active `(k,i)` column:

```cpp
yshear = -qshear_*Omega_0_*x1v(i)*dt
joffset = ceil(yshear/dx2)
eps = fmod(yshear, dx2)/dx2
```

The code calls `RemapFluxPlm()` or `RemapFluxPpm()` depending on reconstruction order, then
stores:

```cpp
roll_buf(k,i,j) = dat(k,i,j) - (pflux(j+1) - pflux(j))
```

This is the fractional part of the old `BlockFFTGravity::RollUnroll()` algorithm.

### IntegerShiftAll

`IntegerShiftAll(dt)` applies the whole-cell y shift. For each x-column the shift is split
into block crossings and overlap cells:

```cpp
joffset = int(yshear/dx2)
Ngrids = abs(joffset)/nx2
joverlap = abs(joffset) - Ngrids*nx2
```

Positive-y shear uses the A.1/A.2 segment decomposition:

- A.1: `[je-(joverlap-1):je] → [js:js+(joverlap-1)]`
- A.2: `[js:je-joverlap] → [js+joverlap:je]`

Negative-y shear uses the B.1/B.2 decomposition:

- B.1: `[js:js+(joverlap-1)] → [je-(joverlap-1):je]`
- B.2: `[js+joverlap:je] → [js:je-joverlap]`

Same-rank segments copy directly from `roll_buf` to `dat`. Remote segments are batched by
rank using the same deterministic packed-message scheme as ghost fill.

`RollUnrollBlock()` was removed after `BlockFFTGravity` was restored to its original
implementation; `ShearingRemapper` only keeps the all-local-MeshBlock path used by
`FFTGravity`.

---

## 10. Profiling

When `<gravity>/profile_gravity = true`, `FFTGravityDriver` writes
`<problem_id>.fft_gravity_time.txt`. The shearing path reports separate timings for:

- particle density accumulation
- source load
- shearing source remap
- forward FFT
- kernel
- backward FFT
- retrieval/unroll
- gravity boundary task list

`BlockFFTGravity` has its own legacy profiling output and is not routed through
`ShearingRemapper`.
