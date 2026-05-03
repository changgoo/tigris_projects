# AthenaFFT — Class Hierarchy, Index Mapping, and Cuboid Decomposition

This document describes the `athena_fft` wrapper stack: class structure, data members,
the axis-permutation / index-mapping machinery, the cuboid multi-meshblock layout,
and how the fftMPI C++ library is called.  It also documents `FFTGravity`, the
multi-meshblock Poisson solver built on top of `FFTBlock`.

---

## 1. Class Overview

```
FFTDriver                        (one per simulation; orchestrates)
  └── FFTBlock  pmy_fb           (one per MPI rank; the FFT domain)
        ├── AthenaFFTIndex  f_in_    forward FFT:  input  layout
        ├── AthenaFFTIndex  f_out_   forward FFT:  output layout
        ├── AthenaFFTIndex  b_in_    backward FFT: input  layout
        ├── AthenaFFTIndex  b_out_   backward FFT: output layout
        ├── AthenaFFTPlan*  fplan_   forward  plan (fftMPI or FFTW)
        └── AthenaFFTPlan*  bplan_   backward plan
              │
              └── FFTGravity  (gravity solver; extends FFTBlock)
```

All base classes are declared in `src/fft/athena_fft.hpp`; implemented in
`src/fft/athena_fft.cpp` and `src/fft/fft_driver.cpp`.  `FFTGravity` lives in
`src/gravity/fft_gravity.hpp` / `fft_gravity.cpp`.

---

## 2. FFTDriver

**File**: `src/fft/fft_driver.cpp`

### Purpose
Manages the global FFT mesh layout: how MPI ranks are mapped to FFT domains, computes the
cuboid grouping of MeshBlocks, and creates the single `FFTBlock` owned by each rank.

### Key public members

| Member | Type | Meaning |
|--------|------|---------|
| `npx1, npx2, npx3` | `int` | Number of FFT blocks (≡ MPI ranks) along x1, x2, x3 |
| `nmb` | `int` | MeshBlocks per rank = `nbx1 * nbx2 * nbx3` |
| `pmy_fb` | `FFTBlock*` | The one FFT block owned by this rank |

### Protected members used by FFTBlock

| Member | Meaning |
|--------|---------|
| `gcnt_` | Total cell count of the FFT mesh |
| `fft_mesh_size_` | `RegionSize` of the full FFT mesh |
| `fft_block_size_` | `RegionSize` of one FFT block (= union of MeshBlocks on this rank) |
| `fft_loclist_[nranks_]` | `LogicalLocation` of each rank's FFT block |
| `decomp_`, `pdim_` | Bitmask and count of decomposed axes |

### Constructor logic (`fft_driver.cpp:32-147`)

1. Copies mesh rank/block lists from `Mesh`.
2. Scans all MeshBlocks on this rank (`ns..ne`) to find the bounding box
   `(lx1min..lx1max, lx2min..lx2max, lx3min..lx3max)` in logical-location coordinates.
3. Computes `nbx1, nbx2, nbx3` = number of MeshBlocks along each axis within the bounding box,
   and `nmb = nbx1 * nbx2 * nbx3`.
4. **Cuboid constraint check** (line 93): `pm->nbtotal / nmb` must equal `nranks_`.
   This enforces a uniform cuboid assignment: every rank owns the same rectangular block of
   MeshBlocks.
5. Computes `fft_loclist_[n]` for each rank by dividing logical location by `(nbx1,nbx2,nbx3)`.
   This is the logical location of the FFT block (treating the cuboid as a single unit).
6. Sets `npx1 = nrbx1 / nbx1` (etc.) — number of FFT blocks along each axis.
7. Sets `fft_block_size_` = `fft_mesh_size_ / (npx1, npx2, npx3)`.
8. Computes `decomp_` bitmask and `pdim_` (number of distributed axes).

### `InitializeFFTBlock(bool set_norm)`

Creates the `FFTBlock` for this rank and optionally sets the normalization factor to `1/gcnt_`.

---

## 3. FFTBlock

**File**: `src/fft/athena_fft.hpp` (declaration), `src/fft/athena_fft.cpp` (implementation)

### Purpose
Represents the FFT domain of one MPI rank. The domain spans the union of all MeshBlocks on that
rank (the cuboid). Owns the complex data buffers and the four `AthenaFFTIndex` objects that
describe how data is laid out at each stage of the FFT pipeline.

### Public data members

| Member | Type | Meaning |
|--------|------|---------|
| `Nx[3]` | `int` | Global FFT mesh size along (x1, x2, x3) |
| `nx[3]` | `int` | Local size of this FFT block in real space |
| `disp[3]` | `int` | Global offset (displacement) of this block in real space |
| `kNx[3]` | `int` | Global size in Fourier space |
| `knx[3]` | `int` | Local size in Fourier space |
| `kdisp[3]` | `int` | Global offset in Fourier space |
| `dkx[3]` | `Real` | Wavenumber spacing: `2π / Lx[i]` |
| `dx1,dx2,dx3` | `Real` | Cell spacing |

These are filled in the constructor from `f_in_` (real space) and `b_in_` (Fourier space) using
`iloc[]` to map from internal permuted order back to physical (x1,x2,x3) order.

### Protected data members

| Member | Type | Meaning |
|--------|------|---------|
| `in_`, `out_` | `complex<Real>*` | Input / output FFT buffers, size `≥ cnt_` |
| `fplan_`, `bplan_` | `AthenaFFTPlan*` | Forward / backward plans |
| `f_in_`, `f_out_` | `AthenaFFTIndex*` | Index descriptors for forward FFT in/out |
| `b_in_`, `b_out_` | `AthenaFFTIndex*` | Index descriptors for backward FFT in/out |
| `norm_factor_` | `Real` | Applied during `RetrieveResult` |
| `loc_`, `msize_`, `bsize_` | — | FFT block logical location and sizes |
| `orig_idx_` | `AthenaFFTIndex` | Canonical unpermuted index (x1 fast, x3 slow) |
| `decomp_`, `pdim_` | `int` | Decomposition bitmask and axis count |
| `permute0_`, `permute1_`, `permute2_` | `int` | Permutation counts at each FFT stage |
| `swap1_`, `swap2_` | `bool` | Whether to swap mid↔slow before fwd / before bwd |

### Buffer size

`cnt_ = fft_block_size_.GetTotalCells()` = the total number of cells in the FFT block
(sum of all MeshBlocks on this rank).  `out_` may grow beyond `cnt_` if any pencil
decomposition in the fftMPI plan is larger than the brick (see §7).

---

## 4. AthenaFFTIndex

**File**: `src/fft/athena_fft.hpp` (declaration), `src/fft/athena_fft.cpp:542-653`

### Purpose
Describes the data layout at one stage of the distributed FFT: which axis is "fast" (contiguous
in memory), how axes are permuted, and what the local index range is on this MPI rank.

### Members

| Member | Type | Meaning |
|--------|------|---------|
| `Lx[3]` | `Real` | Physical domain length along each axis (in internal permuted order) |
| `Nx[3]` | `int` | Global cell count along each axis |
| `np[3]` | `int` | Number of MPI ranks along each axis |
| `ip[3]` | `int` | This rank's position along each axis |
| `nx[3]` | `int` | Local cell count = `Nx[i] / np[i]` |
| `is[3]`, `ie[3]` | `int` | Local start/end indices: `is[i] = ip[i]*nx[i]` |
| `iloc[3]` | `int` | Maps internal axis index → physical axis: `iloc[i] = j` means internal axis `i` corresponds to physical axis `j` |
| `ploc[3]` | `int` | Maps internal process axis → physical axis |

### Construction

`AthenaFFTIndex(int dim, LogicalLocation loc, RegionSize msize, RegionSize bsize)`

Initializes with the canonical (unpermuted) layout:
- `Nx[0] = msize.nx1`, `np[0] = msize.nx1/bsize.nx1`, `ip[0] = loc.lx1`  (x1 = fast)
- `Nx[1] = msize.nx2`, `np[1] = msize.nx2/bsize.nx2`, `ip[1] = loc.lx2`  (x2 = mid)
- `Nx[2] = msize.nx3`, `np[2] = msize.nx3/bsize.nx3`, `ip[2] = loc.lx3`  (x3 = slow)
- `iloc = {0,1,2}`, `ploc = {0,1,2}` (identity)

### Mutations

All mutations operate on the three internal arrays simultaneously:

| Method | Effect |
|--------|--------|
| `PermuteAxis(n)` | Cyclic left-rotate `iloc`, `Nx`, `Lx` by `n` positions |
| `PermuteProc(n)` | Cyclic left-rotate `ploc`, `np`, `ip` by `n` positions |
| `SwapAxis(ref)` | Swap the two axes OTHER than `ref` in `iloc`, `Nx`, `Lx` |
| `SwapProc(ref)` | Swap the two axes OTHER than `ref` in `ploc`, `np`, `ip` |
| `SetLocalIndex()` | Recompute `nx`, `is`, `ie` from `Nx`, `np`, `ip` |

### `GetIndex(i, j, k, pidx)` — index into the FFT buffer

```cpp
new_idx[0] = old_idx[pidx->iloc[0]];
new_idx[1] = old_idx[pidx->iloc[1]];
new_idx[2] = old_idx[pidx->iloc[2]];
return new_idx[0] + pidx->nx[0] * (new_idx[1] + pidx->nx[1] * new_idx[2]);
```

`iloc` permutes the physical (i,j,k) indices into the internal storage order before computing
the flat index. This is the key indirection that makes the same `in_` buffer serve multiple
stages without copying.

---

## 5. InitializeMPI — Building the Four Index States

**File**: `src/fft/athena_fft.cpp:470-536`

`InitializeMPI()` is called from the `FFTBlock` constructor (MPI builds only). It sets
`permute0_`, `permute1_`, `permute2_`, `swap1_`, `swap2_` based on the decomposition, then
constructs the four `AthenaFFTIndex` objects.

### Permutation strategy

For 1D or 2D process decomposition with a 3D FFT, fftMPI needs the
**undecomposed axis** to be the "fast" (innermost) axis. `permute0_` reorders axes so
the long (undecomposed) axis is fast:

| Decomposition | Long axis | `permute0_` | Effect on `(i,j,k)` |
|--------------|-----------|-------------|----------------------|
| `yz_decomp`  | x1        | 0           | `(i,j,k)` unchanged  |
| `xz_decomp`  | x2        | 1           | `(j,k,i)`            |
| `xy_decomp`  | x3        | 2           | `(k,i,j)`            |

`swap1_ = swap2_ = true` and `permute1_ = permute2_ = 2` for all 1D/2D decompositions.

For full 3D decomposition (`pdim_==3`): all permutations are 0, no swaps — two extra remaps
are performed by fftMPI internally.

### The four index states

```
orig_idx_            canonical layout: (x1 fast, x2 mid, x3 slow)
  │
  ├─ PermuteAxis(permute0_)
  ├─ PermuteProc(permute0_)
  ├─ [SwapAxis(0) + SwapProc(0) if swap1_]
  └─ SetLocalIndex()
       ↓
     f_in_          layout expected by fftMPI for forward FFT input

  f_in_ → PermuteAxis(permute1_) → SetLocalIndex()
       ↓
     f_out_         layout of forward FFT output (= Fourier space, permuted)

  f_out_ → [SwapAxis(0) + SwapProc(0) if swap2_] → SetLocalIndex()
       ↓
     b_in_          layout expected by fftMPI for backward FFT input

  b_in_ → PermuteAxis(permute2_) → SetLocalIndex()
       ↓
     b_out_         layout of backward FFT output (= back in real-space, permuted)
```

After the full forward+backward cycle the data is in `b_out_` layout, which has the same
`iloc` order as `f_in_` (the permutation cancels). `RetrieveResult` uses `b_out_` to map
back to physical (i,j,k).

### Worked example: `yz_decomp` (x2 and x3 distributed, x1 undecomposed)

```
orig_idx_:  iloc=(0,1,2)  Nx=(N1,N2,N3)  np=(1,P2,P3)
permute0_=0 → no change
swap1_=true → SwapAxis(0): swap axes 1↔2 → iloc=(0,2,1)  Nx=(N1,N3,N2)  np=(1,P3,P2)
f_in_:  iloc=(0,2,1)  fast=x1(local), mid=x3/P3, slow=x2/P2

permute1_=2 → PermuteAxis(2): shift left twice: (0,2,1)→(2,1,0)→(1,0,2)
f_out_: iloc=(1,0,2)  — Fourier output is in x2-pencil order

swap2_=true → SwapAxis(0) on f_out_: swap axes 1↔2 → iloc=(1,2,0)
b_in_:  iloc=(1,2,0)

permute2_=2 → PermuteAxis(2) on b_in_: shift left twice → iloc=(0,1,2)
b_out_: iloc=(0,1,2) = identity → same as original, ready to copy back
```

---

## 6. Cuboid Multi-MeshBlock Mapping

### LoadSource — filling the FFT buffer from MeshBlock data

```cpp
void FFTBlock::LoadSource(const AthenaArray<Real> &src, bool nu, int ngh,
                          LogicalLocation loc, RegionSize bsize)
```

Called once per MeshBlock on this rank. The MeshBlock's position within the FFT block is:

```cpp
is = loc.lx1 * bsize.nx1 - loc_.lx1 * bsize_.nx1;   // x1 offset within FFTBlock
js = loc.lx2 * bsize.nx2 - loc_.lx2 * bsize_.nx2;   // x2 offset
ks = loc.lx3 * bsize.nx3 - loc_.lx3 * bsize_.nx3;   // x3 offset
```

`loc_` and `bsize_` are the FFT block's own logical location and size (the cuboid origin).
`loc` and `bsize` are the current MeshBlock's location and size.

Then for each cell `(mi, mj, mk)` in the MeshBlock's active zone:
```cpp
std::int64_t idx = GetIndex(mi, mj, mk, f_in_);
dst[idx] = {src(n, k, j, i), 0.0};
```

`GetIndex(..., f_in_)` applies the `iloc` permutation to map the physical cell coordinates
to the storage order expected by fftMPI.

### RetrieveResult — copying FFT output back to MeshBlock

Same offset computation. For each cell:
```cpp
std::int64_t idx = GetIndex(mi, mj, mk, b_out_);
dst(k,j,i) = std::real(src[idx]) * norm_factor_;
```

Uses `b_out_` (backward FFT output layout) to map back to physical order.

---

## 7. CreatePlan — Building fftMPI Plans

**File**: `src/fft/athena_fft.cpp:336-415` (3D case shown)

```cpp
AthenaFFTPlan *FFTBlock::CreatePlan(int nfast, int nmid, int nslow,
                                    std::complex<Real> *data,
                                    AthenaFFTDirection dir)
```

### Forward plan

```cpp
plan->pf3d = new FFTMPI_NS::FFT3d(MPI_COMM_WORLD, 2);
plan->pf3d->scaled = 0;
// Compute output index ranges: apply inverse of permute1_ to f_out_->is/ie
for (int l=0; l<dim_; l++) {
  ois[l] = f_out_->is[(l+(dim_-permute1_)) % dim_];
  oie[l] = f_out_->ie[(l+(dim_-permute1_)) % dim_];
}
plan->pf3d->setup(nfast, nmid, nslow,
                  f_in_->is[0], f_in_->ie[0],   // input local range
                  f_in_->is[1], f_in_->ie[1],
                  f_in_->is[2], f_in_->ie[2],
                  ois[0], oie[0],                // output local range
                  ois[1], oie[1],
                  ois[2], oie[2],
                  permute1_, fftsize, sendsize, recvsize);
```

### Backward plan

Same structure using `b_in_->is/ie` as input and `b_out_` ranges as output, with `permute2_`.

### Buffer sizing

After `setup()`, fftMPI may require intermediate pencil sizes larger than the brick:
```cpp
// fast_cnt, mid_cnt, slow_cnt are the pencil sizes at each stage
int maxcnt = max({cnt_, fast_cnt, mid_cnt, slow_cnt});
if (maxcnt > buf_size_) {
  delete[] out_;
  out_ = new std::complex<Real>[maxcnt];
  buf_size_ = maxcnt;
}
```

This ensures `out_` is always large enough for any intermediate remap stage.

---

## 8. Execute — The fftMPI Remap+FFT Chain

**File**: `src/fft/athena_fft.cpp:418-460`

```cpp
void FFTBlock::Execute(AthenaFFTPlan *plan)
```

For a 3D plan the execution is an explicit five-step pipeline:

```
in_  →  [remap_prefast?]  →  fft_fast  →  remap_fastmid  →  fft_mid
     →  remap_midslow  →  fft_slow  →  [remap_postslow?]  →  out_
```

```cpp
FFTMPI_NS::FFT3d *pf3d = plan->pf3d;
if (pf3d->remap_prefast)
    pf3d->remap(in_, out_, pf3d->remap_prefast);   // optional pre-remap
else
    memcpy(out_, in_, 2*cnt_*sizeof(Real));
pf3d->perform_ffts(out_, plan->dir, pf3d->fft_fast);
pf3d->remap(out_, out_, pf3d->remap_fastmid);
pf3d->perform_ffts(out_, plan->dir, pf3d->fft_mid);
pf3d->remap(out_, out_, pf3d->remap_midslow);
pf3d->perform_ffts(out_, plan->dir, pf3d->fft_slow);
if (pf3d->remap_postslow)
    pf3d->remap(out_, out_, pf3d->remap_postslow); // optional post-remap
```

`plan->dir` is `FFTW_FORWARD` or `FFTW_BACKWARD`; fftMPI uses the same flag for both the
FFT direction and the remap direction.  `remap_prefast` and `remap_postslow` are null when
the input/output layout already matches the fast-pencil layout (identity case).

For a 2D plan the chain reduces to:
```
[remap_prefast?] → fft_fast → [remap_fastslow?] → fft_slow → [remap_postslow?]
```

**Serial fallback** (no MPI): `fftw_execute_dft(plan->plan, in_, out_)`. No fftMPI objects
are created.

---

## 9. AthenaFFTPlan Struct

```cpp
// MPI build (src/fft/athena_fft.hpp)
struct AthenaFFTPlan {
  FFTMPI_NS::FFT2d *pf2d;   // 2D fftMPI plan (MPI)
  FFTMPI_NS::FFT3d *pf3d;   // 3D fftMPI plan (MPI)
  fftw_plan plan;            // FFTW plan (serial)
  int dir;                   // FFTW_FORWARD or FFTW_BACKWARD
  int dim;                   // 1, 2, or 3
};
```

Only one of `pf3d`, `pf2d`, `plan` is non-null at a time. `Execute()` dispatches on `dim`.

---

## 10. Normalization

Set by `SetNormFactor(1.0/gcnt_)` where `gcnt_ = Nx1*Nx2*Nx3` (total global cells).
Applied multiplicatively in `RetrieveResult`. The forward FFT is unnormalized; the user is
responsible for calling `ExecuteForward` then `ExecuteBackward` and the result is scaled by
`1/N` automatically.

---

## 11. Key Design Constraints

1. **Uniform mesh only**: `FFTDriver` checks `use_uniform_meshgen_fn_` — non-uniform spacing
   is rejected.
2. **Cuboid constraint**: All ranks must own exactly `nmb = nbx1*nbx2*nbx3` MeshBlocks
   forming a rectangular cuboid. `nbtotal/nmb` must equal `nranks_`.
3. **One FFTBlock per rank**: `FFTDriver` creates exactly one `FFTBlock` (stored as `pmy_fb`).
4. **`AthenaFFTIndex` is the portability layer**: All axis-permutation and local-range logic
   lives here. The same `is/ie` ranges are passed to `FFTMPI_NS::FFT3d::setup()` for both
   input and output, making the plan self-describing.
5. **`remap_prefast`/`remap_postslow` may be null**: when the brick layout already matches
   the fast-pencil layout no pre/post remap is needed; `Execute` does a `memcpy` instead.

---

## 12. Relationship to BlockFFT

Both `FFTBlock` and `BlockFFT` now use the fftMPI C++ library (`FFTMPI_NS::FFT3d`).
The main remaining differences are in decomposition scope and API exposure:

| Property | `FFTBlock` (AthenaFFT) | `BlockFFT` |
|----------|------------------------|------------|
| MeshBlocks per rank | Many (cuboid, `nmb ≥ 1`) | Exactly 1 |
| Serial support | Yes (FFTW direct) | No |
| Index mapping | `AthenaFFTIndex` (explicit permutation objects) | Implicit in `FFT3d::setup()` ranges |
| Plan type | `AthenaFFTPlan` wrapping `FFT3d*`/`FFT2d*` | `FFTMPI_NS::FFT3d` directly |
| Remap stages exposed | Via `plan->pf3d->remap_*` (friend access) | Via `pf3d->remap_*` directly |
| Gravity solver | `FFTGravity` (all 4 BCs; see §13) | `BlockFFTGravity` (all 4 BCs; 1 MB/rank) |

`BlockFFT` / `BlockFFTGravity` are retained for the shearing-periodic BC regression test
(`swing.py`) and will be removed in a follow-up PR (MPIFFT.md Step 8).

---

## 13. FFTGravity — Multi-Meshblock Poisson Solver

**Files**: `src/gravity/fft_gravity.hpp`, `src/gravity/fft_gravity.cpp`

### Class hierarchy

```
FFTGravityDriver : public FFTDriver   (one per simulation)
  └── FFTGravity  : public FFTBlock   (one per MPI rank; extends FFTBlock)
```

`FFTGravity` inherits all of `FFTBlock`'s cuboid multi-meshblock machinery and overrides
`ApplyKernel`, `ExecuteForward`, and `ExecuteBackward` to dispatch on the gravity BC.

`FFTGravityDriver::Solve(stage, mode, gas_only)` handles source loading, FFT execution,
and result retrieval for all BC modes.

### Boundary conditions

| BC | `gbflag` | Algorithm | Extra data |
|----|----------|-----------|------------|
| Periodic | `GravityBoundaryFlag::periodic` | `ApplyKernel(mode)`: discrete (mode=0) or continuous (mode=1) Poisson kernel applied to `out_`; result in `in_` | — |
| Open | `GravityBoundaryFlag::open` | 8-parity convolution with cell-averaged Green's function `grf_` on 2× domain; 8 forward/backward FFT passes accumulate into `phi` | `grf_` (k-space GRF), `pf3dgrf_` (FFT3d for 2× domain) |
| Disk | `GravityBoundaryFlag::disk` | Even/odd z-pencil split: xy-FFT → phase shift → z-FFT on each half → kernel → iFFT → combine → ixy-FFT; decomposition-independent via dedicated `pf3d_disk_` | `in_e_`, `in_o_` (pencil buffers), `pf3d_disk_` |
| Shearing-periodic | (shear_periodic flag) | Global density assembly via `MPI_Allreduce`; y-roll applied in `LoadShearingSource`/`RetrieveShearingResult`; standard forward/backward FFT with shear kernel | global src/phi vectors |

### Key private members of FFTGravity

| Member | Meaning |
|--------|---------|
| `grf_` | k-space Green's function buffer for open BC (size ≥ max pencil on 2× domain) |
| `in_e_`, `in_o_` | Even/odd z-pencil buffers for disk BC |
| `pf3dgrf_` | `FFT3d*` for the 2× extended domain used in `InitGreen()` |
| `pf3d_disk_` | `FFT3d*` in physical (x,y,z) brick layout for disk BC forward/backward |
| `permute_disk_` | Fixed to 2 (z-pencil output) for `pf3d_disk_` |
| `four_pi_G_`, `Lx3_` | Disk/open BC physical constants |
| `qomt_` | Shear parameter `q*Omega*dt` set each timestep |
| `slow_nx1_`, `slow_nx2_`, `slow_nx3_` | Local dimensions of the z-pencil (disk BC) |
| `slow_ilo_`, `slow_jlo_`, `slow_klo_` | Global offsets of the z-pencil (disk BC, MPI) |

### Open BC: Green's function initialization

`InitGreen(four_pi_G)` (called once at construction):
1. Creates `pf3dgrf_`: an `FFT3d` for the `2Nx1 × 2Nx2 × 2Nx3` domain with the same
   permutation as `fplan_->pf3d` but at 2× scale.
2. Fills `grf_` with the cell-averaged Green's function (via `_GetIGF`) on the 2× domain,
   folded to `[-Nx, Nx-1]` range.
3. Forward-FFTs `grf_` in place so `MultiplyGreen(px, py, pz)` can multiply directly in
   k-space at parity offset `(2*ki+px, 2*kj+py, 2*kk+pz)`.

### Disk BC: decomposition independence

The standard `fplan_->pf3d` permutation depends on the MPI decomposition and may not keep
z as the slow axis.  `pf3d_disk_` is set up independently in the physical `(x,y,z)` layout
(`permute=2`, z-pencil output), so the even/odd z-split always operates on global z-indices
regardless of the process grid.

### Multi-meshblock support

`FFTGravity` inherits `FFTBlock`'s cuboid layout, so `nmb > 1` meshblocks per rank work
for all BC modes.  The shearing-periodic path assembles the full global density via
`MPI_Allreduce` before calling `LoadShearingSource`, which is the standard approach for
that BC and does not benefit from the cuboid structure.
