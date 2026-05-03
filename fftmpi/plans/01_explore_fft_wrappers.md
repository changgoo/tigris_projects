# Plan 1: Code Exploration — FFT MPI Wrapper Structure

## Objective
Understand the existing FFT MPI wrappers (`athena_fft` / `block_fft`) well enough to assess
feasibility of migrating `athena_fft` to use the `fftmpi` C++ backend.

---

## Class Hierarchy

### AthenaFFT Stack (`athena_fft.hpp/.cpp`, `fft_driver.cpp`)
Uses the older `plimpton/` C library.

```
FFTDriver
  └── FFTBlock          (one per MPI rank, spans all MeshBlocks on that rank)
        ├── AthenaFFTIndex  f_in_, f_out_, b_in_, b_out_   (index mapping)
        └── AthenaFFTPlan   fplan_, bplan_                 (plan wrappers)
              └── struct fft_plan_3d* plan3d  (or plan2d, fftw_plan)
```

`AthenaFFTIndex` tracks global mesh size, local extents, axis permutations, and MPI
decomposition geometry. It is the core of the cuboid multi-block aggregation.

### BlockFFT Stack (`block_fft.hpp/.cpp`)
Uses the newer `fftmpi/` C++ library.

```
BlockFFT
  └── FFTMPI_NS::FFT3d* pf3d    (one fftMPI object per MeshBlock)
BlockFFTGravity : public BlockFFT
  └── FFTMPI_NS::FFT3d* pf3dgrf_   (second FFT3d for Green's function, open BC only)
```

---

## Cuboid Multi-Block Mapping (FFTBlock)

`FFTDriver` enforces a cuboid constraint: all MeshBlocks on a rank must tile a rectangular
box. `nmb = nbx1 * nbx2 * nbx3` (e.g., 2×2×2 = 8 MeshBlocks/rank).

`FFTBlock` represents the *union* of all MeshBlocks on the rank as a single FFT domain.

**LoadSource / RetrieveResult** compute the relative offset of each MeshBlock within the
FFTBlock:
```cpp
is = loc.lx1 * bsize.nx1 - loc_.lx1 * bsize_.nx1
```
then iterate over each MeshBlock's cells, mapping through `AthenaFFTIndex` to the global
FFT index.

`BlockFFT` has no equivalent: it works with exactly one MeshBlock per instance.

---

## Plimpton (C) API Surface — used in `athena_fft`

| Call | Purpose |
|------|---------|
| `fft_3d_create_plan(comm, nfast,nmid,nslow, ilo,ihi,jlo,jhi,klo,khi, out_*,permute, &nbuf)` | Create 3D plan |
| `fft_2d_create_plan(comm, nfast,nslow, ilo,ihi,jlo,jhi, out_*,permute, &nbuf)` | Create 2D plan |
| `fft_3d(in, out, dir, plan3d)` | Execute 3D FFT |
| `fft_2d(in, out, dir, plan2d)` | Execute 2D FFT |
| `fft_3d_destroy_plan(plan3d)` / `fft_2d_destroy_plan(plan2d)` | Destroy plan |

---

## fftMPI (C++) API Surface — used in `block_fft`

| Call | Purpose |
|------|---------|
| `new FFTMPI_NS::FFT3d(comm, 2)` | Construct (2 = double precision) |
| `pf3d->setup(Nx,Ny,Nz, in_ilo,ihi,jlo,jhi,klo,khi, out_ilo,...,permute, &fftsize,&sendsize,&recvsize)` | Set up pencil decomposition |
| `pf3d->perform_ffts(data, dir, fft_struct)` | Execute 1D FFT along one axis |
| `pf3d->remap(in, out, remap_struct)` | Transpose between pencil decompositions |
| `delete pf3d` | Destroy |

Public members used by `BlockFFTGravity`:
- `fast_ilo/ihi/jlo/jhi/klo/khi`, `mid_*`, `slow_*` — pencil index ranges
- `remap_prefast`, `remap_fastmid`, `remap_midslow`, `remap_postslow` — remap plans
- `fft_fast`, `fft_mid`, `fft_slow` — 1D FFT structs

---

## Structural Compatibility

Both libraries implement the same algorithmic pipeline:

```
Input (brick) → remap → FFT(x) → remap → FFT(y) → remap → FFT(z) → output (z-pencil)
```

Both use FFTW3 for 1D FFTs and MPI all-to-all for remaps.

| Aspect | Plimpton (C) | fftMPI (C++) |
|--------|-------------|--------------|
| API style | C functions + opaque structs | C++ class |
| Multi-block/rank | Yes (via FFTBlock aggregation) | No (1 block per instance) |
| Remap management | Opaque inside plan | Exposed as public struct pointers |
| Axis permutation | Explicit `permute0/1/2` | Single `permute` flag in `setup()` |
| Index range input | `(ilo,ihi,jlo,jhi,klo,khi)` pair for in+out | Same convention |

**Bottom line**: The pipelines are compatible. Replacing the Plimpton calls with fftMPI calls
is feasible at the call-site level. The main complexity is:
1. Preserving the cuboid multi-block aggregation (unique to `FFTBlock` / `AthenaFFTIndex`).
2. Mapping the axis-permutation logic from `AthenaFFTIndex` to fftMPI's `permute` flag in
   `setup()`.

---

## BlockFFTGravity — How It Uses BlockFFT

`BlockFFTGravity : public BlockFFT` overrides three methods:

| Method | What it does |
|--------|-------------|
| `ExecuteForward()` | For periodic/open BC: delegates to `BlockFFT::ExecuteForward()`. For disk BC: manual remap+FFT steps with even/odd splitting and phase shifts. |
| `ApplyKernel()` | Multiplies Fourier coefficients by wavenumber-dependent kernel (periodic, open, or disk BC). |
| `ExecuteBackward()` | Inverse of above; applies 1/8 normalization for open BC (8× extended domain). |

Gravity solver uses these public `BlockFFT` members directly:
`in_`, `Nx1/Nx2/Nx3`, `nx1/nx2/nx3`, `in_ilo/ihi/jlo/jhi/klo/khi`,
`slow_ilo/ihi/jlo/jhi/klo/khi`, `remap_*`, `fft_fast/mid/slow`.

A replacement `athena_fft_gravity` would need equivalent public access to the same data.

---

## Next Step

→ **Plan 2**: Full documentation of `AthenaFFT` / `FFTBlock` class hierarchy, index mapping
logic, and cuboid decomposition. Output: `AthenaFFT.md`.
