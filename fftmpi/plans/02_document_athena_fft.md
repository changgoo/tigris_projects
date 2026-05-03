# Plan 2: AthenaFFT Documentation

## Objective
Write full documentation of the `athena_fft` wrapper: class hierarchy, index mapping, cuboid
decomposition, and Plimpton API usage. Output: `AthenaFFT.md` in repo root.

## Status: COMPLETE

Output file: `AthenaFFT.md`

## Summary of findings

### Class hierarchy
```
FFTDriver → FFTBlock → {AthenaFFTIndex ×4, AthenaFFTPlan ×2}
```

### Cuboid decomposition
- `FFTDriver` constructor scans all MeshBlocks on this rank to find the bounding box.
- `nmb = nbx1 * nbx2 * nbx3` MeshBlocks per rank; must be uniform across all ranks.
- The FFT block spans the union; `fft_block_size_` = `fft_mesh_size_ / (npx1, npx2, npx3)`.
- `LoadSource` / `RetrieveResult` compute per-MeshBlock offsets:
  `is = loc.lx1*bsize.nx1 - loc_.lx1*bsize_.nx1`

### Index mapping (AthenaFFTIndex)
- Tracks `Nx`, `np`, `ip`, `nx`, `is`, `ie`, and `iloc` (axis permutation map).
- Four states: `f_in_` → `f_out_` → `b_in_` → `b_out_`.
- Built in `InitializeMPI()` via `PermuteAxis`, `PermuteProc`, `SwapAxis`, `SwapProc`.
- `GetIndex(i,j,k, pidx)` uses `iloc` to permute coordinates before computing flat index.

### Plimpton API
- `fft_3d_create_plan(comm, nfast,nmid,nslow, in_ilo,...,in_khi, out_ilo,...,out_khi, 0, permute, &nbuf)`
- `fft_3d(in, out, dir, plan3d)`
- `fft_3d_destroy_plan(plan3d)`
- Input/output ranges come directly from `f_in_->is/ie` and (inverse-permuted) `f_out_->is/ie`.

### Key insight for migration
`AthenaFFTIndex.is[]` and `ie[]` are exactly the local index ranges that `fftMPI::FFT3d::setup()`
also expects. The permutation flags (`permute1_`, `permute2_`) map directly to `fftMPI`'s
`permute` argument. The main work is:
1. Replace `fft_3d_create_plan` + `fft_3d` with `FFT3d::setup` + `FFT3d::perform_ffts`.
2. Preserve `AthenaFFTIndex` for the cuboid multi-block aggregation — fftMPI has no equivalent.

## Next Step
→ **Plan 3**: Feasibility assessment — can Plimpton C calls in `athena_fft` be replaced by
fftMPI C++ equivalents without architectural changes?
