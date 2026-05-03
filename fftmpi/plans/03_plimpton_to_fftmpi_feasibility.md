# Plan 3: Feasibility — Replace Plimpton with fftMPI in `athena_fft`

## Verdict: **Feasible**

Replacing the Plimpton C calls with fftMPI C++ calls is architecturally sound. The two libraries
implement the same remap+FFT pipeline and share nearly identical `setup()` / `create_plan()`
parameter conventions. The changes are contained within `athena_fft.hpp` and `athena_fft.cpp`
(and `fft_driver.cpp` for the serial path), with no changes required to callers
(`TurbulenceDriver`, `PerturbationGenerator`, `FFTAnalysis`, gravity solvers).

---

## API Comparison

### Plan creation

| Plimpton | fftMPI |
|----------|--------|
| `fft_3d_create_plan(comm, nfast,nmid,nslow, in_ilo,in_ihi,...,in_khi, out_ilo,...,out_khi, 0, permute, &nbuf)` → `fft_plan_3d*` | `new FFT3d(comm, 2)` then `pf3d->setup(nfast,nmid,nslow, in_ilo,...,in_khi, out_ilo,...,out_khi, permute, fftsize, sendsize, recvsize)` |

Parameter mapping is one-to-one. `scaled=0` (no scaling) is the fftMPI default. `nbuf` → split
into `fftsize`, `sendsize`, `recvsize` (only sizes matter for buffer allocation).

### Execution

| Plimpton | fftMPI |
|----------|--------|
| `fft_3d(in_, out_, FFTW_FORWARD, plan3d)` | Explicit: `remap(prefast)` → `perform_ffts(FFTW_FORWARD, fft_fast)` → `remap(fastmid)` → ... |
| `fft_3d(in_, out_, FFTW_BACKWARD, plan3d)` | Explicit: `perform_ffts(FFTW_BACKWARD, fft_slow)` → `remap(slowmid)` → ... |

fftMPI exposes the remap+FFT steps explicitly (as `BlockFFT::ExecuteForward/Backward` already
shows). This is more verbose but identical logic.

### Destruction

| Plimpton | fftMPI |
|----------|--------|
| `fft_3d_destroy_plan(plan3d)` | `delete pf3d` |

---

## Required Architectural Changes

### 1. One `FFT3d*` replaces two `AthenaFFTPlan*`

Currently `FFTBlock` has separate `fplan_` and `bplan_` (forward and backward plans), each an
`AthenaFFTPlan` wrapping a `fft_plan_3d*`. With fftMPI, a **single** `FFT3d` instance handles
both directions — the forward remap chain and backward remap chain are both set up by one
`setup()` call. Replace:

```cpp
// OLD
AthenaFFTPlan *fplan_, *bplan_;

// NEW
#ifdef MPI_PARALLEL
FFTMPI_NS::FFT3d *pf3d_;
#endif
```

`AthenaFFTPlan` is still needed for the serial (non-MPI) FFTW path — keep the `fftw_plan`
members; only the `plan3d` / `plan2d` members are replaced.

### 2. In-place single buffer; `out_` eliminated for MPI path

Plimpton writes to a separate `out_` buffer: `fft_3d(in_, out_, dir, plan)`.
fftMPI works in-place (same buffer for input and output), exactly like `BlockFFT`.

- Drop `out_` for the MPI case.
- Resize `in_` to `max(cnt_, fast_cnt, mid_cnt, slow_cnt)` (pencil decompositions can be
  larger than the block — `BlockFFT` already handles this correctly, lines 118-125 of
  `block_fft.cpp`).
- For the serial FFTW path, keep the existing `out_` or use the in-place FFTW form
  `fftw_plan_dft_3d(..., in, in, ...)`.

### 3. `InitializeMPI()` simplifies — four index states reduce to one

Currently `InitializeMPI()` constructs four `AthenaFFTIndex` objects (`f_in_`, `f_out_`,
`b_in_`, `b_out_`) via a sequence of permutations and swaps. This was needed because Plimpton's
`create_plan` requires pre-permuted axis ranges.

With fftMPI, `setup()` takes the unpermuted global index ranges directly and handles internal
axis ordering via the `permute` flag. Only `orig_idx_` (the canonical unpermuted layout) is
needed to provide `is/ie` ranges.

New simplified MPI initialization:

```cpp
void FFTBlock::InitializeMPI() {
    // Use orig_idx_ directly — no pre-permutation needed
    int permute = 2;  // output layout: (slow, mid, fast) after forward FFT
    int fftsize, sendsize, recvsize;

    // Probe setup to get pencil decomposition sizes
    pf3d_ = new FFTMPI_NS::FFT3d(MPI_COMM_FFT, 2);
    pf3d_->setup(orig_idx_.Nx[0], orig_idx_.Nx[1], orig_idx_.Nx[2],
                 orig_idx_.is[0], orig_idx_.ie[0],
                 orig_idx_.is[1], orig_idx_.ie[1],
                 orig_idx_.is[2], orig_idx_.ie[2],
                 orig_idx_.is[0], orig_idx_.ie[0],  // out = same as in (probe)
                 orig_idx_.is[1], orig_idx_.ie[1],
                 orig_idx_.is[2], orig_idx_.ie[2],
                 permute, fftsize, sendsize, recvsize);

    // Read back pencil ranges, then recreate with correct output range
    int slow_ilo = pf3d_->slow_ilo, ...; // etc.
    delete pf3d_;
    pf3d_ = new FFTMPI_NS::FFT3d(MPI_COMM_FFT, 2);
    pf3d_->setup(orig_idx_.Nx[0], ...,
                 orig_idx_.is[0], ...,      // input
                 slow_ilo, slow_ihi, ...,    // output = slow-pencil decomposition
                 permute, fftsize, sendsize, recvsize);

    // f_in_ is still needed (for Nx, nx, disp, kNx, knx, kdisp, dkx public members)
    // but f_out_, b_in_, b_out_ are no longer needed
    f_in_ = new AthenaFFTIndex(&orig_idx_);  // keep for metadata
    // b_in_ now encodes the Fourier-space (slow-pencil) layout:
    //   populate from pf3d_->slow_* members
}
```

`f_out_`, `b_in_`, `b_out_` can be dropped or collapsed into derived metadata members.

### 4. `LoadSource` / `RetrieveResult` simplify — no `iloc` permutation

Currently `LoadSource` calls `GetIndex(mi, mj, mk, f_in_)` which applies the `iloc`
permutation. With fftMPI, the input layout is always the natural brick layout (no
pre-permutation), so the flat index is simply:

```cpp
// OLD: idx = GetIndex(mi, mj, mk, f_in_)   (applies iloc permutation)
// NEW (same as BlockFFT, but with multi-block offset):
idx = mi + bsize_.nx1 * (mj + bsize_.nx2 * mk);
```

The multi-block offset calculation (lines 173-176 in `athena_fft.cpp`) is unchanged.

`RetrieveResult` similarly simplifies — use `orig_idx_` or `b_out_` dimensions without `iloc`.

### 5. `Execute` → explicit remap+perform_ffts

Replace:
```cpp
// OLD
void FFTBlock::Execute(AthenaFFTPlan *plan) {
    fft_3d(reinterpret_cast<fftw_complex*>(in_),
           reinterpret_cast<fftw_complex*>(out_),
           plan->dir, plan->plan3d);
}
```

With the same pattern as `BlockFFT::ExecuteForward/Backward`:
```cpp
// NEW
void FFTBlock::ExecuteForward() {
    FFT_SCALAR *data = reinterpret_cast<FFT_SCALAR*>(in_);
    if (pf3d_->remap_prefast) pf3d_->remap(data, data, pf3d_->remap_prefast);
    pf3d_->perform_ffts(reinterpret_cast<FFT_DATA*>(data), FFTW_FORWARD, pf3d_->fft_fast);
    pf3d_->remap(data, data, pf3d_->remap_fastmid);
    pf3d_->perform_ffts(reinterpret_cast<FFT_DATA*>(data), FFTW_FORWARD, pf3d_->fft_mid);
    pf3d_->remap(data, data, pf3d_->remap_midslow);
    pf3d_->perform_ffts(reinterpret_cast<FFT_DATA*>(data), FFTW_FORWARD, pf3d_->fft_slow);
}
```

Note: `remap_prefast`, `perform_ffts`, `remap`, `fft_fast/mid/slow` are **private** members
of `FFT3d` (unlike `BlockFFT` which gains access via `friend`). We need to either:
- (a) Add `FFTBlock` as a `friend` of `FFT3d` in `fftmpi/fft3d.h` (same pattern as `BlockFFT`)
- (b) Or use `pf3d_->compute(data, data, 1)` / `compute(data, data, -1)` if round-trip layout
  is symmetric (simpler but less flexible)

Option (a) is preferred for consistency with `BlockFFT`.

### 6. `ApplyKernel` simplifies

Currently `ApplyKernel` copies `out_[idx_out]` → `in_[idx_in]`. With in-place transforms,
data is already in `in_` after `ExecuteForward()`. `ApplyKernel` only needs to multiply
in-place (or do nothing for the base class).

### 7. Public metadata (`Nx`, `nx`, `disp`, `kNx`, `knx`, `kdisp`, `dkx`)

These are currently filled from `f_in_` and `b_in_` using `iloc[]`. With fftMPI:
- Real-space metadata (`Nx`, `nx`, `disp`) comes from `orig_idx_`
- Fourier-space metadata (`kNx`, `knx`, `kdisp`) comes from `pf3d_->slow_*` members

`dkx[i] = TWO_PI / Lx[i]` is unchanged.

---

## What Is NOT Changed

| Component | Change needed? |
|-----------|---------------|
| `AthenaFFTIndex` class | **Kept** — still needed for multi-block offset in `LoadSource`/`RetrieveResult`, and for metadata (`Nx`, `nx`, `disp`, etc.) |
| `FFTDriver` constructor | No change — cuboid detection and `nmb` logic is independent of FFT library |
| `LoadSource` / `RetrieveResult` offset math | No change — `is = loc.lx1*bsize.nx1 - loc_.lx1*bsize_.nx1` stays |
| Callers (TurbulenceDriver, gravity, etc.) | No change to public API of `FFTBlock` |
| Serial FFTW path | No change — `fftw_plan_dft_3d` + `fftw_execute_dft` unchanged |
| 2D FFT path | Replace `fft_2d_create_plan` → `FFTMPI_NS::FFT2d::setup()` (same pattern) |

---

## Risk / Complications

| Risk | Severity | Mitigation |
|------|----------|------------|
| `remap_prefast` etc. are `private` in `FFT3d` | Medium | Add `friend class FFTBlock` to `fft3d.h` (one-line change, same as existing `BlockFFT` friend) |
| Pencil decomposition buffer sizing | Low | Follow exact `BlockFFT` pattern: `maxcnt = max(cnt_, fast_cnt, mid_cnt, slow_cnt)` |
| Permute flag semantics | Low | Use `permute=2` (same as `BlockFFT`); output layout = (slow,mid,fast) for both classes |
| Forward/backward `b_in_` encoding Fourier layout | Low | Populate `kNx/knx/kdisp` from `pf3d_->slow_*` after setup |
| 2D FFT callers (TurbulenceDriver?) | Low | Assess separately; `FFT2d` API mirrors `FFT3d` |

---

## Conclusion

The swap is feasible with moderate refactoring of `FFTBlock`/`AthenaFFTPlan` internals.
No changes to the public API or callers. The resulting code is cleaner:
- 4 `AthenaFFTIndex` states → 1 (or 2)
- Separate `in_`/`out_` buffers → single in-place buffer
- Opaque `fft_3d()` → explicit remap+FFT sequence (same as `BlockFFT`, easier to debug)

---

## Next Step

→ **Plan 4**: Implement `athena_fft_gravity` — a counterpart to `block_fft_gravity` that uses
`athena_fft` (with fftMPI backend) for gravity, supporting multiple meshblocks per rank.

Before starting implementation, we need to decide on branch strategy:
- Branch for fftMPI backend swap in `athena_fft` (this plan's changes)
- Branch for `athena_fft_gravity` (Plan 4)
- Or combine into one branch if the gravity work depends on the backend swap
