# Plan 4: Extend `FFTGravity` to Full BC Parity with `BlockFFTGravity`

## Goal

Extend `FFTGravity : public FFTBlock` and `FFTGravityDriver` to support all gravity boundary
conditions currently handled by `BlockFFTGravity`: periodic, open (3D), and disk (horizontal-
periodic + vertical-open). Shearing-periodic BC is documented as Phase 2 and deferred (per
MPIFFT.md Step 6). The two gravity solvers coexist on the same branch (`fftmpi-athena-fft`)
until `BlockFFT`/`BlockFFTGravity` is retired (Plan 8).

---

## Current State

| Class | Backend | BC support | Meshblocks/rank |
|-------|---------|-----------|-----------------|
| `FFTGravityDriver` + `FFTGravity` | `FFTBlock` (fftMPI) | periodic only | multiple (cuboid) |
| `BlockFFTGravityDriver` + `BlockFFTGravity` | `BlockFFT` (fftMPI) | periodic, open, disk, shearing | 1 (fixed) |

`FFTGravity::ApplyKernel(int mode)` currently implements the periodic Poisson kernel only
(modes 0 and 1). `FFTGravityDriver::Solve()` has a single linear path:
`LoadSource → ExecuteForward → ApplyKernel → ExecuteBackward → RetrieveResult`.

---

## Architecture Overview

### Key differences: FFTBlock vs. BlockFFT

| | `FFTBlock` (used by `FFTGravity`) | `BlockFFT` (used by `BlockFFTGravity`) |
|-|-----------------------------------|----------------------------------------|
| MeshBlocks/rank | Multiple (cuboid) | 1 |
| Buffer layout | `in_` (source), `out_` (result), out-of-place | `in_` only, in-place |
| Axis index tracking | 4 `AthenaFFTIndex` objects: `f_in_`, `f_out_`, `b_in_`, `b_out_` | Exposed `fast/mid/slow_*` ranges copied from `pf3d` |
| `pf3d` ownership | Per-plan: `fplan_->pf3d`, `bplan_->pf3d` | Single `pf3d` per `BlockFFT` instance |
| k-space metadata | `knx[3]`, `kdisp[3]`, `kNx[3]`, `dkx[3]` | `slow_nx1/2/3`, `slow_ilo/ihi/...` |

After `pfb->ExecuteForward()` in `FFTGravityDriver::Solve()`:
- k-space data is in `out_` of the `FFTBlock`
- Local k-space extents: `pfb->knx[0,1,2]`
- Global k-space offsets: `pfb->kdisp[0,1,2]`
- Total k-space extents: `pfb->kNx[0,1,2]`
- `ApplyKernel` reads `out_[GetIndex(i,j,k,f_out_)]`, writes `in_[GetIndex(i,j,k,b_in_)]`

### Friend access to `FFT3d` privates

`FFT3d` (in `fftmpi/fft3d.h`) already grants friend access to:
```cpp
friend class ::BlockFFT;
friend class ::BlockFFTGravity;
friend class ::FFTBlock;
```

`FFTGravity` needs the same grant for disk and open BC implementations that access
`fplan_->pf3d->remap_prefast`, `fplan_->pf3d->fft_fast`, `fplan_->pf3d->slow_ilo`, etc.

**Required change**: add `friend class ::FFTGravity;` to `fft3d.h`.

---

## BC-Specific Designs

### 1. Periodic BC (already implemented)

`ApplyKernel(mode=0)` discrete FT, `mode=1` continuous FT. No changes needed.
`Solve()` path: `LoadSource → ExecuteForward → ApplyKernel → ExecuteBackward → RetrieveResult`.

### 2. Open BC (3D fully open)

**Algorithm** (same as `BlockFFTGravity`):
- Solve via convolution: ρ * G where G is the cell-averaged Green's function
- Avoids zero-padding by decomposing into 8 parity components (px,py,pz ∈ {0,1})
- Each parity: load ρ with phase factors → forward FFT → multiply `grf_` → backward FFT → accumulate into `phi`
- `grf_` is computed once at construction time from the 2× extended domain

**Green's function FFT (`pf3dgrf_`)**:
- A separate `FFTMPI_NS::FFT3d` for the 2× extended domain (`2*Nx1 × 2*Nx2 × 2*Nx3`)
- Set up with the same permute flag as `fplan_->pf3d`, input ranges `[2*in_ilo, 2*in_ihi+1] × ...`,
  output ranges `[2*slow_ilo, 2*slow_ihi+1] × ...`
- Output layout of `pf3dgrf_` then spans exactly twice the index range of `fplan_->pf3d` output,
  so `grf_[2*kx+px, 2*ky+py, 2*kz+pz]` maps directly to the parity-offset of k-space element `(kx,ky,kz)`
- Requires friend access to read `fplan_->pf3d->in_ilo/ihi`, `slow_ilo/ihi`, `permute` for setup

**`MultiplyGreen(px, py, pz)`**: for each local k-space element `(i,j,k)`:
```
out_[GetIndex(i,j,k,f_out_)] *= grf_[grf_local_index(i+kdisp[0], j+kdisp[1], k+kdisp[2], px,py,pz)]
```
where `grf_local_index` maps the 2× domain k-space position to a local index within `grf_`.

**`LoadOBCSource(src, px, py, pz)`**: loads density into `in_` with the parity-specific phase shift
`exp(i*π*(px*x/Nx1 + py*y/Nx2 + pz*z/Nx3))` applied cell-by-cell. Same formula as `BlockFFTGravity`.

**`RetrieveOBCResult(dst, px, py, pz)`**: reads real part of `in_` (after backward FFT) and adds
to `dst` with appropriate phase factor, accumulating over the 8 parity passes.

**`grf_` buffer sizing**: `8 × max(knx[0]*knx[1]*knx[2])` across all parity setups
(conservatively: `8 × max(kNx[0]*kNx[1]*kNx[2] / nranks_)`, or follow BlockFFTGravity's
pattern using the max over fast/mid/slow pencil sizes).

**`Solve()` open BC path**:
```cpp
phi.ZeroClear();
for pz in {0,1}: for py in {0,1}: for px in {0,1}:
    LoadOBCSource(rho, px, py, pz)
    pfb->ExecuteForward()          // normal FFT
    pfb_grav->MultiplyGreen(px,py,pz)
    pfb->ExecuteBackward()         // normal iFFT
    pfb->ApplyKernel(mode=-1)      // no-op for open — skip; multiplication already done
    RetrieveOBCResult(phi, px, py, pz)
```
Note: for open BC, `ExecuteForward`/`Backward` are the standard calls — no interception needed.
The normalization factor must be adjusted: divide by 8 after accumulation (or per-pass).

### 3. Disk BC (horizontal periodic, vertical open)

**Algorithm** (same as `BlockFFTGravity`):
- Forward transform in x and y only (not z)
- Split into even (`in_e_`, l=2p) and odd (`in_o_`, l=2p+1) z-modes
- Apply phase shift to odd: `in_o_[k] *= exp(-iπk/Nz)`
- FFT in z separately on `in_e_` and `in_o_`
- Apply kernel to each: `in_e_[idx] *= kernel_e`, `in_o_[idx] *= kernel_o`
- Inverse FFT in z on `in_e_` and `in_o_`
- Combine: `in_[k] = in_e_[k] + exp(iπk/Nz)*in_o_[k]`
- Inverse transform in y and x

**Key challenge**: The standard `FFTBlock::Execute()` performs the full 3-axis transform at once.
Disk BC requires intercepting between the y-FFT and z-FFT stages. This requires:
1. Making `FFTBlock::ExecuteForward()`/`ExecuteBackward()` virtual (mirrors `BlockFFT`)
2. Overriding them in `FFTGravity` to dispatch on `gbflag`

For the z-split step, `FFTGravity::ExecuteForward()` (disk path) must access `fplan_->pf3d`'s
internal remap stages:
```cpp
// disk BC forward in FFTGravity:
FFT_SCALAR *data = reinterpret_cast<FFT_SCALAR*>(in_);
FFT_SCALAR *out  = reinterpret_cast<FFT_SCALAR*>(out_);
FFTMPI_NS::FFT3d *pf3d = fplan_->pf3d;
if (pf3d->remap_prefast) pf3d->remap(data, out, pf3d->remap_prefast);
else std::memcpy(out, data, 2*cnt_*sizeof(Real));
pf3d->perform_ffts(FFT_DATA(out), FFTW_FORWARD, pf3d->fft_fast);   // x-FFT
pf3d->remap(out, out, pf3d->remap_fastmid);
pf3d->perform_ffts(FFT_DATA(out), FFTW_FORWARD, pf3d->fft_mid);    // y-FFT
pf3d->remap(out, out, pf3d->remap_midslow);
// STOP — copy to in_e_, in_o_, apply phase, then z-FFT each separately
std::memcpy(in_e_, out_, sizeof(std::complex<Real>)*knx[0]*knx[1]*knx[2]);
std::memcpy(in_o_, out_, sizeof(std::complex<Real>)*knx[0]*knx[1]*knx[2]);
// phase shift on in_o_, then:
pf3d->perform_ffts(FFT_DATA(in_e_), FFTW_FORWARD, pf3d->fft_slow);
pf3d->perform_ffts(FFT_DATA(in_o_), FFTW_FORWARD, pf3d->fft_slow);
```

**Axis ordering constraint**: The above requires `fplan_->pf3d->fft_slow` to be the z-FFT.
In `FFTBlock`, the axis mapping depends on the AthenaFFTIndex permutation, which in turn depends
on domain decomposition. We must verify (and enforce if needed) that z is the "slow" axis in
the `fplan_->pf3d` setup when gbflag==disk. In practice, disk simulations decompose in x and y
only (npz=1), so z is never remapped and naturally remains as the slow dimension.

**`in_e_`, `in_o_` buffer sizing**: `knx[0]*knx[1]*knx[2]` complex values (k-space local size).
These match `slow_nx1*slow_nx2*slow_nx3` in `BlockFFTGravity`.

**`ApplyKernel` disk path**: inline in `FFTGravity::ApplyKernel(int mode)` dispatching on
`gbflag==disk`. Operates on `in_e_` and `in_o_` instead of `out_` and `in_`. Kernel formula
identical to `BlockFFTGravity::ApplyKernel()` (disk branch). Uses `kdisp[l]` for global k-index
offsets and `kNx[l]` for the full domain sizes.

**`ExecuteBackward()` disk path**: inverse of forward — iFFT_z on `in_e_`/`in_o_`, combine,
then remap backward (slowmid→iFFTy→midfast→iFFTx→postfast). Reads from `fplan_->pf3d` backward
remap chain (or `bplan_->pf3d`). Must verify which plan contains the backward remap chain.

---

## Required Changes

### `fftmpi/fft3d.h`

```cpp
// Add forward declaration and friend grant:
class FFTGravity;

namespace FFTMPI_NS {
class FFT3d {
  friend class ::BlockFFT;
  friend class ::BlockFFTGravity;
  friend class ::FFTBlock;
  friend class ::FFTGravity;       // NEW
  ...
```

### `src/gravity/block_fft_gravity.hpp`

Move `GravityBoundaryFlag` enum and `GetGravityBoundaryFlag()` to a shared header
(e.g., `src/gravity/gravity_bc.hpp`) so both `FFTGravity` and `BlockFFTGravity` can use it.
Alternatively, duplicate the enum in `fft_gravity.hpp` temporarily.

### `src/fft/athena_fft.hpp` — `FFTBlock`

Make `ExecuteForward()` and `ExecuteBackward()` virtual:
```cpp
virtual void ExecuteForward() { Execute(fplan_); }
virtual void ExecuteBackward() { Execute(bplan_); }
```
This mirrors the `BlockFFT` pattern and lets `FFTGravity` override them for disk BC.

### `src/gravity/fft_gravity.hpp` — `FFTGravity`

```cpp
class FFTGravity : public FFTBlock {
 public:
  FFTGravity(FFTDriver *pfd, LogicalLocation iloc, int igid,
             RegionSize msize, RegionSize bsize);
  ~FFTGravity();

  GravityBoundaryFlag gbflag;
  Real time_int;  // time at potential calculation (shearing BC, Phase 2)

  void ApplyKernel(int mode) final;
  void ExecuteForward() override;        // dispatches on gbflag (disk needs override)
  void ExecuteBackward() override;       // dispatches on gbflag (disk needs override)

  void InitGreen();
  void LoadOBCSource(const AthenaArray<Real> &src, int px, int py, int pz);
  void RetrieveOBCResult(AthenaArray<Real> &dst, int px, int py, int pz);
  void MultiplyGreen(int px, int py, int pz);

 private:
  Real dx1_, dx2_, dx3_;
  Real Lx1_, Lx2_, Lx3_;
  const std::complex<Real> I_;
  std::complex<Real> *grf_;       // Green's function k-space buffer (open BC)
  std::complex<Real> *in_e_;      // even z-mode buffer (disk BC)
  std::complex<Real> *in_o_;      // odd z-mode buffer (disk BC)
#ifdef MPI_PARALLEL
#ifdef FFT
  FFTMPI_NS::FFT3d *pf3dgrf_;    // FFT3d for 2× domain (open BC Green's function)
#endif
#endif
};
```

### `src/gravity/fft_gravity.hpp` — `FFTGravityDriver`

```cpp
class FFTGravityDriver : public FFTDriver {
 public:
  FFTGravityDriver(Mesh *pm, ParameterInput *pin);
  ~FFTGravityDriver();
  void Solve(int stage, int mode, bool gas_only=false);

 private:
  Real four_pi_G_;
  GravityBoundaryTaskList *gtlist_;
  // Shearing BC (Phase 2 — declared here, implemented later):
  // Real Omega_0_, qshear_;
  // AthenaArray<Real> roll_var, roll_buf, send_buf, recv_buf, pflux;
  // AthenaArray<Real> send_gbuf, recv_gbuf;
  // void RollUnroll(AthenaArray<Real> &dat, Real dt);
  // Real ShearTimeShift();
};
```

### `src/gravity/fft_gravity.cpp`

New/modified functions:
- `FFTGravity::FFTGravity()` constructor — reads `grav_bc` from input, allocates `grf_`/`in_e_`/`in_o_`, calls `InitGreen()` for open BC
- `FFTGravity::~FFTGravity()` — delete `grf_`, `in_e_`, `in_o_`, `pf3dgrf_`
- `FFTGravity::InitGreen()` — compute cell-averaged GRF on 2× domain, FFT with `pf3dgrf_`
- `FFTGravity::LoadOBCSource()` — phase-shifted density load into `in_`
- `FFTGravity::RetrieveOBCResult()` — phase-shifted accumulation from `in_` into `dst`
- `FFTGravity::MultiplyGreen()` — k-space multiply `out_` by `grf_` (for open BC)
- `FFTGravity::ApplyKernel()` — dispatches on gbflag (periodic: existing; disk: new; open: no-op since MultiplyGreen handles it)
- `FFTGravity::ExecuteForward()` — disk path intercepts z-FFT; others call base
- `FFTGravity::ExecuteBackward()` — disk path intercepts z-FFT; others call base
- `FFTGravityDriver::FFTGravityDriver()` — reads `grav_bc`, passes to `FFTGravity` ctor
- `FFTGravityDriver::Solve()` — dispatches on gbflag:
  - periodic: existing path
  - open: parity loop (ZeroClear → 8×(Load+Fwd+MultiplyGreen+Bwd+Retrieve))
  - disk: periodic-like path but ExecuteForward/Backward handle z-split internally

---

## `_GetIGF` and Green's function initialization

The cell-averaged Green's function formula (`_GetIGF`) in `block_fft_gravity.hpp` is a free
function. Move to a shared translation unit or duplicate in `fft_gravity.cpp`. The computation
is identical — no algorithmic changes needed.

---

## k-space index mapping: `pf3dgrf_` setup

The 2× Green's function FFT3d is set up in `FFTGravity` constructor (after `QuickCreatePlan()`
has been called, so `fplan_->pf3d` is available):

```cpp
// Mirror fplan_->pf3d at 2× scale
FFTMPI_NS::FFT3d *pf = fplan_->pf3d;
int permute = pf->permute;
pf3dgrf_ = new FFTMPI_NS::FFT3d(MPI_COMM_FFT, 2);
int fftsize, sendsize, recvsize;
pf3dgrf_->setup(2*pf->nfast, 2*pf->nmid, 2*pf->nslow,
                2*pf->in_ilo, 2*pf->in_ihi+1,
                2*pf->in_jlo, 2*pf->in_jhi+1,
                2*pf->in_klo, 2*pf->in_khi+1,
                2*pf->slow_ilo, 2*pf->slow_ihi+1,
                2*pf->slow_jlo, 2*pf->slow_jhi+1,
                2*pf->slow_klo, 2*pf->slow_khi+1,
                permute, fftsize, sendsize, recvsize);
// Allocate grf_ at 2× scale
int gcnt2x = (2*pf->slow_ihi-2*pf->slow_ilo+2) * ... ;  // 2× k-space pencil size
grf_ = new std::complex<Real>[8 * gcnt2x];
```

After `pf3dgrf_->compute()` in `InitGreen()`, the k-space element at global position
`(2*kx+px, 2*ky+py, 2*kz+pz)` in the 2× domain (where `(kx,ky,kz)` is the normal-domain
k-space position) lives at a computable local index within `grf_`.

---

## `GravityBoundaryFlag` sharing

`GravityBoundaryFlag` enum and `GetGravityBoundaryFlag()` are currently in
`block_fft_gravity.hpp`. For Phase 1, simplest path: duplicate in `fft_gravity.hpp`.
Consolidation (extract to `gravity.hpp` or a new `gravity_bc.hpp`) is a follow-up cleanup.

---

## Implementation Sequence

1. **Add `friend class ::FFTGravity;`** to `fftmpi/fft3d.h` (forward decl + friend)
2. **Duplicate `GravityBoundaryFlag`** in `fft_gravity.hpp` (or extract to shared header)
3. **Make `ExecuteForward()`/`ExecuteBackward()` virtual** in `FFTBlock` (athena_fft.hpp)
4. **Extend `FFTGravity` class**: new members + constructor/destructor
5. **Implement open BC**: `InitGreen()`, `LoadOBCSource()`, `RetrieveOBCResult()`, `MultiplyGreen()`
6. **Extend `FFTGravity::ApplyKernel()`**: add open (no-op) and disk kernel dispatch
7. **Implement disk BC `ExecuteForward()`/`ExecuteBackward()`** in `FFTGravity`
8. **Extend `FFTGravityDriver::Solve()`**: add open and disk dispatch paths
9. **Compile and run `jeans_3d.py`** regression (periodic) at 1/2/4 ranks to ensure no regression
10. **Test open BC** with an existing test problem (e.g., `poisson_open.py` if it exists, or a new jeans-in-open-domain test)
11. **Test disk BC** with an existing test problem

---

## Deferred: Shearing BC (Phase 2)

Shearing-periodic BC requires `RollUnroll()` and `ShearTimeShift()`, which operate on real-space
data before and after the FFT. In `BlockFFTGravity` these are per-meshblock operations with MPI
communication along y. Porting to `FFTGravity` (multi-meshblock) requires:
- Collecting the full-domain density into a format amenable to rolling
- Managing the fractional y-shift communication across meshblock boundaries
- Using `RemapFluxPlm`/`RemapFluxPpm` from reconstruction

This is architecturally non-trivial and deferred to Phase 2. `swing.py` regression (MPIFFT.md
Step 6) tests only `BlockFFTGravity` shearing; no shearing regression is required for Plan 4.

---

## Open Questions / Risks

| Issue | Severity | Mitigation |
|-------|----------|------------|
| z must be "slow" axis in `fplan_->pf3d` for disk BC | High | Enforce/verify at construction; disk sims typically have npz=1 |
| `MPI_COMM_FFT` scope: `pf3dgrf_` must use same communicator as `fplan_->pf3d` | Medium | Access via `pmy_driver_->MPI_COMM_FFT` (protected in `FFTDriver`) |
| `pf3dgrf_` index ranges: ensure output layout alignment with `fplan_->pf3d` | Medium | Use same permute, double the in/slow ranges — verify with unit test |
| `grf_` local index in `MultiplyGreen` when permute≠0 | Medium | Follow BlockFFTGravity's `grf_` indexing with adjusted offsets |
| Normalization factor for open BC (divide by 8) | Low | Apply same as BlockFFTGravity; verify in test |
| Buffer `out_` large enough for disk BC in_e_/in_o_ | Low | Resize in constructor if needed |

---

## Files Changed

| File | Change |
|------|--------|
| `src/fft/fftmpi/fft3d.h` | Add `friend class ::FFTGravity;` |
| `src/fft/athena_fft.hpp` | Make `ExecuteForward`/`ExecuteBackward` virtual |
| `src/gravity/fft_gravity.hpp` | Extend `FFTGravity`; add `GravityBoundaryFlag`; extend `FFTGravityDriver` |
| `src/gravity/fft_gravity.cpp` | Implement all new methods; extend `Solve()` dispatch |

`BlockFFTGravity` files are untouched in this plan.
