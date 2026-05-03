# Goal

Migrate `athena_fft` (the general-purpose FFT wrapper) to use the newer `fftmpi` C++ library
as its backend, replacing the older `plimpton` C library. Then migrate `block_fft_gravity` to
use `athena_fft` instead of `block_fft`, so that gravity can benefit from multi-meshblock-per-rank
support. If the backend swap is not feasible, extend `block_fft` to support multiple meshblocks
per rank as an alternative, and document which approach is better.

# Background

Two FFT MPI wrappers exist in `src/fft/`:

| Directory | Language | Used by | Notes |
|-----------|----------|---------|-------|
| `plimpton/` | C | `athena_fft.hpp` / `AthenaFFTPlan` | Older version of FFTMPI |
| `fftmpi/`   | C++ | `block_fft.hpp` / `BlockFFT` | Newer version of FFTMPI; no license conflict |

The two libraries are algorithmically similar; `fftmpi` is a C++ rewrite of `plimpton`.

**Key difference in wrapper capabilities:**

- `athena_fft` supports a *cuboid* arrangement of meshblocks per rank — i.e., one MPI rank
  can own multiple meshblocks arranged in a 3D sub-grid (e.g., 2, 4, or 8 meshblocks/rank).
  See `FFTBlock` in `athena_fft.cpp` for how MeshBlock data is mapped into the FFT domain.
- `block_fft` supports only one meshblock per rank and has no serial (non-MPI) support.

The preference is to use `fftmpi` everywhere, while preserving the multi-meshblock capability
of `athena_fft`.

# Steps

1. **Explore code structure** — read `athena_fft.hpp`, `block_fft.hpp`, both library sources
   (`plimpton/`, `fftmpi/`), and `gravity/block_fft_gravity`. Study `FFTBlock` carefully to
   document how MeshBlock data is loaded and how the index mapping works for cuboid configurations.
   Save findings to `plans/plan1.md`.

2. **Document `athena_fft` fully** — write complete documentation of the `AthenaFFT` / `FFTBlock`
   class hierarchy, index mapping logic, and cuboid decomposition. Save to `AthenaFFT.md` in the
   repo root. Save plan notes to `plans/plan2.md`.

3. **Feasibility: replace `plimpton` with `fftmpi` in `athena_fft`** — assess whether the
   `plimpton` C API calls in `athena_fft` can be swapped for `fftmpi` C++ equivalents without
   architectural changes. Document findings in `plans/plan3.md`.
   - **If feasible**: proceed to Steps 4–6 on a new branch.
   - **If not feasible**: assess whether extending `block_fft` to support multiple meshblocks/rank
     is a better path. Document the trade-offs and the recommended approach in `plans/plan3.md`
     before proceeding.

4. **Implement `athena_fft_gravity`** — write a counterpart to `block_fft_gravity` that uses
   `athena_fft` instead of `block_fft`, on a new branch. Save design notes to `plans/plan4.md`.

5. **Design a unit test for Poisson solver** — write a simple test problem (e.g., eigenmode
   density → known analytic solution) that exercises `block_fft_gravity` and `athena_fft_gravity`
   and compares results at 1, 2, 4, and 8 meshblocks/rank. Check `tst/regression/scripts/tests/grav/`
   for existing test patterns to follow. Save to `plans/plan5.md`.

6. **Regression test with `swing.py`** — run `tst/regression/scripts/tests/grav/swing.py`
   (shearing-periodic BCs, currently supported only by `block_fft_gravity`) to verify no
   regression. Note: `athena_fft_gravity` does not need to support shearing-periodic BCs at
   this stage. Save results/notes to `plans/plan6.md`.

7. **Review and update documentation** — review all changes, update `AthenaFFT.md` and `CLAUDE.md`
   as needed. Save to `plans/plan7.md`.

8. **[Deferred] Remove `block_fft`** — once Steps 1–7 are merged and verified in production runs,
   open a follow-up PR to remove `block_fft` and `plimpton/`. Do not include in the current
   PR series.

# Rules

1. Create `plans/planN.md` (N = step number) for each step before writing code. Ask clarifying
   questions at the start of each step if needed.
2. Any significant code modification must be on a new branch (never commit directly to
   `tigris-master`).
3. Prefer `fftmpi` over `plimpton` wherever feasible; document the reason when making a choice.
