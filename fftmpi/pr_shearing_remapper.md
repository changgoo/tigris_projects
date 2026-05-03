## Title

ShearingRemapper: replace global-row exchange with MeshBlock-local remap for FFT gravity

## Base branch

optimize-athenafft-shear-column

## Body

## Summary

- Extracts `BlockFFTGravity::RollUnroll()` into a new `ShearingRemapper` helper
  class (`src/gravity/shearing_remap.{hpp,cpp}`)
- Replaces `FFTGravity`'s global y-row exchange (`ExchangeShearingRows` /
  `ApplyShearingSource` / `RetrieveAppliedShearingResult`) with the same helper,
  eliminating the architectural bottleneck causing the ~14 ms gap between `fft`
  and `blockfft` gravity paths
- Adds rank-aware communication: same-rank MeshBlock copies go direct (no MPI
  self-send); cross-rank uses one `Irecv`/`Send`/`Waitall` per operation with a
  reserved physics tag
- Aggregates per-segment MPI into per-rank buffers to avoid tag pressure at
  high MB/rank counts
- Removes ~500 lines of old global-row remap code from `fft_gravity.cpp`

## Phase-gate validation results

| Phase | Gate | Result |
|-------|------|--------|
| 1 — Extract to ShearingRemapper | Bit-exact `phi` vs pre-refactor `blockfft` at 1 MB/rank | **BIT-EXACT** |
| 2 — Rank-aware communication | `BlockFFTGravity` timing unchanged; `ShearSourceMax` not regressed | **PASS** (20% faster `ShearSourceMax`) |
| 3 — FFTGravity source path | `ShearSourceMax ≤ 24 ms` | **PASS** |
| 4 — FFTGravity result path | `RetrieveMax ≤ 13 ms`; `SelfGravity` within 5% of `blockfft` | **PASS** |
| 5 — Remove old code | Old global-row remap removed; no regressions | **PASS** |

## Test plan

- [ ] Build with `--grav=blockfft` and `--grav=fft`
- [ ] Run 1 MB/rank case: `compare_phi.py` reports BIT-EXACT vs pre-refactor baseline
- [ ] Run timing benchmark: `fft_gravity_timing.py` shows `fft` path within gates above
- [ ] Confirm no regressions in `blockfft` timing baseline
