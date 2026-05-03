# Step 6: Swing Regression

## Goal

Run `tst/regression/scripts/tests/grav/swing.py` to verify that the existing
`block_fft_gravity` shearing-periodic path still works after the AthenaFFT gravity
work. AthenaFFT gravity does not need shearing-periodic support at this stage.

## Test Target

- Regression script: `tst/regression/scripts/tests/grav/swing.py`
- Configure target from script: `mpi`, `fft`, `prob=msa`, `grav=blockfft`
- Runtime shape:
  - 4 MPI ranks
  - resolutions 32 and 64
  - restart comparison from `SwingAmplification_64.00001.rst`

## Notes

- This is intentionally a regression guard for `block_fft_gravity`, not a new
  AthenaFFT gravity feature test.
- Results will be recorded below after the run.

## Results

Pending.
