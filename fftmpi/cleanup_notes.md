# Cleanup Notes

Follow-up ideas from reviewing the FFTGravity open/disk boundary-condition PR.

1. Use one gravity boundary-condition parser.
   `GetFFTGravityBoundaryFlag()` duplicates `GetGravityBoundaryFlag()`. Since
   `GravityBoundaryFlag` now lives in `gravity.hpp`, one shared parser can serve both
   FFT gravity solvers.

2. Factor source-loading code in `FFTGravityDriver::Solve()`.
   The open and periodic/disk paths duplicate density loading and particle-density
   accumulation. A helper such as `LoadGravitySource(...)` would make the solve flow
   easier to read.

3. Share the disk kernel formula.
   MPI and serial disk BC paths duplicate the even/odd kernel expression with only
   indexing differences. A small helper returning `kernel_e` and `kernel_o` from
   global `(kx, ky, kz)` indices would reduce maintenance risk.

4. Skip disk layout repacking when possible.
   The disk path always repacks AthenaFFT layout into physical brick layout. The
   identity-layout case could avoid this copy.

5. Share regression-test boilerplate.
   `poisson_fft.py` and `poisson_fft_disk.py` have nearly identical build/run/analyze
   structure. A local helper parameterized by boundary condition, problem ID, and
   tolerances would make future FFT gravity tests cheaper to add.

6. Clean stale comments.
   `permute_disk_` is fixed to physical-layout `permute=2`; comments should avoid
   suggesting it is derived from `f_in_->ploc[2]`.
