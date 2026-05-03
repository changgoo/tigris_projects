# tigris_projects

Notes, plans, and documentation generated during development of
[TIGRIS](https://github.com/PrincetonUniversity/tigris) (a CRMHD fork of Athena++).

---

## Root — reference documents

| File | Description |
|------|-------------|
| `code_structure.md` | Source layout, execution flow, key physics, coding conventions |
| `task_flow.md` | Full task-dependency graph for CRMHD + shearing-box (RK2 + ops-split) |

---

## fofc/ — First-Order Flux Correction

| File | Description |
|------|-------------|
| `conservation_fix.md` | FOFC conservation violation: root cause, options, and fix |

**Related PRs**

| # | Title | State |
|---|-------|-------|
| [#255](https://github.com/PrincetonUniversity/tigris/pull/255) | FOFC for CR and shear BCs | MERGED |
| [#260](https://github.com/PrincetonUniversity/tigris/pull/260) | Multiple Fixes for FOFC | MERGED |
| [#263](https://github.com/PrincetonUniversity/tigris/pull/263) | Implement implicit solver and optimize fofc array handling | MERGED |
| [#279](https://github.com/PrincetonUniversity/tigris/pull/279) | Guard FOFC test for isothermal EOS | MERGED |

---

## outputs/ — output format documentation

| File | Description |
|------|-------------|
| `zprof_outputs.md` | ZprofOutput expected column definitions |

**Related PRs**

| # | Title | State |
|---|-------|-------|
| [#230](https://github.com/PrincetonUniversity/tigris/pull/230) | Add heating calculation for CR in ZprofOutput | MERGED |
| [#233](https://github.com/PrincetonUniversity/tigris/pull/233) | Additional zprof outputs | MERGED |
| [#251](https://github.com/PrincetonUniversity/tigris/pull/251) | Add the effective CR transport velocity calculation in zprof output | MERGED |
| [#265](https://github.com/PrincetonUniversity/tigris/pull/265) | Add CR flux decomposition along B to zprof | MERGED |
| [#272](https://github.com/PrincetonUniversity/tigris/pull/272) | Added work terms due to magnetic fields (pressure/tension separated) | MERGED |

---

## fftmpi/ — fftMPI migration and FFTGravity project

Goal: migrate `athena_fft` from the Plimpton C backend to fftMPI, implement
`FFTGravity` (multi-meshblock Poisson solver with periodic / open / disk / shearing BCs).

| File | Description |
|------|-------------|
| `project_goals.md` | Project overview: goal, background, step-by-step plan |
| `athena_fft_classes.md` | Full class documentation: FFTDriver, FFTBlock, AthenaFFTIndex, FFTGravity |
| `optimization_context.md` | Handoff context for the shearing-remap optimization work |
| `shearing_remap_plan.md` | Design plan for replacing global row-exchange with MeshBlock-local remap |
| `shearing_remap_review.md` | Review of the shearing remap refactor plan |
| `pr_shearing_remapper.md` | PR title/body for the ShearingRemapper change |
| `cleanup_notes.md` | Follow-up cleanup ideas after the open/disk BC PR |

### fftmpi/plans/ — sequential exploration and design notes

| File | Description |
|------|-------------|
| `01_explore_fft_wrappers.md` | Code exploration: FFT MPI wrapper structures and API comparison |
| `02_document_athena_fft.md` | Plan and summary for writing `athena_fft_classes.md` |
| `03_plimpton_to_fftmpi_feasibility.md` | Feasibility assessment: replacing Plimpton C with fftMPI C++ |
| `04_fft_gravity_bc_parity.md` | Design: extending FFTGravity to open and disk BCs |
| `06_swing_regression.md` | Regression test plan using `swing.py` (shearing-periodic BC guard) |

**Related PRs**

| # | Title | State |
|---|-------|-------|
| [#280](https://github.com/PrincetonUniversity/tigris/pull/280) | Extend AthenaFFT gravity for open and disk boundary conditions | MERGED |
| [#281](https://github.com/PrincetonUniversity/tigris/pull/281) | FFT gravity: shearing-periodic BC + decomp-independent disk BC + streamlining | MERGED |
| [#282](https://github.com/PrincetonUniversity/tigris/pull/282) | Use column-wise shearing remap for AthenaFFT gravity | MERGED |
| [#283](https://github.com/PrincetonUniversity/tigris/pull/283) | ShearingRemapper: replace global-row exchange with MeshBlock-local remap for FFT gravity | MERGED |
| [#277](https://github.com/PrincetonUniversity/tigris/pull/277) | Complete FFT gravity migration: fftMPI backend + full BC support in FFTGravity | OPEN |

---

## particles/ — particle system work

| File | Description |
|------|-------------|
| `overlap_context.md` | Session context: particle overlap logic, flags, CV checks, known issues |
| `accretion_conservation.md` | Accretion mass conservation: diagnosis and fixes for shear-periodic and pid=NEW bugs |

**Related PRs**

| # | Title | State |
|---|-------|-------|
| [#234](https://github.com/PrincetonUniversity/tigris/pull/234) | Complex Particles | MERGED |
| [#235](https://github.com/PrincetonUniversity/tigris/pull/235) | Accretion | MERGED |
| [#238](https://github.com/PrincetonUniversity/tigris/pull/238) | Mass Return | MERGED |
| [#267](https://github.com/PrincetonUniversity/tigris/pull/267) | Add particle momentum and kinetic energy for history | MERGED |
| [#268](https://github.com/PrincetonUniversity/tigris/pull/268) | Ghost particle output in parbin + Python flag helper | MERGED |

---

## particles_p2p/ — Refactor particle return communication

Goal: replace per-MeshBlock collective communication in `ExchangeGhostAccretionDelta()`
and `MassReturn::CollectParticlesInfo()` with task-flow-aware, rank-level exchanges,
removing the O(nranks) scaling path where local communication is sufficient and
supporting multi-MeshBlock/rank particle return paths. Non-uniform MeshBlocks/rank is a
desirable follow-up if it does not fall out naturally.

| File | Description |
|------|-------------|
| `context.md` | Problem background, current code analysis, existing P2P infrastructure |
| `plan.md` | Active implementation design: task-flow hook, rank-level exchange, origin tracking, boundary-safe mass return |
| `review.md` | First review of the original P2P plan |
| `review_taskflow.md` | Second review focused on task flow, multi-MeshBlock/rank, and boundary compatibility |
| `review_final.md` | Review of rev 3: five gaps identified (task flow, exchange manager, REFRESH_MR_GHOSTS, INTERACT split, r_return==0 form) |
| `review_rev4.md` | Review of rev 4: caught task-list ownership issue and neighbor-stencil scope decisions |
| `plan_rev2_obsolete.md` | Superseded rev-2 plan kept only for historical comparison |
| `plan_rev3_obsolete.md` | Superseded rev-3 plan kept only for historical comparison |

**Related issues / PRs**

| # | Title | State |
|---|-------|-------|
| [#269](https://github.com/PrincetonUniversity/tigris/issues/269) | Replace MPI_Allgatherv with point-to-point comm for ghost particle return data | OPEN |
