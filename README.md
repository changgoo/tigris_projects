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

---

## outputs/ — output format documentation

| File | Description |
|------|-------------|
| `zprof_outputs.md` | ZprofOutput expected column definitions |

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

---

## particles/ — particle system work

| File | Description |
|------|-------------|
| `overlap_context.md` | Session context: particle overlap logic, flags, CV checks, known issues |
| `accretion_conservation.md` | Accretion mass conservation: diagnosis and fixes for shear-periodic and pid=NEW bugs |
