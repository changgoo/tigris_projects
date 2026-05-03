# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A documentation and planning companion for [TIGRIS](https://github.com/PrincetonUniversity/tigris), a private CRMHD fork of Athena++. There is no buildable source here — all Markdown files are design notes, plans, reviews, and context documents that inform changes to the upstream TIGRIS checkout at `$HOME/tigris` (branch `tigris-master`).

## Navigation

Search before editing:

```bash
rg "term" .                              # find decisions, PR numbers, function names
git diff -- README.md fftmpi/ particles_p2p/   # review doc changes before committing
markdownlint '**/*.md'                   # catch heading/list/spacing issues (if available)
```

Key reference documents at root level:

| File | What it contains |
|------|-----------------|
| `code_structure.md` | TIGRIS source layout, execution flow, key physics, coding conventions |
| `task_flow.md` | Full per-cycle task-dependency graph (RK2 stages + OperatorSplitTaskList) |
| `README.md` | Index of all project folders with related PR/issue tables |

## Project folders

Each folder owns a focused topic. When adding a new one, create the directory, put one-topic-per-file Markdown inside, and add a short table entry + related PRs to `README.md`.

| Folder | Topic |
|--------|-------|
| `fftmpi/` | fftMPI migration, FFTGravity BCs, shearing remap |
| `fftmpi/plans/` | Sequential numbered exploration/design records (preserve numeric order) |
| `particles/` | Ghost particle boundary logic, accretion conservation |
| `particles_p2p/` | P2P refactor replacing `MPI_Allgatherv` in ghost-return paths |
| `fofc/` | First-order flux correction diagnostics and boundary conservation fix |
| `outputs/` | Output format notes (z-profile columns, etc.) |

## TIGRIS source conventions (for documents that reference upstream code)

- C++11, BSD 3-Clause, Athena++ style: `snake_case`, Doxygen comments
- Style check: `tst/style/check_athena_cpp_style.sh`
- Regression tests live under `tst/regression`; functional changes need coverage
- Python scripts: linted with `flake8`
- Reference issue numbers in commits: `Fixes #42`

When making technical claims about TIGRIS behavior, cite upstream file paths, function names, PR numbers, or issue numbers so the claim is auditable. Do not invent behavior — cross-check existing notes and, when in doubt, consult the upstream source.

## Document style

- One topic per Markdown file, lowercase `snake_case` filenames
- Sequential plans use two-digit numeric prefixes (`01_`, `02_`, …)
- Tables for PR/issue inventories; short paragraphs for prose
- Commits scoped to one documentation topic; subjects use short imperative or descriptive form
