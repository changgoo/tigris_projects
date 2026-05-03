# Repository Guidelines

## Project Structure & Module Organization

This repository is a documentation and planning companion for the private TIGRIS codebase. Keep top-level reference material in `README.md`, `code_structure.md`, and `task_flow.md`. Use project folders for focused work:

- `fftmpi/`: FFTGravity and fftMPI migration plans, reviews, PR text, and follow-up notes.
- `fftmpi/plans/`: numbered exploration and design records; preserve numeric ordering.
- `particles/` and `particles_p2p/`: particle-system context, conservation analysis, and P2P communication plans.
- `fofc/`: first-order flux correction diagnostics and fixes.
- `outputs/`: output-format notes such as z-profile column definitions.

Prefer one topic per Markdown file. When adding a new project area, create a directory and update `README.md` with a short table entry and related PRs or issues.

## Build, Test, and Development Commands

There is no local build for this notes repo. Common checks are documentation-oriented:

- `rg "term" .`: find existing decisions, file names, or PR references before editing.
- `git diff -- README.md fftmpi/ particles_p2p/`: review documentation changes before committing.
- `markdownlint '**/*.md'`: run if available to catch heading, list, and spacing issues.

Implementation validation belongs in the upstream TIGRIS checkout, referenced in `code_structure.md` as `$HOME/tigris`.

## Coding Style & Naming Conventions

Write concise Markdown with descriptive headings, short paragraphs, and tables for inventories or PR lists. Use lowercase snake_case file names such as `shearing_remap_plan.md`. For sequential plans, use two-digit prefixes like `01_explore_fft_wrappers.md`.

When documenting TIGRIS code, follow the project conventions already captured in `code_structure.md`: C++11, Athena++ style, `snake_case`, Doxygen comments, and regression coverage for functional changes.

## Testing Guidelines

For this repository, verify links, issue numbers, PR states, paths, and command examples. If a document describes code behavior, include enough upstream file/function names to make the claim auditable. For TIGRIS changes, note expected regression coverage under `tst/regression` and style checks such as `tst/style/check_athena_cpp_style.sh`.

## Commit & Pull Request Guidelines

Recent commits use short imperative or descriptive subjects, for example `Add particles_p2p/ project: plan to replace MPI_Allgatherv with P2P` and `Update particles_p2p/plan.md based on advisor review`. Keep commits scoped to one documentation topic.

Pull requests should summarize the document purpose, list changed files, link related TIGRIS PRs/issues, and call out whether any upstream code or regression tests were consulted.

## Agent-Specific Instructions

Do not invent TIGRIS behavior. Cross-check existing notes first, then cite upstream paths, functions, PRs, or issues when adding technical claims.
