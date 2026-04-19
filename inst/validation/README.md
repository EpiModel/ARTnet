# Validation Infrastructure

This directory supports the joint g-computation refactor (issues
[#61](https://github.com/EpiModel/ARTnet/issues/61)–[#65](https://github.com/EpiModel/ARTnet/issues/65))
by giving us two things a standard `testthat` suite cannot:

1. A **byte-for-byte reference snapshot** of the output that
   `build_netparams()` and `build_netstats()` produce on the pre-refactor
   `main` branch, captured once and compared against on every subsequent
   commit.
2. A pinned copy of the **downstream consumer** code
   (`EpiModelHIV-Template/R/A-networks/`) so we always know exactly which
   fields of `netstats` must remain stable — no guessing.

## Files

- `epimodelhiv_template_ref/` — verbatim copies of
  `~/git/EpiModelHIV-Template/R/A-networks/{initialize,model_main,model_casl,model_ooff}.R`
  taken on 2026-04-19. These are the ERGM specifications that consume
  `netstats`; they define the backward-compatibility contract. Do not edit
  unless the upstream template changes.
- `netstats_contract.md` — distilled list of exactly which `netstats` fields
  the template scripts read, by layer.
- `validate_backward_compat.R` — `capture_snapshot()` and
  `compare_to_snapshot()` functions. Run `capture_snapshot()` on pre-refactor
  `main`; run `compare_to_snapshot()` on the refactor branch with
  `method = "existing"` (or equivalent default) and expect zero diffs.
- `snapshots/` — created on first capture run. `.gitignore`d by default;
  the captured `.rds` files are large and should not be checked in.

## Workflow

### Step A — Before starting the refactor (pre-capture)

On the pre-refactor `main` branch, with `ARTnetData` installed:

```r
devtools::load_all()  # from the ARTnet repo root
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
capture_snapshot()    # writes inst/validation/snapshots/*.rds
```

This saves one snapshot per parameter set (see the `PARAM_SETS` list in
`validate_backward_compat.R`). Commit the snapshot files only if they are
small enough; otherwise keep them locally and rely on a hash digest that
**is** committed.

### Step B — During/after the refactor (compare)

On the refactor branch, with the new joint-GLM code in place:

```r
devtools::load_all()
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
compare_to_snapshot(method = "existing")
```

The call must report `ALL MATCH` before the PR is considered mergeable.
Any field-level diff is a backward-compatibility regression.

## Why not just `testthat::expect_equal()`?

Two reasons:
1. These runs require `ARTnetData` (private) and take minutes to execute —
   they do not belong in CI.
2. `testthat` snapshots are text-based and don't roundtrip well for deeply
   nested lists containing S3 objects (`glm`, `lm`, `dissolution_coefs`).
   `saveRDS()` + `all.equal()` is the simplest reliable approach here.

Unit tests for individual joint-GLM behaviors (convergence, marginal
recovery, coefficient sanity — see CLAUDE.md §4.5) should still live in
`tests/testthat/` as normal.
