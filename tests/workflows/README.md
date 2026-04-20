# Manual integration workflows

These scripts exercise the full ARTnet → EpiModelHIV estimation pipeline:
`build_epistats()` → `build_netparams()` → `build_netstats()` → `netest()` →
`netdx()`. They live outside `tests/testthat/` because:

- Each run takes tens of seconds to minutes (ERGM estimation + diagnostics).
- They require both `ARTnetData` and `EpiModelHIV` installed.
- They produce interactive `print()` / `plot()` output useful for human
  inspection but noisy in an automated test harness.

`devtools::test()` / `test_local()` (and the `tests/testthat.R` runner
invoked by `R CMD check`) only scan `tests/testthat/`, so these scripts
are not run automatically.

Run them manually after a non-trivial refactor touching `build_netparams()`
or `build_netstats()`:

```r
source("tests/workflows/workflow-standard.R")
source("tests/workflows/workflow-cea.R")
```

Inspect the `netdx` plots and printed diagnostics to confirm the ERGMs
estimate cleanly.
