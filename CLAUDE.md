# CLAUDE.md — ARTnet

Guidance for Claude Code working in this repo. Keep this file the first thing a fresh session reads.

## 1. What ARTnet is

ARTnet (R package, GitHub `EpiModel/ARTnet`) parameterizes epidemic and network models for men who have sex with men (MSM) from the ART-Net 2017–2018 egocentric network survey. It is a production dependency in the **EpiModelHIV modeling ecosystem**.

The package takes data from the private `ARTnetData` package and produces three objects consumed downstream:

1. `build_epistats()` — epidemic parameters (HIV prevalence by race/age, HIV mortality, etc.).
2. `build_netparams()` — estimated per-respondent / per-partnership network statistics (degree, mixing, durations).
3. `build_netstats()` — population-scaled target statistics formatted for ERGM estimation in EpiModelHIV.

Downstream consumers include `~/git/EpiModelHIV-Template/` (template simulation workflow) and `~/git/EpiModelHIV-p/` (model package). Any change to ARTnet's output structure or numerics is load-bearing for an entire simulation ecosystem.

**Critical constraint — backward compatibility:** any refactor must be backward compatible by default. Existing EpiModelHIV users calling `build_netstats()` without new arguments must get byte-identical output. **Adding** fields to the return list is safe. **Removing** or **renaming** fields is not. Verified by the snapshot harness in §5.

## 2. Core workflow files

- [R/NetParams.R](R/NetParams.R) (~1,400 lines after refactor) — fits per-respondent and per-partnership models for each target statistic, layer by layer (main / casl / inst). Both legacy univariate fits and joint g-computation models are produced.
- [R/NetStats.R](R/NetStats.R) (~1,000 lines after refactor) — consumes `netparams` output, builds the synthetic population, and scales per-respondent stats up to population-level ERGM target stats. Branches by `method` argument.
- [R/EpiStats.R](R/EpiStats.R) — act rates, condom use, starting HIV prevalence.
- [R/ARTnet-package.R](R/ARTnet-package.R), [R/globals.R](R/globals.R) — package docs and global-binding declarations.

## 3. Public API

Exported functions (see [NAMESPACE](NAMESPACE)):
`build_epistats`, `build_netparams`, `build_netstats`, `make_race_combo`, `reweight_age_pyr`, `trim_epistats`, `trim_netstats`, `update_asmr`.

### Canonical legacy pipeline (byte-identical to pre-refactor):

```r
epistats  <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                            init.hiv.prev = c(0.33, 0.137, 0.084),
                            race = TRUE, time.unit = 7)
netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
netstats  <- build_netstats(epistats, netparams,
                            expect.mort = 0.000478213, network.size = 5000)
```

### Joint g-computation pipeline (opt-in, post-refactor):

```r
netparams <- build_netparams(epistats, smooth.main.dur = TRUE,
                             method = "joint",                 # joint Poisson + binomial fits
                             duration.method = "joint_lm")     # log-linear duration regression
netstats  <- build_netstats(epistats, netparams,
                            expect.mort = 0.000478213, network.size = 5000,
                            method = "joint",                  # g-computation aggregation
                            target_pop = NULL)                 # optional post-stratification
```

`target_pop` accepts a named list (per-marginal overrides), a data.frame (one row per node, full bypass of attribute sampling), or a character string (reserved for future built-in geography bundles — currently raises a clear not-yet-implemented error). All defaults preserve byte-identical legacy behavior.

Pre-built example objects ship in `inst/` (`epistats-example.rds`, `netstats-example.rds`) for users without `ARTnetData` access.

## 4. The joint g-computation refactor (Phase 1 complete)

### 4.1 Why this happened

The legacy approach fit a separate univariate Poisson / binomial / linear regression for each target statistic, one attribute at a time. Three problems followed:

1. **Marginal-vs-joint bias.** Each per-attribute estimate carried ARTnet's conditional joint distribution of the *other* attributes baked in. When `build_netstats()` applied these to a synthetic target population whose joint attribute distribution differed from ARTnet's (any city-specific MSM target, or any AMIS-projected population), bias propagated systematically.
2. **Internal inconsistency across target stats.** Edges from `md.main * num / 2` did not equal edges from `Σ_r table(race) * nf.race[r] / 2`. The legacy `edges.avg` argument was a tacit acknowledgement.
3. **Patchwork target population.** `build_netstats()` mixed reference sources (NCHS age pyramid + `ARTnetData::race.dist` + ARTnet's own degree/role distributions) with no coherent specification of who the target was.

The refactor replaces the per-attribute univariate fits with joint Poisson and binomial GLMs (and a log-linear regression for durations) that condition on all attributes simultaneously, then aggregates per-respondent or per-dyad predictions against an explicitly-specified synthetic target via g-computation. Direct standardization in the sense of Hernán & Robins (2020).

### 4.2 What got built (PRs and what they did)

| PR | Issue | Contribution |
|---|---|---|
| #66 | — | Backward-compat snapshot harness + pinned `EpiModelHIV-Template/R/A-networks/` reference |
| #67 | #61 | Joint Poisson GLMs per layer in `build_netparams` (formation) |
| #68 | #62 | G-computation in `build_netstats` for edges / nodefactor / concurrent |
| #69 | #63 (phases 1–2) | Joint logistic for nodematch, joint Gaussian for absdiff |
| #70 | — | Package housekeeping: R CMD check 0/0/0, silent tests, GitHub Actions CI |
| #71 | #63 (phase 3) | `duration.method` flag with `empirical` and `joint_lm` options |
| #74 | #73 | Synth-population aggregation of joint_lm durations in `build_netstats` |
| #75 | #59 | Validate `init.hiv.prev` length against `race` flag |
| #76 | #65 | Method-comparison validation suite + four-city comparison report |
| #77 | #64 | `target_pop` argument: list / data.frame / character forms |
| #78 | — | Research report: `inst/validation/method_refactor_report.md` |

### 4.3 What `method = "joint"` actually fits

Per network layer (main, casual, one-time), under `method = "joint"`:

- **Joint Poisson GLM** of partnership count on `age.grp + sqrt(age.grp) + race + cross-layer-degree + hiv2`, with AIC-based selection over `age.grp:race` and `age.grp:cross-layer-degree` interactions. Stored at `netparams$<layer>$joint_model`.
- **Joint binomial GLM** of the concurrency indicator (`deg > 1`) on the same RHS, for main and casual. Stored at `netparams$<layer>$joint_concurrent_model`.
- **Joint binomial GLMs** for `same.age.grp` and `same.race` (when race-stratified), fit on long-form partnership data with ego attributes only on the RHS. Stored at `netparams$<layer>$joint_nm_age_model` and `joint_nm_race_model`.
- **Joint Gaussian regressions** for `ad` (absdiff age) and `ad.sr` (absdiff sqrt-age). Stored at `joint_absdiff_age_model` and `joint_absdiff_sqrtage_model`.
- **Joint log-linear regression** of `log(duration.time)` among ongoing partnerships under `duration.method = "joint_lm"`, with ego + partner + matching terms. Stored at `joint_duration_model`.

Under `method = "joint"` in `build_netstats`, edges become `sum(pred_deg) / 2`, `nodefactor_<attr>[level] = sum(pred_deg | attr == level)`, and dyad-level statistics aggregate `pred_deg * pred_<dyad>` per stratum / 2. Internal consistency `Σ_<level> nodefactor_<attr>[level] = 2 × edges` holds to machine precision.

### 4.4 The duration estimand argument (worth understanding)

TERGM dissolution uses `~offset(edges) + offset(nodematch("age.grp", diff = TRUE))` — a constant-hazard (geometric) per-stratum offset. Under that simulation:
- Mean simulated full partnership duration **equals** mean simulated age of extant ties at cross-section, by the inspection paradox / memoryless property, *but only if the underlying duration distribution is exponential*.
- Under non-exponential durations (Weibull `k ≠ 1`), these two quantities **diverge** by a factor of `(1 + CV²) / 2`.
- Per-stratum Weibull shape parameters from a length-bias-corrected MLE (development-time exploration) showed `k ≈ 0.6` for casual partnerships (decreasing hazard) and `k > 2` for older-matched main partnerships (increasing hazard).

Implication: the geometric simulation can match observed mean age of extant ties exactly, but cannot honor a distinct mean full duration estimate. The `empirical` and `joint_lm` methods both target the cross-sectional age statistic — *not* mean full partnership duration — because that is the quantity TERGM can faithfully reproduce. A proper Weibull-based duration estimate was attempted and discarded for this reason; the attempt also showed that naive `survreg` on ARTnet's elapsed-duration data is catastrophically biased without length-bias correction (matched.5 main partnerships extrapolating to 1.5M weeks).

### 4.5 The geogYN insight (worth understanding)

The legacy pipeline used `geogYN` as a **main-effect covariate** in each per-attribute regression, then predicted at `geogYN = 1` for the city of interest. This was a defensible design choice — filtering ARTnet to Atlanta would have left only n = 206 respondents — but it did not accomplish demographic rebalancing, because:

| Population | % Black |
|---|---:|
| ARTnet full sample | 5.5% |
| ARTnet Atlanta sub-sample | 12.6% |
| Atlanta MSM target population (`ARTnetData::race.dist`) | 51.5% |

The Atlanta sub-sample is more Black than the full sample (so `geogYN = 1` does shift toward Atlanta-relevant rates) but is still far from Atlanta's actual MSM composition. Online-recruitment selection bias on race operates within geographic strata, not just at the national level. The per-respondent rate (`md.main = 0.398`) thus reflects the Atlanta sub-sample's race composition, while the synthetic population is drawn from a *different* source (`ARTnetData::race.dist`). The two halves of the legacy pipeline implicitly assumed different race compositions and never reconciled. The joint g-computation closes the loop by aggregating per-respondent predictions across the synth's actual joint distribution.

### 4.6 Headline empirical findings

From `inst/validation/method_comparison.md` (4 city scenarios × ~96 cells each = 384 target-stat cells total):

- **60% of cells (232 / 384) shift > 5%** between joint and existing methods.
- **Cross-city ordering tracks demographic distance from ARTnet sample**: Seattle (closest, 87.4% W/Other vs ARTnet's 80.7%) shifts 46% of cells; Atlanta / Boston / Chicago all shift 64–67%.
- **Largest single shifts**: dissolution durations in matched-and-old strata (−47% to −66%), one-time-partnership nodematch in oldest age groups (−51%), high-deg.main casual nodefactor (+40%).
- **Atlanta main-edges decomposition**: legacy `0.398 × 5000 / 2 = 995`; joint `0.338 × 5000 / 2 = 845`. The −15% shift comes entirely from race-composition aggregation: ARTnet 80.7% W/O × 0.420 main degree → 0.398; Atlanta 51.5% Black × 0.272 → 0.338.
- **Coefficient strengthening under joint adjustment** (main Poisson): `deg.casl` slope −0.24 → −0.55; `hiv2` slope +0.09 → +0.25; AIC selects `age.grp:deg.casl` interaction (+0.10).
- **End-to-end ERGM convergence**: clean for both methods on `EpiModelHIV-Template`-equivalent estimation. Static `netdx` matches all formation targets within `|Z| ≤ 2.05` and `|% diff| ≤ 4.2%` over 1000 sims.

## 5. Validation infrastructure (`inst/validation/`)

- `validate_backward_compat.R` — `capture_snapshot()` and `compare_to_snapshot()`. Iterates over `PARAM_SETS` (Atlanta+race, national no-geog, Atlanta no-race) and full-object-diffs `netparams` and `netstats` at machine precision. Snapshot `.rds` files at `inst/validation/snapshots/` are gitignored (~12 MB each). Strips `joint_*_model` additive fields before comparison.
- `method_comparison.R` — `compare_methods()`, `summarize_comparison()`, `render_comparison_report()`. Runs both methods across four city scenarios (Atlanta, Boston, Chicago, Seattle), produces long-format comparison data frame with all per-stat shifts, renders the markdown table.
- `method_comparison.md` — committed canonical comparison report (regenerated by the harness).
- `method_refactor_report.md` — full research-style writeup (~4,500 words, 6 tables, 3 figures). Suitable for export to Word for collaborator review.
- `figures.R` + `figures/` — ggplot2 rendering script for the report's figures (joint-vs-existing scatter, dissolution-duration comparison, race-composition bar chart). PNGs gitignored (regenerate via `render_all_figures()`).
- `epimodelhiv_template_ref/` — verbatim pinned copies of `EpiModelHIV-Template/R/A-networks/{initialize,model_main,model_casl,model_ooff}.R` documenting the public contract.
- `netstats_contract.md` — distilled list of `netstats` fields the template ERGM specs read.

### Snapshot regression workflow

```r
# Step A: capture once on pre-refactor main (already done)
devtools::load_all()
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
capture_snapshot()

# Step B: verify on any subsequent change
compare_to_snapshot()                      # default args path
compare_to_snapshot(method = "existing")   # explicit legacy combo
# Both must report ALL MATCH 3/3 before merging anything.
```

### Method-comparison workflow

```r
source(system.file("validation/method_comparison.R", package = "ARTnet"))
res <- compare_methods()           # ~30s on 4 scenarios × 2 methods
summarize_comparison(res)
render_comparison_report(res)      # writes inst/validation/method_comparison.md
```

The `expect.mort` value used in the comparison harness is `0.00025`, slightly lower than the EpiModelHIV-Template default `0.000478213`. This is because Seattle's empirical matched.5 main duration is 1381 weeks (small ARTnet sub-sample noise) and the standard mortality rate trips `dissolution_coefs()`'s competing-risk-of-departure check. The lower rate accommodates Seattle uniformly across cities for the comparison.

## 6. Open work

- **Issue #72 — formation-stat sampling bias.** ARTnet's cross-sectional sampling design exposes formation-stat estimates to length-biased sampling (on partnership-pair targets like `nodematch` and `absdiff`) and to partnership-count truncation at 2 main / 3 casual partners (on degree-based targets). Neither is corrected by the joint g-computation. The first concern is partially addressable by restricting partnership-pair fits to ongoing partnerships only; the second requires a truncated-Poisson likelihood. Has not been started.
- **Geography-named `target_pop` bundles.** The character form of `target_pop` (`"atlanta"`, `"us_msm_male"`) raises a not-yet-implemented error. Implementation is a lookup table from name → list using NCHS age pyramid + `ARTnetData::race.dist` (no new data needed, pure code change). Convenience for users.
- **Methods paper.** `inst/validation/method_refactor_report.md` is currently a research report. To become a publishable methods paper it would need: a simulation study with known ground truth quantifying joint-vs-marginal bias as a function of population divergence; deeper engagement with prior literature (Krivitsky & Morris, ergm.ego); formal statistical claims (variance, robustness); and a downstream EpiModelHIV-p simulation showing differential incidence over a multi-year horizon. Path A (applications note, ~3 weeks) and Path B (full methods paper, 3–6 months) discussed in chat.

## 7. Repo conventions

- **R style:** package-standard roxygen2 docs; 2-space indentation; base R + `dplyr` mix.
- **Function naming:** snake_case, verb-first (`build_*`, `trim_*`, `update_*`).
- **Lint:** [.lintr](.lintr) — 120-char lines, cyclocomp and object-name linters disabled.
- **Roxygen:** markdown mode (`Roxygen: list(markdown = TRUE)` in DESCRIPTION).
- **Dependencies:** `Imports` is `dplyr`; `Depends` is `EpiModel`. `Suggests` is `ARTnetData`, `knitr`, `rmarkdown`, `testthat`. `EpiModelHIV` was removed from Suggests in #70 because the package has no runtime dependency on it (only docstring references).
- **Testing:** `tests/testthat/` for unit tests in proper `test_that()` form (currently 571 assertions across 9 test files). `tests/workflows/` for end-to-end ERGM estimation scripts that require `EpiModelHIV` and produce console output — not picked up by `testthat::test_local()`.
- **ARTnetData access:** the data package is private and requires a `GITHUB_PAT`. Functions check for it with `system.file(package = "ARTnetData") == ""`.

## 8. Running checks

```r
devtools::document()                       # regenerate man/ from roxygen
devtools::test()                           # run testthat suite (571 assertions)
devtools::check(error_on = "never")        # full R CMD check (target: 0/0/0)
```

## 9. CI

GitHub Actions workflow at `.github/workflows/R-CMD-check.yaml` runs `R CMD check` on every push and PR.

- **Single OS / R combination**: `ubuntu-latest` × R release. The package is not CRAN-published; multi-OS portability matrix would be overkill.
- **`error-on: "warning"`** explicitly set so any new WARNING fails the build.
- **Private-repo access**: requires a fine-grained PAT with read access to `EpiModel/ARTnetData`, configured as repository secret `EPIMODEL_PAT`. Falls back to `GITHUB_TOKEN` if absent (so the workflow file is valid on forks). Without the secret, dependency resolution fails because `ARTnetData` is private.
- **`ARTnetData` is installed** in CI so that joint-fit functionality is exercised, not skipped.

## 10. Git / PR workflow

- Branch off `main` with descriptive name (`feature/...`, `fix/...`, `chore/...`, `docs/...`).
- Commits per logical unit. Match recent commit message style.
- PR description: summary, design decisions, validation results, test plan checklist.
- **Auto-close keywords are tricky.** GitHub parses `close[sd]?`, `fix(e[sd])?`, `resolve[sd]?` followed by `#N` *anywhere in the PR body or any commit message*, including in relative-clause references like "the PR that closes #N". Two issues (#63, #69) were inadvertently closed early by such phrasing during the refactor. Going forward: grep PR bodies for `(close|closes|closed|fix|fixed|fixes|resolve|resolves|resolved)\s+#\d+` before submitting; verify each match references the PR's own scope. Use "lands", "addresses", or "ships" for forward-references to other PRs.
- Don't push failing CI.
- Don't auto-commit during interactive sessions unless asked.

## 11. External references

### Sister project: ARTnetPredict
The detailed methodological analysis that motivated this refactor lives in `~/git/ARTnetPredict/`:
- `Plan/Step08-BuildLog.md` — covariate-shift analysis, raw AMIS vs model.
- `Plan/Step09-Plan.md` §2 — most directly relevant.
- `CLAUDE.md` — sister-project context.

### Methodological references
- Hernán MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC, 2020. (Chapter 13: g-computation / direct standardization.)
- Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, Jenness SM. Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study. *Epidemics* 2020; 30: 100386. — canonical ARTnet methodology paper.
- Krivitsky PN, Morris M. Inference for social network models from egocentrically sampled data. *Annals of Applied Statistics* 2017; 11(1): 427–455.

## 12. Things to flag to the PI if encountered

- Joint GLM convergence failures (may need regularization or model simplification).
- Snapshot regression failure under default args or `method = "existing"` — that's a backward-compat break. Stop and investigate.
- Any change to the legacy output structure (field names, types, byte values) — even additive fields should be vetted.
- Ambiguity about which interaction terms to include in joint fits — the AIC-based selection in the existing code is the default, but specific applications may want hand-specified terms.
- Issues that touch sampling design (length-bias, truncation) — those go to #72.
