# CLAUDE.md — ARTnet

Guidance for Claude Code working in this repo. Keep this file the first thing a fresh session reads.

## 1. What ARTnet is

ARTnet (R package, GitHub `EpiModel/ARTnet`) parameterizes epidemic and network models for men who have sex with men (MSM) from the ART-Net 2017–2018 egocentric network survey. It is a production dependency in the **EpiModelHIV modeling ecosystem**.

The package takes data from the private `ARTnetData` package and produces three objects consumed downstream:

1. `build_epistats()` — epidemic parameters (HIV prevalence by race/age, HIV mortality, etc.).
2. `build_netparams()` — estimated per-respondent network statistics (mean degrees, match/factor rates, durations).
3. `build_netstats()` — population-scaled target statistics formatted for ERGM estimation in EpiModelHIV.

Downstream consumers include `~/git/EpiModelHIV-Template/` (template simulation workflow) and `~/git/EpiModelHIV-p/` (model package). Any change to ARTnet's output structure or numerics is load-bearing for an entire simulation ecosystem.

**Critical constraint — backward compatibility:** any refactor must be backward compatible by default. Existing EpiModelHIV users calling `build_netstats()` without new arguments must get byte-identical output. **Adding** fields to the return list is safe. **Removing** or **renaming** fields is not.

## 2. Core workflow files (where most work happens)

This session's refactor work — and most ongoing work in this repo — concentrates in two files:

- [R/NetParams.R](R/NetParams.R) (~1,168 lines) — fits per-respondent Poisson / binomial / linear models for each target statistic, layer by layer (main / casl / inst). Currently **one univariate model per target stat**. This is the main surface of the joint g-comp refactor (issue [#61](https://github.com/EpiModel/ARTnet/issues/61)).
- [R/NetStats.R](R/NetStats.R) (~723 lines) — consumes `netparams` output, builds the synthetic population, and scales per-respondent stats up to population-level ERGM target stats. This is the main surface of [#62](https://github.com/EpiModel/ARTnet/issues/62).

Other files:
- [R/EpiStats.R](R/EpiStats.R) — act rates, condom use, starting HIV prevalence.
- [R/ARTnet-package.R](R/ARTnet-package.R), [R/globals.R](R/globals.R) — package docs and global-binding declarations.

## 3. Public API

Exported functions (see [NAMESPACE](NAMESPACE)):
`build_epistats`, `build_netparams`, `build_netstats`, `make_race_combo`, `reweight_age_pyr`, `trim_epistats`, `trim_netstats`, `update_asmr`.

Canonical pipeline:

```r
epistats  <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta", race = TRUE, time.unit = 7)
netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
netstats  <- build_netstats(epistats, netparams, expect.mort = 0.000478213, network.size = 5000)
```

Pre-built example objects ship in `inst/` (`epistats-example.rds`, `netstats-example.rds`) for users without `ARTnetData` access.

## 4. The current refactor — joint g-computation

### 4.1 Motivation (why we are doing this)

For each ERGM target stat, `R/NetParams.R` currently fits a **separate univariate Poisson GLM** (see [R/NetParams.R:220-450](R/NetParams.R) for the main layer; casl + inst mirror it):

```r
mod_md    <- glm(deg.main ~ 1,                           family = poisson())
mod_age   <- glm(deg.main ~ age.grp + sqrt(age.grp),     family = poisson())
mod_race  <- glm(deg.main ~ as.factor(race.cat.num),     family = poisson())
mod_dcas  <- glm(deg.main ~ deg.casl,                    family = poisson())
mod_hiv   <- glm(deg.main ~ hiv2,                        family = poisson())
```

Each target stat (`md.main`, `nf.age.grp`, `nf.race`, `nf.deg.casl`, `nf.diag.status`) is estimated independently, each marginalizing over ARTnet's conditional joint distribution of the other attributes. This causes three problems:

**Problem 1 — Marginal means don't carry to a different population.** When `build_netstats()` applies `nf.age.grp[k]` (estimated under ARTnet's `{race, deg.casl, hiv2, ...}` distribution given age = k) to a synthetic population with a *different* joint distribution, ARTnet's baked-in conditional distribution creates systematic bias whenever interactions on degree exist.

**Problem 2 — Internal inconsistency across target stats.** In `build_netstats()`:

```r
edges_a = (md.main * n) / 2                        # option 1
edges_b = sum(table(race) * nf.race) / 2           # option 2
```

These disagree whenever the target population's race distribution ≠ ARTnet's. The `edges.avg` argument (default `FALSE`) is a tell — the package already knows the math doesn't close.

**Problem 3 — Frankenstein reference distribution.** `build_netstats()` samples attributes from a mix of sources (NCHS general-population age pyramid; `ARTnetData::race.dist`; ARTnet's own `deg.casl` / `deg.main` dists; uniform `risk.grp`). There is no coherent single target population.

### 4.2 The proposed fix — joint Poisson GLM + g-computation

Fit **one joint Poisson GLM per layer** with relevant attributes + interactions:

```r
m_joint_main <- glm(deg.main ~ age.grp + sqrt(age.grp) +
                      as.factor(race.cat.num) + deg.casl + hiv2 +
                      age.grp:as.factor(race.cat.num),
                    data = d, family = poisson())
```

In a later refactor (#62), `build_netstats()` will build a synthetic population from a single reference distribution, predict per-node expected degree, and aggregate to get target stats that are **internally consistent by construction**:

```r
pred_deg <- predict(m_joint_main, newdata = synth, type = "response")
edges_target      <- sum(pred_deg) / 2
nf.age.grp[k]     <- sum(pred_deg[synth$age.grp == k])
nf.race[r]        <- sum(pred_deg[synth$race.cat.num == r])
# ...
```

### 4.3 Empirical evidence motivating this

From the sibling project `EpiModel/ARTnetPredict`:

- Age-standardized covariate shifts in AMIS MSM 2017 → 2024 are substantial *within* strata: PrEP ever 11% → 48% (+350%); STI any past 12 mo 10% → 15% (+49%); NH Black share 7.7% → 20.8% (+169% age-std); UAS 67% → 84% (+24%).
- Raw AMIS self-reported partners/yr: +38% (OANUM) 2017 → 2024, age-std. Current ARTnet marginal method was the baseline; the PI flagged that the marginal-vs-joint gap may be non-trivial given these compositional shifts.

If attribute interactions on degree are meaningful, the current marginal approach has been introducing quietly biased target stats into every EpiModelHIV simulation. The bias is small when target ≈ ARTnet, but may be substantial for NHBS MSM demographics or 2022–24 AMIS projections.

### 4.4 Phase roadmap

```
Phase 1 (ARTnet): methodology fix
  #61 (joint GLM fit)                  ← ACTIVE
  #62 (build_netstats g-comp refactor)
  #63 (nodematch joint modeling)
  #64 (post-stratification target population API)
  #65 (validation suite: marginal vs joint)

Phase 2+ (ARTnetPredict repo): application
  regenerate 2017-18 baseline, re-run 2022-24 projection,
  post-stratify AMIS to NHBS, gonorrhea simulation, AJE paper
```

Each arrow is a hard dependency. #61 unblocks #62 and #65.

### 4.5 Active task — issue #61

[Refactor `build_netparams()` to fit joint Poisson GLM across attributes](https://github.com/EpiModel/ARTnet/issues/61).

**In scope for #61:**
- Fit one joint Poisson GLM per layer (main, casl, inst) predicting degree from `age.grp`, `race.cat.num`, `deg.{casl|main}`, `hiv2` + interactions.
- Use AIC / LRT to select interaction terms; document choices.
- Store fitted `glm` objects in the return list at `netparams$main$joint_model`, `netparams$casl$joint_model`, `netparams$inst$joint_model`.
- **Keep all existing univariate marginal output fields unchanged**. Joint models are an additive output.

**Out of scope for #61** (separate issues):
- `build_netstats()` refactor — #62.
- Nodematch joint modeling — #63.
- User-supplied target population API — #64.
- Marginal-vs-joint validation suite — #65.

**Design decisions to make in #61:**
1. **Which interactions?** At minimum `age.grp × race.cat.num`. Consider also `age.grp × deg.casl`. Select via AIC, document reasoning.
2. **GLM vs GAM vs ML?** Start with GLM + explicit interactions (interpretable, matches existing ARTnet idiom). Escalate to `mgcv::gam()` only if interaction structure is complex. Avoid ML here — interpretability matters.
3. **Regularization?** Not needed by default. If unregularized GLM fails to converge for sparse cells, consider `glmnet::glmnet(family = "poisson")` as fallback.
4. **Geography.** Include `geogYN` as a main effect (matches existing pattern); no interactions with `geogYN` in the first pass.

**Validation before marking #61 done:**
1. Joint model converges on ARTnet 2017–18 data without warnings, for each layer.
2. Marginal recovery: `mean(predict(joint, type = 'response'))` on the ARTnet sample is within 1% of `mean(observed deg)`.
3. Coefficient sanity: age slope negative (older → lower degree); race coefficients in expected direction.
4. Regression safety: all existing univariate output fields are byte-identical to before the refactor — enforced by the snapshot harness in [`inst/validation/`](inst/validation/) (see §4.7).

### 4.7 Backward-compatibility snapshot harness

[`inst/validation/`](inst/validation/) contains a pre/post regression harness plus a pinned copy of the downstream consumer ([`EpiModelHIV-Template/R/A-networks/`](inst/validation/epimodelhiv_template_ref/)) and the [field-level contract](inst/validation/netstats_contract.md) it reads.

Workflow:

```r
# Step A — run once on pre-refactor `main` (captures golden snapshots)
devtools::load_all()
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
capture_snapshot()

# Step B — run on refactor branch; must report ALL MATCH before merge
devtools::load_all()
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
compare_to_snapshot(method = "existing")   # or whatever the legacy flag is named
```

The harness iterates over `PARAM_SETS` (Atlanta+race, national no-geog, Atlanta no-race at minimum — add coverage if the refactor risks touching more paths) and full-object-diffs both `netparams` and `netstats`. New additive fields like `$joint_model` are stripped before comparison so they don't cause spurious diffs. Snapshot `.rds` files live under [`inst/validation/snapshots/`](inst/validation/snapshots/) and are gitignored (large, local).

**Naming:** when introducing the `method` argument, default should be the legacy behavior so downstream EpiModelHIV users are unaffected until they opt in. A clean progression is `method = c("existing", "joint")` with `"existing"` the default in the refactor release, transitioning the default to `"joint"` in a later minor/major version after validation work is done.

### 4.6 Things to flag to the PI if encountered
- Joint GLM convergence failures (may need regularization or model simplification).
- >10% difference in predicted mean(deg) between marginal and joint on ARTnet-self — would indicate the basic sanity check fails.
- Any accidental change to existing `netparams` output structure (field names, types, values).
- Ambiguity about which interaction terms to include — worth a quick discussion rather than guessing.

## 5. Repo conventions

- **R style:** package-standard roxygen2 docs; 2-space indentation; base R + `dplyr` mix. Match surrounding style.
- **Function naming:** snake_case, verb-first (`build_*`, `trim_*`, `update_*`).
- **Lint:** see [.lintr](.lintr) — 120-char lines, cyclocomp and object-name linters disabled.
- **Roxygen:** markdown mode (`Roxygen: list(markdown = TRUE)` in DESCRIPTION).
- **Dependencies:** prefer `stats::glm` over `glmnet` / `mgcv` unless functionally needed. DESCRIPTION currently imports only `dplyr` and depends on `EpiModel`.
- **Testing:** files under `tests/testthat/`. Note that the existing `test-workflow-*.R` files are end-to-end workflow scripts, not `test_that()` blocks — they require `ARTnetData` + `EpiModelHIV` and execute `netest` / `netdx`. Add new unit tests in proper `test_that()` form alongside them.
- **ARTnetData access:** the data package is private and requires a `GITHUB_PAT`. Functions check for it with `system.file(package = "ARTnetData") == ""`. Document any new direct dependency on specific columns of `ARTnetData::ARTnet.wide` / `ARTnet.long`.

## 6. Running checks

```r
devtools::document()   # regenerate man/ from roxygen
devtools::test()       # run testthat suite
devtools::check()      # full R CMD check
```

Do not push failing CI.

## 7. Git / PR workflow

- Branch off `main` with descriptive name (e.g., `feature/joint-gcomp-netparams`).
- Commits per logical unit (model fitting, tests, docs).
- PR description: summary, design decisions made, validation results (the four checks in §4.5).
- PR should reference the issue (`Closes #61` or `Addresses #61`).
- Request review from PI + package maintainers.

## 8. External references

### Sister project: ARTnetPredict
The detailed methodological analysis that motivated issues #61–#65 lives in `~/git/ARTnetPredict/`:
- `Plan/Step08-BuildLog.md` — §13 age-adjusted partners/yr; §17 raw AMIS vs model; §18 covariate-shift analysis; §20 target-year comparison; §21 3-year pooled recommendation + methodology limitations.
- `Plan/Step09-Plan.md` §2 — most directly relevant to #61 work.
- `CLAUDE.md` — sister-project context.

### Methodological references
- Hernán & Robins, *Causal Inference: What If*, ch. 13 (g-computation / direct standardization).
- Weiss KM, Goodreau SM, Jenness SM et al. (2020), "Egocentric sexual networks of men who have sex with men in the United States: Results from the ARTnet study." *Epidemics* — canonical ARTnet methodology paper.
