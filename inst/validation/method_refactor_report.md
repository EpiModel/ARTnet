# Joint G-Computation for ARTnet-Based MSM Network Parameterization: A Methodological Refactor

_ARTnet 2.9.0_

## Introduction

The ARTnet study (Weiss et al., *Epidemics* 2020) is an anonymous cross-sectional web-based survey of men who have sex with men (MSM) in the United States, conducted in 2017–2018, that has become a canonical input for parameterizing temporal exponential random graph models (TERGMs) of sexual partnership networks in the EpiModelHIV ecosystem. The `ARTnet` R package transforms the survey microdata into three objects consumed downstream by `EpiModelHIV-p`: epidemic parameters (`build_epistats`), per-respondent network statistics (`build_netparams`), and population-scaled ERGM target statistics (`build_netstats`). The current pipeline has carried the same analytic conventions for several years, with refinements layered into geographic stratification, age-range flexibility, and race-stratification options, but the underlying estimation strategy for each ERGM target statistic — a separate univariate Poisson, binomial, or linear regression of that statistic on a single attribute at a time — remained unchanged.

The motivation for the present work arose from an ongoing project — ARTnetPredict — that is building machine-learning-based forward projections of ARTnet target statistics to incorporate the 2022–2024 American Men's Internet Survey (AMIS) waves. Specifying that forward projection surfaced a methodological concern applicable to the ARTnet baseline itself. Under the legacy "univariate marginal" approach, each per-attribute estimate carries ARTnet's *conditional joint distribution of the other attributes* baked in. When `build_netstats()` then applies these per-attribute estimates to a synthetic target population whose attribute joint distribution differs from ARTnet's — any city-specific MSM target, or any forward-projected population — three problems emerge. First, attribute interactions on degree (e.g., age × race effects on partnership formation) are structurally invisible to univariate fits, so the marginal estimate is biased in any direction in which such interactions exist. Second, the resulting target statistics are mutually inconsistent: edges computed from `md.main * num / 2` need not equal edges computed as `Σ_r table(race) * nf.race[r] / 2`, and the long-standing `edges.avg` argument in `build_netstats()` was a tacit acknowledgement that the math did not close. Third, the target population sampled in `build_netstats()` was a patchwork of reference sources (NCHS general-population age pyramid, `ARTnetData::race.dist`, ARTnet's own degree and role distributions), with no coherent specification of who exactly the population was.

The within-ARTnet baseline analytics needed to be defensible on their own terms before any ML-based forward projection could rest on them. The methodology presented here replaces the legacy per-attribute univariate fits with joint Poisson and binomial generalized linear models that condition on all attributes simultaneously, applies g-computation to aggregate per-respondent or per-dyad predictions against an explicitly-specified synthetic target population, and provides a unified post-stratification interface for that target. Conceptually, the move is from estimating marginal per-attribute effects (which transport poorly across populations) to estimating attribute-conditional effects and predicting at the target's joint distribution — direct standardization in the sense of Hernán and Robins (2020). Strict backward compatibility with the legacy analytic pipeline is preserved as the default behavior, verified at every step by a snapshot-based regression harness.

## Methods

The new methodology is exposed through three arguments on the public ARTnet API: `method = c("existing", "joint")` on both `build_netparams()` and `build_netstats()`; `duration.method = c("empirical", "joint_lm")` on `build_netparams()`; and `target_pop` (`NULL`, list, data.frame, or character) on `build_netstats()`. All defaults preserve byte-identical legacy behavior.

**Joint estimation in `build_netparams`.** Under `method = "joint"`, seven model objects are fit per network layer (main, casual, one-time partnerships):

- A **joint Poisson GLM** of partnership count on `age.grp + sqrt(age.grp) + race + cross-layer-degree + hiv2`, with AIC-based selection over `age.grp:race` and `age.grp:cross-layer-degree` interactions. This replaces the five univariate Poisson fits the legacy method ran for `md`, `nf.age.grp`, `nf.race`, `nf.deg.<x>`, and `nf.diag.status`.
- A **joint binomial GLM** of the concurrency indicator (`deg > 1`) on the same RHS, for the main and casual layers. The binomial form was chosen over a Poisson-derived `P(deg > 1) = 1 - exp(-λ)(1+λ)` because `deg.main` and `deg.casl` are truncated in the training data at 2 and 3 respectively; the Poisson-implied tail probability is biased upward by orders of magnitude.
- **Joint binomial GLMs** for `nodematch` on `same.age.grp` and `same.race`, fit on long-form partnership data with ego attributes only on the right-hand side. Partnership pairs are treated as marginally generated in this first pass; the ego side comes from the corrected target population.
- **Joint Gaussian regressions** for `absdiff_age` and `absdiff_sqrt.age`, with the same ego-attribute RHS structure.
- A **joint log-linear regression** of `log(duration.time)` among ongoing partnerships under `duration.method = "joint_lm"`, with ego attributes plus `same.age.grp` and `same.race` matching terms.

**G-computation aggregation in `build_netstats`.** Under `method = "joint"`, target statistics are computed by aggregating per-synthetic-node predictions against the synthetic population. Edges become `sum(pred_deg) / 2`, and `nodefactor_<attr>[level] = sum(pred_deg | attr == level)`. By construction, `Σ_<level> nodefactor_<attr>[level] = 2 × edges` to machine precision — the long-standing `edges.avg` inconsistency is resolved structurally. For dyad-level targets, per-ego predicted partnership properties are weighted by per-ego predicted degree: `nodematch_<attr>[level] = sum(pred_deg * pred_<dyad>)[ego in level] / 2`. For dissolution offsets, partner-race uncertainty is marginalized via the `joint_nm_race_model` predictions before per-stratum medians are computed and run through the existing geometric transformation `mean.dur.adj = 1 / (1 - 2^(-1 / median))`.

An important methodological consideration concerns the choice of estimand for the duration target. Under TERGM's geometric (constant-hazard) dissolution offset, the simulated network's mean age of extant ties at cross-section coincides with mean full partnership duration when partnerships are exponentially distributed, and diverges otherwise. The simulation can therefore match the former exactly but not the latter. A length-bias-corrected Weibull MLE explored during development confirmed that ARTnet's per-stratum partnership-duration shape parameter is far from 1 (`k ≈ 0.6` for casual partnerships, `k > 2` for older-matched main partnerships), so the direct mean-duration estimate from a Weibull AFT fit could not be honored by the geometric simulation in any case. Both the empirical and joint_lm methods therefore target the cross-sectional age statistic that the geometric simulation *can* match — mean age of extant ties — rather than mean full partnership duration.

**Post-stratification API.** The `target_pop` argument unifies what was previously several arguments. As a named list, it overrides any subset of `{age.pyramid, race.prop, deg.casl, deg.main, deg.tot, role.class, risk.grp}` against the legacy default sources. As a data.frame, it supplies the synthetic population in full — required columns `age`, `deg.casl`, `deg.main`, `role.class`, `risk.grp` (plus `race` when race-stratified) — bypassing attribute sampling entirely. A character form is reserved for a future bundle of geography-specific defaults (e.g., `"atlanta"`, `"us_msm_male"`) constructed from data already in `ARTnetData::race.dist` plus the NCHS age pyramid.

**Validation and reproducibility infrastructure.** A snapshot-based regression harness captures golden-reference output of the legacy pipeline on three reference parameter sets and compares against any subsequent code state at machine precision; legacy behavior is preserved exactly across all releases of the new methodology. A separate method-comparison harness runs both methods across the four city scenarios reported below and produces the side-by-side comparison tables. Continuous integration runs the standard R package check on every change, with `ARTnetData` installed so that the joint-fit functionality is exercised in CI rather than skipped. The unit-test suite contains 571 assertions covering the joint formation, dyad-level, and duration models, plus the post-stratification API and parameterization edge cases.

## Results

The legacy approach (`method = "existing"`) was compared against the joint approach (`method = "joint"` paired with `duration.method = "joint_lm"`) across four city scenarios drawn from settings in which ARTnet has been used for EpiModelHIV-based simulation: Atlanta, Boston, Chicago, and Seattle. Each scenario stratifies by race and draws the synthetic population's race composition from the city-specific entry in `ARTnetData::race.dist`. The cities span the range of demographic distance from ARTnet's overall sample composition (80.7% White/Other, 13.8% Hispanic, 5.5% Black): Seattle is closest, Atlanta is farthest (Table 1). All four scenarios use a common mortality rate (`expect.mort = 0.00025`, slightly lower than the EpiModelHIV-Template baseline of 0.000478213) to accommodate Seattle's small-sample noise on the matched-and-old main-partnership duration stratum; this choice shifts absolute dissolution coefficient values uniformly across cities and preserves the within-comparison consistency.

**Table 1.** Comparison scenarios. Race composition from `ARTnetData::race.dist`.

| Scenario | White / Other | Black | Hispanic | Distance from ARTnet sample |
|---|---:|---:|---:|---|
| `seattle_default` | 87.4% | 6.1% | 6.5% | very close |
| `boston_default` | 59.7% | 20.9% | 19.4% | medium |
| `chicago_default` | 41.8% | 29.2% | 29.0% | medium |
| `atlanta_default` | 43.9% | 51.5% | 4.6% | far |
| _ARTnet sample (reference)_ | 80.7% | 5.5% | 13.8% | — |

Across **384 target-statistic cells** (every numeric entry in the per-layer netstats namespace, summed across the four scenarios), **232 cells (60%) shift by more than 5%** under joint-vs-existing (Table 2; Figure 1). The materially-shifted fraction tracks the demographic-distance ordering: Seattle (closest to ARTnet sample) shows the fewest material shifts at 46%; Atlanta, Boston, and Chicago all sit at 64–67%. This is consistent with the marginal-vs-joint hypothesis: the larger the gap between target-population composition and ARTnet sample composition, the more the correction matters.

**Table 2.** Material-shift summary per scenario.

| Scenario | Cells | Shifted (\|%Δ\| > 5%) | % shifted |
|---|---:|---:|---:|
| `seattle_default` | 96 | 44 | 46% |
| `chicago_default` | 96 | 61 | 64% |
| `atlanta_default` | 96 | 63 | 66% |
| `boston_default` | 96 | 64 | 67% |
| **Total** | **384** | **232** | **60%** |

![Figure 1: Joint vs existing target statistics across all 384 cells. Each point is one (scenario, layer, stat, level) cell, plotted on log–log axes; the dashed line is y = x (no correction). Cells off the diagonal are bias corrections, with distance proportional to the magnitude of the shift. Coloured by network layer; faceted by city. Seattle (top right; demographically closest to ARTnet sample) sits visibly tighter to the diagonal than the others.](figures/fig1_joint_vs_existing.png)

The largest shifts cluster on three target-stat families (Table 3): dissolution durations in matched-and-old strata, one-time-partnership nodematch in older age groups, and high-degree casual nodefactor. The first two move downward (the joint method shrinks small-cell estimates toward the population), the third upward (joint adjustment finds a weaker negative cross-degree substitution than the marginal fit).

**Table 3.** Top 12 shifts by |%Δ|, atlanta_default scenario. Positive %Δ means joint > existing.

| Layer | Stat | Level | Existing | Joint | %Δ |
|---|---|---:|---:|---:|---:|
| inst | nodematch_age.grp | 5 | 7.90 | 3.86 | −51.1% |
| main | dissolution_duration | 6 | 934.45 | 491.48 | −47.4% |
| inst | nodefactor_diag.status | 2 | 150.36 | 89.69 | −40.3% |
| main | dissolution_duration | 4 | 539.20 | 323.74 | −40.0% |
| casl | nodefactor_deg.main | 3 | 11.74 | 16.38 | +39.5% |
| main | dissolution_duration | 5 | 682.59 | 428.11 | −37.3% |
| casl | dissolution_duration | 2 | 50.44 | 68.93 | +36.7% |
| inst | nodefactor_age.grp | 1 | 46.94 | 63.36 | +35.0% |
| inst | nodefactor_age.grp | 5 | 85.81 | 57.25 | −33.3% |
| inst | nodematch_age.grp | 4 | 11.22 | 7.99 | −28.8% |
| casl | dissolution_duration | 4 | 113.16 | 81.76 | −27.7% |
| casl | dissolution_duration | 3 | 73.17 | 92.29 | +26.1% |

**Dissolution durations on matched-and-old strata** are where the methods diverge most prominently (Figure 2). The dissolution-duration target for main partnerships matched within the oldest age group moves from 934.4 to 491.5 weeks (−47%) on default Atlanta, and the same direction holds across all four city scenarios — with even larger drops for the smaller-sample cities (Seattle: −66%, Chicago: −52%, Boston: −50%). The empirical method takes a stratum-level median over a small per-cell sample whose tail is dominated by a few long ongoing partnerships; the joint log-linear regression smooths across the full sample via the multivariate fit. The legacy `smooth.main.dur` heuristic (which averages the oldest two age-matched strata) reduces the unsmoothed empirical estimate from 1186 to 934 weeks, but the joint method's estimate for the same stratum is already 476 before smoothing — the residual ≈2× gap between the smoothed empirical and joint methods is the difference between two-neighbor averaging and full multivariate borrowing of strength. **One-time-partnership nodematch in older age groups** shifts in proportion: the count of within-age-group matches for one-time partnerships in the oldest age group moves from 7.9 to 3.9 (−51%) on default Atlanta. The dyad-level joint binomial picked up an age × race interaction that was structurally invisible to the marginal fit, and the older one-time-partnership age strata are sparse enough that the marginal estimate was dominated by sampling noise that the multivariate fit shrinks toward the population mean. **High-degree casual nodefactor** shifts in the opposite direction: the casual-layer endpoint count among respondents reporting two main partners rises from 11.7 to 16.4 (+40%). The joint Poisson on casual degree finds a weaker negative cross-degree (main → casual) substitution effect once age and HIV status are controlled for, which raises the predicted casual-layer edge contribution from respondents at the high end of the main-degree distribution.

![Figure 2: Stratum-level dissolution duration (mean.dur.adj, weeks) across the six main-layer strata. Empirical without smoothing (open circles) shows the matched.5 small-sample noise; the legacy smooth.main.dur option averages matched.5 with matched.4 (filled circles). joint_lm (red triangles) borrows strength from the full sample via the multivariate fit and produces a smoother stratum-level profile.](figures/fig2_dissolution_durations.png)

A decomposition of the **−15% main-edges shift** on Atlanta default illustrates the correction concretely (Table 4; Figure 3). ARTnet's sample is 80.7% White/Other, 13.8% Hispanic, and 5.5% Black; Atlanta's MSM population (per `ARTnetData::race.dist`) is 51.5% Black, 4.6% Hispanic, and 43.9% White/Other. The legacy pipeline's edges target was `md.main * num / 2 = 0.398 * 5000 / 2 = 995`, where `md.main = 0.398` is ARTnet's overall marginal mean main degree. The joint method aggregates per-race predictions across Atlanta's actual race composition: `0.515 × 0.272 + 0.046 × 0.408 + 0.439 × 0.408 = 0.338`, giving `0.338 * 5000 / 2 = 845` edges — exactly the −15% shift. The legacy pipeline's edges target implicitly assumed Atlanta MSM had ARTnet's sample race mix.

**Table 4.** Race composition mismatch and per-race main degree, decomposing the −15% main-edges shift.

| Race | ARTnet sample | Atlanta MSM target | Joint per-race E[deg.main] |
|---|---:|---:|---:|
| Black (race = 1) | 5.5% | 51.5% | 0.272 |
| Hispanic (race = 2) | 13.8% | 4.6% | 0.408 |
| White / Other (race = 3) | 80.7% | 43.9% | 0.408 |
| **Marginal mean** (legacy) | **0.398** | — | — |
| **Joint mean at Atlanta target** (this work) | — | **0.338** | — |

![Figure 3: Race composition mismatch between the ARTnet sample and the Atlanta MSM target population. ARTnet over-represents White / Other respondents; Atlanta's actual MSM population is roughly half Black. The legacy pipeline's main-edges target multiplied ARTnet's overall mean main degree by network size, implicitly assuming the target population has ARTnet's race mix.](figures/fig3_race_composition.png)

Coefficient comparisons on the joint Poisson for main degree show several attribute effects strengthen substantially after adjustment (Table 5). The `deg.casl` slope moves from −0.24 (marginal) to −0.55 (joint) — the apparent within-respondent substitution between casual and main partnerships is more than twice as strong once age is controlled for, consistent with the marginal estimate being attenuated by the positive correlation between `deg.casl` and youth (and youth's positive direct effect on main degree). The `hiv2` slope moves from +0.09 to +0.25 in the same direction. The age profile (`age.grp + sqrt(age.grp)`) sharpens by approximately 20% — the steeper slope post-adjustment reflects that older respondents in ARTnet differ from younger respondents on more than age alone. AIC selects an `age.grp:deg.casl` interaction (+0.10) that is structurally invisible to the univariate approach.

**Table 5.** Coefficient comparison: marginal univariate vs joint multivariate Poisson, main layer (log-link). Marginal coefficients are from separate `glm(deg.main ~ <single attribute>, family = poisson)` fits, mirroring the legacy method.

| Term | Marginal | Joint | Δ |
|---|---:|---:|---:|
| `age.grp` | −1.001 | −1.185 | −0.185 |
| `sqrt(age.grp)` | +2.968 | +3.475 | +0.507 |
| `as.factor(race.cat.num)2` (Hispanic) | +0.411 | +0.450 | +0.039 |
| `as.factor(race.cat.num)3` (W/Other) | +0.403 | +0.467 | +0.064 |
| `deg.casl` | −0.237 | −0.553 | −0.316 |
| `hiv2` | +0.094 | +0.246 | +0.152 |
| `age.grp:deg.casl` (AIC-selected) | — | +0.103 | — |

End-to-end ERGM estimation under both methods (Atlanta defaults, network.size = 5000, Stochastic-Approximation) converges cleanly across all six target ERGMs (3 layers × 2 methods). Static `netdx` diagnostics on the default-method main model match all formation target stats within `|Z| ≤ 2.05` and `|% diff| ≤ 4.2%` over 1000 simulations. The joint method's `coef.form[edges]` is consistent with its lower edge target, and no degeneracy or sampling pathology appears in either path.

The complete per-scenario, per-target-statistic comparison is provided as a supplementary table. Backward-compatibility regression checks confirm that the legacy pipeline output is reproduced exactly under default arguments and under explicit selection of the legacy method (across all three reference parameter sets), at machine precision; the unit-test suite (571 assertions) and the standard R package check both pass cleanly.

## Discussion

The 60% material-shift rate is large but methodologically expected: the legacy univariate approach computes target stats one attribute at a time from a sample whose joint attribute distribution does not match the synthetic populations the package targets. The size and pattern of the shifts confirm that this confounding is not a small effect, particularly for derived quantities (nodefactor counts in sparse strata, dissolution durations in small per-stratum cells) where the joint method's effective sample size advantage is largest. The cross-city pattern is particularly informative: Seattle's racial composition is closest to ARTnet's sample, and Seattle has the smallest material-shift rate (46%); the three more-different cities (Atlanta, Boston, Chicago) all sit at 64–67%. Shift directions are stable across scenarios — the joint method is not introducing noise but applying a systematic correction whose magnitude scales with the demographic distance between target and sample.

A practical question worth addressing: the legacy pipeline used `geogYN` as a main-effect covariate in each per-attribute regression, rather than filtering ARTnet to the geographic stratum of interest. That was a defensible design choice — filtering to Atlanta would have reduced the analytic sample to just 206 respondents (out of the full n = 4863), with per-race subsets dropping to a few dozen, and the resulting estimates would have been dominated by sampling noise. The `geogYN` covariate let the regression borrow strength across the whole sample while still specializing to Atlanta at prediction time. And it did partial demographic rebalancing: ARTnet's Atlanta sub-sample is 12.6% Black versus the full sample's 5.5%, so `predict(..., geogYN = 1)` shifts toward Atlanta-relevant rates — just not enough. The unaddressed problem is that the Atlanta sub-sample remains 78.6% White/Other (the online-recruitment selection bias on race operates within geographic strata, not just at the national level) while Atlanta's actual MSM population is 43.9% White/Other. The per-respondent mean (`md.main = 0.398`) thus reflects the Atlanta sub-sample's race composition, while the synthetic population gets drawn separately with the *target* Atlanta race composition from `ARTnetData::race.dist`. The two halves of the pipeline implicitly assumed different race compositions and never reconciled. The joint g-computation closes that loop by aggregating per-respondent predictions across the synth's actual joint distribution rather than the sample's. The legacy `geogYN` adjustment was the right move statistically but couldn't fix a problem that lived in the synth-construction step, not the regression step.

For ongoing EpiModelHIV-p simulations, the immediate implication is that any model parameterized via the legacy ARTnet pipeline carries baseline target statistics that approximately reflect ARTnet's race and (to a lesser extent) age composition, not the modeled target population's composition. For Atlanta-specific models — the dominant use case — the legacy method over-targets main edges by approximately 15%, under-targets casual edges by approximately 8%, and shifts the `nodematch` and dissolution stats heterogeneously across age strata. Whether this propagates to materially different epidemic dynamics depends on the simulation question, but the corrections move in directions one would expect to affect cumulative HIV/STI incidence over a multi-year horizon. Re-parameterizing existing EpiModelHIV-p workflows with `method = "joint"` is now possible without code changes downstream — the netstats object retains its full legacy field set; only the *values* shift.

Several methodological limitations remain. First, TERGM dissolution is constant-hazard. Per-stratum Weibull shape parameters fit during development indicated that the geometric assumption is reasonable for casual partnerships (`k ≈ 0.6` across strata, decreasing hazard) but mis-specified for older-matched main partnerships (`k > 2`, increasing hazard). A duration-dependent dissolution offset that would allow the simulation to honor these shape differences is outside the present scope and would require extension of the underlying network simulation framework. Second, ARTnet's cross-sectional sampling design exposes formation-statistic estimates to length-biased sampling on partnership-pair targets (`nodematch`, `absdiff`) and to partnership-count truncation at two main / three casual partners on degree-based targets. The first concern is partially addressable by restricting partnership-pair fits to ongoing partnerships only; the second requires a truncated-Poisson likelihood and is more substantial. Neither correction is implemented in the present work; both are open extensions. Third, the joint duration model uses ongoing partnerships only, consistent with the legacy convention but inheriting its length-bias in observed elapsed durations. As discussed above, the estimand here is the cross-sectional mean age of extant ties — a quantity the geometric simulation can target exactly — rather than mean full partnership duration.

The corrected within-ARTnet baseline opens three concrete directions for the ARTnetPredict project. First, any ML-based forward projection should now rest on the joint-corrected 2017–2018 baseline rather than the marginal-method baseline. Second, post-stratification to a 2022–2024 AMIS-derived demographic posterior can be specified directly through the new `target_pop` interface, without dependency on package-shipped reference data. Third, post-stratification to any city's overall male age × race distribution — the actual production use case — becomes a one-line argument rather than a manual reconstruction of the synthetic population.

A natural future direction is a methodological paper formalizing the joint-versus-marginal comparison with simulation-based validation: fit synthetic data with known joint-attribute structure, demonstrate recovery under both methods, and quantify the bias as a function of target-versus-sample population divergence. The survival-versus-network-dynamics framing of the duration estimand — particularly the inspection-paradox argument that mean age of extant ties is the appropriate inferential target given geometric simulation — would sharpen the formal exposition. Re-analysis of one or two published HIV-simulation findings with the corrected baseline would give concrete substantive stakes to the methodological argument, and would naturally connect the formation-statistic sampling-bias corrections (length-bias and partnership-count truncation) discussed above to the broader literature on egocentric network inference.

The joint-corrected analytic pipeline is now available as opt-in arguments to the package, with legacy behavior preserved by default and verified by a snapshot-based regression harness. Downstream simulations adopting the corrected target statistics require no change to consumer code beyond setting the new arguments at the parameterization step.

## References

- Hernán MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC, 2020. (Chapter 13: g-computation / direct standardization.)
- Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, Jenness SM. Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study. *Epidemics* 2020; 30: 100386.
- Krivitsky PN, Morris M. Inference for social network models from egocentrically sampled data, with application to understanding persistent racial disparities in HIV prevalence in the US. *Annals of Applied Statistics* 2017; 11(1): 427–455.

## Reproducibility

```r
# Backward-compat regression
source(system.file("validation/validate_backward_compat.R", package = "ARTnet"))
compare_to_snapshot(method = "existing")    # ALL MATCH 3/3

# Method comparison across scenarios
source(system.file("validation/method_comparison.R", package = "ARTnet"))
res <- compare_methods()
summarize_comparison(res)
render_comparison_report(res)               # writes inst/validation/method_comparison.md
```

The snapshot files used by the regression harness are gitignored at ~12 MB each; capture them once via `capture_snapshot()` on a clean main checkout.
