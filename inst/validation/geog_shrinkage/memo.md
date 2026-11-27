# Partial pooling of ARTnet city effects: Detroit, Miami, and the four calibrated cities

November 26, 2026. ARTnet 2.10.0 (main at e6b7654), `method = "joint"`, `duration.method = "joint_lm"`, DoxyPEP-Resist settings (seed 20261107, `expect.mort = 0.000478213`, `duration.cap = 941`, network size = `msm.pop`, race composition from each city's `city.csv`). The current targets reproduce DoxyPEP-Resist's fits to the reported precision.

## Question

Under the joint path, each city enters every model as a single unpooled indicator (`geogYN`, city vs the rest of ARTnet). Detroit (69 respondents) and Miami (86) are the smallest cities DoxyPEP-Resist models. We asked whether shrinking those city effects toward division, region, and national values changes the network targets enough to matter (more than about 10% on an edges target or a main or casual duration).

## Methods

For each joint model that sets a target (degree per layer, concurrency, nodematch age and race, absdiff, the log-duration regressions, and the partner-race multinomial), we dropped `geogYN` and refit the same fixed effects (including the AIC-selected interactions) with random intercepts over ARTnet's 24 `city2` groups (15 named cities plus 9 division-level "other" groups), in two forms: flat `(1 | city2)` and nested `(1 | region) + (1 | division) + (1 | city2)`, with glmmTMB (ML; REML for Gaussian models). Count models used a negative binomial likelihood, because the Poisson SEs understate uncertainty badly for the one-off count (Detroit: Poisson SE 0.10, quasi-Poisson SE 0.50). The partner-race multinomial was handled as two baseline-category binary logits. The shrunk city effect is the city's predicted random-effect total minus its mean over the rest of ARTnet, which is the quantity `geogYN` estimates. We substituted it for the `geogYN` coefficient in each model and reran `build_netstats()` with everything else fixed, including the seed. The city's raw cross-layer degree distributions, which enter the synthetic population unpooled, were shrunk toward the national distribution by a Dirichlet-multinomial weight estimated from between-city variation. A brms fit of the nested model (half-normal(0, 0.5) priors on the SDs) checked the three degree models.

## Is there between-city variation to estimate?

Likelihood-ratio tests of a 23-df city factor, given age, race, cross-layer degree, and HIV status:

| Model | Chi-square (23 df) | p |
|---|---:|---:|
| Main degree | 10.4 | 0.99 |
| Main concurrency | 6.4 | 1.00 |
| Main duration | 27.9 | 0.22 |
| Casual degree (NB) | 44.2 | 0.005 |
| Casual duration | 50.4 | 0.0008 |
| One-off count (NB) | 72.6 | < 0.001 |
| Same-race, main / casual / one-off | 38.3 / 74.9 / 57.5 | 0.02 / < 0.001 / < 0.001 |

Once composition is accounted for, ARTnet shows no detectable city variation in main partnerships, so the between-city SD is estimated at zero and every city's main degree, concurrency, and duration effect is pooled to the national value (brms agrees: city SD 0.03, 95% CrI 0.001 to 0.08). Casual degree, casual duration, one-off counts, and race homophily do vary by city, and those effects are partially pooled (city SD about 0.09 to 0.11 for casual degree and 0.20 to 0.24 for one-off counts on the log scale).

## Current city effects for Detroit and Miami

Ratios vs the rest of ARTnet, adjusted, with 95% CIs (quasi-Poisson SE for counts), and the nested shrunk ratio:

| City | Model | Current (95% CI) | Shrunk |
|---|---|---|---:|
| Detroit | Main degree | 0.92 (0.67 to 1.27) | 1.00 |
| Detroit | Casual degree | 0.78 (0.51 to 1.18) | 0.96 |
| Detroit | One-off count | 0.36 (0.14 to 0.96) | 0.74 |
| Detroit | Main duration | 0.88 (0.53 to 1.45) | 1.00 |
| Detroit | Casual duration | 1.34 (0.83 to 2.15) | 1.14 |
| Miami | Main degree | 1.02 (0.77 to 1.37) | 1.00 |
| Miami | Casual degree | 1.30 (0.97 to 1.73) | 1.07 |
| Miami | One-off count | 1.15 (0.73 to 1.82) | 1.05 |
| Miami | Main duration | 0.94 (0.60 to 1.49) | 1.00 |
| Miami | Casual duration | 0.76 (0.54 to 1.06) | 0.88 |

Detroit's low one-off rate comes from the tail, not from more men reporting none: 48% of Detroit respondents report no one-off partners, against 46% in the rest of ARTnet, but the largest Detroit count is 13, against a mean of 4.4 and a maximum of 250 elsewhere. A mean driven by a few high-count men is what a sample of 69 estimates worst.

## Targets, current vs shrunk (nested, with degree-distribution pooling)

Edges per 1,000 nodes, durations in weeks (edge-weighted mean over the age-match strata), percent change in parentheses:

| City | Main edges | Casual edges | One-off edges | Main duration | Casual duration |
|---|---|---|---|---|---|
| Detroit | 168.6 to 180.2 (+7%) | 234.7 to 286.7 (+22%) | 14.3 to 29.2 (+105%) | 200 to 227 (+13%) | 117 to 100 (-14%) |
| Miami | 188.7 to 191.6 (+2%) | 415.3 to 340.2 (-18%) | 56.6 to 51.6 (-9%) | 207 to 220 (+6%) | 62 to 72 (+16%) |
| Atlanta | 179.1 to 177.6 (-1%) | 289.7 to 301.1 (+4%) | 35.5 to 36.0 (+2%) | 248 to 230 (-7%) | 85 to 85 (+1%) |
| Los Angeles | 173.8 to 198.6 (+14%) | 316.9 to 325.5 (+3%) | 42.0 to 45.3 (+8%) | 189 to 217 (+15%) | 65 to 71 (+9%) |
| Seattle | 206.0 to 196.0 (-5%) | 325.4 to 315.3 (-3%) | 44.0 to 52.3 (+19%) | 235 to 229 (-3%) | 64 to 71 (+10%) |
| San Francisco | 224.6 to 193.5 (-14%) | 440.1 to 350.7 (-20%) | 92.5 to 76.8 (-17%) | 198 to 226 (+14%) | 83 to 78 (-5%) |

The flat and nested forms give the same edges targets to within 1 percentage point and durations to within 5. The brms nested fit pools Detroit's one-off effect slightly less (ratio 0.67, 95% CrI 0.42 to 1.01), which would put Detroit's one-off target near 26.5 (+85%) rather than 29.2.

Race homophily: the same-race share of edges changes by at most 1 percentage point in Detroit; in Miami it falls from 71.6% to 65.0% in the main layer and by 1 to 2 points in the others. Per-race nodematch counts move mostly with the edges they belong to (Detroit one-off nodematch roughly doubles in every race; Miami's White/Other main nodematch, a small cell, falls from 519 to 303). Concurrency targets move more in relative terms because they rest on few events (main concurrency in Miami 16.1 to 29.1 per 1,000; casual concurrency in Detroit 109 to 156).

## Reading

1. Both new cities cross the 10% threshold: Detroit on casual and one-off edges and both durations, Miami on casual edges and casual duration. Detroit's one-off rate in particular is not supported by the data at its current value; under any of the three pooled fits it is 1.85 to 2 times higher.
2. The four calibrated cities move by comparable amounts. San Francisco loses 14% to 20% of its edges in every layer and Los Angeles gains 14% on main edges and 15% on main duration. These are the same mechanism, not a Detroit-specific problem: unpooled city effects in main partnerships are noise around a common value (p = 0.99), and the current targets carry that noise into every city.
3. The current approach overstates between-city differences. Main edges per 1,000 range from 169 to 225 now and from 178 to 199 under pooling; one-off edges range from 14 to 92 now and from 29 to 77. Since the shared tetM fitness-cost identification relies on real differences between cities, part of the spread it now sees is sampling noise in ARTnet.

## Corroboration from calibration prep

DoxyPEP-Resist's calibration prep (reported November 26, 2026) found that, at Atlanta's per-act STI probabilities, STIs die out on Detroit's current network, and Detroit's starting probabilities had to be 1.60, 1.63, and 2.05 times Atlanta's for gonorrhea, chlamydia, and syphilis. Miami went the other way: at Los Angeles's probabilities its STIs ran 3 to 7 times its targets, and its starting probabilities are 0.72 to 0.78 times Los Angeles's. San Francisco's calibrated probabilities are already 30% to 40% below the other three calibrated cities'. Each of these goes in the direction pooling would correct: Detroit's one-off target roughly doubles, Miami's casual target falls 18%, and San Francisco's targets fall 14% to 20% in every layer. This is consistent with the calibrated per-act probabilities absorbing sampling noise in the city network targets, which would put that noise into the cross-city differences used to identify the shared fitness cost.

## Recommendation

Hold the Detroit and Miami calibrations, and move all six cities to partially pooled targets in one refit, rather than applying pooling to the two new cities only. Pooling two cities and not the others would put the method difference inside the cross-city contrasts the project uses for identification. Doing this properly needs the pooling built into ARTnet as an opt-in `build_netparams()` option (for example `geog.pooling = "hierarchical"`, default unchanged so legacy output stays byte-identical), not the coefficient substitution used here. If refitting and recalibrating the four calibrated cities is not affordable now, the fallback is to proceed on current targets for all six and report the pooled targets as a sensitivity analysis, with Detroit's one-off rate named as the least reliable input.

## Limitations

The random effects are intercepts only, so city does not modify the age, race, or cross-layer slopes (as with `geogYN` now). Between-city variances rest on 24 groups, and the ML estimates for the main-layer models sit at zero, while the Bayesian fit puts small positive values on them; the conclusions are the same. The "other" groups pool many smaller metros within each division. The pooled targets do not address length-biased sampling or partnership-count truncation (issue #72).

## Implementation in ARTnet 2.11.0

`build_netparams(geog.pooling = "hierarchical")` implements the nested form above, with one change found while building it. For overdispersed count models (Pearson dispersion 1.4 for casual degree and 27 for the one-off count; main degree is underdispersed at 0.7 and is left as fit), the AIC step that chooses the age interactions now uses quasi-AIC, and the pooled refit uses a negative binomial; the point estimates stay Poisson, which keeps predicted counts summing to the observed total over the sample. Quasi-AIC drops both interactions (age by race, age by total degree) from the one-off model that Poisson AIC kept, and keeps the casual model's, which moves the one-off targets of every city (Atlanta +5% with no change in its geography effect). The analysis above used the Poisson-AIC formulas. Main and casual targets and durations match the analysis to within 0.1%. Targets from the package, current vs pooled:

| City | Main edges | Casual edges | One-off edges | Main duration | Casual duration |
|---|---|---|---|---|---|
| Detroit | 168.6 to 180.2 (+7%) | 234.7 to 286.7 (+22%) | 14.3 to 30.2 (+112%) | 200 to 227 (+13%) | 117 to 100 (-14%) |
| Miami | 188.7 to 191.6 (+2%) | 415.3 to 340.2 (-18%) | 56.6 to 46.2 (-18%) | 207 to 220 (+6%) | 62 to 72 (+16%) |
| Atlanta | 179.1 to 177.6 (-1%) | 289.7 to 301.1 (+4%) | 35.5 to 37.2 (+5%) | 248 to 230 (-7%) | 85 to 85 (+1%) |
| Los Angeles | 173.8 to 198.6 (+14%) | 316.9 to 325.5 (+3%) | 42.0 to 45.0 (+7%) | 189 to 217 (+15%) | 65 to 71 (+9%) |
| Seattle | 206.0 to 196.0 (-5%) | 325.4 to 315.3 (-3%) | 44.0 to 48.8 (+11%) | 235 to 229 (-3%) | 64 to 71 (+10%) |
| San Francisco | 224.6 to 193.5 (-14%) | 440.1 to 350.7 (-20%) | 92.5 to 66.4 (-28%) | 198 to 226 (+14%) | 83 to 78 (-5%) |

## Files

Scripts are in `inst/validation/geog_shrinkage/` on branch `feature/geog-pooling`: `01_build.R` (current targets), `03_models.R` (random-effect fits), `04_targets.R` (rebuilt targets), `05_het.R` (heterogeneity tests), `06_brms.R` (Bayesian check), `07_tables.R` (tables). They read and write `.rds` files in the working directory and read the city inputs from DoxyPEP-Resist.
