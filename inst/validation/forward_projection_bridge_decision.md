# Forward-projection bridge: the estimand decision (ARTnetPredict to ARTnet)

Status: decision memo for the PI. Written 2026-06-08 during the deep-dive review of the
joint g-computation refactor and its use as the bridge from ARTnetPredict's forward-projected
populations to updated `netstats`.

This memo exists because one question on the ARTnetPredict path is a modeling decision, not a
bug, and it determines the architecture of everything downstream. The package-side bugs found
during the review have been fixed (see the last section); this is the part that needs a call.

## The question in one sentence

When we push ARTnetPredict's forward-projected 2022-24 population into ARTnet to get updated
network targets, should the temporal change enter through the population's **attribute
distribution** (Option A) or through ARTnetPredict's **projected degrees** (Option B)? They
give different answers and only one of them can ride through `build_netstats(target_pop=)`.

## What the refactor gives us, and the catch

`build_netstats(method = "joint", target_pop = df)` takes a one-row-per-node population and
runs g-computation (direct standardization): it predicts each node's expected degree from the
joint Poisson GLMs **fit on ARTnet 2017-18**, then aggregates `edges = sum(pred_deg)/2`,
`nodefactor[level] = sum(pred_deg | level)`, internally consistent by construction.

The catch, verified live against current HEAD: the `deg.main` / `deg.casl` columns you supply
are used **only as covariates**; the edge counts are re-predicted, not taken from the supplied
degrees.

```
supplied mean deg.casl = 1.50  ->  implied casl edges 3753  ;  g-comp casl edges 1182
supplied mean deg.main = 1.01  ->  implied main edges 2528  ;  g-comp main edges  680
```

There is no fixed point. Forcing everyone to a given casual degree changes the predicted main
degree but never reproduces the input:

```
supplied deg.casl = 0 (all)  ->  predicted mean main degree 0.459
supplied deg.casl = 3 (all)  ->  predicted mean main degree 0.222
```

This is correct behavior for g-computation. It is also the crux of the problem.

## Why it matters for ARTnetPredict

ARTnetPredict's headline result is a **conditional behavioral** change: after standardizing to
the 2017-18 age distribution, casual degree still grows +20 to +70% (Step 8b), and total
partners/year is +31% pooled 2022-24. That change lives in the conditional behavior, not in the
demographic composition.

`target_pop` holds conditional behavior fixed at ARTnet 2017-18 and only swaps the attribute
distribution. So feeding ARTnetPredict's projected population through `target_pop` would
transport the demographic shift (age, race, HIV mix) but reproduce roughly 2017-18 partner-count
behavior on top of it. It cannot, on its own, carry the +31% finding. The HIV mix does move
edges (`diag.status` is a degree predictor: all-negative gives main degree 0.346, all-positive
0.442, about a 28% swing), so composition is not nothing, but it is not the behavioral signal.

## The three options

**Option A: attribute standardization (what target_pop does today).**
Supply the projected attribute distribution; let ARTnet's frozen 2017-18 GLMs predict degrees.
- Answers: "what would 2017-18 ARTnet network structure look like on a 2022-24 demographic
  population?"
- Loses: the projected behavioral change. The +31% does not appear.
- Cost: near zero. Already built, validated, snapshot-backed.

**Option B: transport ARTnetPredict's projected degrees.**
Make the projected degrees drive the targets.
- Answers: "what are the 2022-24 network targets given projected behavior?"
- Requirement: `target_pop` is the wrong hook (it discards supplied degrees). Needs either a
  custom netstats assembler that uses ARTnetPredict's degrees directly for the degree-based
  targets, or re-fitting ARTnet's joint GLMs on projected node-level (attribute, degree) data.
- Cost: real engineering; this is essentially the original Step 9 Part B plan.

**Option C: refit on projected data, reuse ARTnet's machinery (recommended).**
Re-fit the joint Poisson/binomial GLMs on ARTnetPredict's projected per-node (attribute,
degree) pairs, then reuse ARTnet's g-computation aggregation for internal consistency and
borrow ARTnet's joint dyad / mixing / duration models only where AMIS lacks the data
(partner attributes, role.class). This gets both the projected behavior and the
internal-consistency guarantee, without a second from-scratch g-comp implementation.
- This reframes the earlier "Step 9 Part B is ~80-90% subsumed by the refactor" read: the
  refactor subsumes the **aggregation and post-stratification machinery**, but **not** the
  fit-on-projected-data step, which is exactly what carries the behavioral change.

## The go/no-go experiment (do this before committing to an architecture)

Build the realized 2022-24 node data.frame, run it through Option A, and compare the resulting
`edges` / `nodefactor` to ARTnetPredict's own projected `md.*` / `nf.*`. If Option A does not
reproduce the +31% casual shift (expected, since attributes alone should not carry it), that is
the definitive evidence that Option A is insufficient and the project needs Option B/C. Cheap,
and it settles the fork with data rather than argument.

## Secondary decisions that ride along

- **HIV prevalence for the projected population.** `diag.status` is a degree predictor (about
  28% swing), so the chosen prevalence scales edges, not just `nodefactor_diag.status`. AMIS/SL
  do not project it. Decide: supply `df$diag.status` from a chosen 2022-24 prevalence, or fall
  back to `epistats$init.hiv.prev`.
- **role.class.** AMIS has no sex-role variable; carry it from the ARTnet canonical
  distribution. Impact is confined to node composition (`nodematch_role.class` is structurally
  c(0,0)), so carrying it forward is defensible.
- **Durations.** ARTnetPredict projects scalar mean durations (D_main -36%, D_casl +38%);
  ARTnet aggregates `joint_lm` log-duration predictions into the `diss.byage` offset object the
  template consumes, targeting a different estimand (cross-sectional extant-tie age, not full
  duration). Decide whether to trust ARTnet's `diss.byage` over the projected population
  (internally consistent) or splice the projected scalars (reintroduces the inconsistency the
  refactor fixed).
- **Population size N.** Independent of the numerics (counts scale linearly with N; rates and
  mixing are N-invariant), but the data.frame form forces `network.size = nrow(df)`, so N is
  chosen at df-assembly time. Can proceed in parallel with the popsize2026 work.

## Dependency notes for whoever wires this

- ARTnetPredict's `renv.lock` pins ARTnet to the pre-refactor SHA `a19ffe98`; the joint API
  does not exist in the pinned dependency. Bump the pin to a post-refactor commit and
  `renv::restore` before any of this is callable in the project environment.
- ARTnetPredict currently emits aggregate targets plus a per-node table of soft predictions
  (`resp_preds` / `slot_preds`), not realized integer node attributes. A realization adapter
  (sample/round the slot decomposition into integer `deg.main`/`deg.casl`/`deg.tot`, carry
  `role.class`/`diag.status`, integer-code race to ARTnet's level order, clamp to fitted
  support, drop/recode age>65) is the one genuinely non-subsumed engineering piece.

## Package fixes already made (so this memo stands alone)

These were unambiguous bugs found during the review and are fixed in `R/NetStats.R`, with
backward compatibility verified byte-identical (snapshot 3/3) and regression tests added in
`tests/testthat/test-target-pop-robustness.R`:

1. `method="joint"` + `duration.method="joint_lm"` + a `target_pop` with an empty age stratum
   (e.g. a population restricted below 55) crashed in `dissolution_coefs`. Now falls back to the
   within-ARTnet byage durations per empty stratum.
2. Degree / role / risk values outside the fitted support are clamped with a warning so the
   `sum(nodefactor) = 2 * edges` invariant holds.
3. Nodefactor tabulations now span the full level set, so an unobserved level
   (e.g. an absent risk quintile) no longer silently recycles a per-level vector.
4. Optional `seed` argument to `build_netstats` for reproducibility of the sampling path.
