# Smoke-test the build_* pipeline across parameterization modes that
# aren't already covered by test-joint-{model,netstats,dyad}.R (which
# all use Atlanta + race = TRUE). Each test_that() block runs the full
# pipeline under a different parameterization and asserts that the
# downstream-facing fields of netstats have the right shape and no NAs
# in places where they'd break EpiModelHIV-Template's ERGM formulas.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

# Check the public contract that model_{main,casl,ooff}.R in
# EpiModelHIV-Template actually reads from netstats. Keep this in sync
# with inst/validation/netstats_contract.md.
expect_netstats_contract <- function(ns, race = TRUE) {
  # Attribute vectors used by initialize.R
  for (a in c("age", "sqrt.age", "age.grp", "active.sex",
              "race", "deg.casl", "deg.main", "deg.tot",
              "risk.grp", "role.class", "diag.status")) {
    expect_true(!is.null(ns$attr[[a]]),
                info = paste("missing attr:", a))
  }
  # Per-layer fields used by model_*.R
  for (layer in c("main", "casl", "inst")) {
    expect_true(is.numeric(ns[[layer]]$edges) && length(ns[[layer]]$edges) == 1)
    expect_true(ns[[layer]]$edges > 0)
    expect_true(is.numeric(ns[[layer]]$nodefactor_age.grp))
    expect_false(any(is.na(ns[[layer]]$nodefactor_age.grp)))
    if (isTRUE(race)) {
      expect_true(is.numeric(ns[[layer]]$nodefactor_race))
      expect_false(any(is.na(ns[[layer]]$nodefactor_race)))
      expect_true(is.numeric(ns[[layer]]$nodematch_race_diffF))
      expect_true(length(ns[[layer]]$nodematch_race_diffF) == 1)
    }
  }
  # concurrent exists for main/casl (not inst)
  for (layer in c("main", "casl")) {
    expect_true(is.numeric(ns[[layer]]$concurrent))
    expect_true(ns[[layer]]$concurrent >= 0)
  }
  # diss objects for main/casl
  expect_s3_class(ns$main$diss.byage, "disscoef")
  expect_s3_class(ns$casl$diss.byage, "disscoef")
}


# ---------------------------------------------------------------------------
# National (no geographic stratification)
# ---------------------------------------------------------------------------
test_that("no-geog (national) parameterization: build_* pipeline works", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(race = TRUE, time.unit = 7,
                             init.hiv.prev = c(0.33, 0.137, 0.084))
  expect_null(epistats$geog.lvl)

  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "existing")
  expect_netstats_contract(
    build_netstats(epistats, np, expect.mort = 0.000478213,
                   network.size = 3000, method = "existing")
  )

  # Joint path must also run under no-geog (no geogYN term in formulas).
  set.seed(20260419L)
  np_j <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  ns_j <- build_netstats(epistats, np_j, expect.mort = 0.000478213,
                         network.size = 3000, method = "joint")
  expect_netstats_contract(ns_j)
  # Internal consistency under joint carries through without geog too.
  expect_equal(sum(ns_j$main$nodefactor_race), 2 * ns_j$main$edges,
               tolerance = 1e-9)
})


# ---------------------------------------------------------------------------
# Sexual cessation (CEA-style: ages 15-100, sexual cessation at 65)
# ---------------------------------------------------------------------------
test_that("sex.cess.mod parameterization: build_* pipeline works", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7,
    age.limits = c(15, 100),
    age.sexual.cessation = 65
  )
  expect_true(epistats$sex.cess.mod)

  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "existing")
  set.seed(20260419L)
  ns <- build_netstats(epistats, np, expect.mort = 0.000478213,
                       network.size = 3000, young.prop = 0.99,
                       method = "existing")
  expect_netstats_contract(ns)

  # active.sex is 0 for the post-cessation age group, 1 otherwise.
  expect_true(all(ns$attr$active.sex %in% c(0L, 1L)))
  expect_gt(sum(ns$attr$active.sex == 0L), 0)  # some inactive nodes exist
  expect_gt(sum(ns$attr$active.sex == 1L), 0)  # some active nodes exist

  # nodefactor_age.grp[last] should be 0 under sex.cess.mod -- inactive
  # nodes contribute no edge endpoints.
  last <- length(ns$main$nodefactor_age.grp)
  expect_equal(ns$main$nodefactor_age.grp[last], 0)
  expect_equal(ns$casl$nodefactor_age.grp[last], 0)

  # Joint path under sex.cess.mod
  set.seed(20260419L)
  np_j <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  set.seed(20260419L)
  ns_j <- build_netstats(epistats, np_j, expect.mort = 0.000478213,
                         network.size = 3000, young.prop = 0.99,
                         method = "joint")
  expect_netstats_contract(ns_j)
  # Internal consistency still holds: inactive egos have pred_deg zeroed,
  # so their contribution to nodefactor sums is 0 and the identity holds.
  expect_equal(sum(ns_j$main$nodefactor_age.grp), 2 * ns_j$main$edges,
               tolerance = 1e-9)
})


# ---------------------------------------------------------------------------
# Non-Atlanta city (smoke test for other geog.cat values)
# ---------------------------------------------------------------------------
test_that("non-Atlanta city parameterization works (New York City)", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "New York City",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  expect_equal(epistats$geog.lvl, "city")
  expect_equal(epistats$geog.cat, "New York City")

  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "existing")
  set.seed(20260419L)
  ns <- build_netstats(epistats, np, expect.mort = 0.000478213,
                       network.size = 3000, method = "existing")
  expect_netstats_contract(ns)
})


# ---------------------------------------------------------------------------
# Custom age.breaks and age.limits
# ---------------------------------------------------------------------------
test_that("custom age.breaks / age.limits parameterization works", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "state", geog.cat = "GA",
    race = TRUE, time.unit = 7,
    age.limits = c(20, 50),
    age.breaks = c(30, 40)
  )
  expect_equal(epistats$age.limits, c(20, 50))
  # age.grps = length(age.breaks) + 1 = 3 groups: (20,30], (30,40], (40,50]
  expect_equal(epistats$age.grps, 3)

  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "existing")
  set.seed(20260419L)
  ns <- build_netstats(epistats, np, expect.mort = 0.000478213,
                       network.size = 3000, method = "existing")
  expect_netstats_contract(ns)

  # nodefactor_age.grp length matches the configured age groups
  expect_length(ns$main$nodefactor_age.grp, 3)
  expect_length(ns$casl$nodefactor_age.grp, 3)
  expect_length(ns$inst$nodefactor_age.grp, 3)

  # Sampled ages fall within the configured limits
  expect_true(all(ns$attr$age >= 20 & ns$attr$age < 50))
})


# ---------------------------------------------------------------------------
# Non-default time.unit (monthly instead of weekly)
# ---------------------------------------------------------------------------
test_that("non-default time.unit parameterization works", {
  skip_without_artnetdata()
  set.seed(20260419L)
  ep_weekly <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  ep_monthly <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 30
  )
  expect_equal(ep_weekly$time.unit, 7)
  expect_equal(ep_monthly$time.unit, 30)

  set.seed(20260419L)
  np_w <- build_netparams(ep_weekly, smooth.main.dur = TRUE, method = "existing")
  set.seed(20260419L)
  np_m <- build_netparams(ep_monthly, smooth.main.dur = TRUE, method = "existing")

  # md.main is per-respondent mean degree: unit-invariant (just a count)
  expect_equal(np_w$main$md.main, np_m$main$md.main, tolerance = 1e-9)

  # md.inst is per-time-unit (annual count / (364 / time.unit)).
  # Monthly should be larger than weekly by factor (30 / 7).
  expect_equal(np_m$inst$md.inst / np_w$inst$md.inst, 30 / 7,
               tolerance = 1e-6)
})
