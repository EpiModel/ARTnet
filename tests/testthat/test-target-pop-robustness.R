# Robustness tests for build_netstats() with target_pop, covering bug fixes
# surfaced during the ARTnetPredict forward-projection bridge review:
#  1. joint_lm + target_pop must not crash on empty/sparse age strata
#     (per-stratum fallback to netparams byage durations).
#  2. Out-of-support degree/role/risk values are clamped (with a warning) so
#     the sum(nodefactor) == 2 * edges invariant is preserved.
#  3. nodefactor over a fixed level set never recycles when a level is
#     unobserved in the target population.
#  4. The `seed` argument makes the stochastic sampling path reproducible.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

have_data <- system.file(package = "ARTnetData") != ""
if (have_data) {
  set.seed(20260420L)
  .ep <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                        init.hiv.prev = c(0.33, 0.137, 0.084),
                        race = TRUE, time.unit = 7)
  set.seed(20260420L)
  .np_joint <- build_netparams(.ep, smooth.main.dur = TRUE,
                               method = "joint", duration.method = "joint_lm")
}

mk_df <- function(n, age_lo, age_hi, risk_levels = 1:5) {
  data.frame(
    age        = stats::runif(n, age_lo, age_hi),
    race       = sample(1:3, n, replace = TRUE, prob = c(0.5, 0.05, 0.45)),
    deg.casl   = sample(0:3, n, replace = TRUE),
    deg.main   = sample(0:2, n, replace = TRUE),
    role.class = sample(0:2, n, replace = TRUE),
    risk.grp   = sample(risk_levels, n, replace = TRUE)
  )
}


# ---- 1. joint_lm + empty age stratum must not crash ------------------------

test_that("joint_lm + target_pop with an empty top age stratum does not crash", {
  skip_without_artnetdata()
  set.seed(11L)
  # age in [20, 54] leaves the 55-64 age group (group 5) empty under the
  # default breaks 15/25/35/45/55/65. Before the fix this aborted in
  # dissolution_coefs with "missing value where TRUE/FALSE needed".
  df <- mk_df(3000, 20, 54)
  expect_true(all(df$age < 55))
  ns <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                       method = "joint", target_pop = df)
  expect_true(is.finite(ns$main$edges))
  expect_true(all(is.finite(ns$main$diss.byage$coef.diss)))
  expect_true(all(is.finite(ns$casl$diss.byage$coef.diss)))
})

test_that("joint_lm empty-stratum durations fall back to netparams byage", {
  skip_without_artnetdata()
  set.seed(12L)
  df <- mk_df(3000, 20, 54)
  ns <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                       method = "joint", target_pop = df)
  # diss.byage offset must be fully populated (one edges coef + one per
  # age-group nodematch offset), all finite, no NA leaking through.
  expect_false(anyNA(ns$main$diss.byage$coef.diss))
  expect_false(anyNA(ns$casl$diss.byage$coef.diss))
})


# ---- 2. out-of-support degrees are clamped, invariant preserved ------------

test_that("out-of-support deg.casl is clamped with a warning", {
  skip_without_artnetdata()
  set.seed(21L)
  df <- mk_df(2000, 15, 64)
  df$deg.casl[1:50] <- 4L  # 4 is outside the fitted support 0:3
  expect_warning(
    build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                   method = "joint", target_pop = df),
    regexp = "deg.casl.*outside the fitted support"
  )
})

test_that("clamping preserves the deg.casl nodefactor invariant under joint", {
  skip_without_artnetdata()
  set.seed(22L)
  df <- mk_df(2000, 15, 64)
  df$deg.casl[1:50] <- 4L
  ns <- suppressWarnings(
    build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                   method = "joint", target_pop = df))
  # After clamping, sum over the deg.casl-keyed nodefactor must equal 2*edges.
  expect_equal(sum(ns$main$nodefactor_deg.casl), 2 * ns$main$edges,
               tolerance = 1e-8)
})

test_that("in-support degrees produce no clamp warning", {
  skip_without_artnetdata()
  set.seed(23L)
  df <- mk_df(1000, 15, 64)
  expect_no_warning(
    build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                   method = "joint", target_pop = df)
  )
})


# ---- 3. nodefactor over a fixed level set does not recycle -----------------

test_that("nodefactor_risk.grp spans all quintiles when a level is absent", {
  skip_without_artnetdata()
  set.seed(31L)
  nquants <- length(.np_joint$inst$nf.risk.grp)
  # supply risk.grp only in {1, 2}; levels 3..nquants are unobserved
  df <- mk_df(1500, 15, 64, risk_levels = 1:2)
  ns <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                       method = "joint", target_pop = df)
  expect_length(ns$inst$nodefactor_risk.grp, nquants)
  # unobserved quintiles contribute exactly zero (no recycling)
  expect_equal(as.numeric(ns$inst$nodefactor_risk.grp[3:nquants]),
               rep(0, nquants - 2))
})


# ---- 4. seed reproducibility for the sampling path -------------------------

test_that("seed makes the sampling path reproducible under method = 'joint'", {
  skip_without_artnetdata()
  ns1 <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                        network.size = 2000, method = "joint", seed = 123L)
  ns2 <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                        network.size = 2000, method = "joint", seed = 123L)
  expect_equal(ns1$attr$age, ns2$attr$age)
  expect_equal(ns1$main$edges, ns2$main$edges)
  expect_equal(ns1$casl$edges, ns2$casl$edges)
})

test_that("different seeds give different draws", {
  skip_without_artnetdata()
  ns1 <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                        network.size = 2000, method = "joint", seed = 123L)
  ns3 <- build_netstats(.ep, .np_joint, expect.mort = 0.000478213,
                        network.size = 2000, method = "joint", seed = 456L)
  expect_false(isTRUE(all.equal(ns1$main$edges, ns3$main$edges)))
})
