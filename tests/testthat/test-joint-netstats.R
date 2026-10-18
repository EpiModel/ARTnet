# Tests for the method = "joint" path added to build_netstats() in issue #62.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

cache_env <- new.env(parent = emptyenv())

get_stats <- function(netparams_method, netstats_method) {
  key <- paste0(netparams_method, "_", netstats_method)
  if (!is.null(cache_env[[key]])) return(cache_env[[key]])

  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = netparams_method)
  set.seed(20260419L)
  ns <- build_netstats(epistats, np,
                      expect.mort = 0.000478213, network.size = 5000,
                      method = netstats_method)
  out <- list(epistats = epistats, netparams = np, netstats = ns)
  cache_env[[key]] <- out
  out
}

test_that("method = 'joint' requires netparams built with method = 'joint'", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np_existing <- build_netparams(epistats, smooth.main.dur = TRUE,
                                 method = "existing")
  expect_error(
    build_netstats(epistats, np_existing,
                   expect.mort = 0.000478213, network.size = 5000,
                   method = "joint"),
    regexp = "method = 'joint' requires netparams"
  )
})

test_that("joint netstats are internally consistent: 2 * edges == sum(nodefactor_*)", {
  skip_without_artnetdata()
  ns <- get_stats("joint", "joint")$netstats
  for (layer in c("main", "casl", "inst")) {
    e <- ns[[layer]]$edges
    expect_equal(sum(ns[[layer]]$nodefactor_race),      2 * e,
                 tolerance = 1e-9,
                 info = paste0(layer, ": 2*edges vs sum(nf_race)"))
    expect_equal(sum(ns[[layer]]$nodefactor_age.grp),   2 * e,
                 tolerance = 1e-9,
                 info = paste0(layer, ": 2*edges vs sum(nf_age.grp)"))
    expect_equal(sum(ns[[layer]]$nodefactor_diag.status), 2 * e,
                 tolerance = 1e-9,
                 info = paste0(layer, ": 2*edges vs sum(nf_diag.status)"))
  }
  # Layer-specific deg-attr nodefactor sums
  expect_equal(sum(ns$main$nodefactor_deg.casl), 2 * ns$main$edges,
               tolerance = 1e-9)
  expect_equal(sum(ns$casl$nodefactor_deg.main), 2 * ns$casl$edges,
               tolerance = 1e-9)
  expect_equal(sum(ns$inst$nodefactor_deg.tot),  2 * ns$inst$edges,
               tolerance = 1e-9)
})

test_that("joint netstats preserve existing dissolution coefs (numeric content)", {
  skip_without_artnetdata()
  existing <- get_stats("existing", "existing")$netstats
  joint    <- get_stats("joint", "joint")$netstats
  # dissolution_coefs() returns a list containing a formula, which captures
  # the caller's scope environment — so full-object expect_equal() differs
  # between the existing and joint code paths even though coefficients,
  # durations, and rates are identical. Compare numeric content only.
  diss_numeric <- function(x) {
    x$dissolution <- NULL  # drop formula/env
    x
  }
  expect_equal(diss_numeric(existing$main$diss.homog), diss_numeric(joint$main$diss.homog))
  expect_equal(diss_numeric(existing$main$diss.byage), diss_numeric(joint$main$diss.byage))
  expect_equal(diss_numeric(existing$casl$diss.homog), diss_numeric(joint$casl$diss.homog))
  expect_equal(diss_numeric(existing$casl$diss.byage), diss_numeric(joint$casl$diss.byage))
})

test_that("joint netstats preserve non-refactored fields exactly", {
  skip_without_artnetdata()
  existing <- get_stats("existing", "existing")$netstats
  joint    <- get_stats("joint", "joint")$netstats
  # risk.grp nodefactor (inst) is NOT refactored — keeps univariate.
  expect_equal(existing$inst$nodefactor_risk.grp, joint$inst$nodefactor_risk.grp)
  # attributes are sampled identically when seeded, so they should match
  expect_equal(existing$attr, joint$attr)
  # demography likewise unchanged
  expect_equal(existing$demog$num, joint$demog$num)
})

test_that("joint concurrent is in a reasonable range (not inflated by Poisson)", {
  skip_without_artnetdata()
  existing <- get_stats("existing", "existing")$netstats
  joint    <- get_stats("joint", "joint")$netstats
  # Joint concurrent should stay in the same order of magnitude as existing.
  # A Poisson-derived P(deg>1) would inflate main concurrent ~5x because of
  # the deg.main truncation at 2 in the training data; the binomial joint
  # model avoids that.
  ratio_main <- joint$main$concurrent / existing$main$concurrent
  expect_gt(ratio_main, 0.5)
  expect_lt(ratio_main, 2.0)
  ratio_casl <- joint$casl$concurrent / existing$casl$concurrent
  expect_gt(ratio_casl, 0.5)
  expect_lt(ratio_casl, 2.0)
})

test_that("joint edges differ materially from existing when population shifts", {
  skip_without_artnetdata()
  # Run under a race distribution that diverges from Atlanta defaults.
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  set.seed(20260419L)
  ns_ex <- build_netstats(epistats, np,
                          expect.mort = 0.000478213, network.size = 5000,
                          race.prop = c(0.35, 0.25, 0.40),
                          method = "existing")
  set.seed(20260419L)
  ns_jt <- build_netstats(epistats, np,
                          expect.mort = 0.000478213, network.size = 5000,
                          race.prop = c(0.35, 0.25, 0.40),
                          method = "joint")
  # Expect >1% divergence between methods on at least one layer under a
  # shifted population — this is the empirical point of the refactor.
  divergence <- abs(ns_jt$main$edges - ns_ex$main$edges) / ns_ex$main$edges
  expect_gt(divergence, 0.01)
  # Joint stays internally consistent even with shifted race.prop.
  expect_equal(sum(ns_jt$main$nodefactor_race), 2 * ns_jt$main$edges,
               tolerance = 1e-9)
})

test_that("joint netstats work with race = FALSE", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = FALSE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  set.seed(20260419L)
  ns <- build_netstats(epistats, np,
                       expect.mort = 0.000478213, network.size = 5000,
                       method = "joint")
  # race nodefactor/match not populated when race = FALSE
  expect_null(ns$main$nodefactor_race)
  expect_null(ns$main$nodematch_race)
  # but edges and age.grp still internally consistent
  expect_equal(sum(ns$main$nodefactor_age.grp), 2 * ns$main$edges,
               tolerance = 1e-9)
})
