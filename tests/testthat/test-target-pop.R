# Tests for the target_pop API on build_netstats (#64).
# Three input forms: list of marginal distributions, data.frame of nodes,
# character (not yet implemented). NULL preserves legacy behavior.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

setup_pipeline <- function() {
  set.seed(20260420L)
  ep <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                       init.hiv.prev = c(0.33, 0.137, 0.084),
                       race = TRUE, time.unit = 7)
  set.seed(20260420L)
  np <- build_netparams(ep, smooth.main.dur = TRUE)
  list(epistats = ep, netparams = np)
}

test_that("target_pop = NULL is byte-identical to no target_pop arg", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  ns_default <- build_netstats(s$epistats, s$netparams,
                               expect.mort = 0.000478213, network.size = 2000)
  set.seed(20260420L)
  ns_explicit <- build_netstats(s$epistats, s$netparams,
                                expect.mort = 0.000478213, network.size = 2000,
                                target_pop = NULL)
  expect_equal(ns_default$attr, ns_explicit$attr)
  expect_equal(ns_default$main$edges, ns_explicit$main$edges)
})


# ---- list form -------------------------------------------------------------

test_that("list form: race.prop override produces matching race composition", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213, network.size = 5000,
                       target_pop = list(race.prop = c(0.4, 0.2, 0.4)))
  obs <- prop.table(table(ns$attr$race))
  expect_equal(as.numeric(obs), c(0.4, 0.2, 0.4), tolerance = 0.01)
})

test_that("list form: deg.casl override produces matching distribution", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213, network.size = 5000,
                       target_pop = list(deg.casl = c(0.5, 0.3, 0.15, 0.05)))
  obs <- prop.table(table(ns$attr$deg.casl))
  expect_equal(as.numeric(obs), c(0.5, 0.3, 0.15, 0.05), tolerance = 0.005)
})

test_that("list form: race.props alias is normalized to race.prop", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213, network.size = 2000,
                       target_pop = list(race.props = c(0.5, 0.25, 0.25)))
  obs <- prop.table(table(ns$attr$race))
  expect_equal(as.numeric(obs), c(0.5, 0.25, 0.25), tolerance = 0.02)
})

test_that("list form: unknown elements raise an informative error", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  expect_error(
    build_netstats(s$epistats, s$netparams,
                   expect.mort = 0.000478213, network.size = 1000,
                   target_pop = list(foo = 1, age.pyramid = NULL)),
    regexp = "unknown elements: foo"
  )
})


# ---- data.frame form -------------------------------------------------------

test_that("data.frame form bypasses sampling and respects user attrs", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  df <- data.frame(
    age        = sample(15:64, 1500, replace = TRUE),
    race       = sample(1:3, 1500, replace = TRUE, prob = c(0.4, 0.2, 0.4)),
    deg.casl   = sample(0:3, 1500, replace = TRUE),
    deg.main   = sample(0:2, 1500, replace = TRUE),
    role.class = sample(0:2, 1500, replace = TRUE),
    risk.grp   = sample(1:5, 1500, replace = TRUE),
    diag.status = rbinom(1500, 1, 0.15)
  )
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213,
                       network.size = 99999,  # must be ignored
                       target_pop = df)
  expect_equal(ns$demog$num, 1500)
  expect_length(ns$attr$age, 1500)
  expect_equal(ns$attr$age, df$age)
  expect_equal(ns$attr$race, df$race)
  expect_equal(ns$attr$deg.casl, df$deg.casl)
  expect_equal(ns$attr$diag.status, as.integer(df$diag.status))
  # Derived attrs filled in
  expect_length(ns$attr$sqrt.age, 1500)
  expect_length(ns$attr$age.grp, 1500)
  expect_length(ns$attr$deg.tot, 1500)
  expect_length(ns$attr$active.sex, 1500)
})

test_that("data.frame form: derived deg.tot caps at 3", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  n <- 500
  df <- data.frame(
    age        = sample(15:64, n, replace = TRUE),
    race       = sample(1:3, n, replace = TRUE),
    role.class = sample(0:2, n, replace = TRUE),
    risk.grp   = sample(1:5, n, replace = TRUE),
    deg.casl   = c(rep(3L, 5), sample(0:3, n - 5, replace = TRUE)),
    deg.main   = c(rep(2L, 5), sample(0:2, n - 5, replace = TRUE))
  )
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213, network.size = 100,
                       target_pop = df)
  # First 5 rows have deg.casl = 3, deg.main = 2: raw sum 5, capped to 3.
  expect_equal(ns$attr$deg.tot[1:5], rep(3L, 5))
  # All values must satisfy the cap.
  expect_true(all(ns$attr$deg.tot <= 3L))
  expect_true(all(ns$attr$deg.tot >= 0L))
})

test_that("data.frame form: missing required columns raise informative error", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  expect_error(
    build_netstats(s$epistats, s$netparams,
                   expect.mort = 0.000478213, network.size = 100,
                   target_pop = data.frame(age = 1:5)),
    regexp = "data.frame missing required columns"
  )
})

test_that("data.frame form: diag.status falls back to epistats when absent", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  set.seed(20260420L)
  df <- data.frame(
    age        = sample(15:64, 500, replace = TRUE),
    race       = sample(1:3, 500, replace = TRUE),
    deg.casl   = sample(0:3, 500, replace = TRUE),
    deg.main   = sample(0:2, 500, replace = TRUE),
    role.class = sample(0:2, 500, replace = TRUE),
    risk.grp   = sample(1:5, 500, replace = TRUE)
  )
  ns <- build_netstats(s$epistats, s$netparams,
                       expect.mort = 0.000478213, network.size = 100,
                       target_pop = df)
  expect_length(ns$attr$diag.status, 500)
  expect_true(all(ns$attr$diag.status %in% c(0L, 1L)))
})

test_that("data.frame form: composes with method = 'joint'", {
  skip_without_artnetdata()
  set.seed(20260420L)
  ep <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                       init.hiv.prev = c(0.33, 0.137, 0.084),
                       race = TRUE, time.unit = 7)
  set.seed(20260420L)
  np <- build_netparams(ep, smooth.main.dur = TRUE,
                        method = "joint", duration.method = "joint_lm")
  df <- data.frame(
    age        = sample(15:64, 800, replace = TRUE),
    race       = sample(1:3, 800, replace = TRUE),
    deg.casl   = sample(0:3, 800, replace = TRUE),
    deg.main   = sample(0:2, 800, replace = TRUE),
    role.class = sample(0:2, 800, replace = TRUE),
    risk.grp   = sample(1:5, 800, replace = TRUE),
    diag.status = rbinom(800, 1, 0.15)
  )
  ns <- build_netstats(ep, np,
                       expect.mort = 0.000478213, network.size = 100,
                       target_pop = df, method = "joint")
  # Internal consistency under joint must still hold
  expect_equal(sum(ns$main$nodefactor_race), 2 * ns$main$edges,
               tolerance = 1e-9)
})


# ---- character form --------------------------------------------------------

test_that("character form raises informative not-yet-implemented error", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  expect_error(
    build_netstats(s$epistats, s$netparams,
                   expect.mort = 0.000478213, network.size = 1000,
                   target_pop = "nhbs_msm_2022"),
    regexp = "not yet implemented"
  )
})


# ---- bad input -------------------------------------------------------------

test_that("non-list non-data.frame non-character input raises error", {
  skip_without_artnetdata()
  s <- setup_pipeline()
  expect_error(
    build_netstats(s$epistats, s$netparams,
                   expect.mort = 0.000478213, network.size = 1000,
                   target_pop = 42),
    regexp = "must be NULL, a list, a data.frame"
  )
})
