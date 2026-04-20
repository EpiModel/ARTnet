# Tests for the dyad-level (nodematch + absdiff) joint GLMs added for
# issue #63 phases 1 and 2. See test-joint-model.R and test-joint-netstats.R
# for the ego-level (#61/#62) coverage.

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
  np <- build_netparams(epistats, smooth.main.dur = TRUE,
                       method = netparams_method)
  set.seed(20260419L)
  ns <- build_netstats(epistats, np,
                      expect.mort = 0.000478213, network.size = 5000,
                      method = netstats_method)
  out <- list(epistats = epistats, netparams = np, netstats = ns)
  cache_env[[key]] <- out
  out
}

test_that("method='joint' produces all dyad models per layer", {
  skip_without_artnetdata()
  np <- get_stats("joint", "joint")$netparams
  for (layer in c("main", "casl", "inst")) {
    expect_s3_class(np[[layer]]$joint_nm_age_model, "glm")
    expect_s3_class(np[[layer]]$joint_nm_race_model, "glm")
    expect_s3_class(np[[layer]]$joint_absdiff_age_model, "glm")
    expect_s3_class(np[[layer]]$joint_absdiff_sqrtage_model, "glm")
  }
})

test_that("method='existing' does not produce dyad models", {
  skip_without_artnetdata()
  np <- get_stats("existing", "existing")$netparams
  for (layer in c("main", "casl", "inst")) {
    expect_null(np[[layer]]$joint_nm_age_model)
    expect_null(np[[layer]]$joint_nm_race_model)
    expect_null(np[[layer]]$joint_absdiff_age_model)
    expect_null(np[[layer]]$joint_absdiff_sqrtage_model)
  }
})

test_that("dyad models converge with expected families", {
  skip_without_artnetdata()
  np <- get_stats("joint", "joint")$netparams
  for (layer in c("main", "casl", "inst")) {
    expect_true(isTRUE(np[[layer]]$joint_nm_age_model$converged))
    expect_true(isTRUE(np[[layer]]$joint_nm_race_model$converged))
    expect_true(isTRUE(np[[layer]]$joint_absdiff_age_model$converged))
    expect_true(isTRUE(np[[layer]]$joint_absdiff_sqrtage_model$converged))
    expect_identical(np[[layer]]$joint_nm_age_model$family$family,  "binomial")
    expect_identical(np[[layer]]$joint_nm_race_model$family$family, "binomial")
    expect_identical(np[[layer]]$joint_absdiff_age_model$family$family,     "gaussian")
    expect_identical(np[[layer]]$joint_absdiff_sqrtage_model$family$family, "gaussian")
  }
})

test_that("dyad models recover training-data marginal means within 1%", {
  skip_without_artnetdata()
  np <- get_stats("joint", "joint")$netparams
  responses <- list(
    joint_nm_age_model  = "same.age.grp",
    joint_nm_race_model = "same.race",
    joint_absdiff_age_model     = "ad",
    joint_absdiff_sqrtage_model = "ad.sr"
  )
  for (layer in c("main", "casl", "inst")) {
    for (modname in names(responses)) {
      m <- np[[layer]][[modname]]
      obs <- mean(m$data[[responses[[modname]]]], na.rm = TRUE)
      pred <- mean(predict(m, type = "response"))
      expect_lt(abs(pred - obs) / obs, 0.01,
                label = sprintf("[%s/%s] marginal recovery", layer, modname))
    }
  }
})

test_that("joint nodematch_race_diffF equals sum of nodematch_race (diff)", {
  skip_without_artnetdata()
  ns <- get_stats("joint", "joint")$netstats
  for (layer in c("main", "casl", "inst")) {
    expect_equal(sum(ns[[layer]]$nodematch_race),
                 ns[[layer]]$nodematch_race_diffF,
                 tolerance = 1e-9,
                 info = paste("layer =", layer))
  }
})

test_that("joint nodematch counts do not exceed edge counts", {
  skip_without_artnetdata()
  ns <- get_stats("joint", "joint")$netstats
  # Edges where both endpoints have attribute r can be at most the total
  # edges. sum(nodematch_<attr>) <= edges.
  for (layer in c("main", "casl", "inst")) {
    expect_lte(sum(ns[[layer]]$nodematch_age.grp), ns[[layer]]$edges + 1e-6)
    expect_lte(sum(ns[[layer]]$nodematch_race),    ns[[layer]]$edges + 1e-6)
  }
})

test_that("joint nodematch/absdiff differ from existing method", {
  skip_without_artnetdata()
  # Concrete evidence this PR actually changes something under method='joint'
  ex <- get_stats("existing", "existing")$netstats
  jt <- get_stats("joint", "joint")$netstats
  # The absdiff_age target should differ by at least 1% on the main layer —
  # this is the effect the refactor is supposed to create.
  rel_diff <- abs(jt$main$absdiff_age - ex$main$absdiff_age) / ex$main$absdiff_age
  expect_gt(rel_diff, 0.01)
  rel_diff2 <- abs(jt$main$nodematch_race_diffF -
                   ex$main$nodematch_race_diffF) / ex$main$nodematch_race_diffF
  expect_gt(rel_diff2, 0.01)
})

test_that("race = FALSE skips nm_race model gracefully", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = FALSE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  # nm_age, absdiff models should be present; nm_race should not
  expect_s3_class(np$main$joint_nm_age_model, "glm")
  expect_null(np$main$joint_nm_race_model)
  expect_s3_class(np$main$joint_absdiff_age_model, "glm")

  set.seed(20260419L)
  ns <- build_netstats(epistats, np,
                      expect.mort = 0.000478213, network.size = 5000,
                      method = "joint")
  # nodematch_race shouldn't exist under race = FALSE
  expect_null(ns$main$nodematch_race)
  expect_null(ns$main$nodematch_race_diffF)
  # nodematch_age.grp still valid
  expect_lte(sum(ns$main$nodematch_age.grp), ns$main$edges + 1e-6)
})
