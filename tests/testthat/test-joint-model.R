# Tests for the method = "joint" path added in issue #61.
# Require ARTnetData; skipped silently otherwise so CI stays green without it.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

cache_env <- new.env(parent = emptyenv())

get_netparams <- function(method) {
  key <- paste0("netparams_", method)
  if (!is.null(cache_env[[key]])) return(cache_env[[key]])
  # Seed because build_netparams() does a stochastic NA-imputation of
  # partner races via sample(); without a seed, nm.race values drift
  # between back-to-back calls.
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = method)
  cache_env[[key]] <- np
  np
}

test_that("default method is 'existing' and omits joint_model fields", {
  skip_without_artnetdata()
  np_default <- get_netparams("existing")
  expect_null(np_default$main$joint_model)
  expect_null(np_default$casl$joint_model)
  expect_null(np_default$inst$joint_model)
})

test_that("method='joint' adds a converged glm per layer", {
  skip_without_artnetdata()
  np <- get_netparams("joint")
  for (layer in c("main", "casl", "inst")) {
    m <- np[[layer]]$joint_model
    expect_s3_class(m, "glm")
    expect_true(isTRUE(m$converged),
                info = paste("joint_model did not converge for layer", layer))
  }
})

test_that("joint models recover observed marginal means within 1%", {
  skip_without_artnetdata()
  np <- get_netparams("joint")
  responses <- c(main = "deg.main", casl = "deg.casl",
                 inst = "count.oo.part")
  for (layer in names(responses)) {
    m <- np[[layer]]$joint_model
    obs <- mean(m$data[[responses[[layer]]]], na.rm = TRUE)
    pred <- mean(predict(m, type = "response"))
    expect_lt(abs(pred - obs) / obs, 0.01,
              label = paste0("[", layer, "] relative diff"))
  }
})

test_that("joint model coefficients match expected signs (sanity)", {
  skip_without_artnetdata()
  np <- get_netparams("joint")
  co <- coef(np$main$joint_model)
  # The age profile is peaked-then-declining, which in this parameterization
  # means age.grp slope is negative and sqrt(age.grp) slope is positive.
  expect_lt(co[["age.grp"]], 0)
  expect_gt(co[["sqrt(age.grp)"]], 0)
  # Race ref level (1 = NH Black) should be lower than 2 and 3 — so coefs
  # for the non-reference levels should be positive.
  expect_gt(co[["as.factor(race.cat.num)2"]], 0)
  expect_gt(co[["as.factor(race.cat.num)3"]], 0)
})

test_that("AIC-based interaction selection records what it considered", {
  skip_without_artnetdata()
  np <- get_netparams("joint")
  for (layer in c("main", "casl", "inst")) {
    m <- np[[layer]]$joint_model
    # Every selected interaction must come from the candidate pool.
    sel <- attr(m, "selected_interactions")
    cand <- attr(m, "candidate_interactions")
    expect_type(sel, "character")
    expect_type(cand, "character")
    expect_true(all(sel %in% cand))
  }
})

test_that("method='joint' leaves existing marginal fields untouched", {
  skip_without_artnetdata()
  np_existing <- get_netparams("existing")
  np_joint <- get_netparams("joint")
  for (layer in c("main", "casl", "inst", "all")) {
    if (is.null(np_existing[[layer]])) next
    fields <- setdiff(names(np_existing[[layer]]), "joint_model")
    for (f in fields) {
      expect_equal(np_joint[[layer]][[f]], np_existing[[layer]][[f]],
                   info = paste0("[", layer, "$", f, "] diverged under joint"))
    }
  }
})

test_that("method='joint' works with race = FALSE (skips race terms)", {
  skip_without_artnetdata()
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = FALSE, time.unit = 7
  )
  np <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
  m <- np$main$joint_model
  expect_s3_class(m, "glm")
  expect_true(m$converged)
  # No race.cat.num term should appear
  expect_false(any(grepl("race.cat.num", names(coef(m)))))
})
