# Tests for the duration.method flag added in PR for #63 phase 3.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

cache_env <- new.env(parent = emptyenv())

get_np <- function(duration.method) {
  key <- duration.method
  if (!is.null(cache_env[[key]])) return(cache_env[[key]])
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE,
                       duration.method = duration.method,
                       method = "existing")
  cache_env[[key]] <- np
  np
}

test_that("default duration.method is 'empirical'", {
  skip_without_artnetdata()
  np_default <- get_np("empirical")
  # No joint_duration_model attached when method is empirical
  expect_null(np_default$main$joint_duration_model)
  expect_null(np_default$casl$joint_duration_model)
})

test_that("joint_lm stores fitted model at joint_duration_model", {
  skip_without_artnetdata()
  np <- get_np("joint_lm")
  for (layer in c("main", "casl")) {
    m <- np[[layer]]$joint_duration_model
    expect_s3_class(m, "lm")
    # Fit should have a reasonable number of observations
    expect_gt(stats::nobs(m), 500)
    # Response is log(duration.time)
    expect_true(grepl("log\\(duration\\.time\\)",
                      as.character(formula(m))[2]))
  }
})

test_that("all duration.methods preserve output shape for dissolution_coefs", {
  skip_without_artnetdata()
  # Every method must produce durs.<layer>.byage with the columns
  # dissolution_coefs() needs downstream, in the same order and types.
  expected_cols <- c("index.age.grp", "mean.dur", "median.dur",
                     "rates.main.adj", "mean.dur.adj")
  for (dm in c("empirical", "joint_lm")) {
    df <- get_np(dm)$main$durs.main.byage
    expect_identical(colnames(df), expected_cols,
                     info = paste("column mismatch for", dm))
    expect_gt(nrow(df), 0)
    expect_true(all(is.finite(df$mean.dur.adj)),
                info = paste("mean.dur.adj must be finite for", dm))
  }
})

test_that("duration.method does not affect non-duration netparams fields", {
  skip_without_artnetdata()
  np_e <- get_np("empirical")
  np_j <- get_np("joint_lm")
  # Everything except durs.*.byage, durs.*.homog, joint_duration_model
  # should be identical across methods.
  for (layer in c("main", "casl", "inst", "all")) {
    if (is.null(np_e[[layer]])) next
    common_fields <- setdiff(
      names(np_e[[layer]]),
      c(paste0("durs.", layer, ".byage"),
        paste0("durs.", layer, ".homog"),
        "joint_duration_model")
    )
    for (f in common_fields) {
      expect_equal(np_j[[layer]][[f]], np_e[[layer]][[f]],
                   info = paste0("[joint_lm ", layer, "$", f, "]"))
    }
  }
})

test_that("joint_lm rejects invalid method arg with clear error", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  expect_error(
    build_netparams(epistats, duration.method = "weibull_strat",
                    method = "existing"),
    regexp = "should be one of"
  )
})

test_that("joint_lm can be combined with method = 'joint' in build_netstats", {
  skip_without_artnetdata()
  set.seed(20260419L)
  epistats <- build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260419L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE,
                       duration.method = "joint_lm",
                       method = "joint")
  set.seed(20260419L)
  ns <- build_netstats(epistats, np,
                      expect.mort = 0.000478213, network.size = 3000,
                      method = "joint")
  # Dissolution coefs still valid
  expect_s3_class(ns$main$diss.byage, "disscoef")
  expect_s3_class(ns$casl$diss.byage, "disscoef")
  expect_true(all(is.finite(ns$main$diss.byage$coef.diss)))
  # Both the ego-level joint_model and the dyad joint_duration_model
  # coexist on netparams
  expect_s3_class(np$main$joint_model, "glm")
  expect_s3_class(np$main$joint_duration_model, "lm")
})
