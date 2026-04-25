# Tests for the method-comparison validation harness in
# inst/validation/method_comparison.R (Phase 1.5 / issue #65).
#
# These tests exercise the helper structure on a single small scenario
# rather than the full suite (full suite is intentionally slow and
# produces inst/validation/method_comparison.md as its real output).

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

source_helper <- function() {
  path <- system.file("validation/method_comparison.R", package = "ARTnet")
  if (!nzchar(path)) {
    path <- "inst/validation/method_comparison.R"
  }
  testthat::skip_if(!file.exists(path), "method_comparison.R not found")
  source(path, local = parent.frame())
}

mini_scenario <- list(
  list(
    name = "test_atlanta",
    epistats = list(geog.lvl = "city", geog.cat = "Atlanta",
                    init.hiv.prev = c(0.33, 0.137, 0.084),
                    race = TRUE, time.unit = 7),
    netparams = list(smooth.main.dur = TRUE),
    netstats  = list(expect.mort = 0.000478213, network.size = 2000)
  )
)

test_that("compare_methods returns expected long-format structure", {
  skip_without_artnetdata()
  source_helper()
  res <- suppressMessages(compare_methods(mini_scenario))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("scenario", "layer", "stat", "level",
                      "existing", "joint", "abs_diff", "pct_diff"),
               ignore.order = TRUE)
  expect_true(all(res$layer %in% c("main", "casl", "inst")))
  expect_true(all(res$scenario == "test_atlanta"))
  expect_true(all(c("edges", "nodefactor_race", "nodefactor_age.grp",
                    "dissolution_duration") %in% res$stat))
  # All stats produced both methods
  expect_false(any(is.na(res$existing)))
  expect_false(any(is.na(res$joint)))
})

test_that("abs_diff and pct_diff are computed consistently", {
  skip_without_artnetdata()
  source_helper()
  res <- suppressMessages(compare_methods(mini_scenario))
  expect_equal(res$abs_diff, res$joint - res$existing, tolerance = 1e-9)
  # When existing is non-zero, pct_diff matches the formula
  ok <- abs(res$existing) > 1e-12 & !is.na(res$pct_diff)
  expect_equal(res$pct_diff[ok],
               100 * res$abs_diff[ok] / res$existing[ok],
               tolerance = 1e-9)
})

test_that("comparison includes dissolution_duration for main and casl, not inst", {
  skip_without_artnetdata()
  source_helper()
  res <- suppressMessages(compare_methods(mini_scenario))
  dur <- res[res$stat == "dissolution_duration", , drop = FALSE]
  expect_true(all(dur$layer %in% c("main", "casl")))
  expect_true(any(dur$layer == "main"))
  expect_true(any(dur$layer == "casl"))
})

test_that("at least one cell shifts > 5% between methods on Atlanta default", {
  skip_without_artnetdata()
  source_helper()
  res <- suppressMessages(compare_methods(mini_scenario))
  # The whole point of the refactor is that joint differs from existing
  # in at least some places. If this isn't true, something is broken.
  expect_true(sum(abs(res$pct_diff) > 5, na.rm = TRUE) > 5)
})
