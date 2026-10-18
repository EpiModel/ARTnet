# Saved epistats and netstats objects should not carry the ARTnet data.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

test_that("trim_epistats strips hiv.mod and keeps its predictions", {
  skip_without_artnetdata()
  epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                             race = TRUE, time.unit = 7)
  trimmed <- trim_epistats(epistats)
  nd <- epistats$hiv.mod$data
  expect_null(trimmed$hiv.mod$data)
  expect_null(trimmed$hiv.mod$model)
  expect_null(attr(trimmed$hiv.mod$terms, ".Environment"))
  expect_equal(predict(trimmed$hiv.mod, newdata = nd, type = "response"),
               predict(epistats$hiv.mod, newdata = nd, type = "response"))
})

test_that("trim_epistats works when there is no hiv.mod", {
  skip_without_artnetdata()
  epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                             init.hiv.prev = c(0.33, 0.137, 0.084),
                             race = TRUE, time.unit = 7)
  trimmed <- trim_epistats(epistats)
  expect_null(trimmed$hiv.mod)
  expect_null(trimmed$acts.mod$data)
})

test_that("netstats dissolution formulas do not keep the build environment", {
  skip_without_artnetdata()
  epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                             race = TRUE, time.unit = 7)
  netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
  netstats <- build_netstats(epistats, netparams, network.size = 1000)
  for (layer in c("main", "casl")) {
    for (diss in c("diss.homog", "diss.byage")) {
      expect_identical(environment(netstats[[layer]][[diss]]$dissolution),
                       globalenv())
    }
  }
})

test_that("the example objects carry no captured environments", {
  netstats <- readRDS(system.file("netstats-example.rds", package = "ARTnet"))
  for (layer in c("main", "casl")) {
    for (diss in c("diss.homog", "diss.byage")) {
      expect_identical(environment(netstats[[layer]][[diss]]$dissolution),
                       globalenv())
    }
  }
  epistats <- readRDS(system.file("epistats-example.rds", package = "ARTnet"))
  for (m in c("acts.mod", "cond.mc.mod", "cond.oo.mod")) {
    expect_null(epistats[[m]]$data)
    expect_null(attr(epistats[[m]]$terms, ".Environment"))
  }
})
