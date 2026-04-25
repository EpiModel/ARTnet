# Tests for the init.hiv.prev length-validation contract (#59).
# Length 1 is fine when race = FALSE (the user's request).
# Length matching race.level is required when race = TRUE.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

test_that("length-1 init.hiv.prev works when race = FALSE (#59)", {
  skip_without_artnetdata()
  expect_silent(
    ep <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                         init.hiv.prev = 0.33, race = FALSE, time.unit = 7)
  )
  expect_equal(ep$init.hiv.prev, 0.33)
  expect_false(ep$race)
})

test_that("length-3 init.hiv.prev still works when race = TRUE", {
  skip_without_artnetdata()
  expect_silent(
    ep <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                         init.hiv.prev = c(0.33, 0.137, 0.084),
                         race = TRUE, time.unit = 7)
  )
  expect_equal(ep$init.hiv.prev, c(0.33, 0.137, 0.084))
})

test_that("length-1 init.hiv.prev with race = TRUE raises a clear error", {
  skip_without_artnetdata()
  expect_error(
    build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                   init.hiv.prev = 0.33, race = TRUE, time.unit = 7),
    regexp = "init.hiv.prev must have length 3"
  )
})

test_that("length-2 init.hiv.prev with race = TRUE raises a clear error", {
  skip_without_artnetdata()
  expect_error(
    build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                   init.hiv.prev = c(0.33, 0.1),
                   race = TRUE, time.unit = 7),
    regexp = "init.hiv.prev must have length 3"
  )
})

test_that("out-of-range init.hiv.prev still rejected with original message", {
  skip_without_artnetdata()
  expect_error(
    build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                   init.hiv.prev = c(0.33, 1.5, 0.084),
                   race = TRUE, time.unit = 7),
    regexp = "between 0 and 1 non-inclusive"
  )
})
