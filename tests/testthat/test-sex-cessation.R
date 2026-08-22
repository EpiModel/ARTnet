# Tests for the sexual-cessation configuration of build_netstats(), where
# `age.limits[2]` sits above `age.sexual.cessation` so the population keeps
# aging past the point where it stops forming partnerships.
#
# Two defects motivated these:
#  1. `active.sex` was zeroed on `max(attr_age.grp)` as realized rather than on
#     the age group the breaks declare. Whenever the post-cessation group came
#     out empty, that marked the oldest *sexually active* group inactive and
#     zeroed its degree, silently and without any error.
#  2. `young.prop = 1` empties the post-cessation group by construction, which
#     triggered (1) and then failed much later in `netest()` on a target.stats
#     length mismatch, because an age group with no members carries no `ergm`
#     term.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

have_data <- system.file(package = "ARTnetData") != ""
if (have_data) {
  # Sexual cessation at 65, population to 100: six age groups, the sixth being
  # the post-cessation band.
  set.seed(20260821L)
  .ep_cess <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                             init.hiv.prev = c(0.33, 0.137, 0.084),
                             race = TRUE, time.unit = 7,
                             age.limits = c(15, 100),
                             age.sexual.cessation = 65)
  .np_cess <- build_netparams(.ep_cess, smooth.main.dur = TRUE)

  # The default configuration, where the population limit and the cessation age
  # coincide and there is no post-cessation band at all.
  set.seed(20260821L)
  .ep_plain <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
                              init.hiv.prev = c(0.33, 0.137, 0.084),
                              race = TRUE, time.unit = 7)
  .np_plain <- build_netparams(.ep_plain, smooth.main.dur = TRUE)
}

mk_pop <- function(n, age_lo, age_hi) {
  data.frame(
    age        = stats::runif(n, age_lo, age_hi),
    race       = sample(1:3, n, replace = TRUE, prob = c(0.5, 0.05, 0.45)),
    deg.casl   = sample(0:3, n, replace = TRUE),
    deg.main   = sample(0:2, n, replace = TRUE),
    role.class = sample(0:2, n, replace = TRUE),
    risk.grp   = sample(1:5, n, replace = TRUE)
  )
}


# ---- The configuration itself ------------------------------------------------

test_that("sexual cessation adds a post-cessation age group and real mortality", {
  skip_without_artnetdata()
  expect_true(.ep_cess$sex.cess.mod)
  expect_equal(.ep_cess$age.breaks, c(15, 25, 35, 45, 55, 65, 100))
  expect_equal(.ep_cess$age.grps, 6)

  set.seed(101L)
  ns <- build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                       network.size = 3000)
  asmr <- ns$demog$asmr
  rates <- as.matrix(asmr[asmr$age >= 65 & asmr$age <= 99, -1])
  # Under the default 15-to-65 configuration mortality is forced to 1 at 65.
  # Here the cut has to move to the upper age limit instead.
  expect_true(all(rates > 0 & rates < 1))
  expect_true(all(asmr[asmr$age == 100, -1] == 1))
  expect_equal(min(asmr$age[apply(asmr[, -1] >= 1, 1, any)]), 100)
})


# ---- active.sex keys on the declared band, not the realized maximum ----------

test_that("active.sex is 0 exactly for nodes at or above the cessation age", {
  skip_without_artnetdata()
  set.seed(102L)
  ns <- build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                       network.size = 3000)
  expect_gt(sum(ns$attr$age >= 65), 0)
  expect_equal(ns$attr$active.sex == 0, ns$attr$age >= 65)
  expect_equal(ns$attr$active.sex == 0, ns$attr$age.grp == 6)
  # Inactive nodes cannot carry degree on any layer.
  inactive <- ns$attr$active.sex == 0
  expect_true(all(ns$attr$deg.main[inactive] == 0))
  expect_true(all(ns$attr$deg.casl[inactive] == 0))
  expect_true(all(ns$attr$deg.tot[inactive] == 0))
})

test_that("an empty post-cessation band does not retire the 55-64 group", {
  skip_without_artnetdata()
  # This is the regression test for defect (1). Every member of this target
  # population is sexually active, so the post-cessation group is empty and the
  # realized maximum age group is 5. Reading the band off the data marked all
  # 55-to-64 year olds inactive and zeroed their degree.
  set.seed(103L)
  df <- mk_pop(2000, 20, 64)
  expect_warning(
    ns <- build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                         target_pop = df),
    "post-cessation age group"
  )
  expect_equal(sum(ns$attr$age.grp == 6), 0)
  expect_gt(sum(ns$attr$age.grp == 5), 0)
  expect_true(all(ns$attr$active.sex == 1))
  expect_true(any(ns$attr$deg.main > 0))
})

test_that("a populated band is unaffected by the fix", {
  skip_without_artnetdata()
  set.seed(104L)
  df <- mk_pop(2000, 20, 90)
  ns <- build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                       target_pop = df)
  expect_gt(sum(ns$attr$age.grp == 6), 0)
  expect_equal(ns$attr$active.sex == 0, ns$attr$age >= 65)
})


# ---- young.prop endpoints -----------------------------------------------------

test_that("young.prop = 1 is rejected with an actionable message", {
  skip_without_artnetdata()
  expect_error(
    build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                   network.size = 3000, young.prop = 1),
    "must be less than 1"
  )
  expect_error(
    build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                   network.size = 3000, young.prop = 1),
    "0.999"
  )
})

test_that("young.prop = 0 is rejected", {
  skip_without_artnetdata()
  expect_error(
    build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                   network.size = 3000, young.prop = 0),
    "must be greater than 0"
  )
})

test_that("young.prop just below 1 seeds the band and nothing else changes", {
  skip_without_artnetdata()
  set.seed(105L)
  ns <- build_netstats(.ep_cess, .np_cess, expect.mort = 0.000478213,
                       network.size = 25000, young.prop = 0.999)
  n_retired <- sum(ns$attr$age.grp == 6)
  expect_gt(n_retired, 0)
  expect_lt(n_retired / ns$demog$num, 0.01)
  expect_equal(ns$attr$active.sex == 0, ns$attr$age >= 65)
  # The band's target statistics are zero on every layer, which is what makes
  # `ergm` pin those terms off and keep retired nodes tie-free.
  for (layer in c("main", "casl", "inst")) {
    nf <- ns[[layer]]$nodefactor_age.grp
    expect_length(nf, 6)
    expect_equal(unname(nf[6]), 0)
    expect_gt(sum(nf[1:5]), 0)
  }
})


# ---- The default configuration is untouched ----------------------------------

test_that("young.prop is inert without a sexual-cessation band", {
  skip_without_artnetdata()
  # With `age.limits[2]` equal to the cessation age there is no split to
  # reweight, so no value of young.prop should be rejected and none should
  # change the result.
  set.seed(106L)
  a <- build_netstats(.ep_plain, .np_plain, expect.mort = 0.000478213,
                      network.size = 3000, young.prop = 1)
  set.seed(106L)
  b <- build_netstats(.ep_plain, .np_plain, expect.mort = 0.000478213,
                      network.size = 3000)
  expect_equal(a$attr, b$attr)
  expect_true(all(a$attr$active.sex == 1))
  expect_equal(length(.ep_plain$age.breaks) - 1L, 5L)
})
