# Tests for the synth-aggregated duration g-computation (#73). Under
# method = "joint" + duration.method = "joint_lm", build_netstats overrides
# the within-ARTnet stratum medians from build_netparams with synth-aggregated
# medians: predict joint_lm log-duration per synth ego, marginalize over
# partner-race uncertainty via joint_nm_race_model, take median per stratum.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

build_setup <- function(race = TRUE, dur_method = "joint_lm",
                        netparams_method = "joint", race.prop = NULL,
                        sex.cess = FALSE) {
  set.seed(20260420L)
  age_args <- if (sex.cess) {
    list(age.limits = c(15, 100), age.sexual.cessation = 65)
  } else list()
  epistats <- do.call(build_epistats, c(list(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = race, time.unit = 7
  ), age_args))
  set.seed(20260420L)
  np <- build_netparams(epistats, smooth.main.dur = TRUE,
                        method = netparams_method,
                        duration.method = dur_method)
  set.seed(20260420L)
  young.prop <- if (sex.cess) 0.99 else NULL
  ns <- build_netstats(epistats, np,
                       expect.mort = 0.000478213, network.size = 5000,
                       race.prop = race.prop, young.prop = young.prop,
                       method = if (netparams_method == "joint") "joint" else "existing")
  list(epistats = epistats, netparams = np, netstats = ns)
}


test_that("synth override fires under joint + joint_lm and differs from netparams", {
  skip_without_artnetdata()
  obj <- build_setup()
  netparams_main <- obj$netparams$main$durs.main.byage$mean.dur.adj
  diss_main      <- obj$netstats$main$diss.byage$duration
  # Shapes match
  expect_length(diss_main, length(netparams_main))
  # At least one stratum's value diverges (synth attribute distribution
  # differs from ARTnet's), confirming the override fires.
  expect_true(any(abs(diss_main - netparams_main) > 0.5))
})


test_that("synth override does not fire under duration.method = 'empirical'", {
  skip_without_artnetdata()
  obj <- build_setup(dur_method = "empirical")
  expect_null(obj$netparams$main$joint_duration_model)
  netparams_main <- obj$netparams$main$durs.main.byage$mean.dur.adj
  diss_main      <- obj$netstats$main$diss.byage$duration
  # Without joint_duration_model, dissolution_coefs falls back to the
  # netparams values exactly.
  expect_equal(diss_main, netparams_main)
})


test_that("synth override does not fire under method = 'existing'", {
  skip_without_artnetdata()
  obj <- build_setup(netparams_method = "existing", dur_method = "joint_lm")
  netparams_main <- obj$netparams$main$durs.main.byage$mean.dur.adj
  diss_main      <- obj$netstats$main$diss.byage$duration
  # Even though joint_duration_model exists in netparams, build_netstats
  # under method = "existing" doesn't construct synth predictions and so
  # uses the within-ARTnet aggregation directly.
  expect_equal(diss_main, netparams_main)
})


test_that("synth-aggregated durations diverge under shifted race.prop", {
  skip_without_artnetdata()
  default_run <- build_setup()  # ARTnetData::race.dist Atlanta
  shifted_run <- build_setup(race.prop = c(0.35, 0.25, 0.40))
  d_default <- default_run$netstats$casl$diss.byage$duration
  d_shifted <- shifted_run$netstats$casl$diss.byage$duration
  # The casl joint_lm has stronger race-related effects than main; we
  # expect at least one stratum to diverge by > 1% under the population
  # shift.
  rel_diff <- abs(d_default - d_shifted) / d_default
  expect_true(any(rel_diff > 0.01),
              info = paste("max relative diff =",
                           sprintf("%.4f", max(rel_diff, na.rm = TRUE))))
})


test_that("sex.cess.mod preserves the deterministic post-cessation row", {
  skip_without_artnetdata()
  obj <- build_setup(sex.cess = TRUE)
  d_main <- obj$netstats$main$diss.byage$duration
  d_casl <- obj$netstats$casl$diss.byage$duration
  # Last row should equal 1 (deterministic dissolution after sexual cessation)
  expect_equal(d_main[length(d_main)], 1)
  expect_equal(d_casl[length(d_casl)], 1)
  # Length matches netparams shape (1 nonmatch + N age-grps + 1 dead row)
  expect_equal(length(d_main),
               nrow(obj$netparams$main$durs.main.byage))
})


test_that("dissolution_coefs object is well-formed under override", {
  skip_without_artnetdata()
  obj <- build_setup()
  for (layer in c("main", "casl")) {
    diss <- obj$netstats[[layer]]$diss.byage
    expect_s3_class(diss, "disscoef")
    expect_true(all(is.finite(diss$coef.diss)))
    # No NaN / Inf in d.rate adjustment
    expect_true(all(is.finite(diss$d.rate)))
  }
})


test_that("race = FALSE skips the partner-race marginalization gracefully", {
  skip_without_artnetdata()
  obj <- build_setup(race = FALSE)
  # joint_nm_race_model is not fit when race = FALSE; the helper should
  # treat partner-race probability as 0 and predict only at same.race = 0.
  diss <- obj$netstats$main$diss.byage$duration
  expect_true(all(is.finite(diss)))
  expect_true(all(diss > 0))
})
