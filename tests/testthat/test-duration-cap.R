# Tests for build_netstats(duration.cap) and the warning that a by-age
# duration edited in netparams is ignored under duration.method = "joint_lm".

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

cap_env <- new.env(parent = emptyenv())

get_inputs <- function() {
  if (is.null(cap_env$epi)) {
    set.seed(20260419L)
    cap_env$epi <- build_epistats(
      geog.lvl = "city", geog.cat = "Atlanta",
      init.hiv.prev = c(0.33, 0.137, 0.084),
      race = TRUE, time.unit = 7
    )
    cap_env$np_e <- build_netparams(cap_env$epi, smooth.main.dur = TRUE, seed = 1L)
    cap_env$np_j <- build_netparams(cap_env$epi, smooth.main.dur = TRUE, seed = 1L,
                                    method = "joint", duration.method = "joint_lm")
  }
  cap_env
}

build <- function(np, method, ...) {
  build_netstats(get_inputs()$epi, np, network.size = 2000, expect.mort = 0.00025,
                 method = method, seed = 1L, check.feasibility = FALSE, ...)
}

test_that("duration.cap caps the by-age and homogeneous durations under both methods", {
  skip_without_artnetdata()
  inp <- get_inputs()
  cap <- 100
  for (method in c("existing", "joint")) {
    np <- if (method == "joint") inp$np_j else inp$np_e
    ns0 <- suppressWarnings(build(np, method))
    expect_true(any(ns0$main$diss.byage$duration > cap))   # the cap binds
    expect_message(ns1 <- suppressWarnings(build(np, method, duration.cap = cap)),
                   "duration.cap: capping")
    for (layer in c("main", "casl")) {
      for (diss in c("diss.byage", "diss.homog")) {
        expect_true(all(ns1[[layer]][[diss]]$duration <= cap + 1e-9),
                    info = paste(method, layer, diss))
        expect_equal(ns1[[layer]][[diss]]$duration,
                     pmin(ns0[[layer]][[diss]]$duration, cap),
                     info = paste(method, layer, diss))
      }
    }
    # a cap that does not bind changes nothing
    ns2 <- suppressWarnings(build(np, method, duration.cap = 1e6))
    for (layer in c("main", "casl")) {
      expect_equal(ns2[[layer]]$diss.byage$duration, ns0[[layer]]$diss.byage$duration)
      expect_equal(ns2[[layer]]$diss.byage$coef.adj, ns0[[layer]]$diss.byage$coef.adj)
      expect_equal(ns2[[layer]]$diss.homog$duration, ns0[[layer]]$diss.homog$duration)
    }
  }
})

test_that("duration.cap is validated", {
  skip_without_artnetdata()
  inp <- get_inputs()
  expect_error(build(inp$np_e, "existing", duration.cap = 0), "duration.cap")
  expect_error(build(inp$np_e, "existing", duration.cap = "a"), "duration.cap")
  expect_error(build(inp$np_e, "existing", duration.cap = c(10, 20)), "duration.cap")
})

test_that("build_netparams records the joint_lm by-age durations it produced", {
  skip_without_artnetdata()
  inp <- get_inputs()
  for (layer in c("main", "casl")) {
    df <- inp$np_j[[layer]][[paste0("durs.", layer, ".byage")]]
    expect_equal(attr(df, "mean.dur.adj.built"), df$mean.dur.adj)
    expect_null(attr(inp$np_e[[layer]][[paste0("durs.", layer, ".byage")]],
                     "mean.dur.adj.built"))
  }
})

test_that("an edited by-age duration is honored under empirical and warned about under joint_lm", {
  skip_without_artnetdata()
  inp <- get_inputs()
  # legacy: the edit flows through to the dissolution coefficients
  np_e2 <- inp$np_e
  np_e2$main$durs.main.byage$mean.dur.adj[1] <- 50
  expect_no_warning(ns_e <- build(np_e2, "existing"))
  expect_equal(ns_e$main$diss.byage$duration[1], 50)
  # joint_lm: the synth aggregation replaces the column, so the edit is lost;
  # build_netstats says so
  ns_j0 <- expect_no_warning(build(inp$np_j, "joint"))
  np_j2 <- inp$np_j
  np_j2$main$durs.main.byage$mean.dur.adj[1] <- 50
  expect_warning(ns_j <- build(np_j2, "joint"),
                 "durs.main.byage\\$mean.dur.adj differs from what build_netparams")
  expect_equal(ns_j$main$diss.byage$duration, ns_j0$main$diss.byage$duration)
  expect_false(isTRUE(all.equal(ns_j$main$diss.byage$duration[1], 50)))
})
