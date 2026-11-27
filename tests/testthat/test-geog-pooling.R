# Tests for build_netparams(geog.pooling = "hierarchical").

skip_without_pooling_deps <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "", "ARTnetData not installed")
  testthat::skip_if_not_installed("glmmTMB")
}

pool_cache <- new.env(parent = emptyenv())

get_pooled <- function() {
  if (!is.null(pool_cache$np)) return(pool_cache)
  ep <- build_epistats(geog.lvl = "city", geog.cat = "Detroit",
                       init.hiv.prev = c(0.33, 0.137, 0.084),
                       race = TRUE, time.unit = 7)
  np_none <- build_netparams(ep, smooth.main.dur = TRUE, method = "joint",
                             duration.method = "joint_lm", seed = 1)
  np <- build_netparams(ep, smooth.main.dur = TRUE, method = "joint",
                        duration.method = "joint_lm", seed = 1,
                        geog.pooling = "hierarchical")
  pool_cache$ep <- ep
  pool_cache$np_none <- np_none
  pool_cache$np <- np
  pool_cache
}

test_that("geog.pooling validates its requirements", {
  skip_without_pooling_deps()
  x <- get_pooled()
  expect_error(build_netparams(x$ep, method = "existing", geog.pooling = "hierarchical"),
               "requires method = \"joint\"")
  expect_error(build_netparams(x$ep, method = "joint", duration.method = "empirical",
                               geog.pooling = "hierarchical"),
               "requires duration.method = \"joint_lm\"")
  ep_nat <- build_epistats(init.hiv.prev = c(0.33, 0.137, 0.084), race = TRUE, time.unit = 7)
  expect_error(build_netparams(ep_nat, method = "joint", duration.method = "joint_lm",
                               geog.pooling = "hierarchical"),
               "requires epistats built with a geog.lvl")
  expect_error(build_netparams(x$ep, method = "joint", geog.pooling = "bayes"))
})

test_that("default geog.pooling leaves netparams unchanged", {
  skip_without_pooling_deps()
  x <- get_pooled()
  np_explicit <- build_netparams(x$ep, smooth.main.dur = TRUE, method = "joint",
                                 duration.method = "joint_lm", seed = 1,
                                 geog.pooling = "none")
  expect_null(x$np_none$geog.pooling)
  expect_equal(np_explicit, x$np_none, ignore_formula_env = TRUE,
               ignore_function_env = TRUE)
})

test_that("hierarchical pooling replaces every geogYN coefficient with the pooled effect", {
  skip_without_pooling_deps()
  x <- get_pooled()
  gp <- x$np$geog.pooling
  expect_identical(gp$method, "hierarchical")
  expect_identical(gp$groups, c("REGCODE", "DIVCODE", "city2"))
  eff <- gp$effects
  expect_false(anyNA(eff$pooled))
  expect_true(all(eff$n_geog > 0))

  glm_rows <- eff[eff$term == "geogYN", ]
  for (i in seq_len(nrow(glm_rows))) {
    o <- x$np[[glm_rows$layer[i]]][[glm_rows$model[i]]]
    expect_equal(unname(coef(o)["geogYN"]), glm_rows$pooled[i])
    expect_false(is.null(attr(o, "geog_pooling")))
  }
  # Unpooled values match the default fit where the formula is unchanged.
  dur <- glm_rows[glm_rows$model == "joint_duration_model", ]
  for (i in seq_len(nrow(dur))) {
    expect_equal(dur$unpooled[i],
                 unname(coef(x$np_none[[dur$layer[i]]]$joint_duration_model)["geogYN"]))
  }

  pr <- eff[eff$model == "joint_prace_model", ]
  for (layer in unique(pr$layer)) {
    cf <- coef(x$np[[layer]]$joint_prace_model)
    sub <- pr[pr$layer == layer, ]
    expect_equal(unname(cf[sub("geogYN:", "", sub$term), "geogYN"]), sub$pooled)
  }
})

test_that("pooled multinomial predictions match the softmax of its coefficients", {
  skip_without_pooling_deps()
  x <- get_pooled()
  m <- x$np$casl$joint_prace_model
  nd <- x$np$casl$joint_nm_age_model$data
  nd <- nd[complete.cases(nd[, c("age.grp", "hiv2", "race.cat.num", "geogYN")]), ]
  X <- model.matrix(delete.response(m$terms), nd)
  eta <- cbind(0, X %*% t(coef(m)))
  p <- exp(eta) / rowSums(exp(eta))
  expect_equal(unname(p), unname(predict(m, nd, type = "probs")), tolerance = 1e-10)
})

test_that("main-layer city effects pool to the national value", {
  skip_without_pooling_deps()
  x <- get_pooled()
  eff <- x$np$geog.pooling$effects
  main_deg <- eff[eff$layer == "main" & eff$model == "joint_model", ]
  expect_lt(abs(main_deg$pooled), 0.01)
  oo <- eff[eff$layer == "inst" & eff$model == "joint_model", ]
  expect_lt(abs(oo$pooled), abs(oo$unpooled))
})

test_that("raw geography inputs are pooled", {
  skip_without_pooling_deps()
  x <- get_pooled()
  np <- x$np
  expect_length(np$main$deg.casl.dist, 4)
  expect_length(np$casl$deg.main.dist, 3)
  expect_length(np$inst$deg.tot.dist, 4)
  for (p in list(np$main$deg.casl.dist, np$casl$deg.main.dist, np$inst$deg.tot.dist)) {
    expect_equal(sum(p), 1)
  }
  w <- np$geog.pooling$dist_weights
  expect_true(all(w >= 0 & w < 1))
  r <- np$geog.pooling$risk.grp_ratio
  expect_equal(np$inst$nf.risk.grp, x$np_none$inst$nf.risk.grp * r[["pooled"]] / r[["raw"]])
})

test_that("build_netstats uses pooled netparams under joint and refuses existing", {
  skip_without_pooling_deps()
  x <- get_pooled()
  ns <- build_netstats(x$ep, x$np, expect.mort = 0.000478213, network.size = 5000,
                       method = "joint", seed = 1)
  for (layer in c("main", "casl", "inst")) {
    expect_equal(sum(ns[[layer]]$nodefactor_race), 2 * ns[[layer]]$edges)
  }
  ns_none <- build_netstats(x$ep, x$np_none, expect.mort = 0.000478213, network.size = 5000,
                            method = "joint", seed = 1)
  expect_gt(ns$inst$edges, ns_none$inst$edges)
  expect_error(build_netstats(x$ep, x$np, expect.mort = 0.000478213, network.size = 5000,
                              method = "existing"),
               "use method = \"joint\"")
})

test_that(".dm_shrink_dist returns the national distribution without between-group variation", {
  x <- rep(0:2, times = 40)
  grp <- rep(c("a", "b", "c", "d"), each = 30)
  inc <- grp == "a"
  res <- ARTnet:::.dm_shrink_dist(x, grp, inc, 0:2)
  expect_equal(res$weight, 0)
  expect_equal(res$p, rep(1 / 3, 3))

  x2 <- c(rep(0, 30), rep(1, 30), rep(2, 30), rep(0:2, 10))
  res2 <- ARTnet:::.dm_shrink_dist(x2, grp, inc, 0:2)
  expect_gt(res2$weight, 0)
  expect_lt(res2$weight, 1)
  expect_equal(sum(res2$p), 1)
})
