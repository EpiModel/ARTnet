# Tests for the race-mixing helpers (issues #9, #10, #11): the balanced
# partner-race mixing matrix under method = "joint", the opt-in random-mixing
# floor, the feasibility check, and the build_netparams seed.

skip_without_artnetdata <- function() {
  testthat::skip_if(system.file(package = "ARTnetData") == "",
                    "ARTnetData not installed")
}

# Seattle main-layer reported endpoint matrix from the planning prototype
# (joint, 10,000 nodes, default race.prop). Rows are the ego race (Black,
# Hispanic, White/Other), columns the partner race. White egos report more
# Black partners (231) than Black egos have endpoints in total (174).
seattle_main <- matrix(c(91.7,   7.3,   74.5,
                         17.4,  75.7,  177.4,
                         231.1, 277.4, 3181.7),
                       nrow = 3, byrow = TRUE)

slack_of <- function(B) {
  nf <- rowSums(B)
  nf - sum(nf) / 2 - sum(diag(B)) / 2
}

# ---- .balance_race_mixing ---------------------------------------------------

test_that("a symmetric matrix is a fixed point of the balance", {
  A <- matrix(c(10, 2, 1,
                2, 20, 3,
                1, 3, 30), nrow = 3)
  B <- .balance_race_mixing(A)
  expect_equal(unclass(B)[, ], A, tolerance = 1e-9)
  expect_equal(attr(B, "scaling"), rep(1, 3), tolerance = 1e-9)
})

test_that("the balance holds the margins, is symmetric, and is feasible on Seattle", {
  A <- seattle_main
  expect_gt(max(abs(A - t(A))), 100)  # the reports do not balance
  B <- .balance_race_mixing(A)
  expect_equal(rowSums(B), rowSums(A), tolerance = 1e-8)
  expect_equal(unclass(B)[, ], t(unclass(B)[, ]), tolerance = 1e-9)
  expect_true(all(B >= 0))
  expect_true(all(diag(B) > 0))
  expect_true(all(slack_of(B) < 0))
  # the homophily odds ratios of the symmetrized input are preserved
  S <- (A + t(A)) / 2
  or <- function(M, r, s) M[r, r] * M[s, s] / M[r, s]^2
  for (r in 1:2) for (s in (r + 1):3) {
    expect_equal(or(B, r, s), or(S, r, s), tolerance = 1e-8)
  }
})

test_that("plain averaging with the margins held goes negative where the balance does not", {
  # Seattle casual from the prototype: averaging drives the Hispanic
  # diagonal negative; the symmetric scaling keeps every cell positive.
  A <- matrix(c(191.3,  51.5,  175.3,
                46.6,  181.6,  249.7,
                443.8, 961.6, 4131.8), nrow = 3, byrow = TRUE)
  expect_gt(max(slack_of(A)), 0)  # the unbalanced reports are infeasible as targets
  S <- (A + t(A)) / 2
  off <- S
  diag(off) <- 0
  avg_diag <- rowSums(A) - rowSums(off)
  expect_lt(min(avg_diag), 0)
  B <- .balance_race_mixing(A)
  expect_true(all(diag(B) > 0))
  expect_true(all(slack_of(B) < 0))
})

test_that("a race with no endpoints gets an empty row and column", {
  A <- seattle_main
  A[2, ] <- 0
  B <- .balance_race_mixing(A)
  expect_equal(B[2, ], rep(0, 3))
  expect_equal(B[, 2], rep(0, 3))
  expect_equal(rowSums(B), rowSums(A), tolerance = 1e-8)
  expect_equal(attr(B, "scaling")[2], 0)
})

test_that("the balance rejects bad input", {
  expect_error(.balance_race_mixing(matrix(1, 2, 3)), "square")
  expect_error(.balance_race_mixing(matrix(c(1, -1, 0, 1), 2)), "non-negative")
})

# ---- .reported_race_mixing / .predict_prace --------------------------------

test_that("the reported matrix sums predicted degree by ego race and partner race", {
  pred_deg <- c(1, 2, 3)
  probs <- rbind(c(0.5, 0.5), c(0.2, 0.8), c(0.9, 0.1))
  race <- c(1L, 2L, 2L)
  A <- .reported_race_mixing(pred_deg, probs, race, 2)
  expect_equal(A[1, ], c(0.5, 0.5))
  expect_equal(A[2, ], c(2 * 0.2 + 3 * 0.9, 2 * 0.8 + 3 * 0.1))
  expect_equal(rowSums(A), c(1, 5))
})

test_that(".predict_prace returns an n x K matrix for two and three levels", {
  set.seed(1)
  d <- data.frame(x = rnorm(300))
  d$y3 <- factor(sample(1:3, 300, TRUE), levels = 1:3)
  d$y2 <- factor(sample(1:2, 300, TRUE), levels = 1:2)
  m3 <- nnet::multinom(y3 ~ x, data = d, trace = FALSE)
  m2 <- nnet::multinom(y2 ~ x, data = d, trace = FALSE)
  nd <- data.frame(x = c(-1, 0, 1))
  p3 <- .predict_prace(m3, nd, 3)
  p2 <- .predict_prace(m2, nd, 2)
  expect_equal(dim(p3), c(3, 3))
  expect_equal(dim(p2), c(3, 2))
  expect_equal(rowSums(p3), rep(1, 3), tolerance = 1e-9)
  expect_equal(rowSums(p2), rep(1, 3), tolerance = 1e-9)
  # single-row newdata
  expect_equal(dim(.predict_prace(m3, nd[1, , drop = FALSE], 3)), c(1, 3))
  expect_equal(dim(.predict_prace(m2, nd[1, , drop = FALSE], 2)), c(1, 2))
})

# ---- check_race_feasibility / floor ----------------------------------------

fake_netstats <- function(nf, edges, M) {
  list(race = TRUE,
       demog = list(props = data.frame(Black = 0.1, Hispanic = 0.1, White.Other = 0.8)),
       main = list(edges = edges, nodefactor_race = nf, nodematch_race_diffF = M),
       casl = list(edges = edges, nodefactor_race = nf, nodematch_race_diffF = M),
       inst = list(edges = edges, nodefactor_race = nf, nodematch_race_diffF = M))
}

test_that("check_race_feasibility reports slack per layer and race", {
  ns <- fake_netstats(nf = c(100, 100, 800), edges = 500, M = 400)
  res <- check_race_feasibility(ns)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 9)
  expect_named(res, c("layer", "race", "label", "nodefactor", "edges",
                      "nodematch_race_diffF", "slack", "feasible"))
  expect_equal(res$slack[res$layer == "main"], c(100, 100, 800) - 500 - 400)
  expect_true(all(res$feasible))
  # labels: integer codes by default (demog$props columns are not in code
  # order), explicit when supplied, from mixmat_race dimnames when present
  expect_equal(unique(res$label), c("1", "2", "3"))
  lab <- c("Black", "Hispanic", "White.Other")
  expect_equal(unique(check_race_feasibility(ns, labels = lab)$label), lab)
  ns$main$mixmat_race <- matrix(0, 3, 3, dimnames = list(lab, lab))
  res2 <- check_race_feasibility(ns)
  expect_equal(res2$label[res2$layer == "main"], lab)
  expect_equal(res2$label[res2$layer == "casl"], c("1", "2", "3"))
  # M too small for the large race: nodefactor 800 > 500 + 250
  bad <- check_race_feasibility(fake_netstats(c(100, 100, 800), 500, 250))
  expect_false(all(bad$feasible))
  expect_true(all(bad$race[!bad$feasible] == 3))
  expect_equal(bad$slack[bad$race == 3], rep(50, 3))
})

test_that("check_race_feasibility is empty without race targets", {
  expect_message(res <- check_race_feasibility(list(race = FALSE, main = list(edges = 1))),
                 "nothing to check")
  expect_equal(nrow(res), 0)
})

test_that("the random-mixing floor binds only below random mixing", {
  layer <- list(edges = 500, nodefactor_race = c(100, 100, 800),
                nodematch_race_diffF = 250)
  p <- c(0.1, 0.1, 0.8)
  rnd <- 500 * sum(p^2)
  expect_message(out <- .apply_race_floor(layer, "main"), "random-mixing floor")
  expect_equal(out$nodematch_race_diffF, rnd)
  # feasible after flooring
  expect_true(all(c(100, 100, 800) - 500 - rnd <= 0))
  above <- layer
  above$nodematch_race_diffF <- 400
  expect_silent(out2 <- .apply_race_floor(above, "main"))
  expect_equal(out2$nodematch_race_diffF, 400)
  # layers without race targets pass through
  expect_identical(.apply_race_floor(list(edges = 1), "inst"), list(edges = 1))
})

# ---- end to end on Seattle (ARTnetData) -------------------------------------

seattle_env <- new.env(parent = emptyenv())

seattle_epistats <- function() {
  if (is.null(seattle_env$epi)) {
    set.seed(20260420L)
    seattle_env$epi <- build_epistats(
      geog.lvl = "city", geog.cat = "Seattle",
      init.hiv.prev = c(0.33, 0.137, 0.084),
      race = TRUE, time.unit = 7
    )
  }
  seattle_env$epi
}

seattle_netparams <- function(method) {
  key <- paste0("np_", method)
  if (is.null(seattle_env[[key]])) {
    epi <- seattle_epistats()
    dur <- if (method == "joint") "joint_lm" else "empirical"
    seattle_env[[key]] <- build_netparams(epi, smooth.main.dur = TRUE,
                                          method = method, duration.method = dur,
                                          seed = 20260420L)
  }
  seattle_env[[key]]
}

test_that("build_netstats stops on Seattle's infeasible legacy race targets", {
  skip_without_artnetdata()
  epi <- seattle_epistats()
  np <- seattle_netparams("existing")
  expect_error(
    build_netstats(epi, np, network.size = 5000, expect.mort = 0.00025,
                   seed = 20260420L),
    regexp = "infeasible.*nodefactor_race\\[3\\] \\(White.Other\\)"
  )
  ns <- build_netstats(epi, np, network.size = 5000, expect.mort = 0.00025,
                       seed = 20260420L, check.feasibility = FALSE)
  feas <- check_race_feasibility(ns)
  expect_false(all(feas$feasible))
  expect_true(all(feas$feasible[feas$race != 3]))
  expect_true(all(!feas$feasible[feas$race == 3]))
})

test_that("the floor makes Seattle's legacy race targets feasible", {
  skip_without_artnetdata()
  epi <- seattle_epistats()
  np <- seattle_netparams("existing")
  expect_message(
    ns <- build_netstats(epi, np, network.size = 5000, expect.mort = 0.00025,
                         seed = 20260420L, nodematch.race.floor = TRUE),
    "random-mixing floor"
  )
  feas <- check_race_feasibility(ns)
  expect_true(all(feas$feasible))
  for (layer in c("main", "casl", "inst")) {
    p <- ns[[layer]]$nodefactor_race / sum(ns[[layer]]$nodefactor_race)
    expect_equal(ns[[layer]]$nodematch_race_diffF,
                 ns[[layer]]$edges * sum(p^2), tolerance = 1e-9,
                 info = layer)
  }
})

test_that("the joint method balances Seattle's race mixing and passes the check", {
  skip_without_artnetdata()
  epi <- seattle_epistats()
  np <- seattle_netparams("joint")
  ns <- build_netstats(epi, np, network.size = 5000, expect.mort = 0.00025,
                       method = "joint", seed = 20260420L)
  feas <- check_race_feasibility(ns)
  expect_true(all(feas$feasible))
  for (layer in c("main", "casl", "inst")) {
    L <- ns[[layer]]
    expect_s3_class(np[[layer]]$joint_prace_model, "multinom")
    expect_equal(sum(L$nodematch_race), L$nodematch_race_diffF, tolerance = 1e-9)
    expect_true(all(L$nodematch_race <= L$nodefactor_race / 2 + 1e-9))
    expect_true(all(L$nodematch_race > 0))
    # balanced matrix: symmetric edge counts; each race's endpoints
    # (twice the diagonal plus the off-diagonal row) equal nodefactor_race,
    # and the distinct cells sum to edges
    M <- L$mixmat_race
    expect_equal(dim(M), c(3, 3))
    expect_equal(M, t(M), tolerance = 1e-9)
    expect_equal(unname(2 * diag(M) + (rowSums(M) - diag(M))), L$nodefactor_race,
                 tolerance = 1e-8, info = layer)
    expect_equal(sum(M[upper.tri(M, diag = TRUE)]), L$edges, tolerance = 1e-8)
    # reported matrix: rows sum to nodefactor_race, and it is not symmetric
    R <- L$mixmat_race_reported
    expect_equal(unname(rowSums(R)), L$nodefactor_race, tolerance = 1e-8)
    expect_gt(max(abs(R - t(R))), 1)
    expect_equal(rownames(M), c("Black", "Hispanic", "White.Other"))
    # the checker picks the labels up from the matrix
    expect_equal(feas$label[feas$layer == layer], rownames(M))
  }
  # same-race shares are assortative relative to random mixing in every layer
  for (layer in c("main", "casl", "inst")) {
    L <- ns[[layer]]
    p <- L$nodefactor_race / sum(L$nodefactor_race)
    expect_gt(L$nodematch_race_diffF, L$edges * sum(p^2))
  }
})

test_that("build_netparams(seed) makes the partner-race imputation reproducible", {
  skip_without_artnetdata()
  epi <- seattle_epistats()
  a <- build_netparams(epi, smooth.main.dur = TRUE, seed = 1L)
  b <- build_netparams(epi, smooth.main.dur = TRUE, seed = 1L)
  c2 <- build_netparams(epi, smooth.main.dur = TRUE, seed = 2L)
  for (layer in c("main", "casl", "inst")) {
    expect_identical(a[[layer]]$nm.race, b[[layer]]$nm.race)
    expect_identical(a[[layer]]$nm.race_diffF, b[[layer]]$nm.race_diffF)
    expect_identical(a[[layer]]$md.main, c2[[layer]]$md.main)
  }
  moved <- vapply(c("main", "casl", "inst"),
                  function(l) !isTRUE(all.equal(a[[l]]$nm.race, c2[[l]]$nm.race)),
                  logical(1))
  expect_true(any(moved))
})
