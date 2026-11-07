# Race-mixing helpers for build_netstats(): the balanced partner-race mixing
# matrix under method = "joint" (issue #9), the opt-in random-mixing floor
# (#9), and the feasibility check on the race targets (#10).
#
# Race levels are integer-coded 1..K in the order of epistats$race.level
# (default 1 = Black, 2 = Hispanic, 3 = White/Other).


# Predicted partner-race probabilities from an nnet::multinom fit, always
# returned as an n x K matrix with columns in race-level order. nnet returns
# a bare vector of P(level 2) when there are two levels, and a bare vector of
# length K when newdata has a single row.
.predict_prace <- function(model, newdata, n_levels) {
  p <- predict(model, newdata = newdata, type = "probs")
  if (is.null(dim(p))) {
    p <- if (n_levels == 2) cbind(1 - p, p) else matrix(p, nrow = 1)
  }
  p <- unname(as.matrix(p))
  if (ncol(p) != n_levels || nrow(p) != nrow(newdata)) {
    stop("Partner-race prediction has dimension ", nrow(p), " x ", ncol(p),
         "; expected ", nrow(newdata), " x ", n_levels, ".", call. = FALSE)
  }
  p
}


# Reported endpoint matrix at the synthetic population. A[r, s] is the
# expected number of partnership endpoints held by egos of race r whose
# partner is of race s: sum over egos of race r of predicted degree times
# the predicted probability that a partner is of race s. Rows are the ego
# side and sum to nodefactor_race by construction. In a closed population
# A would be symmetric; egocentric reports standardized on the ego side
# only are not.
.reported_race_mixing <- function(pred_deg, probs, race, n_levels) {
  A <- matrix(0, n_levels, n_levels)
  for (r in seq_len(n_levels)) {
    sel <- race == r
    if (any(sel)) {
      A[r, ] <- colSums(pred_deg[sel] * probs[sel, , drop = FALSE])
    }
  }
  A
}


# Balance a reported endpoint matrix so that it describes a closed
# population, holding each race's total endpoints (nodefactor_race) fixed.
#
# With S = (A + t(A)) / 2, find the diagonal scaling d such that
# B = diag(d) S diag(d) has rowSums(B) == rowSums(A). B is symmetric, so it
# is the endpoint matrix of a real network: B[r, r] / 2 same-race edges and
# B[r, s] cross-race edges. The scaling keeps the homophily odds ratios of
# S (B[r, r] B[s, s] / B[r, s]^2 == S[r, r] S[s, s] / S[r, s]^2) and refits
# the margins to the target population, so it is the log-linear transport
# of ARTnet's selective mixing to a population with a different race
# composition. Plain averaging of the two reports of each cross-race cell
# is not used because, with the margins held, it can drive a minority
# group's same-race cell negative when the majority's reports of that
# group exceed what the group's own degree can supply.
#
# Solved by the symmetric Sinkhorn iteration d <- d * sqrt(m / (S d)).
# A race with no endpoints (m == 0) gets d == 0 and an all-zero row and
# column. Returns B with the scaling factors and the iteration count as
# attributes.
.balance_race_mixing <- function(A, tol = 1e-10, maxit = 5000L) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("A must be a square matrix.", call. = FALSE)
  }
  if (any(A < 0)) stop("A must be non-negative.", call. = FALSE)
  m <- rowSums(A)
  S <- (A + t(A)) / 2
  active <- m > 0
  d <- as.numeric(active)
  it <- 0L
  converged <- !any(active)
  while (!converged && it < maxit) {
    it <- it + 1L
    Sd <- as.vector(S %*% d)
    if (any(Sd[active] <= 0)) {
      stop("Race mixing matrix cannot be balanced: a race with partnership ",
           "endpoints has no partners of any race present in the population.",
           call. = FALSE)
    }
    # rowSums(B) = d * (S d); ratio is the target margin over the current one
    ratio <- rep(1, length(m))
    ratio[active] <- m[active] / (d[active] * Sd[active])
    d[active] <- d[active] * sqrt(ratio[active])
    converged <- max(abs(ratio - 1)) < tol
  }
  if (!converged) {
    stop("Race mixing balance did not converge in ", maxit, " iterations.",
         call. = FALSE)
  }
  B <- outer(d, d) * S
  dimnames(B) <- dimnames(A)
  attr(B, "scaling") <- d
  attr(B, "iterations") <- it
  B
}


# Full joint-method race-mixing step for one layer. Predicts partner race
# for every synthetic node, forms the reported endpoint matrix, balances it,
# and returns the nodematch targets plus the two matrices that
# build_netstats stores as additive fields under method = "joint".
.joint_race_mixing <- function(prace_model, synth, pred_deg, race, n_levels,
                               race_labels = NULL, layer = "") {
  if (is.null(prace_model)) {
    stop("method = 'joint' needs the partner-race model ",
         "netparams$", layer, "$joint_prace_model, which build_netparams() ",
         "fits under method = 'joint' from ARTnet 2.10.0. Re-run ",
         "build_netparams(..., method = 'joint').", call. = FALSE)
  }
  probs <- .predict_prace(prace_model, synth, n_levels)
  A <- .reported_race_mixing(pred_deg, probs, race, n_levels)
  if (!is.null(race_labels) && length(race_labels) == n_levels) {
    dimnames(A) <- list(race_labels, race_labels)
  }
  B <- .balance_race_mixing(A)
  nodematch_race <- diag(B) / 2
  mixmat <- B
  diag(mixmat) <- diag(B) / 2
  attr(mixmat, "scaling") <- NULL
  attr(mixmat, "iterations") <- NULL
  list(
    nodematch_race        = unname(nodematch_race),
    nodematch_race_diffF  = sum(nodematch_race),
    mixmat_race           = mixmat,
    mixmat_race_reported  = A
  )
}


# Opt-in random-mixing floor on nodematch_race_diffF for one layer (#9).
# With endpoint shares p = nodefactor_race / sum(nodefactor_race), random
# mixing gives edges * sum(p^2) same-race edges, which every network with
# those nodefactor targets can realize (the feasibility condition
# nodefactor[r] <= edges + nodematch at that value reduces to
# (1 - p[r])^2 >= 0). The floor raises the target to that value when it
# falls below and says so. Only nodematch_race_diffF is touched, the
# target that the EpiModelHIV template formulas consume.
.apply_race_floor <- function(layer_stats, layer_name) {
  E  <- layer_stats$edges
  nf <- layer_stats$nodefactor_race
  M  <- layer_stats$nodematch_race_diffF
  if (is.null(E) || is.null(nf) || is.null(M)) return(layer_stats)
  p <- as.numeric(nf) / sum(nf)
  floor_val <- E * sum(p^2)
  if (M < floor_val) {
    message(sprintf(
      paste0("nodematch.race.floor: raising the %s same-race target from ",
             "%.1f to %.1f edges (share %.3f -> %.3f), the random-mixing floor."),
      layer_name, M, floor_val, M / E, floor_val / E))
    layer_stats$nodematch_race_diffF <- floor_val
  }
  layer_stats
}


#' Check Feasibility of Race-Mixing Target Statistics
#'
#' @description Checks, for each network layer, whether the race target
#'   statistics in a `netstats` object can be realized by any network.
#'
#' @param netstats Output from [`build_netstats`].
#' @param tol Numeric slack above zero that is still treated as feasible,
#'   to absorb floating-point noise.
#' @param labels Optional character vector of race labels in race-code order,
#'   for the `label` column. When `NULL`, the labels are taken from the
#'   dimnames of `mixmat_race` when present (`method = "joint"`), and are the
#'   integer codes otherwise.
#'
#' @details
#' With a single same-race count `M` (the `nodematch("race", diff = FALSE)`
#' target, `nodematch_race_diffF`) and `edges` `E`, every race `r` needs
#' `nodefactor_race[r] <= E + M`: a race cannot hold more partnership
#' endpoints than the same-race edges can contribute twice plus the
#' cross-race edges can contribute once. The function reports
#' `slack = nodefactor_race - edges - nodematch_race_diffF` for every layer
#' and race; a positive slack means no network satisfies the targets.
#'
#' Infeasible targets fail in `ergm` MCMLE estimation with "The MPLE does
#' not exist!" and then the simulated-density guard, usually after a
#' warning that `nodefactor.race.<r>` is linearly dependent on
#' `edges + nodematch.race`, which is this same boundary seen from the
#' other side. Stochastic-Approximation estimation settles on the nearest
#' reachable network without complaint, so the problem can go unnoticed
#' until a model is refit with MCMLE.
#'
#' `build_netstats()` runs this check by default (`check.feasibility = TRUE`)
#' and stops on a positive slack. The condition covers the
#' `nodematch("race", diff = FALSE)` specification used by the EpiModelHIV
#' template formulas; the per-race `nodematch_race` vector is not checked.
#'
#' @return A data frame with one row per layer and race: `layer`, `race`
#'   (integer code), `label`, `nodefactor`, `edges`, `nodematch_race_diffF`,
#'   `slack`, and `feasible`. Zero rows when `netstats` has no race targets.
#'
#' @export
#'
#' @examples
#' netstats <- readRDS(system.file("netstats-example.rds", package = "ARTnet"))
#' check_race_feasibility(netstats)
#'
check_race_feasibility <- function(netstats, tol = 1e-6, labels = NULL) {
  rows <- list()
  for (l in c("main", "casl", "inst")) {
    x <- netstats[[l]]
    if (is.null(x$edges) || is.null(x$nodefactor_race) ||
        is.null(x$nodematch_race_diffF)) {
      next
    }
    nf <- as.numeric(x$nodefactor_race)
    # demog$props columns follow ARTnetData::race.dist, not the race codes,
    # so they are not used for labels.
    lab <- labels
    if (is.null(lab) && !is.null(x$mixmat_race)) lab <- rownames(x$mixmat_race)
    if (is.null(lab) || length(lab) != length(nf)) lab <- as.character(seq_along(nf))
    slack <- nf - x$edges - x$nodematch_race_diffF
    rows[[l]] <- data.frame(
      layer = l,
      race = seq_along(nf),
      label = lab,
      nodefactor = nf,
      edges = x$edges,
      nodematch_race_diffF = x$nodematch_race_diffF,
      slack = slack,
      feasible = slack <= tol,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    message("No race targets in netstats (race = FALSE); nothing to check.")
    return(data.frame(layer = character(0), race = integer(0),
                      label = character(0), nodefactor = numeric(0),
                      edges = numeric(0), nodematch_race_diffF = numeric(0),
                      slack = numeric(0), feasible = logical(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# Stop with a layer-by-race account of the excess when the race targets
# are infeasible. Called at the end of build_netstats().
.stop_if_race_infeasible <- function(netstats, labels = NULL) {
  feas <- check_race_feasibility(netstats, labels = labels)
  bad <- feas[!feas$feasible, , drop = FALSE]
  if (nrow(bad) == 0) return(invisible(TRUE))
  lines <- sprintf(
    paste0("  %s: nodefactor_race[%d] (%s) exceeds edges + nodematch_race_diffF ",
           "by %.1f endpoints (%.1f%% of nodefactor)"),
    bad$layer, bad$race, bad$label, bad$slack, 100 * bad$slack / bad$nodefactor)
  stop("Race-mixing targets are infeasible: no network can carry these ",
       "nodefactor_race and nodematch_race_diffF values together, and MCMLE ",
       "estimation will fail on them.\n",
       paste(lines, collapse = "\n"), "\n",
       "Options: method = 'joint' (partner-race mixing balanced to the target ",
       "population), nodematch.race.floor = TRUE (raise the same-race target ",
       "to random mixing), or check.feasibility = FALSE to return the targets ",
       "anyway. See ?check_race_feasibility.", call. = FALSE)
}
