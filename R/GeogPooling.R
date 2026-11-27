# Partial pooling of the geography effect for build_netparams(geog.pooling =
# "hierarchical").
#
# Under method = "joint" every joint model carries the target geography as a
# single main-effect indicator (geogYN: in the target geography vs the rest of
# ARTnet), estimated without pooling. For a geography with few respondents
# that coefficient is noisy, and build_netstats() carries the noise straight
# into the targets. The hierarchical option replaces each geogYN coefficient
# with the geography's partially pooled effect from a refit of the same fixed
# effects with crossed random intercepts over the nested geographic levels of
# ARTnet (region, division, then state / county / city), so a small geography
# borrows from its division and region before the nation. The effect is the
# geography's predicted random-effect total minus its mean over the rest of
# ARTnet, which is the contrast geogYN estimates. The model objects keep their
# class, so build_netstats() predicts from them unchanged.
#
# The same option pools the raw geography-specific inputs that bypass the
# models: the cross-layer degree distributions used to build the synthetic
# population (Dirichlet-multinomial shrinkage toward the national
# distribution) and the one-off risk-group rates (rescaled by the pooled
# one-off rate ratio instead of the raw one). For overdispersed count models
# the interaction selection uses quasi-AIC and the pooled refit uses a negative
# binomial, since Poisson standard errors and AIC are far too optimistic for
# partner counts; the point estimates stay Poisson.


# Random-effect grouping levels for each geog.lvl, coarsest first. city2 is not
# strictly nested in division (a few metros span divisions), which is why the
# intercepts are specified as crossed rather than nested.
.geog_pool_groups <- function(geog.lvl) {
  switch(geog.lvl,
    region   = "REGCODE",
    division = c("REGCODE", "DIVCODE"),
    state    = c("REGCODE", "DIVCODE", "State"),
    county   = c("REGCODE", "DIVCODE", "State", "COUNTYFIPS"),
    city     = c("REGCODE", "DIVCODE", "city2"),
    stop("geog.pooling = \"hierarchical\" does not support geog.lvl = \"",
         geog.lvl, "\".", call. = FALSE)
  )
}


# Attach the grouping columns (as .pool_<level> factors) by respondent ID.
.with_geog_groups <- function(dat, groups) {
  w <- ARTnetData::ARTnet.wide
  idx <- match(dat$AMIS_ID, w$AMIS_ID)
  if (anyNA(idx)) {
    stop("Geographic pooling could not match every respondent to ARTnet.wide.",
         call. = FALSE)
  }
  for (g in groups) {
    v <- as.character(w[[g]][idx])
    v[is.na(v)] <- "(missing)"
    dat[[paste0(".pool_", g)]] <- factor(v)
  }
  dat
}


# Fit `fml` (which contains geogYN) without geogYN and with crossed random
# intercepts over `groups`, and return the pooled geography contrast. NULL if
# the geography has no rows in the model data; list(error = msg) if the fit
# fails.
.pooled_geog_contrast <- function(fml, dat, family, groups) {
  dat <- .with_geog_groups(dat, groups)
  vars <- intersect(all.vars(fml), names(dat))
  dat <- dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
  inc <- dat$geogYN == 1
  if (!any(inc) || all(inc)) return(NULL)
  re <- paste0("(1 | .pool_", groups, ")", collapse = " + ")
  f_re <- stats::update(fml, stats::as.formula(paste(". ~ . - geogYN +", re)))
  fit <- tryCatch(
    suppressWarnings(glmmTMB::glmmTMB(
      f_re, data = dat, family = family,
      REML = identical(family$family, "gaussian"),
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  re_list <- glmmTMB::ranef(fit)$cond
  u <- numeric(nrow(dat))
  for (g in names(re_list)) {
    v <- stats::setNames(re_list[[g]][, 1], rownames(re_list[[g]]))
    u <- u + unname(v[as.character(dat[[g]])])
  }
  sds <- vapply(glmmTMB::VarCorr(fit)$cond,
                function(x) unname(attr(x, "stddev")[1]), numeric(1))
  names(sds) <- sub("^\\.pool_", "", names(sds))
  list(effect = mean(u[inc]) - mean(u[!inc]), sd = sds, n_geog = sum(inc))
}


# Overdispersed count models: rerun fit_joint_glm()'s forward AIC selection
# over the same candidate interactions with quasi-AIC, taking the dispersion
# from the Poisson fit with every candidate included. The Poisson likelihood
# understates uncertainty for overdispersed partner counts, so plain AIC keeps
# interactions the data do not support. The point estimates stay Poisson
# (quasi-likelihood), which keeps the property that predictions sum to the
# observed total over the fitting sample. A model with dispersion <= 1 is
# returned as fit. The dispersion is recorded so the pooled refit uses a
# negative binomial for the same models.
.reselect_joint_count_qaic <- function(o) {
  selected <- attr(o, "selected_interactions")
  candidates <- attr(o, "candidate_interactions")
  base <- stats::formula(o)
  for (ix in selected) base <- stats::update(base, paste("~ . -", ix))
  full <- base
  for (ix in candidates) full <- stats::update(full, paste("~ . +", ix))
  d <- o$data
  fit_glm <- function(f) {
    tryCatch(suppressWarnings(stats::glm(f, data = d, family = stats::poisson())),
             error = function(e) NULL)
  }
  fit_full <- fit_glm(full)
  phi <- if (is.null(fit_full)) NA_real_ else
    sum(stats::residuals(fit_full, type = "pearson")^2) / fit_full$df.residual
  attr(o, "dispersion") <- phi
  if (is.na(phi) || phi <= 1) return(o)
  qaic <- function(f) -2 * as.numeric(stats::logLik(f)) / phi + 2 * (length(stats::coef(f)) + 1)
  best <- fit_glm(base)
  if (is.null(best)) return(o)
  new_sel <- character(0)
  for (ix in candidates) {
    cand <- fit_glm(stats::update(stats::formula(best), paste("~ . +", ix)))
    if (!is.null(cand) && isTRUE(cand$converged) && qaic(cand) < qaic(best)) {
      best <- cand
      new_sel <- c(new_sel, ix)
    }
  }
  if (setequal(new_sel, selected)) return(o)
  attr(best, "selected_interactions") <- new_sel
  attr(best, "candidate_interactions") <- candidates
  attr(best, "dispersion") <- phi
  best
}


.tmb_family <- function(o) {
  fam <- stats::family(o)$family
  phi <- attr(o, "dispersion")
  if (fam == "poisson" && !is.null(phi) && isTRUE(phi > 1)) return(glmmTMB::nbinom2())
  switch(fam,
    poisson  = stats::poisson(),
    binomial = stats::binomial(),
    gaussian = stats::gaussian(),
    stop("No pooled refit for family ", fam, ".", call. = FALSE)
  )
}


# Replace the geogYN coefficient of a glm / lm with the pooled contrast.
# Returns list(model, row) where row records the unpooled and pooled values.
.pool_glm_geog <- function(o, dat, family, groups, layer, model_name) {
  unpooled <- unname(stats::coef(o)["geogYN"])
  res <- .pooled_geog_contrast(stats::formula(o), dat, family, groups)
  row <- .pool_row(layer, model_name, "geogYN", unpooled, res, groups)
  if (!is.null(res) && is.null(res$error) && !is.na(unpooled)) {
    o$coefficients[["geogYN"]] <- res$effect
    attr(o, "geog_pooling") <- row
  }
  list(model = o, row = row)
}


.pool_row <- function(layer, model_name, term, unpooled, res, groups) {
  sd <- stats::setNames(rep(NA_real_, length(groups)), paste0("sd_", groups))
  pooled <- NA_real_
  n_geog <- NA_integer_
  if (!is.null(res) && !is.null(res$error)) {
    warning("Geographic pooling failed for ", layer, " ", model_name, " (", res$error,
            "); the unpooled geogYN coefficient is kept.", call. = FALSE)
  } else if (!is.null(res)) {
    pooled <- res$effect
    n_geog <- res$n_geog
    sd[paste0("sd_", names(res$sd))] <- res$sd
  }
  cbind(data.frame(layer = layer, model = model_name, term = term,
                   unpooled = unpooled, pooled = pooled, n_geog = n_geog,
                   stringsAsFactors = FALSE),
        as.data.frame(as.list(sd)))
}


# Pool the partner-race multinomial: each non-reference category's geogYN
# coefficient is replaced by the pooled contrast from the corresponding
# baseline-category binary logit (category k vs the reference) with the same
# right-hand side. The coefficients live in the nnet weight vector, so the
# edit carries through predict() and coef().
.pool_multinom_geog <- function(o, ldat, groups, layer) {
  lev <- o$lev
  ldat$.prace <- factor(ldat$p_race.cat.num, levels = lev)
  j <- match("geogYN", o$vcoefnames)
  cf <- stats::coef(o)
  rows <- list()
  for (k in seq_along(lev)[-1]) {
    sub <- ldat[!is.na(ldat$.prace) & ldat$.prace %in% lev[c(1, k)], , drop = FALSE]
    sub$.yk <- as.integer(sub$.prace == lev[k])
    fk <- stats::update(stats::formula(o$terms), .yk ~ .)
    unpooled <- if (is.matrix(cf)) cf[lev[k], "geogYN"] else cf[["geogYN"]]
    res <- .pooled_geog_contrast(fk, sub, stats::binomial(), groups)
    row <- .pool_row(layer, "joint_prace_model", paste0("geogYN:", lev[k]),
                     unname(unpooled), res, groups)
    if (!is.null(res) && is.null(res$error)) {
      idx <- if (length(lev) == 2L) 1L + j else (k - 1L) * (length(o$wts) / o$n[3L]) + 1L + j
      o$wts[idx] <- res$effect
    }
    rows[[length(rows) + 1]] <- row
  }
  rows <- do.call(rbind, rows)
  attr(o, "geog_pooling") <- rows
  list(model = o, rows = rows)
}


# Dirichlet-multinomial shrinkage of the target geography's distribution of x
# toward the national one. The intraclass correlation is the one-way ANOVA
# estimator over `grp`, averaged over the category indicators; the geography's
# weight is n / (n + (1 - icc) / icc).
.dm_shrink_dist <- function(x, grp, inc, levels) {
  keep <- !is.na(x) & !is.na(grp)
  x <- x[keep]
  grp <- as.character(grp[keep])
  inc <- inc[keep]
  n_g <- table(grp)
  n_grp <- length(n_g)
  n <- length(x)
  n0 <- (n - sum(n_g^2) / n) / (n_grp - 1)
  icc <- vapply(levels, function(k) {
    y <- as.numeric(x == k)
    m <- tapply(y, grp, mean)[names(n_g)]
    msb <- sum(n_g * (m - mean(y))^2) / (n_grp - 1)
    msw <- sum(tapply(y, grp, function(v) sum((v - mean(v))^2))) / (n - n_grp)
    den <- msb + (n0 - 1) * msw
    if (den <= 0) 0 else max(0, (msb - msw) / den)
  }, numeric(1))
  rho <- mean(icc)
  n_geog <- sum(inc)
  w <- if (rho <= 0 || n_geog == 0) 0 else n_geog / (n_geog + (1 - rho) / rho)
  p_geog <- if (n_geog > 0) tabulate(match(x[inc], levels), length(levels)) / n_geog else 0
  p_nat <- tabulate(match(x, levels), length(levels)) / n
  list(p = w * p_geog + (1 - w) * p_nat, weight = w, icc = rho)
}


# Post-pass over a build_netparams() output under geog.pooling =
# "hierarchical": NB refits of the count models, pooled geogYN coefficients in
# every joint formation model, pooled cross-layer degree distributions and
# one-off risk-group rates, and a summary at out$geog.pooling. The duration
# models are pooled inside fit_joint_lm_dur(), where the by-age durations are
# computed; their rows are collected here from the model attributes.
.pool_netparams <- function(out, geog.lvl) {
  groups <- .geog_pool_groups(geog.lvl)
  rows <- list()
  glm_models <- c("joint_model", "joint_concurrent_model", "joint_nm_age_model",
                  "joint_nm_race_model", "joint_absdiff_age_model",
                  "joint_absdiff_sqrtage_model")
  for (layer in c("main", "casl", "inst")) {
    np <- out[[layer]]
    np$joint_model <- .reselect_joint_count_qaic(np$joint_model)
    for (m in glm_models) {
      o <- np[[m]]
      if (is.null(o)) next
      res <- .pool_glm_geog(o, o$data, .tmb_family(o), groups, layer, m)
      np[[m]] <- res$model
      rows[[length(rows) + 1]] <- res$row
    }
    if (!is.null(np$joint_prace_model)) {
      res <- .pool_multinom_geog(np$joint_prace_model, np$joint_nm_age_model$data,
                                 groups, layer)
      np$joint_prace_model <- res$model
      rows[[length(rows) + 1]] <- res$rows
    }
    dur_row <- attr(np$joint_duration_model, "geog_pooling")
    if (!is.null(dur_row)) rows[[length(rows) + 1]] <- dur_row
    out[[layer]] <- np
  }

  # Raw geography-specific distributions that seed the synthetic population.
  finest <- paste0(".pool_", groups[length(groups)])
  shrink <- function(dat, var, levels) {
    dat <- .with_geog_groups(dat, groups)
    .dm_shrink_dist(dat[[var]], dat[[finest]], dat$geogYN == 1, levels)
  }
  dd_casl <- shrink(out$main$joint_model$data, "deg.casl", 0:3)
  dd_main <- shrink(out$casl$joint_model$data, "deg.main", 0:2)
  dd_tot  <- shrink(out$inst$joint_model$data, "deg.tot3", 0:3)
  out$main$deg.casl.dist <- dd_casl$p
  out$casl$deg.main.dist <- dd_main$p
  out$inst$deg.tot.dist  <- dd_tot$p

  # One-off risk-group rates: rescale from the raw geography ratio of mean
  # one-off rates to the pooled one-off rate ratio.
  di <- out$inst$joint_model$data
  raw_ratio <- mean(di$rate.oo.part[di$geogYN == 1], na.rm = TRUE) /
    mean(di$rate.oo.part, na.rm = TRUE)
  pooled_ratio <- exp(unname(stats::coef(out$inst$joint_model)["geogYN"]))
  if (is.finite(raw_ratio) && raw_ratio > 0 && is.finite(pooled_ratio)) {
    out$inst$nf.risk.grp <- out$inst$nf.risk.grp * pooled_ratio / raw_ratio
  }

  effects <- do.call(rbind, rows)
  rownames(effects) <- NULL
  out$geog.pooling <- list(
    method = "hierarchical",
    groups = groups,
    effects = effects,
    count_dispersion = vapply(c("main", "casl", "inst"), function(l) {
      unname(attr(out[[l]]$joint_model, "dispersion"))
    }, numeric(1)),
    dist_weights = c(deg.casl = dd_casl$weight, deg.main = dd_main$weight,
                     deg.tot = dd_tot$weight),
    risk.grp_ratio = c(raw = raw_ratio, pooled = pooled_ratio)
  )
  out
}
