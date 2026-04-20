# Weibull AFT fit for partnership durations with right-censoring.
# In survreg, `event = 1` means observed / completed (ONGOING = 0) and
# `event = 0` means censored / still ongoing (ONGOING = 1). Uses both
# completed and ongoing partnerships. Returns a small list with the
# stratum-level summary quantities plus the shape parameter `k` as a
# diagnostic for the geometric-distribution assumption: k ~= 1 means
# the constant-hazard (geometric) assumption holds; k far from 1 means
# the observed hazard is increasing (k > 1) or decreasing (k < 1) in
# partnership age. Returns NULL when the fit is infeasible (too few
# events / too few censored / convergence failure) so callers can fall
# back to empirical.
fit_weibull_dur <- function(data) {
  if (nrow(data) < 10 || !requireNamespace("survival", quietly = TRUE)) {
    return(NULL)
  }
  event <- as.integer(data$ongoing2 == 0)
  n_events <- sum(event)
  n_censored <- nrow(data) - n_events
  # survreg needs some of each; otherwise scale is unidentified.
  if (n_events < 3 || n_censored < 3) return(NULL)
  if (any(data$duration.time <= 0, na.rm = TRUE)) {
    data <- data[data$duration.time > 0 & !is.na(data$duration.time), , drop = FALSE]
    event <- as.integer(data$ongoing2 == 0)
    if (nrow(data) < 10) return(NULL)
  }
  fit <- tryCatch(
    suppressWarnings(survival::survreg(
      survival::Surv(data$duration.time, event) ~ 1,
      dist = "weibull"
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  # survreg parameterization: weibull_shape = 1 / fit$scale;
  # weibull_scale = exp((Intercept)).
  shape <- 1 / fit$scale
  scale_weibull <- exp(coef(fit)[["(Intercept)"]])
  med <- scale_weibull * log(2)^(1 / shape)            # closed-form median
  mean_val <- scale_weibull * gamma(1 + 1 / shape)     # closed-form mean
  list(mean.dur = mean_val, median.dur = med, weibull_shape = shape)
}


# Joint log-linear lm on log(duration.time) for ongoing partnerships.
# Same length-bias convention as the "empirical" default (uses only
# ongoing, relies on the memoryless property to read stratum-specific
# predicted durations as the median full-duration quantity that feeds
# the geometric transformation in build_netstats). Returns the fitted
# model plus a vector of per-row predicted durations on the training
# set so callers can stratify however they like.
fit_joint_lm_dur <- function(l_layer, race, geog.lvl) {
  ongoing <- l_layer[l_layer$ongoing2 == 1 &
                       !is.na(l_layer$duration.time) &
                       l_layer$duration.time > 0, , drop = FALSE]
  if (nrow(ongoing) < 20) return(NULL)
  terms <- c("index.age.grp", "sqrt(index.age.grp)", "hiv2", "same.age.grp")
  if (isTRUE(race)) {
    terms <- c(terms, "as.factor(race.cat.num)", "same.race")
  }
  if (!is.null(geog.lvl)) {
    terms <- c("geogYN", terms)
  }
  fml <- reformulate(terms, response = "log(duration.time)")
  fit <- tryCatch(
    suppressWarnings(lm(fml, data = ongoing)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  fitted_dur <- exp(predict(fit, newdata = ongoing))
  list(model = fit, ongoing = ongoing, fitted_dur = fitted_dur)
}


# Given a partnership-level layer and a non-default duration.method,
# return a list(homog, byage, joint_duration_model) of duration
# summaries matching the shape of the empirical output (so that the
# geometric transformation, smoothing, and sex.cess.mod logic in the
# calling layer block apply uniformly). `byage_strata` is the vector
# of stratum keys expected in the output (first = 0 for non-matched,
# subsequent = index.age.grp values for matched-within-age-group).
compute_alt_durs <- function(l_layer, duration.method, byage_strata,
                             race, geog.lvl) {
  stopifnot(duration.method %in% c("weibull_strat", "joint_lm"))

  # Subset each stratum the same way the empirical block does.
  strata_data <- lapply(byage_strata, function(k) {
    if (k == 0) {
      l_layer[l_layer$same.age.grp == 0, , drop = FALSE]
    } else {
      l_layer[l_layer$same.age.grp == 1 & l_layer$index.age.grp == k,
              , drop = FALSE]
    }
  })
  names(strata_data) <- as.character(byage_strata)

  if (duration.method == "weibull_strat") {
    # Per-stratum Weibull fits.
    per_stratum <- lapply(strata_data, fit_weibull_dur)
    byage <- data.frame(
      index.age.grp = byage_strata,
      mean.dur = vapply(per_stratum,
        function(x) if (is.null(x)) NA_real_ else x$mean.dur, numeric(1)),
      median.dur = vapply(per_stratum,
        function(x) if (is.null(x)) NA_real_ else x$median.dur, numeric(1))
    )
    shapes <- vapply(per_stratum,
      function(x) if (is.null(x)) NA_real_ else x$weibull_shape, numeric(1))

    # Overall (homog) fit on the whole subset used for durations.
    homog_fit <- fit_weibull_dur(l_layer)
    homog <- if (is.null(homog_fit)) {
      data.frame(mean.dur = NA_real_, median.dur = NA_real_)
    } else {
      data.frame(mean.dur = homog_fit$mean.dur,
                 median.dur = homog_fit$median.dur)
    }
    return(list(homog = homog, byage = byage, weibull_shapes = shapes,
                weibull_shape_overall = if (is.null(homog_fit)) NA_real_
                                         else homog_fit$weibull_shape))
  }

  # duration.method == "joint_lm"
  res <- fit_joint_lm_dur(l_layer, race = race, geog.lvl = geog.lvl)
  if (is.null(res)) {
    return(list(homog = data.frame(mean.dur = NA_real_, median.dur = NA_real_),
                byage = data.frame(index.age.grp = byage_strata,
                                   mean.dur = NA_real_, median.dur = NA_real_),
                model = NULL))
  }
  # Stratum medians from predicted durations on the training set.
  byage <- do.call(rbind, lapply(byage_strata, function(k) {
    sub <- if (k == 0) {
      res$ongoing[res$ongoing$same.age.grp == 0, , drop = FALSE]
    } else {
      res$ongoing[res$ongoing$same.age.grp == 1 &
                    res$ongoing$index.age.grp == k, , drop = FALSE]
    }
    if (nrow(sub) == 0) {
      data.frame(index.age.grp = k, mean.dur = NA_real_, median.dur = NA_real_)
    } else {
      pred_sub <- exp(predict(res$model, newdata = sub))
      data.frame(index.age.grp = k,
                 mean.dur = mean(pred_sub, na.rm = TRUE),
                 median.dur = median(pred_sub, na.rm = TRUE))
    }
  }))
  homog <- data.frame(mean.dur = mean(res$fitted_dur, na.rm = TRUE),
                      median.dur = median(res$fitted_dur, na.rm = TRUE))
  list(homog = homog, byage = byage, model = res$model)
}


# Fit a joint GLM predicting `response` from all available individual
# attributes (age, race, the concurrent-layer degree, HIV status,
# geography). Candidate interactions are considered one at a time and
# kept only when they reduce AIC. `family` selects between the Poisson
# model for degree / one-off count and the binomial model for the
# concurrency indicator (1 = deg > 1). Internal helper for the
# `method = "joint"` path of `build_netparams()`; see issues #61 / #62.
fit_joint_glm <- function(d, response, main_terms,
                          race, geog.lvl,
                          interaction_cross_deg = NULL,
                          family = poisson()) {
  if (isTRUE(race)) {
    main_terms <- c(main_terms, "as.factor(race.cat.num)")
  }
  if (!is.null(geog.lvl)) {
    main_terms <- c("geogYN", main_terms)
  }

  base_fml <- reformulate(main_terms, response = response)
  best <- glm(base_fml, data = d, family = family)
  selected <- character(0)

  candidates <- character(0)
  if (isTRUE(race)) {
    candidates <- c(candidates, "age.grp:as.factor(race.cat.num)")
  }
  if (!is.null(interaction_cross_deg)) {
    candidates <- c(candidates, paste0("age.grp:", interaction_cross_deg))
  }

  for (ix in candidates) {
    ix_fml <- update(formula(best), paste("~ . +", ix))
    cand <- tryCatch(
      suppressWarnings(glm(ix_fml, data = d, family = family)),
      error = function(e) NULL
    )
    if (!is.null(cand) && isTRUE(cand$converged) && AIC(cand) < AIC(best)) {
      best <- cand
      selected <- c(selected, ix)
    }
  }

  attr(best, "selected_interactions") <- selected
  attr(best, "candidate_interactions") <- candidates
  best
}


#' Calculate Individual-Level Network Parameters
#'
#' @description Builds statistical models predicting mean degree, mixing, and duration of sexual
#'              partnerships, for use in the EpiModelHIV workflow.
#'
#' @param epistats Output from [`build_epistats`].
#' @param smooth.main.dur If `TRUE`, function averages the main sexual partnership durations for
#'        oldest and second oldest age groups.
#' @param oo.nquants Number of quantiles to split the one-off partnership risk distribution (count
#'        of one-off partners per unit time).
#' @param duration.method Character. Controls how partnership durations for dissolution-coef
#'        estimation are computed. One of:
#'        \itemize{
#'          \item `"empirical"` (default) — mean / median of observed durations among ongoing
#'            partnerships, stratified by (age-match x index.age.grp). Relies on the
#'            geometric / memoryless assumption that median elapsed time in ongoing
#'            partnerships equals median full-partnership duration. Byte-identical to the
#'            pre-refactor behavior.
#'          \item `"weibull_strat"` — Weibull AFT fits per stratum via
#'            `survival::survreg(Surv(duration.time, 1 - ongoing2) ~ 1, ...)`. Uses both
#'            ongoing (right-censored) and completed partnerships. The Weibull shape
#'            parameter `k` is attached to `durs.<layer>.byage` as a `"weibull_shape"`
#'            attribute — diagnostic for the constant-hazard assumption embedded in the
#'            TERGM dissolution (k ~= 1 means geometric is correct; k far from 1 flags
#'            mis-specification). **Caveat**: on the ARTnet data the fitted `k` comes
#'            out well below 1 (decreasing hazard) in every stratum, and the implied
#'            Weibull median can extrapolate to implausibly large values in heavily
#'            censored strata (e.g., the oldest matched age groups). Treat
#'            `"weibull_strat"` as diagnostic — look at the `k` values to decide
#'            whether the geometric assumption is defensible — rather than as a
#'            drop-in production method. `"joint_lm"` is the more production-safe
#'            non-default option.
#'          \item `"joint_lm"` — log-linear `lm(log(duration.time) ~ <joint ego + partner + match
#'            terms>)` on ongoing partnerships; per-stratum medians computed from model
#'            predictions. Fitted model is stored at `netparams$<layer>$joint_duration_model`
#'            for optional per-dyad use by future `build_netstats` refactors.
#'        }
#'        Under all three methods, the stratum-level medians flow through the same
#'        geometric transformation (`mean.dur.adj = 1 / (1 - 2^(-1 / median))`) that
#'        `build_netstats` passes to `dissolution_coefs()`, so the TERGM offset structure
#'        is identical across methods. The one-off layer has no durations, so this flag
#'        only affects main and casual-layer output.
#' @param method Character. Either `"existing"` (default) or `"joint"`. `"existing"` reproduces
#'        the pre-refactor behavior byte-for-byte: a separate univariate Poisson/binomial/linear
#'        fit for each ERGM target statistic. `"joint"` leaves all of those outputs intact **and**
#'        additionally fits joint GLMs per layer with AIC-based interaction selection:
#'        \itemize{
#'          \item \strong{Ego-level}: a Poisson model for degree / one-off count at
#'            `$<layer>$joint_model`, and (main/casual only) a binomial for the
#'            concurrency indicator at `$<layer>$joint_concurrent_model`.
#'          \item \strong{Dyad-level} (partnership data with ego attrs on RHS):
#'            binomial for same-age-group at `$<layer>$joint_nm_age_model`, binomial
#'            for same-race (when `race = TRUE`) at `$<layer>$joint_nm_race_model`,
#'            and Gaussian fits for the age absolute-difference terms at
#'            `$<layer>$joint_absdiff_age_model` and
#'            `$<layer>$joint_absdiff_sqrtage_model`.
#'        }
#'        These models are consumed by [`build_netstats`] under `method = "joint"` to produce
#'        internally-consistent g-computation target statistics.
#' @param browser If `TRUE`, run `build_netparams` in interactive browser mode.
#'
#' @details
#' `build_netparams` is a helper function that constructs the necessary network parameters for use
#' in building network models with [`build_netstats`], building on models estimated using
#' [`build_epistats`].
#'
#' The parameter `smooth.main.dur` is used when partnership duration and mortality compete in the
#' eldest age group; in such a case mean duration is averaged over the oldest and second oldest age
#' groups (as specified by `age.breaks` in [`build_epistats`]). Subsequently, this smoothing is only
#' done if there are three or more age categories specified. Note, this does not affect calculations
#' if an age group after sexual cessation is included; durations averaged only in the oldest two
#' age groups within the bounds of the sexual cessation age.
#'
#' @export
#' @examples
#' # Standard model with default age stratification
#' epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta")
#' netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
#'
#' # Restricted age stratification
#' epistats2 <- build_epistats(geog.lvl = "state", geog.cat = "GA",
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(20, 30, 40))
#' netparams2 <- build_netparams(epistats2, smooth.main.dur = TRUE)
#'
#' # Model with sexual cessation age < age limit
#' epistats3 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams3 <- build_netparams(epistats3, smooth.main.dur = TRUE)
#'
#' # Fit joint Poisson GLMs in addition to the univariate marginals
#' netparams4 <- build_netparams(epistats, smooth.main.dur = TRUE, method = "joint")
#' summary(netparams4$main$joint_model)
#'
build_netparams <- function(epistats,
                            smooth.main.dur = FALSE,
                            oo.nquants = 5,
                            duration.method = c("empirical", "weibull_strat", "joint_lm"),
                            method = c("existing", "joint"),
                            browser = FALSE) {
  method <- match.arg(method)
  duration.method <- match.arg(duration.method)

  if (duration.method == "weibull_strat" &&
      !requireNamespace("survival", quietly = TRUE)) {
    stop("duration.method = 'weibull_strat' requires the 'survival' package. ",
         "Install it with install.packages('survival').")
  }

  # Ensures that ARTnetData is installed
  if (system.file(package = "ARTnetData") == "") stop(missing_data_msg)

  if (browser == TRUE) {
    browser()
  }

  ## Inputs ##
  geog.lvl <- epistats$geog.lvl
  race <- epistats$race
  race.level <- epistats$race.level
  age.limits <- epistats$age.limits
  age.breaks <- epistats$age.breaks
  age.sexual.cessation <- epistats$age.sexual.cessation
  sex.cess.mod <- epistats$sex.cess.mod
  age.grps <- epistats$age.grps
  time.unit <- epistats$time.unit

  # Fix global binding check error
  duration.time <- NULL

  # 0. Data Processing ------------------------------------------------------

  ## Age Processing ##

  ## Data ##
  d <- ARTnetData::ARTnet.wide
  l <- ARTnetData::ARTnet.long

  # p_age_imp initialization for lintr
  p_age_imp <- NULL

  # Subset datasets by lower age limit and age.sexual.cessation
  # Now applies to both index (respondents) and partners for long dataset
  l <- subset(l, age >= age.limits[1] & age < age.sexual.cessation &
                p_age_imp >= age.limits[1] & p_age_imp < age.sexual.cessation)
  d <- subset(d, age >= age.limits[1] & age < age.sexual.cessation)

  l$comb.age <- l$age + l$p_age_imp
  l$diff.age <- abs(l$age - l$p_age_imp)

  l$duration.time <- l$duration * 7 / time.unit

  #Append Data when geog.lvl is defined
  if (!is.null(geog.lvl)) {
    d$geog <- epistats$geog.d
    d$geogYN <- epistats$geogYN.d
    l$geog <- epistats$geog.l
    l$geogYN <- epistats$geogYN.l
  }

  ## Degree calculations ##

  l$ONGOING <- as.numeric(l$ONGOING)
  l$ongoing2 <- ifelse(is.na(l$ONGOING), 0, l$ONGOING)
  l$ONGOING <- NULL

  d <- l %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ptype == 1) %>%
    group_by(AMIS_ID) %>%
    summarise(deg.main = sum(ongoing2)) %>%
    right_join(d, by = "AMIS_ID")

  d <- l %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ptype == 2) %>%
    group_by(AMIS_ID) %>%
    summarise(deg.casl = sum(ongoing2)) %>%
    right_join(d, by = "AMIS_ID")

  # If missing degree, then set to 0
  d$deg.main <- ifelse(is.na(d$deg.main), 0, d$deg.main)
  d$deg.casl <- ifelse(is.na(d$deg.casl), 0, d$deg.casl)

  # recoding to truncate degree
  d$deg.casl <- ifelse(d$deg.casl > 3, 3, d$deg.casl)
  d$deg.main <- ifelse(d$deg.main > 2, 2, d$deg.main)

  d$deg.tot <- d$deg.main + d$deg.casl

  # Concurrency
  d$deg.main.conc <- ifelse(d$deg.main > 1, 1, 0)
  d$deg.casl.conc <- ifelse(d$deg.casl > 1, 1, 0)

  ## one-off calcs ##

  # Total MC anal sex partner count
  d <- l %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ptype %in% 1:2) %>%
    group_by(AMIS_ID) %>%
    count() %>%
    rename(count.mc.part = n) %>%
    right_join(d, by = "AMIS_ID")
  d$count.mc.part <- ifelse(is.na(d$count.mc.part), 0, d$count.mc.part)

  d$count.oo.part <- d$ai.part - d$count.mc.part
  d$count.oo.part <- pmax(0, d$count.oo.part)

  # Truncated OO part
  d$count.oo.part.trunc <- ifelse(d$count.oo.part > 100, 100, d$count.oo.part)


  ## Race ethnicity ##
  if (race == TRUE) {
    mult_race_cat <- c("asian", "ai/an", "mult", "nh/pi")
    flat_race.level <- unlist(race.level)

    # Determine which variables to use in ARTnet
    if (any(flat_race.level %in% mult_race_cat)) {
      l <- merge(l, d[, c("AMIS_ID", "race")], by = "AMIS_ID", all.x = TRUE)
      p_race_var <- "p_race2"
      race_var <- "race"
    } else {
      p_race_var <- "p_race.cat"
      race_var <- "race.cat"
    }

    # Assign race categories based on race.level
    race.categories <- seq_along(race.level)

    d$race.cat.num <- rep(NA, nrow(d))
    l$race.cat.num <- rep(NA, nrow(l))
    l$p_race.cat.num <- rep(NA, nrow(l))

    for (i in seq_along(race.level)) {
      d$race.cat.num[d[[race_var]] %in% race.level[[i]]] <- race.categories[i]
      l$race.cat.num[l[[race_var]] %in% race.level[[i]]] <- race.categories[i]
      l$p_race.cat.num[l[[p_race_var]] %in% race.level[[i]]] <- race.categories[i]
    }

    # Redistribute NAs in proportion to non-missing partner races
    probs <- prop.table(table(l$race.cat.num, l$p_race.cat.num), 1)

    for (i in race.categories) {
      imp_indices <- which(is.na(l$p_race.cat.num) & l$race.cat.num == i)
      if (length(imp_indices) > 0) {
        l$p_race.cat.num[imp_indices] <- sample(race.categories, length(imp_indices), TRUE, probs[i, ])
      }
    }

    # Assign race.combo
    l$race.combo <- make_race_combo(l$race.cat.num, l$p_race.cat.num)

  }

  ## HIV status

  l$p_hiv2 <- as.integer(l$p_hiv == 1)
  table(l$p_hiv, l$p_hiv2, useNA = "always")

  hiv.combo <- rep(NA, nrow(l))
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5

  l$hiv.concord.pos <- as.integer(hiv.combo == 2)

  ## Setup output list ##

  out <- list()


  # 1. Main Model -----------------------------------------------------------

  out$main <- list()
  lmain <- l[l$ptype == 1, ]


  ## edges ----
  if (is.null(geog.lvl)) {
    mod <- glm(deg.main ~ 1, data = d, family = poisson())

    pred <- exp(coef(mod)[[1]])

    out$main$md.main <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geogYN, data = d, family = poisson())

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$md.main <- as.numeric(pred)
  }


  ## nodematch("age.grp") ----

  lmain$index.age.grp <- cut(lmain$age, age.breaks, labels = FALSE,
                             right = FALSE, include.lowest = FALSE)
  lmain$part.age.grp <- cut(as.numeric(lmain$p_age_imp), age.breaks, labels = FALSE,
                            right = FALSE, include.lowest = FALSE)

  lmain$same.age.grp <- as.integer(lmain$index.age.grp == lmain$part.age.grp)

  if (is.null(geog.lvl)) {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = lmain, family = binomial())

    dat <- data.frame(index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geogYN + index.age.grp,
               data = lmain, family = binomial())

    dat <- data.frame(geogYN = 1, index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  lmain$ad <- abs(lmain$age - lmain$p_age_imp)
  lmain$ad.sr <- abs(sqrt(lmain$age) - sqrt(lmain$p_age_imp))

  if (is.null(geog.lvl)) {
    mod <- lm(ad ~ 1, data = lmain)

    pred <- coef(mod)[[1]]

    out$main$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geogYN, data = lmain)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (is.null(geog.lvl)) {
    mod <- lm(ad.sr ~ 1, data = lmain)

    pred <- coef(mod)[[1]]

    out$main$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geogYN, data = lmain)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE,
                   right  = FALSE, include.lowest = FALSE)

  if (is.null(geog.lvl)) {
    mod <- glm(deg.main ~ + age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geogYN + age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.age.grp <- as.numeric(pred)
  }


  ## nodematch("race", diff = TRUE) ----

  if (race == TRUE) {
    lmain$same.race <- as.integer(lmain$race.cat.num == lmain$p_race.cat.num)

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ as.factor(race.cat.num),
                 data = lmain, family = binomial())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN + as.factor(race.cat.num),
                 data = lmain, family = binomial())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ 1, data = lmain, family = binomial())

      pred <- exp(coef(mod)[[1]]) / (1 + exp(coef(mod)[[1]]))

      out$main$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN,
                 data = lmain, family = binomial())

      dat <- data.frame(geogYN = 1)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race_diffF <- as.numeric(pred)
    }

    if (is.null(geog.lvl)) {
      mod <- glm(deg.main ~ as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(deg.main ~ geogYN + as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("deg.casl") ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.main ~ deg.casl,
               data = d, family = poisson())

    dat <- data.frame(deg.casl = sort(unique(d$deg.casl)))
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.deg.casl <- as.numeric(pred)

    deg.casl.dist <- prop.table(table(d$deg.casl))
    out$main$deg.casl.dist <- as.numeric(deg.casl.dist)
  } else {
    mod <- glm(deg.main ~ geogYN + deg.casl,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, deg.casl = sort(unique(d$deg.casl)))
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.deg.casl <- as.numeric(pred)

    deg.casl.dist <- prop.table(table(d$deg.casl[d$geogYN == 1]))
    out$main$deg.casl.dist <- as.numeric(deg.casl.dist)
  }


  ## concurrent ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.main.conc ~ 1, data = d, family = binomial())

    pred <- exp(coef(mod)[[1]]) / (1 + exp(coef(mod)[[1]]))

    out$main$concurrent <- as.numeric(pred)
  } else {
    mod <- glm(deg.main.conc ~ geogYN,
               data = d, family = binomial())

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$concurrent <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.main ~ hiv2, data = d, family = poisson())

    dat <- data.frame(hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geogYN + hiv2,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.diag.status <- as.numeric(pred)
  }


  ## Durations ----

  # overall
  durs.main <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ongoing2 == 1) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()

  # create city weights
  if (!is.null(geog.lvl)) {
    durs.main.geo <- lmain %>%
      filter(RAI == 1 | IAI == 1) %>%
      filter(ongoing2 == 1) %>%
      filter(geogYN == 1) %>%
      summarise(mean.dur = mean(duration.time, na.rm = TRUE),
                median.dur = median(duration.time, na.rm = TRUE)) %>%
      as.data.frame()

    # city-specific weight based on ratio of medians
    wt <- durs.main.geo$median.dur / durs.main$median.dur
  } else {
    wt <- 1
  }

  # The dissolution rate is function of the mean of the geometric distribution
  # which relates to the median as:
  durs.main$rates.main.adj <- 1 - (2^(-1 / (wt * durs.main$median.dur)))

  # Mean duration associated with a geometric distribution that median:
  durs.main$mean.dur.adj <- 1 / (1 - (2^(-1 / (wt * durs.main$median.dur))))
  out$main$durs.main.homog <- durs.main

  # stratified by age

  # first, non-matched by age group
  durs.main.nonmatch <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 0) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()
  durs.main.nonmatch$index.age.grp <- 0

  # then, matched within age-groups
  durs.main.matched <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 1) %>%
    group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()
  durs.main.matched

  durs.main.all <- rbind(durs.main.nonmatch, durs.main.matched)

  durs.main.all$rates.main.adj <- 1 - (2^(-1 / (wt * durs.main.all$median.dur)))
  durs.main.all$mean.dur.adj <- 1 / (1 - (2^(-1 / (wt * durs.main.all$median.dur))))

  durs.main.all <- durs.main.all[, c(3, 1, 2, 4, 5)]
  out$main$durs.main.byage <- durs.main.all


  ## duration.method override (see #63 phase 3) ----
  # Replace the empirical stratum medians with model-based medians and
  # re-apply the geometric transformation, preserving the downstream
  # data.frame shape so smoothing and sex.cess.mod below work unchanged.
  if (duration.method != "empirical") {
    alt <- compute_alt_durs(
      l_layer = lmain[lmain$RAI == 1 | lmain$IAI == 1, , drop = FALSE],
      duration.method = duration.method,
      byage_strata = out$main$durs.main.byage$index.age.grp,
      race = race, geog.lvl = geog.lvl
    )
    # homog: replace mean/median, recompute adj from new median
    if (!is.na(alt$homog$median.dur)) {
      out$main$durs.main.homog$mean.dur   <- alt$homog$mean.dur
      out$main$durs.main.homog$median.dur <- alt$homog$median.dur
      out$main$durs.main.homog$rates.main.adj <-
        1 - 2^(-1 / (wt * out$main$durs.main.homog$median.dur))
      out$main$durs.main.homog$mean.dur.adj <-
        1 / (1 - 2^(-1 / (wt * out$main$durs.main.homog$median.dur)))
    }
    # byage: per-stratum replacement, with empirical fallback if a stratum
    # fit failed (keeps the row populated so dissolution_coefs doesn't NA).
    for (i in seq_len(nrow(out$main$durs.main.byage))) {
      k <- out$main$durs.main.byage$index.age.grp[i]
      j <- which(alt$byage$index.age.grp == k)
      if (length(j) == 1 && !is.na(alt$byage$median.dur[j])) {
        out$main$durs.main.byage$mean.dur[i]        <- alt$byage$mean.dur[j]
        out$main$durs.main.byage$median.dur[i]      <- alt$byage$median.dur[j]
        out$main$durs.main.byage$rates.main.adj[i]  <-
          1 - 2^(-1 / (wt * alt$byage$median.dur[j]))
        out$main$durs.main.byage$mean.dur.adj[i] <-
          1 / (1 - 2^(-1 / (wt * alt$byage$median.dur[j])))
      }
    }
    if (duration.method == "weibull_strat") {
      attr(out$main$durs.main.byage, "weibull_shape") <- alt$weibull_shapes
      attr(out$main$durs.main.homog, "weibull_shape") <- alt$weibull_shape_overall
    }
    if (duration.method == "joint_lm" && !is.null(alt$model)) {
      out$main$joint_duration_model <- alt$model
    }
  }


  if (smooth.main.dur == TRUE) {
    n2 <- nrow(durs.main.all)
    n1 <- n2 - 1
    if (n2 > 3) {
      out$main$durs.main.byage$mean.dur.adj[n2] <-
        mean(out$main$durs.main.byage$mean.dur.adj[n1:n2])
    }
  }

  # If sexual cessation model, then set diss coef for age grp above boundary to 1
  if (sex.cess.mod == TRUE) {
    index.age.grp <- max(out$main$durs.main.byage$index.age.grp) + 1
    df <- data.frame(index.age.grp = index.age.grp, mean.dur = 1, median.dur = 1,
                     rates.main.adj = 1, mean.dur.adj = 1)
    out$main$durs.main.byage <- rbind(out$main$durs.main.byage, df)
  }


  ## joint g-computation models (additive outputs; see issues #61/#62/#63) ----
  if (method == "joint") {
    # Ego-level Poisson + binomial (edges / nodefactor / concurrent)
    out$main$joint_model <- fit_joint_glm(
      d,
      response = "deg.main",
      main_terms = c("age.grp", "sqrt(age.grp)", "deg.casl", "hiv2"),
      race = race,
      geog.lvl = geog.lvl,
      interaction_cross_deg = "deg.casl",
      family = poisson()
    )
    out$main$joint_concurrent_model <- fit_joint_glm(
      d,
      response = "deg.main.conc",
      main_terms = c("age.grp", "sqrt(age.grp)", "deg.casl", "hiv2"),
      race = race,
      geog.lvl = geog.lvl,
      interaction_cross_deg = "deg.casl",
      family = binomial()
    )

    # Dyad-level (nodematch / absdiff) -- fit on the partnership-level
    # data lmain with ego attributes only on the RHS (Option A per #63).
    # Alias index.age.grp -> age.grp so fit_joint_glm's hardcoded
    # interaction-term naming works unchanged.
    lmain$age.grp <- lmain$index.age.grp

    out$main$joint_nm_age_model <- fit_joint_glm(
      lmain,
      response = "same.age.grp",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = binomial()
    )
    if (isTRUE(race)) {
      out$main$joint_nm_race_model <- fit_joint_glm(
        lmain,
        response = "same.race",
        main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
        race = race, geog.lvl = geog.lvl,
        family = binomial()
      )
    }
    out$main$joint_absdiff_age_model <- fit_joint_glm(
      lmain,
      response = "ad",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
    out$main$joint_absdiff_sqrtage_model <- fit_joint_glm(
      lmain,
      response = "ad.sr",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
  }


  # 2. Casual Model ---------------------------------------------------------

  out$casl <- list()
  lcasl <- l[l$ptype == 2, ]


  ## edges ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.casl ~ 1, data = d, family = poisson())

    pred <- exp(coef(mod)[[1]])

    out$casl$md.casl <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geogYN, data = d, family = poisson())

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$md.casl <- as.numeric(pred)
  }

  ## nodematch("age.grp") ----

  lcasl$index.age.grp <- cut(lcasl$age, age.breaks, labels = FALSE, right = FALSE,
                             include.lowest = FALSE)
  lcasl$part.age.grp <- cut(as.numeric(lcasl$p_age_imp), age.breaks,
                            right = FALSE, labels = FALSE, include.lowest = FALSE)

  lcasl$same.age.grp <- as.integer(lcasl$index.age.grp == lcasl$part.age.grp)

  if (is.null(geog.lvl)) {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = lcasl, family = binomial())

    dat <- data.frame(index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geogYN + index.age.grp,
               data = lcasl, family = binomial())

    dat <- data.frame(geogYN = 1, index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  lcasl$ad <- abs(lcasl$age - lcasl$p_age_imp)
  lcasl$ad.sr <- abs(sqrt(lcasl$age) - sqrt(lcasl$p_age_imp))

  if (is.null(geog.lvl)) {
    mod <- lm(ad ~ 1, data = lcasl)

    pred <- coef(mod)[[1]]

    out$casl$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geogYN, data = lcasl)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (is.null(geog.lvl)) {
    mod <- lm(ad.sr ~ 1, data = lcasl)

    pred <- coef(mod)[[1]]

    out$casl$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geogYN, data = lcasl)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE,
                   right = FALSE, include.lowest = FALSE)

  if (is.null(geog.lvl)) {
    mod <- glm(deg.casl ~ age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geogYN + age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.age.grp <- as.numeric(pred)
  }


  if (race == TRUE) {

    ## nodematch("race") ----

    lcasl$same.race <- as.integer(lcasl$race.cat.num == lcasl$p_race.cat.num)

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ as.factor(race.cat.num),
                 data = lcasl, family = binomial())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN + as.factor(race.cat.num),
                 data = lcasl, family = binomial())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ 1,
                 data = lcasl, family = binomial())

      pred <- exp(coef(mod)[[1]]) / (1 + exp(coef(mod)[[1]]))

      out$casl$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN,
                 data = lcasl, family = binomial())

      dat <- data.frame(geogYN = 1)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race_diffF <- as.numeric(pred)
    }

    ## nodefactor("race") ----

    if (is.null(geog.lvl)) {

      mod <- glm(deg.casl ~ as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(deg.casl ~ geogYN + as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("deg.main") ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.casl ~ deg.main,
               data = d, family = poisson())

    dat <- data.frame(deg.main = 0:2)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.deg.main <- as.numeric(pred)

    deg.main.dist <- prop.table(table(d$deg.main))
    out$casl$deg.main.dist <- as.numeric(deg.main.dist)
  } else {
    mod <- glm(deg.casl ~ geogYN + deg.main,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, deg.main = 0:2)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.deg.main <- as.numeric(pred)

    deg.main.dist <- prop.table(table(d$deg.main[d$geogYN == 1]))
    out$casl$deg.main.dist <- as.numeric(deg.main.dist)
  }


  ## concurrent ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.casl.conc ~ 1,
               data = d, family = binomial())

    pred <- exp(coef(mod)[[1]]) / (1 + exp(coef(mod)[[1]]))

    out$casl$concurrent <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl.conc ~ geogYN,
               data = d, family = binomial())

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$concurrent <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (is.null(geog.lvl)) {
    mod <- glm(deg.casl ~ hiv2,
               data = d, family = poisson())

    dat <- data.frame(hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geogYN + hiv2,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.diag.status <- as.numeric(pred)
  }


  ## Durations ----

  # overall
  durs.casl <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ongoing2 == 1) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()

  # create city weights
  if (!is.null(geog.lvl)) {
    durs.casl.geo <- lcasl %>%
      filter(RAI == 1 | IAI == 1) %>%
      filter(ongoing2 == 1) %>%
      filter(geogYN == 1) %>%
      summarise(mean.dur = mean(duration.time, na.rm = TRUE),
                median.dur = median(duration.time, na.rm = TRUE)) %>%
      as.data.frame()

    # city-specific weight based on ratio of medians
    wt <- durs.casl.geo$median.dur / durs.casl$median.dur
  } else {
    wt <- 1
  }

  # The dissolution rate is function of the mean of the geometric distribution
  # which relates to the median as:
  durs.casl$rates.casl.adj <- 1 - (2^(-1 / (wt * durs.casl$median.dur)))

  # Mean duration associated with a geometric distribution that median:
  durs.casl$mean.dur.adj <- 1 / (1 - (2^(-1 / (wt * durs.casl$median.dur))))
  out$casl$durs.casl.homog <- durs.casl

  # stratified by age

  # first, non-matched by age group
  durs.casl.nonmatch <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 0) %>%
    # group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()
  durs.casl.nonmatch$index.age.grp <- 0

  # then, matched within age-groups
  durs.casl.matched <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(same.age.grp == 1) %>%
    filter(ongoing2 == 1) %>%
    group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration.time, na.rm = TRUE),
              median.dur = median(duration.time, na.rm = TRUE)) %>%
    as.data.frame()

  durs.casl.all <- rbind(durs.casl.nonmatch, durs.casl.matched)

  durs.casl.all$rates.casl.adj <- 1 - (2^(-1 / (wt * durs.casl.all$median.dur)))
  durs.casl.all$mean.dur.adj <- 1 / (1 - (2^(-1 / (wt * durs.casl.all$median.dur))))

  durs.casl.all <- durs.casl.all[, c(3, 1, 2, 4, 5)]
  out$casl$durs.casl.byage <- durs.casl.all


  ## duration.method override (see #63 phase 3) ----
  if (duration.method != "empirical") {
    alt <- compute_alt_durs(
      l_layer = lcasl[lcasl$RAI == 1 | lcasl$IAI == 1, , drop = FALSE],
      duration.method = duration.method,
      byage_strata = out$casl$durs.casl.byage$index.age.grp,
      race = race, geog.lvl = geog.lvl
    )
    if (!is.na(alt$homog$median.dur)) {
      out$casl$durs.casl.homog$mean.dur   <- alt$homog$mean.dur
      out$casl$durs.casl.homog$median.dur <- alt$homog$median.dur
      out$casl$durs.casl.homog$rates.casl.adj <-
        1 - 2^(-1 / (wt * out$casl$durs.casl.homog$median.dur))
      out$casl$durs.casl.homog$mean.dur.adj <-
        1 / (1 - 2^(-1 / (wt * out$casl$durs.casl.homog$median.dur)))
    }
    for (i in seq_len(nrow(out$casl$durs.casl.byage))) {
      k <- out$casl$durs.casl.byage$index.age.grp[i]
      j <- which(alt$byage$index.age.grp == k)
      if (length(j) == 1 && !is.na(alt$byage$median.dur[j])) {
        out$casl$durs.casl.byage$mean.dur[i]        <- alt$byage$mean.dur[j]
        out$casl$durs.casl.byage$median.dur[i]      <- alt$byage$median.dur[j]
        out$casl$durs.casl.byage$rates.casl.adj[i]  <-
          1 - 2^(-1 / (wt * alt$byage$median.dur[j]))
        out$casl$durs.casl.byage$mean.dur.adj[i] <-
          1 / (1 - 2^(-1 / (wt * alt$byage$median.dur[j])))
      }
    }
    if (duration.method == "weibull_strat") {
      attr(out$casl$durs.casl.byage, "weibull_shape") <- alt$weibull_shapes
      attr(out$casl$durs.casl.homog, "weibull_shape") <- alt$weibull_shape_overall
    }
    if (duration.method == "joint_lm" && !is.null(alt$model)) {
      out$casl$joint_duration_model <- alt$model
    }
  }


  # If sexual cessation model, then set diss coef for age grp above boundary to 1
  if (sex.cess.mod == TRUE) {
    index.age.grp <- max(out$casl$durs.casl.byage$index.age.grp) + 1
    df <- data.frame(index.age.grp = index.age.grp, mean.dur = 1, median.dur = 1,
                     rates.casl.adj = 1, mean.dur.adj = 1)
    out$casl$durs.casl.byage <- rbind(out$casl$durs.casl.byage, df)
  }


  ## joint g-computation models (additive outputs; see issues #61/#62/#63) ----
  if (method == "joint") {
    out$casl$joint_model <- fit_joint_glm(
      d,
      response = "deg.casl",
      main_terms = c("age.grp", "sqrt(age.grp)", "deg.main", "hiv2"),
      race = race,
      geog.lvl = geog.lvl,
      interaction_cross_deg = "deg.main",
      family = poisson()
    )
    out$casl$joint_concurrent_model <- fit_joint_glm(
      d,
      response = "deg.casl.conc",
      main_terms = c("age.grp", "sqrt(age.grp)", "deg.main", "hiv2"),
      race = race,
      geog.lvl = geog.lvl,
      interaction_cross_deg = "deg.main",
      family = binomial()
    )

    # Dyad-level (nodematch / absdiff) -- see #63.
    lcasl$age.grp <- lcasl$index.age.grp
    out$casl$joint_nm_age_model <- fit_joint_glm(
      lcasl, response = "same.age.grp",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = binomial()
    )
    if (isTRUE(race)) {
      out$casl$joint_nm_race_model <- fit_joint_glm(
        lcasl, response = "same.race",
        main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
        race = race, geog.lvl = geog.lvl,
        family = binomial()
      )
    }
    out$casl$joint_absdiff_age_model <- fit_joint_glm(
      lcasl, response = "ad",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
    out$casl$joint_absdiff_sqrtage_model <- fit_joint_glm(
      lcasl, response = "ad.sr",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
  }


  # 3. One-off Model --------------------------------------------------------

  out$inst <- list()
  linst <- l[l$ptype == 3, ]

  ## edges ----

  head(d$count.oo.part, 25)

  # rate by time unit
  d$rate.oo.part <- d$count.oo.part / (364 / time.unit)

  if (is.null(geog.lvl)) {
    mod <- glm(count.oo.part ~ 1,
               data = d, family = poisson())

    pred <- exp(coef(mod)[[1]]) / (364 / time.unit)

    out$inst$md.inst <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geogYN,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$md.inst <- as.numeric(pred)
  }


  ## nodematch("age.grp") ----

  linst$index.age.grp <- cut(linst$age, age.breaks, labels = FALSE,
                             right = FALSE, include.lowest = FALSE)
  linst$part.age.grp <- cut(as.numeric(linst$p_age_imp), age.breaks, labels = FALSE,
                            right = FALSE, include.lowest = FALSE)

  linst$same.age.grp <- as.integer(linst$index.age.grp == linst$part.age.grp)

  if (is.null(geog.lvl)) {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = linst, family = binomial())

    dat <- data.frame(index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geogYN + index.age.grp,
               data = linst, family = binomial())

    dat <- data.frame(geogYN = 1, index.age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  linst$ad <- abs(linst$age - linst$p_age_imp)
  linst$ad.sr <- abs(sqrt(linst$age) - sqrt(linst$p_age_imp))

  if (is.null(geog.lvl)) {
    mod <- lm(ad ~ 1, data = linst)

    pred <- coef(mod)[[1]]

    out$inst$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geogYN, data = linst)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (is.null(geog.lvl)) {
    mod <- lm(ad.sr ~ 1, data = linst)

    pred <- coef(mod)[[1]]

    out$inst$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geogYN, data = linst)

    dat <- data.frame(geogYN = 1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE,
                   right = FALSE, include.lowest = FALSE)

  if (is.null(geog.lvl)) {
    mod <- glm(count.oo.part ~ age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geogYN + age.grp + sqrt(age.grp),
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, age.grp = 1:age.grps)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.age.grp <- as.numeric(pred)
  }


  if (race == TRUE) {

    ## nodematch("race", diff = TRUE) ----

    linst$same.race <- as.integer(linst$race.cat.num == linst$p_race.cat.num)

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ as.factor(race.cat.num),
                 data = linst, family = binomial())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN + as.factor(race.cat.num),
                 data = linst, family = binomial())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (is.null(geog.lvl)) {
      mod <- glm(same.race ~ 1,
                 data = linst, family = binomial())

      pred <- exp(coef(mod)[[1]]) / (1 + exp(coef(mod)[[1]]))

      out$inst$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geogYN,
                 data = linst, family = binomial())

      dat <- data.frame(geogYN = 1)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race_diffF <- as.numeric(pred)
    }


    ## nodefactor("race") ----

    if (is.null(geog.lvl)) {
      mod <- glm(count.oo.part ~ as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

      out$inst$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(count.oo.part ~ geogYN + as.factor(race.cat.num),
                 data = d, family = poisson())

      dat <- data.frame(geogYN = 1, race.cat.num = race.categories)
      pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

      out$inst$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("risk.grp") ----

  # geography-specific wts

  if (!is.null(geog.lvl)) {
    wt <- mean(d$rate.oo.part[d$geogYN == 1], na.rm = TRUE) /
      mean(d$rate.oo.part, na.rm = TRUE)
  } else {
    wt <- 1
  }
  wt.rate <- d$rate.oo.part * wt

  nquants <- oo.nquants
  oo.quants <- rep(NA, nquants)
  sr <- sort(wt.rate)
  qsize <- floor(length(sr) / nquants)
  for (i in 1:nquants) {
    if (i == 1) {
      oo.quants[i] <- mean(sr[1:qsize])
    } else if (i > 1 && i < nquants) {
      oo.quants[i] <- mean(sr[(((i - 1) * qsize) + 1):(i * qsize)])
    } else if (i == nquants) {
      oo.quants[i] <- mean(sr[(((i - 1) * qsize) + 1):length(sr)])
    }
  }

  # Save it
  out$inst$nf.risk.grp <- oo.quants


  ## nodefactor("deg.tot") ----

  d$deg.tot3 <- ifelse(d$deg.tot >= 3, 3, d$deg.tot)

  if (!is.null(geog.lvl)) {
    deg.tot.dist <- prop.table(table(d$deg.tot3[d$geogYN == 1]))
    out$inst$deg.tot.dist <- as.numeric(deg.tot.dist)
  } else {
    deg.tot.dist <- prop.table(table(d$deg.tot3))
    out$inst$deg.tot.dist <- as.numeric(deg.tot.dist)
  }

  if (is.null(geog.lvl)) {
    mod <- glm(count.oo.part ~ deg.tot3 + sqrt(deg.tot3),
               data = d, family = poisson())

    dat <- data.frame(deg.tot3 = 0:3)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.deg.tot <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geogYN + deg.tot3 + sqrt(deg.tot3),
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, deg.tot3 = 0:3)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.deg.tot <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (is.null(geog.lvl)) {
    mod <- glm(count.oo.part ~ hiv2,
               data = d, family = poisson())

    dat <- data.frame(hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geogYN + hiv2,
               data = d, family = poisson())

    dat <- data.frame(geogYN = 1, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response") / (364 / time.unit)

    out$inst$nf.diag.status <- as.numeric(pred)
  }


  ## joint g-computation models (additive outputs; see issues #61/#62/#63) ----
  # No concurrent target for the one-off layer, so no binomial concurrent
  # model here. Dyad-level nodematch/absdiff still apply.
  if (method == "joint") {
    out$inst$joint_model <- fit_joint_glm(
      d,
      response = "count.oo.part",
      main_terms = c("age.grp", "sqrt(age.grp)", "deg.tot3", "sqrt(deg.tot3)",
                     "hiv2"),
      race = race,
      geog.lvl = geog.lvl,
      interaction_cross_deg = "deg.tot3",
      family = poisson()
    )

    # Dyad-level (nodematch / absdiff) -- see #63.
    linst$age.grp <- linst$index.age.grp
    out$inst$joint_nm_age_model <- fit_joint_glm(
      linst, response = "same.age.grp",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = binomial()
    )
    if (isTRUE(race)) {
      out$inst$joint_nm_race_model <- fit_joint_glm(
        linst, response = "same.race",
        main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
        race = race, geog.lvl = geog.lvl,
        family = binomial()
      )
    }
    out$inst$joint_absdiff_age_model <- fit_joint_glm(
      linst, response = "ad",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
    out$inst$joint_absdiff_sqrtage_model <- fit_joint_glm(
      linst, response = "ad.sr",
      main_terms = c("age.grp", "sqrt(age.grp)", "hiv2"),
      race = race, geog.lvl = geog.lvl,
      family = gaussian()
    )
  }



  # 4. Other Parameters -----------------------------------------------------

  ## Sexual Role ----

  d <- l %>%
    filter(RAI == 1) %>%
    group_by(AMIS_ID) %>%
    count() %>%
    rename(nRAIpart = n) %>%
    right_join(d, by = "AMIS_ID") %>%
    as.data.frame()
  d$nRAIpart <- ifelse(is.na(d$nRAIpart), 0, d$nRAIpart)

  d <- l %>%
    filter(IAI == 1) %>%
    group_by(AMIS_ID) %>%
    count() %>%
    rename(nIAIpart = n) %>%
    right_join(d, by = "AMIS_ID") %>%
    as.data.frame()
  d$nIAIpart <- ifelse(is.na(d$nIAIpart), 0, d$nIAIpart)

  # default NA for no AI
  roletype <- rep(NA, nrow(d))
  roletype[d$nRAIpart == 0 & d$nIAIpart > 0] <- 0
  roletype[d$nIAIpart == 0 & d$nRAIpart > 0] <- 1
  roletype[d$nIAIpart > 0 & d$nRAIpart > 0] <- 2

  out$all$role.type <- prop.table(table(roletype))



  # SAVE --------------------------------------------------------------------

  return(out)
}
