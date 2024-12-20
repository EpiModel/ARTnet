
#' @title Statistical Models for Act Rates, Condom Use, and Starting HIV Prevalence
#'
#' @description Builds statistical models governing act rates and probability of condom use among
#' main, casual and one-time sexual partnerships, and the probability of diagnosed HIV infection,
#' for use in the EpiModelHIV workflow.
#'
#' @param geog.lvl Specifies geographic level for ARTnet statistics.
#' @param geog.cat Specifies one or more geographic strata within the level to base ARTnet
#'        statistics on. If the vector is of length 2+, data from the strata will be combined into
#'        one analysis.
#' @param race If `TRUE`, stratify model estimates by race/ethnic grouping.
#' @param age.limits Lower and upper limit of age range to include in model. Minimum of 15 and
#'        maximum of 100 allowed. Lower limit is inclusive boundary and upper boundary is
#'        exclusive boundary.
#' @param age.breaks Ages that define the upper closed boundary of the age categories. Default is
#'        `c(25, 35, 45, 55)`, which corresponds to `(0, 25], (25, 35], (35, 45], (45, 55], (55, 65]`
#'        with `age.limits = c(15, 65)`.
#' @param age.sexual.cessation Age of cessation of sexual activity, while aging process continues
#'        through the upper age limit. Maximum allowed value of 66.
#' @param init.hiv.prev Initial HIV prevalence to be used in epidemic model estimated model, with a
#'        numerical vector of size 3 corresponding to starting prevalence in three race/ethnic
#'        groups (Black, Hispanic, and White/Other, respectively). If `init.hiv.prev = NULL`,
#'        `build_epistats` will estimate a logistic regression model to predict starting prevalence
#'        as a function of estimated prevalence in ARTnet as a function of race/ethnicity and age.
#' @param time.unit Specifies time unit for time-dependent ARTnet statistics. Default is 7,
#'        corresponding to a weekly time unit. Allowed inputs range from 1 for a daily time unit to
#'        30 for a monthly time unit.
#' @param browser If `TRUE`, run `build_epistats` in interactive browser mode.
#'
#' @details
#' The `build_epistats` function provides a way to input of geographic, age, and racial parameters
#' necessary to build statistical models predicting sexual activity and condom use in sexual
#' partnerships among men who have sex with men (MSM). Estimation of these models is performed using
#' data from the `ARTnetData` package, containing the results of the ARTnet study.
#'
#' ## Explanation of Parameter Values
#' * `geog.lvl`: level of geographic stratification desired. Acceptable values are `"city"`,
#'   `"county"`, `"state"`, `"region"`, and `"division"` corresponding to the metropolitan
#'   statistical area, county, state, census region, and census division, respectively. Default
#'   value is `NULL`, indicating no geographic stratification.
#' * `geog.cat`: given a geographic level above, `"geog.cat"` is a vector comprising the desired
#'   feature(s) of interest. Acceptable values are based on the chosen geographic level:
#'     - `city`: `"Atlanta"`, `"Boston"`, `"Chicago"`, `"Dallas"`, `"Denver"`, `"Detroit"`,
#'       `"Houston"`, `"Los Angeles"`, `"Miami"`, `"New York City"`, `"Philadelphia"`, `"San Diego"`,
#'       `"San Francisco"`, `"Seattle"`, `"Washington DC"`.
#'     - `county`: FIPS codes for the county or county equivalents to be included. Selecting a
#'       single county or a set of smaller counties may lead to insufficient data used to estimate
#'       the models; an error may result.
#'     - `state`: Two-letter postal code for each state.
#'     - `division`: `"1"` (New England), `"2"` (Middle Atlantic), `"3"` (East North Central),
#'       `"4"` (West North Central), `"5"` (South Atlantic), `"6"` (East South Central),
#'       `"7"` (West South Central), `"8"` (Mountain) `"9"` (Pacific).
#'     - `region`: `"1"` (Northeast), `"2"` (Midwest), `"3"` (South), `"4"` (North)
#' * `race`: whether to introduce modeling by racial stratification. `TRUE` or `FALSE`.
#' * `age.limits`: a vector giving the lower and upper limit for the age of interest. Set to
#'   `c(15, 65)` by default. The lower boundary is inclusive, meaning persons may be initialized into
#'   the model at age 15.0; the upper boundary is exclusive, meaning exit from the population will
#'   occur on turning 65.0 years. Although the ARTnet data include respondents from age 15 to 65, this
#'   may be set to restricted values within that range (for example, `c(25, 40)`) for a subsetted
#'   data analysis. The upper boundary may be set up to age 100 (that is, `c(15, 100)`), which may
#'   be used in models where sexual activity ceases before mortality.
#' * `age.breaks`: a vector giving the upper age breaks to categorize data by age. Must be within
#'   the bounds specified by `age.limits`. These should be the interior age breaks only (that is,
#'   there is no need to include the age limit boundaries). If an age of sexual cessation is added,
#'   then this age is also added to the age breaks if it is not explicitly specified.
#' * `age.sexual.cessation`: a numerical value for the age of cessation of sexual activity. This may
#'   be by assumption or given data constraints (ARTnet eligibility were through age 65, so the
#'   maximum value here is 66). This specification is useful for models in which the HIV transmission
#'   process stops at a certain age but aging and other demographic features should continue through
#'   natural mortality.
#' * `time.unit`: a number between 1 and 30 that specifies time units for ARTnet statistics. Set to
#'   `7` by default.
#' * `race.level`:
#'
#' @examples
#' # Age and geographic stratification, for the Atlanta metropolitan statistical area
#' epistats1 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(24, 34, 44))
#'
#' # Default age stratification
#' epistats2 <- build_epistats(geog.lvl = "state", geog.cat = "WA")
#'
#' # Default age stratification, multiple states
#' epistats3 <- build_epistats(geog.lvl = "state", geog.cat = c("ME", "NH", "VT"))
#'
#' # No race stratification
#' epistats4 <- build_epistats(geog.lvl = "state", geog.cat = "GA", race = FALSE)
#'
#' # Age and race stratification, for the municipality (not metro) of New York City
#' # geog.cat values are FIPS codes for the 5 boroughs of NYC
#' epistats5 <- build_epistats(geog.lvl = "county",
#'                             geog.cat = c(36005, 36047, 36061, 36081, 36085),
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(24, 34, 44))
#'
#' # Use broader age range (to age 100) but with sexual cessation at age 66
#' epistats6 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 66)
#'
#' @export
#'
build_epistats <- function(geog.lvl = NULL,
                           geog.cat = NULL,
                           race = TRUE,
                           race.level = list("black", "hispanic", c("white", "other")),
                           age.limits = c(15, 65),
                           age.breaks = c(25, 35, 45, 55),
                           age.sexual.cessation = NULL,
                           init.hiv.prev = NULL,
                           time.unit = 7,
                           browser = FALSE) {
  # Ensures that ARTnetData is installed
  if (system.file(package = "ARTnetData") == "") stop(missing_data_msg)

  # Fix global binding check errors
  duration.time <- anal.acts.time <- anal.acts.time.cp <- NULL

  if (browser == TRUE) {
    browser()
  }


  ## Data ##
  d <- ARTnetData::ARTnet.wide
  l <- ARTnetData::ARTnet.long

  out <- list()

  geog_names <- c("city", "county", "state", "region", "division")
  if (!is.null(geog.lvl)) {
    if (!(geog.lvl %in% geog_names)) {
      stop("Selected geographic level must be one of: city, county, state, region or division")
    }
  }

  # Data Processing ---------------------------------------------------------

  # Geography
  if (length(geog.lvl) > 1) {
    stop("Only one geographical level may be chosen at a time.")
  }

  if (!is.null(geog.lvl)) {
    if (geog.lvl == "city") {
      if (sum(geog.cat %in% unique(d$city)) == 0) {
        stop("None of the city names found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "city2")]))
      l$geogYN <- ifelse(l[, "city2"] %in% geog.cat, 1, 0)
      l$geog <- l$city2
      d$geogYN <- ifelse(d[, "city2"] %in% geog.cat, 1, 0)
      d$geog <- d$city2
    }

    if (geog.lvl == "county") {
      if (sum(geog.cat %in% unique(d$COUNTYFIPS)) == 0) {
         stop("None of the county FIPS codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "COUNTYFIPS")]))
      l$geogYN <- ifelse(l[, "COUNTYFIPS"] %in% geog.cat, 1, 0)
      l$geog <- l$COUNTYFIPS
      d$geogYN <- ifelse(d[, "COUNTYFIPS"] %in% geog.cat, 1, 0)
      d$geog <- d$COUNTYFIPS
    }

    if (geog.lvl == "state") {
      if (sum(geog.cat %in% unique(d$State)) == 0) {
        stop("None of the states found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "State")]))
      l$geogYN <- ifelse(l[, "State"] %in% geog.cat, 1, 0)
      l$geog <- l$State
      d$geogYN <- ifelse(d[, "State"] %in% geog.cat, 1, 0)
      d$geog <- d$State
    }

    if (geog.lvl == "division") {
      if (sum(geog.cat %in% unique(d$DIVCODE)) == 0) {
        stop("None of the census division codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "DIVCODE")]))
      l$geogYN <- ifelse(l[, "DIVCODE"] %in% geog.cat, 1, 0)
      l$geog <- l$DIVCODE
      d$geogYN <- ifelse(d[, "DIVCODE"] %in% geog.cat, 1, 0)
      d$geog <- d$DIVCODE
    }

    if (geog.lvl == "region") {
      if (sum(geog.cat %in% unique(d$REGCODE)) == 0) {
        stop("None of the census region codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "REGCODE")]))
      l$geogYN <- ifelse(l[, "REGCODE"] %in% geog.cat, 1, 0)
      l$geog <- l$REGCODE
      d$geogYN <- ifelse(d[, "REGCODE"] %in% geog.cat, 1, 0)
      d$geog <- d$REGCODE
    }
  }


  ## Age Processing ##
  if (length(age.limits) != 2 || age.limits[1] > age.limits[2]) {
    stop("age.limits must be a vector of length 2, where age.limits[2] > age.limits[1]")
  }

  # Warning if age range is out of allowed range
  flag.ll <- age.limits[1] >= 15 & age.limits[1] <= 100
  flag.ul <- age.limits[2] >= 15 & age.limits[2] <= 100
  flag.lim <- flag.ll * flag.ul
  if (flag.lim == 0) {
    stop("Age range specified in `age.limits` must be >= 15 and <= 100")
  }

  # Warning if age breaks fall outside age limits
  flag.bks <- prod(age.breaks < age.limits[2] & age.breaks >= age.limits[1])
  if (flag.bks == 0) {
    stop("Age breaks must be between specified age limits")
  }

  # Set default age.sexual.cessation and error if > ARTnet data
  if (is.null(age.sexual.cessation)) {
    age.sexual.cessation <- age.limits[2]
  }
  if (age.sexual.cessation > 66) {
    stop("Maximum allowed age of sexual cessation is 66, corresponding to the upper age eligilibity
         criteria of 65 (inclusive) in ARTnet")
  }

  # Composite age.breaks are now union of age.limits, age.breaks, and age.sexual.cessation
  age.breaks <- unique(sort(c(age.limits[1], age.breaks, age.sexual.cessation, age.limits[2])))

  # p_age_imp initialization for lintr
  p_age_imp <- NULL

  # Subset datasets by lower age limit and age.sexual.cessation
  # Now applies to both index (respondents) and partners for long dataset
  l <- subset(l, age >= age.limits[1] & age < age.sexual.cessation &
                p_age_imp >= age.limits[1] & p_age_imp < age.sexual.cessation)
  d <- subset(d, age >= age.limits[1] & age < age.sexual.cessation)

  if (age.limits[2] > age.sexual.cessation) {
    sex.cess.mod <- TRUE
  } else {
    sex.cess.mod <- FALSE
  }

  # Calculate combine age of index and partners
  l$comb.age <- l$age + l$p_age_imp
  l$diff.age <- abs(l$age - l$p_age_imp)

  ## Race ethnicity ##

  if (race == TRUE) {
    mult_race_cat <- c("asian", "ai/an", "mult", "nh/pi")
    flat_race.level <- unlist(race.level)

    # Determine which variables to use in ARTnet
    if (any(flat_race.level %in% mult_race_cat)) {

      d$race.eth <- ifelse(d$hispan == 1, "hispanic", d$race)
      l <- merge(l, d[, c("AMIS_ID", "race.eth")], by = "AMIS_ID", all.x = TRUE)
      l$p_race.eth <- ifelse(l$p_hispan == 1, "hispanic", l$p_race2)

      p_race_var <- "p_race.eth"
      race_var <- "race.eth"
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

    # Initialize race.combo and assign combinations dynamically
    l$race.combo <- rep(NA, nrow(l))
    combo_index <- 1
    for (i in race.categories) {
      # Case 1: Same race as one combination
      l$race.combo[l$race.cat.num == i & l$p_race.cat.num == i] <- combo_index
      combo_index <- combo_index + 1

      # Case 2: Race compared with all other race groups
      l$race.combo[l$race.cat.num == i & l$p_race.cat.num %in% setdiff(race.categories, i)] <- combo_index
      combo_index <- combo_index + 1
    }
  }

  ## HIV diagnosed status of index and partners ##
  l$p_hiv2 <- ifelse(l$p_hiv == 1, 1, 0)

  hiv.combo <- rep(NA, nrow(l))
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5
  l$hiv.concord.pos <- ifelse(hiv.combo == 2, 1, 0)

  ## PrEP ##
  d$prep <- ifelse(d$artnetPREP_CURRENT == 0 | is.na(d$artnetPREP_CURRENT), 0, 1)

  dlim <- select(d, c(AMIS_ID, survey.year, prep))
  l <- left_join(l, dlim, by = "AMIS_ID")

  ## Time unit processing ##

  # Set time.unit limits from 1 to 30
  if (time.unit < 1 || time.unit > 30) {
    stop("time.unit must be between 1 and 30")
  }

  # Scale time-based ARTnet data by time.unit
  l$duration.time <- l$duration * 7 / time.unit
  l$anal.acts.time <- l$anal.acts.week * time.unit / 7
  l$anal.acts.time.cp <- l$anal.acts.week.cp * time.unit / 7


  # Act Rates ---------------------------------------------------------------

  # acts/per week/per partnership for main and casual partnerships

  # Pull Data
  if (race == TRUE) {
    if (is.null(geog.lvl)) {
      la <- select(l, ptype, duration.time, comb.age,
                   race.combo, RAI, IAI, hiv.concord.pos, prep,
                   acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
        filter(ptype %in% 1:2) %>%
        filter(RAI == 1 | IAI == 1)
      la <- select(la, -c(RAI, IAI))
    } else {
      la <- select(l, ptype, duration.time, comb.age, geogYN = geogYN,
                   race.combo, RAI, IAI, hiv.concord.pos, prep,
                   acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
        filter(ptype %in% 1:2) %>%
        filter(RAI == 1 | IAI == 1)
      la <- select(la, -c(RAI, IAI))
    }
  }  else {
    if (is.null(geog.lvl)) {
      la <- select(l, ptype, duration.time, comb.age,
                   RAI, IAI, hiv.concord.pos, prep,
                   acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
        filter(ptype %in% 1:2) %>%
        filter(RAI == 1 | IAI == 1)
      la <- select(la, -c(RAI, IAI))
    } else {
      la <- select(l, ptype, duration.time, comb.age, geogYN = geogYN,
                   RAI, IAI, hiv.concord.pos, prep,
                   acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
        filter(ptype %in% 1:2) %>%
        filter(RAI == 1 | IAI == 1)
      la <- select(la, -c(RAI, IAI))
    }
  }

  # Poisson Model
  if (race == TRUE) {
    if (is.null(geog.lvl)) {
      acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
                        as.factor(race.combo) + as.factor(ptype) + duration.time *
                        as.factor(ptype) + comb.age + I(comb.age^2) + hiv.concord.pos,
                      family = poisson(), data = la)
    } else {
      acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
                        as.factor(race.combo) + as.factor(ptype) +
                        duration.time * as.factor(ptype) + comb.age + I(comb.age^2) +
                        hiv.concord.pos + geogYN,
                      family = poisson(), data = la)
    }
  }  else {
    if (is.null(geog.lvl)) {
      acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
                        as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                        I(comb.age^2) + hiv.concord.pos,
                      family = poisson(), data = la)
    } else {
      acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
                        as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                        I(comb.age^2) + hiv.concord.pos + geogYN,
                      family = poisson(), data = la)
    }
  }

  # Condom Use // Main Casual -----------------------------------------------

  la$prob.cond <- la$cp.acts / la$acts
  la$any.cond <- ifelse(la$prob.cond > 0, 1, 0)
  la$never.cond <- ifelse(la$prob.cond == 0, 1, 0)

  if (race == TRUE) {
    if (is.null(geog.lvl)) {
      cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) + as.factor(race.combo) +
                           as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                           I(comb.age^2) + hiv.concord.pos + prep,
                         family = binomial(), data = la)
    } else {
      cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) + as.factor(race.combo) +
                           as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                           I(comb.age^2) + hiv.concord.pos + prep + geogYN,
                         family = binomial(), data = la)
    }
  }  else {
    if (is.null(geog.lvl)) {
      cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) +
                           as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                           I(comb.age^2) + hiv.concord.pos + prep,
                         family = binomial(), data = la)
    } else {
      cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) +
                           as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
                           I(comb.age ^ 2) + hiv.concord.pos + prep + geogYN,
                         family = binomial(), data = la)
    }
  }

  # Condom Use // Inst ------------------------------------------------------
  if (race == TRUE) {
    if (is.null(geog.lvl)) {
      lb <- select(l, ptype, comb.age,
                   race.combo, hiv.concord.pos, prep,
                   RAI, IAI, RECUAI, INSUAI) %>%
        filter(ptype == 3) %>%
        filter(RAI == 1 | IAI == 1)
    } else {
      lb <- select(l, ptype, comb.age, geogYN = geogYN,
                   race.combo, hiv.concord.pos, prep,
                   RAI, IAI, RECUAI, INSUAI) %>%
        filter(ptype == 3) %>%
        filter(RAI == 1 | IAI == 1)
    }
  } else {
    if (is.null(geog.lvl)) {
      lb <- select(l, ptype, comb.age,
                   hiv.concord.pos, prep,
                   RAI, IAI, RECUAI, INSUAI) %>%
        filter(ptype == 3) %>%
        filter(RAI == 1 | IAI == 1)
    } else {
      lb <- select(l, ptype, comb.age, geogYN = geogYN,
                   hiv.concord.pos, prep,
                   RAI, IAI, RECUAI, INSUAI) %>%
        filter(ptype == 3) %>%
        filter(RAI == 1 | IAI == 1)
    }
  }

  lb$prob.cond <- rep(NA, nrow(lb))
  lb$prob.cond[lb$RAI == 1 & lb$IAI == 0] <- lb$RECUAI[lb$RAI == 1 & lb$IAI == 0] /
    lb$RAI[lb$RAI == 1 & lb$IAI == 0]
  lb$prob.cond[lb$RAI == 0 & lb$IAI == 1] <- lb$INSUAI[lb$RAI == 0 & lb$IAI == 1] /
    lb$IAI[lb$RAI == 0 & lb$IAI == 1]
  lb$prob.cond[lb$RAI == 1 & lb$IAI == 1] <- (lb$RECUAI[lb$RAI == 1 & lb$IAI == 1] +
                                                lb$INSUAI[lb$RAI == 1 & lb$IAI == 1]) /
    (lb$RAI[lb$RAI == 1 & lb$IAI == 1] + lb$IAI[lb$RAI == 1 & lb$IAI == 1])
  lb$prob.cond[which(lb$prob.cond == 0.5)] <- 0
  lb$prob.cond[which(lb$prob.cond %in% c(88, 99, 44))] <- NA
  lb <- select(lb, -c(RAI, IAI, RECUAI, INSUAI))

  if (race == TRUE) {
    if (is.null(geog.lvl)) {
      cond.oo.mod <- glm(prob.cond ~ as.factor(race.combo) +
                           comb.age + I(comb.age^2) +
                           hiv.concord.pos + prep,
                         family = binomial(), data = lb)
    } else {
      cond.oo.mod <- glm(prob.cond ~ as.factor(race.combo) +
                           comb.age + I(comb.age^2) +
                           hiv.concord.pos + prep + geogYN,
                         family = binomial(), data = lb)
    }
  } else {
    if (is.null(geog.lvl)) {
      cond.oo.mod <- glm(prob.cond ~ comb.age + I(comb.age^2) +
                           hiv.concord.pos + prep,
                         family = binomial(), data = lb)
    } else {
      cond.oo.mod <- glm(prob.cond ~ comb.age + I(comb.age^2) +
                           hiv.concord.pos + prep + geogYN,
                         family = binomial(), data = lb)
    }
  }

  # Init HIV Status ---------------------------------------------------------
  if (is.null(init.hiv.prev)) {
    if (race == TRUE) {
      if (is.null(geog.lvl)) {
        d1 <- select(d, race.cat.num, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age + as.factor(race.cat.num),
                       data = d1, family = binomial())
      } else {
        d1 <- select(d, race.cat.num, geogYN, age, hiv2)
        hiv.mod <- glm(hiv2 ~ age + geogYN + as.factor(race.cat.num) + geogYN * as.factor(race.cat.num),
                       data = d1, family = binomial())
      }
    } else {
      if (is.null(geog.lvl)) {
        d1 <- select(d, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age,
                       data = d1, family = binomial())
      } else {
        d1 <- select(d, geogYN, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age + geogYN,
                       data = d1, family = binomial())
      }
    }
    # Output
    out$hiv.mod <- hiv.mod
  } else {
    #if (length(init.hiv.prev) != 3) {
    #  stop("Input parameter init.prev.hiv must be a vector of size three")
    #}
    if (prod(init.hiv.prev < 1) == 0  || prod(init.hiv.prev > 0) == 0) {
      stop("All elements of init.hiv.prev must be between 0 and 1 non-inclusive")
    }
  }

  # Save Out File -----------------------------------------------------------

  if (!is.null(geog.lvl)) {
    out$geogYN.l <- l$geogYN
    out$geogYN.d <- d$geogYN
    out$geog.cat  <- geog.cat
  }

  out$geog.lvl <- geog.lvl
  out$race <- race
  out$race.level <- race.level
  out$acts.mod <- acts.mod
  out$cond.mc.mod <- cond.mc.mod
  out$cond.oo.mod <- cond.oo.mod
  out$geog.l <- as.character(l$geog)
  out$geog.d <- as.character(d$geog)
  out$age.limits <- age.limits
  out$age.breaks <- age.breaks
  out$age.grps <- length(age.breaks) - 1
  out$age.sexual.cessation <- age.sexual.cessation
  out$sex.cess.mod <- sex.cess.mod
  out$init.hiv.prev <- init.hiv.prev
  out$time.unit <- time.unit
  return(out)
}

# strip a `glm` object from all its components, leaving only what is required
# for predicting from new data
strip_glm <- function(cm) {
  root_elts <- c("y", "model", "residuals", "fitted.values", "effects",
                 "linear.predictors", "weights", "prior.weights", "data")
  for (elt in root_elts) cm[[elt]] <- c()

  family_elts <- c("variance", "dev.resids", "aic", "validmu", "simulate")
  for (elt in family_elts) cm$family[[elt]] <- c()

  cm$qr$qr <- c()
  attr(cm$terms, ".Environment") <- c()
  attr(cm$formula, ".Environment") <- c()

  return(cm)
}

#' Reduces the Size of the Epistats Object by Trimming its Models
#'
#' The `epistats` object contains 3 GLM models. These R object take up a lot of
#' space and only a small subset of their functionalities is used by
#' EpiModelHIV. This function trims these models to save up space without
#' compromising EpiModelHIV functionalities.
#'
#' @param epistats the `epistats` object to be trimmed
#'
#' @return a trimmed `epistats`
#'
#' @export
#'
trim_epistats <- function(epistats) {
  epistats$acts.mod <- strip_glm(epistats$acts.mod)
  epistats$cond.mc.mod <- strip_glm(epistats$cond.mc.mod)
  epistats$cond.oo.mod <- strip_glm(epistats$cond.oo.mod)
  return(epistats)
}

