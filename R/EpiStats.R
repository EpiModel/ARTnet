
#' @title Statistical Models for Act Rates, Condom Use, and Starting HIV Prevalence
#'
#' @description Builds statistical models governing act rates and probability of
#'              condom use among main, casual and one-time sexual partnerships, and
#'              the probability of diagnosed HIV infection, for use in the EpiModelHIV
#'              workflow.
#'
#' @param geog.lvl Specifies geographic level for ARTnet statistics.
#' @param geog.cat Specifies one or more geographic strata within the level to base
#'        ARTnet statistics on. If the vector is of length 2+, data from the strata
#'        will be combined into one analysis.
#' @param age.limits Upper and lower limit. Age range to subset ARTnet data by.
#'        Default is 15 to 65.
#' @param age.breaks Ages that define the upper age categories. Default is
#'        \code{c(25, 35, 45, 55, 65)}, which corresponds to (0, 25], (25, 35],
#'        (35, 45], (45, 55], (55, 65], (65, 100].
#' @param race Whether to stratify by racial status. Default is TRUE.
#' @param browser Run function in interactive browser mode. Default is FALSE.
#' @param init.hiv.prev Initial HIV prevalence of estimated model, vector of size
#'         three pertaining to prevalence among three racial classes (black,
#'         Hispanic and white respectively). If \code{init.hiv.prev = NULL},
#'         ARTnet will handle calculation of prevalence through ARTnet data.
#'         Note: if `\code{race = FALSE}, prevalence vector must still be supplied.
#' @param time.unit Specifies time unit for ARTnet statistics. Default is 7 for
#'         weekly time units. Inputs are by days with a maximum of 30 days.
#'         Use 1 for daily time units and 30 for monthly time units.
#'
#' @details
#' \code{build_epistats}, through input of geographic, age and racial
#' parameters, builds the necessary statistical models predicting sexual activity
#' and condom use during sexual activity in main, casual and one-off partnerships
#' among men-who-have-sex-with-men (MSM). Estimation of these linear models is
#' done using data from the ARTnetData package, a package containing the
#' results of the online ARTnet survey of HIV-related risk behaviors, testing
#' and use of preventive services among MSM in the United States. Accepted
#' values for each input parameter are provided below.
#'
#' @section Parameter Values:
#' \itemize{
#'   \item \code{geog.lvl}: level of geographic stratification desired. Acceptable
#'          values are \code{"city"}, \code{"county"}, \code{"state"}, \code{"region"}, and
#'          \code{"division"} corresponding to the metropolitan statistical area,
#'          county, state, census region, and census
#'          division, respectively. Default value is "NULL",
#'          indicating no geographic stratification.
#'    \item \code{geog.cat}: given a geographic level above, \code{"geog.cat"}
#'          is a vector comprising the desired feature(s) of interest. Acceptable values
#'          are based on the chosen geographic level:
#'    \itemize{
#'      \item \code{city}: \code{"Atlanta"}, \code{"Boston"}, \code{"Chicago"},
#'            \code{"Dallas"}, \code{"Denver"}, \code{"Detroit"}, \code{"Houston"},
#'            \code{"Los Angeles"}, \code{"Miami"}, \code{"New York City"},
#'            \code{"Philadelphia"}, \code{"San Diego"}, \code{"San Franciso"},
#'            \code{"Seattle"}, \code{"Washington DC"}
#'      \item \code{county}: FIPS codes for the county or county equivalents to
#'      be included. Be aware that selecting a single county or a set of smaller
#'      counties may lead to insufficient data being used to estimate the model;
#'      if the sample lacks sufficient behavioral diversity (e.g. if nobody
#'      reports 2+ ongoing casual partners), an error may result.
#'      \item \code{state}: \code{"AK"}, \code{"AL"}, \code{"AR"}, \code{"AZ"},
#'            \code{"CA"}, \code{"CO"}, \code{"CT"}, \code{"DC"}, \code{"DE"},
#'            \code{"FL"}, \code{"GA"}, \code{"HI"}, \code{"IA"}, \code{"ID"},
#'            \code{"IL"}, \code{"IN"}, \code{"KS"}, \code{"KY"}, \code{"LA"},
#'            \code{"MA"}, \code{"MD"}, \code{"ME"}, \code{"MI"}, \code{"MN"},
#'            \code{"MO"}, \code{"MS"}, \code{"MT"}, \code{"NC"}, \code{"ND"},
#'            \code{"NE"}, \code{"NH"}, \code{"NJ"}, \code{"NM"}, \code{"NV"},
#'            \code{"NY"}, \code{"OH"}, \code{"OK"}, \code{"OR"}, \code{"PA"},
#'            \code{"RI"}, \code{"SC"}, \code{"SD"}, \code{"TN"}, \code{"TX"},
#'            \code{"UT"}, \code{"VA"}, \code{"VT"}, \code{"WA"}, \code{"WI"},
#'            \code{"WV"}, \code{"WY"}
#'      \item \code{division}: \code{"1"} (New England),
#'            \code{"2"} (Middle Atlantic), \code{"3"} (East North Central),
#'            \code{"4"} (West North Central), \code{"5"} (South Atlantic),
#'            \code{"6"} (East South Central), \code{"7"} (West South Central),
#'            \code{"8"} (Mountain) \code{"9"} (Pacific)
#'      \item \code{region}: \code{"1"} (Northeast), \code{"2"} (Midwest),
#'            \code{"3"} (South), \code{"4"} (North)
#'    }
#'    \item \code{race}: whether to introduce modeling by racial stratification.
#'          \code{TRUE} or \code{FALSE}.
#'    \item \code{age.limits}: a vector giving the lower and upper limit for
#'          the age of interest. Set to \code{c(15, 65)} by default.
#'    \item \code{age.breaks}: a vector giving the upper age breaks to categorize
#'          data by age. Must be within the bounds specified by \code{age.limits}.
#'    \item \code{time.unit}: a number between 1 and 30 that specifies time units
#'          for ARTnet statistics. Set to \code{7} by default.
#' }
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

#' # No race stratification
#' epistats4 <- build_epistats(geog.lvl = "state", geog.cat = "GA",
#'                             race = FALSE)
#'
#' # Age and race stratification, for the municipality (not metro) of New York City
#' # geog.cat values are FIPS codes for the 5 boroughs of NYC
#' epistats5 <- build_epistats(geog.lvl = "county",
#'                             geog.cat = c(36005, 36047, 36061, 36081, 36085),
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(24, 34, 44))
#'
#' @export
#'
build_epistats <- function(geog.lvl = NULL, geog.cat = NULL, race = FALSE,
                           age.limits = c(15, 65), age.breaks = c(25, 35, 45, 55),
                           init.hiv.prev = NULL, time.unit = 7, browser = FALSE) {

  # Fix global binding check errors
  duration.time <- anal.acts.time <- anal.acts.time.cp <- NULL

  if (browser == TRUE) {
    browser()
  }


  ## Data ##
  d <- ARTnet.wide
  l <- ARTnet.long

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



  # Age Processing

  # Subset data by selected age range


  # Warning if age range is out of allowed range
  flag.ll <- age.limits[1] >= 15 & age.limits[1] <= 65
  flag.ul <- age.limits[2] >= 15 & age.limits[2] <= 65
  flag.lim <- flag.ll * flag.ul

  if (flag.lim == 0) {
    stop("Age range must be between 15 and 65")
  }

  age.limits <- c(min(age.limits), max(age.limits))

  flag.bks <- prod(age.breaks < age.limits[2] & age.breaks >= age.limits[1])

  if (flag.bks == 0) {
    stop("Age breaks must be between specified age limits")
  }

  age.breaks <- unique(sort(c(age.limits[1], age.breaks, 100)))

  l <- subset(l, age >= age.limits[1] & age <= age.limits[2])
  d <- subset(d, age >= age.limits[1] & age <= age.limits[2])

  l$comb.age <- l$age + l$p_age_imp
  l$diff.age <- abs(l$age - l$p_age_imp)

  if (race == TRUE) {
    # Race
    d$race.cat3 <- rep(NA, nrow(d))
    d$race.cat3[d$race.cat == "black"] <- 1
    d$race.cat3[d$race.cat == "hispanic"] <- 2
    d$race.cat3[d$race.cat %in% c("white", "other")] <- 3

    l$race.cat3 <- rep(NA, nrow(l))
    l$race.cat3[l$race.cat == "black"] <- 1
    l$race.cat3[l$race.cat == "hispanic"] <- 2
    l$race.cat3[l$race.cat %in% c("white", "other")] <- 3

    l$p_race.cat3 <- rep(NA, nrow(l))
    l$p_race.cat3[l$p_race.cat == "black"] <- 1
    l$p_race.cat3[l$p_race.cat == "hispanic"] <- 2
    l$p_race.cat3[l$p_race.cat %in% c("white", "other")] <- 3

    # redistribute NAs in proportion to non-missing partner races
    probs <- prop.table(table(l$race.cat3, l$p_race.cat3), 1)

    imp_black <- which(is.na(l$p_race.cat3) & l$race.cat3 == 1)
    l$p_race.cat3[imp_black] <- sample(1:3, length(imp_black), TRUE, probs[1, ])

    imp_hisp <- which(is.na(l$p_race.cat3) & l$race.cat3 == 2)
    l$p_race.cat3[imp_hisp] <- sample(1:3, length(imp_hisp), TRUE, probs[2, ])

    imp_white <- which(is.na(l$p_race.cat3) & l$race.cat3 == 3)
    l$p_race.cat3[imp_white] <- sample(1:3, length(imp_white), TRUE, probs[3, ])

    l$race.combo <- rep(NA, nrow(l))
    l$race.combo[l$race.cat3 == 1 & l$p_race.cat3 == 1] <- 1
    l$race.combo[l$race.cat3 == 1 & l$p_race.cat3 %in% 2:3] <- 2
    l$race.combo[l$race.cat3 == 2 & l$p_race.cat3 %in% c(1, 3)] <- 3
    l$race.combo[l$race.cat3 == 2 & l$p_race.cat3 == 2] <- 4
    l$race.combo[l$race.cat3 == 3 & l$p_race.cat3 %in% 1:2] <- 5
    l$race.combo[l$race.cat3 == 3 & l$p_race.cat3 == 3] <- 6

    l <- select(l, -c(race.cat3, p_race.cat3))
  }


  # HIV
  l$p_hiv2 <- ifelse(l$p_hiv == 1, 1, 0)

  hiv.combo <- rep(NA, nrow(l))
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5
  l$hiv.concord.pos <- ifelse(hiv.combo == 2, 1, 0)

  # PrEP
  d$prep <- ifelse(d$artnetPREP_CURRENT == 0 | is.na(d$artnetPREP_CURRENT), 0, 1)

  dlim <- select(d, c(AMIS_ID, survey.year, prep))
  l <- left_join(l, dlim, by = "AMIS_ID")

  # Time unit processing

  ## Set time.unit limits from 1 to 30

  if (time.unit < 1 || time.unit > 30) {
    stop("Time unit must be between 1 and 30")
  }

  ## Create time-based ARTnet statistics using time.unit

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
        d1 <- select(d, race.cat3, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age + as.factor(race.cat3),
                       data = d1, family = binomial())
      } else {
        d1 <- select(d, race.cat3, geogYN, age, hiv2)
        hiv.mod <- glm(hiv2 ~ age + geogYN + as.factor(race.cat3) + geogYN * as.factor(race.cat3),
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
    if (length(init.hiv.prev) != 3) {
      stop("Input parameter init.prev.hiv must be a vector of size three")
    }
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
  out$acts.mod <- acts.mod
  out$cond.mc.mod <- cond.mc.mod
  out$cond.oo.mod <- cond.oo.mod
  out$geog.l <- as.character(l$geog)
  out$geog.d <- as.character(d$geog)
  out$age.limits <- age.limits
  out$age.breaks <- age.breaks
  out$age.grps <- length(age.breaks) - 1
  out$init.hiv.prev <- init.hiv.prev
  out$time.unit <- time.unit
  return(out)
}
