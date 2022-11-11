
#' Calculate Network Target Statistics
#'
#' @description Calculates the final target statistics for the network models by applying
#'              individual-level network statistics against the population size and structure, for
#'              use in the EpiModelHIV workflow.
#'
#' @param epistats Output from [`build_epistats`].
#' @param netparams Output from [`build_netparams`].
#' @param network.size Size of the starting network.
#' @param expect.mort Expected average mortality level to pass into [`dissolution_coefs`] function.
#' @param age.pyramid Numerical vector of length equal to the length of the age range specified in
#'        [`build_epistats`], containing probability distribution of each year of age, summing to
#'        one. If `NULL`, then a uniform distribution is used.
#' @param edges.avg If `TRUE`, calculates the overall edges target statistics as a weighted average
#'        of the statistics for edges by race/ethnicity group; if `FALSE`, takes the raw average.
#' @param race.prop A numerical vector of length 3, containing the proportion of the population with
#'        each of the three values for the nodal attribute "race" in order: White/Other, Black,
#'        and Hispanic).
#' @param browser If `TRUE`, run `build_netparams` in interactive browser mode.
#'
#' @details
#' This function takes output from [`build_epistats`] and [`build_netparams`] to build the relevant
#' population-level network statistics that will be used in network estimation using the [EpiModel]
#' package.
#'
#' The parameter `edge.avg` allows a user specify how the network edges statistics as estimated from
#' ARTnet with [`build_netparams`] should be calculated. With `edges.avg = FALSE`, the edges count
#' is the overall mean degree (divided by 2) times the `network.size`. However, this may create model
#' fitting issues if sample proportions by race/ethnicity in ARTnet do not match the population
#' proportions (i.e., the overall average may be inconsistent with the weighted average due to
#' different race denominator sizes). The `edges.avg = TRUE` setting calculates the overall edges as
#' the sum of mean degrees in each race/ethnicity group (divided by 2) times the size of each group.
#'
#' The `race.prop` argument only needs to be specified if the built-in race/ethnicity distribution
#' data in [`ARTnetData::race.dist`] is not used. This may be the case `geog.lvl = "county"` or if
#' `geog.cat` has length >1; otherwise, the values can be obtained automatically. If `race.prop` is
#' not supplied in either of these cases, national US values will be used..
#'
#' @export
#'
#' @examples
#' # Standard model with default age stratification
#' epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta")
#' netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
#' netstats <- build_netstats(epistats, netparams)
#'
#' # Restricted age stratification
#' epistats2 <- build_epistats(geog.lvl = "state", geog.cat = "GA",
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(20, 30, 40))
#' netparams2 <- build_netparams(epistats2, smooth.main.dur = TRUE)
#' netstats2 <- build_netstats(epistats2, netparams2)
#'
#' # Model with sexual cessation age < age limit
#' epistats3 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams3 <- build_netparams(epistats3, smooth.main.dur = TRUE)
#' netstats3 <- build_netstats(epistats3, netparams3)
#'
build_netstats <- function(epistats, netparams,
                           network.size = 10000,
                           expect.mort = 0.0001,
                           age.pyramid = NULL,
                           edges.avg = FALSE,
                           race.prop = NULL,
                           browser = FALSE) {

  if (browser == TRUE) {
    browser()
  }

  ## Data ##
  race.dist <- ARTnetData::race.dist

  ## Inputs ##
  geog.cat <- epistats$geog.cat
  geog.lvl <- epistats$geog.lvl
  sex.cess.mod <- epistats$sex.cess.mod
  race <- epistats$race
  age.limits <- epistats$age.limits

  time.unit <- epistats$time.unit


  # Demographic Initialization ----------------------------------------------

  out <- list()
  out$demog <- list()
  out$geog.lvl <- geog.lvl
  out$race <- race
  out$time.unit <- time.unit

  # Overall network size
  num <- out$demog$num <- network.size

  # Population size by race group
  # race.dist.3cat

  if (!is.null(race.prop)) {
    props <- as.data.frame(t(race.prop))
    colnames(props) <- c("White.Other", "Black", "Hispanic")
  } else {
    if (!is.null(geog.lvl) && geog.lvl != "county" && length(geog.cat) == 1) {
      props <- race.dist[[geog.lvl]][which(race.dist[[geog.lvl]]$Geog == geog.cat), -c(1, 2)] / 100
    } else {
      props <- race.dist[["national"]][, -c(1, 2)] / 100
    }
  }
  num.B <- out$demog$num.B <- round(num * props$Black)
  num.H <- out$demog$num.H <- round(num * props$Hispanic)
  num.W <- out$demog$num.W <- num - num.B - num.H

  ## Age-sex-specific mortality rates (B, H, W)
  #  in 1-year age decrements starting with age 1
  #  from CDC NCHS Underlying Cause of Death database (for 2020)
  asmr.B <- c(0.00079, 0.00046, 0.00030, 0.00025, 0.00024, 0.00025, 0.00019,
              0.00019, 0.00021, 0.00020, 0.00026, 0.00026, 0.00038, 0.00056,
              0.00077, 0.00100, 0.00151, 0.00227, 0.00271, 0.00264, 0.00297,
              0.00302, 0.00315, 0.00319, 0.00322, 0.00319, 0.00336, 0.00337,
              0.00330, 0.00363, 0.00396, 0.00392, 0.00407, 0.00428, 0.00411,
              0.00453, 0.00485, 0.00486, 0.00533, 0.00513, 0.00575, 0.00580,
              0.00628, 0.00671, 0.00669, 0.00750, 0.00773, 0.00858, 0.00934,
              0.00947, 0.00999, 0.01141, 0.01216, 0.01360, 0.01432, 0.01517,
              0.01699, 0.01853, 0.02021, 0.02099, 0.02366, 0.02547, 0.02877,
              0.02979, 0.03104, 0.03467, 0.03653, 0.03941, 0.04114, 0.04320,
              0.04487, 0.04879, 0.05100, 0.05678, 0.05611, 0.06384, 0.06891,
              0.07399, 0.07682, 0.08209, 0.08938, 0.09737, 0.10400, 0.11336,
              0.16336)
  asmr.H <- c(0.00032, 0.00021, 0.00018, 0.00011, 0.00011, 0.00009, 0.00010,
              0.00009, 0.00009, 0.00012, 0.00013, 0.00015, 0.00016, 0.00025,
              0.00036, 0.00058, 0.00076, 0.00106, 0.00125, 0.00134, 0.00145,
              0.00156, 0.00164, 0.00166, 0.00164, 0.00159, 0.00176, 0.00172,
              0.00201, 0.00198, 0.00192, 0.00191, 0.00202, 0.00204, 0.00219,
              0.00223, 0.00251, 0.00246, 0.00272, 0.00272, 0.00298, 0.00307,
              0.00321, 0.00351, 0.00367, 0.00391, 0.00442, 0.00484, 0.00512,
              0.00521, 0.00616, 0.00649, 0.00714, 0.00790, 0.00863, 0.00938,
              0.00992, 0.01094, 0.01222, 0.01217, 0.01464, 0.01483, 0.01630,
              0.01731, 0.01850, 0.02054, 0.02269, 0.02321, 0.02515, 0.02734,
              0.02937, 0.03064, 0.03349, 0.03670, 0.03980, 0.04387, 0.04724,
              0.05151, 0.05591, 0.05902, 0.06345, 0.07317, 0.07849, 0.08617,
              0.13436)
  asmr.W <- c(0.00034, 0.00023, 0.00019, 0.00014, 0.00014, 0.00010, 0.00010,
              0.00009, 0.00009, 0.00012, 0.00014, 0.00015, 0.00022, 0.00028,
              0.00036, 0.00050, 0.00059, 0.00082, 0.00096, 0.00104, 0.00126,
              0.00128, 0.00134, 0.00144, 0.00153, 0.00163, 0.00172, 0.00186,
              0.00194, 0.00205, 0.00220, 0.00225, 0.00238, 0.00245, 0.00247,
              0.00264, 0.00274, 0.00280, 0.00306, 0.00312, 0.00324, 0.00329,
              0.00344, 0.00354, 0.00371, 0.00405, 0.00442, 0.00479, 0.00511,
              0.00547, 0.00599, 0.00653, 0.00706, 0.00768, 0.00827, 0.00922,
              0.00978, 0.01065, 0.01151, 0.01235, 0.01349, 0.01437, 0.01548,
              0.01664, 0.01730, 0.01879, 0.01986, 0.02140, 0.02263, 0.02419,
              0.02646, 0.02895, 0.03031, 0.03625, 0.03753, 0.04268, 0.04631,
              0.05235, 0.05724, 0.06251, 0.06934, 0.07589, 0.08669, 0.09582,
              0.16601)

  if (race == TRUE) {
    # transformed to rates by time unit
    trans.asmr.B <- 1 - (1 - asmr.B)^(1 / (364 / time.unit))
    trans.asmr.H <- 1 - (1 - asmr.H)^(1 / (364 / time.unit))
    trans.asmr.W <- 1 - (1 - asmr.W)^(1 / (364 / time.unit))

    # Transformed rates, 85+ rate for ages 85 - 99, total rate for 100
    vec.asmr.B <- c(trans.asmr.B, rep(tail(trans.asmr.B, n=1), 14), 1)
    vec.asmr.H <- c(trans.asmr.H, rep(tail(trans.asmr.H, n=1), 14), 1)
    vec.asmr.W <- c(trans.asmr.W, rep(tail(trans.asmr.W, n=1), 14), 1)
    asmr <- data.frame(age = 1:100, vec.asmr.B, vec.asmr.H, vec.asmr.W)

    out$demog$asmr <- asmr
  } else {
    asmr.O <- rbind(asmr.B, asmr.H, asmr.W)
    asmr.O <- colMeans(asmr.O)

    # transformed to rates by time unit
    trans.asmr <- 1 - (1 - asmr.O)^(1 / (364 / time.unit))

    # Transformed rates, 85+ rate for ages 85 - 99, total rate for 100
    vec.asmr <- c(trans.asmr, rep(tail(trans.asmr, n = 1), 14), 1)
    asmr <- data.frame(age = 1:100, vec.asmr, vec.asmr, vec.asmr)
    out$demog$asmr <- asmr
  }

  # Nodal Attribute Initialization ------------------------------------------

  out$attr <- list()

  # age attributes
  # Currently uniform
  nAges <- age.limits[2] - age.limits[1]
  age.vals <- age.limits[1]:(age.limits[2] - 1)
  if (!is.null(age.pyramid)) {
    if (length(age.pyramid) != nAges) {
      stop("Length of age.pyramid vector must be equal to length of unique age values: ", nAges)
    }
  } else {
    age.pyramid <- rep(1/nAges, nAges)
  }

  attr_age <- sample(x = age.vals, size = num, prob = age.pyramid, replace = TRUE)
  age_noise <- runif(num)
  attr_age <- attr_age + age_noise
  out$attr$age <- attr_age

  attr_sqrt.age <- sqrt(attr_age)
  out$attr$sqrt.age <- attr_sqrt.age

  age.breaks <- out$demog$age.breaks <- epistats$age.breaks
  attr_age.grp <- cut(attr_age, age.breaks, labels = FALSE, right = FALSE, include.lowest = FALSE)
  out$attr$age.grp <- attr_age.grp

  # sexually active attribute
  attr_active.sex <- rep(1, num)
  if (sex.cess.mod == TRUE) {
    attr_active.sex[attr_age.grp == max(attr_age.grp)] <- 0
  }
  out$attr$attr_active.sex <- attr_active.sex

  # race attribute
  attr_race <- apportion_lr(num, 1:3, c(num.B / num, num.H / num, num.W / num), shuffled = TRUE)
  out$attr$race <- attr_race

  # deg.casl attribute
  attr_deg.casl <- apportion_lr(num, 0:3, netparams$main$deg.casl.dist, shuffled = TRUE)
  if (sex.cess.mod == TRUE) {
    attr_deg.casl[attr_active.sex == 0] <- 0
  }
  out$attr$deg.casl <- attr_deg.casl

  # deg main attribute
  attr_deg.main <- apportion_lr(num, 0:2, netparams$casl$deg.main.dist, shuffled = TRUE)
  if (sex.cess.mod == TRUE) {
    attr_deg.main[attr_active.sex == 0] <- 0
  }
  out$attr$deg.main <- attr_deg.main

  # deg tot 3 attribute
  attr_deg.tot <- apportion_lr(num, 0:3, netparams$inst$deg.tot.dist, shuffled = TRUE)
  if (sex.cess.mod == TRUE) {
    attr_deg.tot[attr_active.sex == 0] <- 0
  }
  out$attr$deg.tot <- attr_deg.tot

  # risk group
  nquants <- length(netparams$inst$nf.risk.grp)
  attr_risk.grp <- apportion_lr(num, 1:nquants, rep(1/nquants, nquants), shuffled = TRUE)
  out$attr$risk.grp <- attr_risk.grp

  # role class
  attr_role.class <- apportion_lr(num, 0:2, netparams$all$role.type, shuffled = TRUE)
  out$attr$role.class <- attr_role.class

  # diag status
  if (is.null(epistats$init.hiv.prev)) {
    if (race == TRUE) {
      xs <- data.frame(age = attr_age, race.cat3 = attr_race, geogYN = 1)
      preds <- predict(epistats$hiv.mod, newdata = xs, type = "response")
      attr_diag.status <- rbinom(num, 1, preds)
      out$attr$diag.status <- attr_diag.status
    }  else {
      xs <- data.frame(age = attr_age, geogYN = 1)
      preds <- predict(epistats$hiv.mod, newdata = xs, type = "response")
      attr_diag.status <- rbinom(num, 1, preds)
      out$attr$diag.status <- attr_diag.status
    }
  } else {
    if (race == TRUE) {
      init.hiv.prev <- epistats$init.hiv.prev
      samp.B <- which(attr_race == 1)
      exp.B <- ceiling(length(samp.B) * init.hiv.prev[1])
      samp.H <- which(attr_race == 2)
      exp.H <- ceiling(length(samp.H) * init.hiv.prev[2])
      samp.W <- which(attr_race == 3)
      exp.W <- ceiling(length(samp.W) * init.hiv.prev[3])

      attr_diag.status <- rep(0, network.size)

      attr_diag.status[sample(samp.B, exp.B)] <- 1
      attr_diag.status[sample(samp.H, exp.H)] <- 1
      attr_diag.status[sample(samp.W, exp.W)] <- 1

      out$attr$diag.status <- attr_diag.status

    } else {
      init.hiv.prev <- epistats$init.hiv.prev[1]
      samp.size <- ceiling(network.size * init.hiv.prev)
      attr_diag.status <- sample(1:network.size, samp.size)
      out$attr$diag.status <- rep(0, network.size)
      out$attr$diag.status[attr_diag.status] <- 1
    }
  }


  # Main Model -----------------------------------------------------------

  out$main <- list()

  # edges ---
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$main$edges <- (netparams$main$md.main * num) / 2
    } else {
      out$main$edges <- sum(unname(table(out$attr$race)) * netparams$main$nf.race) / 2
    }
    # nodefactor("race") ---
    nodefactor_race <- table(out$attr$race) * netparams$main$nf.race
    out$main$nodefactor_race <- unname(nodefactor_race)

    # nodematch("race") ---
    nodematch_race <- nodefactor_race / 2 * netparams$main$nm.race
    out$main$nodematch_race <- unname(nodematch_race)

    # nodematch("race", diff = FALSE) ---
    nodematch_race <- out$main$edges * netparams$main$nm.race_diffF
    out$main$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$main$edges <- (netparams$main$md.main * num) / 2
  }

  # nodefactor("age.grp") ---
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$main$nf.age.grp
  if (sex.cess.mod == TRUE) {
    nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
  }
  out$main$nodefactor_age.grp <- unname(nodefactor_age.grp)

  # nodematch("age.grp") ---
  nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$main$nm.age.grp
  out$main$nodematch_age.grp <- unname(nodematch_age.grp)

  # absdiff("age") ---
  absdiff_age <- out$main$edges * netparams$main$absdiff.age
  out$main$absdiff_age <- absdiff_age

  # absdiff("sqrt.age") ---
  absdiff_sqrt.age <- out$main$edges * netparams$main$absdiff.sqrt.age
  out$main$absdiff_sqrt.age <- absdiff_sqrt.age

  # nodefactor("deg.casl") ---
  out$main$nodefactor_deg.casl <- num * netparams$main$deg.casl.dist * netparams$main$nf.deg.casl

  # concurrent ---
  out$main$concurrent <- num * netparams$main$concurrent

  # nodefactor("diag.status") ---
  nodefactor_diag.status <- table(out$attr$diag.status) * netparams$main$nf.diag.status
  out$main$nodefactor_diag.status <- unname(nodefactor_diag.status)

  # Dissolution ---
  out$main$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$main$durs.main.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$main$diss.byage <- dissolution_coefs(dissolution = ~offset(edges) +
                                             offset(nodematch("age.grp", diff = TRUE)),
                                           duration = netparams$main$durs.main.byage$mean.dur.adj,
                                           d.rate = expect.mort)



  # Casual Model ------------------------------------------------------------

  out$casl <- list()

  # edges ---
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$casl$edges <- (netparams$casl$md.casl * num) / 2
    } else {
      out$casl$edges <- sum(unname(table(out$attr$race)) * netparams$casl$nf.race) / 2
    }
    # nodefactor("race") ---
    nodefactor_race <- table(out$attr$race) * netparams$casl$nf.race
    out$casl$nodefactor_race <- unname(nodefactor_race)

    # nodematch("race") ---
    nodematch_race <- nodefactor_race / 2 * netparams$casl$nm.race
    out$casl$nodematch_race <- unname(nodematch_race)

    # nodematch("race", diff = FALSE) ---
    nodematch_race <- out$casl$edges * netparams$casl$nm.race_diffF
    out$casl$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$casl$edges <- (netparams$casl$md.casl * num) / 2
  }

  # nodefactor("age.grp") ---
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$casl$nf.age.grp
  if (sex.cess.mod == TRUE) {
    nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
  }
  out$casl$nodefactor_age.grp <- unname(nodefactor_age.grp)

  # nodematch("age.grp") ---
  nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$casl$nm.age.grp
  out$casl$nodematch_age.grp <- unname(nodematch_age.grp)

  # absdiff("age") ---
  absdiff_age <- out$casl$edges * netparams$casl$absdiff.age
  out$casl$absdiff_age <- absdiff_age

  # absdiff("sqrt.age") ---
  absdiff_sqrt.age <- out$casl$edges * netparams$casl$absdiff.sqrt.age
  out$casl$absdiff_sqrt.age <- absdiff_sqrt.age

  # nodefactor("deg.main") ---
  out$casl$nodefactor_deg.main <- num * netparams$casl$deg.main.dist * netparams$casl$nf.deg.main

  # concurrent ---
  out$casl$concurrent <- num * netparams$casl$concurrent

  # nodefactor("diag.status") ---
  nodefactor_diag.status <- table(out$attr$diag.status) * netparams$casl$nf.diag.status
  out$casl$nodefactor_diag.status <- unname(nodefactor_diag.status)

  # Dissolution
  out$casl$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$casl$durs.casl.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$casl$diss.byage <- dissolution_coefs(dissolution = ~offset(edges) +
                                             offset(nodematch("age.grp", diff = TRUE)),
                                           duration = netparams$casl$durs.casl.byage$mean.dur.adj,
                                           d.rate = expect.mort)


  # One-Time Model ----------------------------------------------------------

  out$inst <- list()

  # edges ---
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$inst$edges <- (netparams$inst$md.inst * num) / 2
    } else {
      out$inst$edges <- sum(unname(table(out$attr$race)) * netparams$inst$nf.race) / 2
    }
    # nodefactor("race") ---
    nodefactor_race <- table(out$attr$race) * netparams$inst$nf.race
    out$inst$nodefactor_race <- unname(nodefactor_race)

    # nodematch("race") ---
    nodematch_race <- nodefactor_race / 2 * netparams$inst$nm.race
    out$inst$nodematch_race <- unname(nodematch_race)

    # nodematch("race", diff = FALSE) ---
    nodematch_race <- out$inst$edges * netparams$inst$nm.race_diffF
    out$inst$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$inst$edges <- (netparams$inst$md.inst * num) / 2
  }

  # nodefactor("age.grp") ---
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$inst$nf.age.grp
  if (sex.cess.mod == TRUE) {
    nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
  }
  out$inst$nodefactor_age.grp <- unname(nodefactor_age.grp)

  # nodematch("age.grp") ---
  nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$inst$nm.age.grp
  out$inst$nodematch_age.grp <- unname(nodematch_age.grp)

  # absdiff("age") ---
  absdiff_age <- out$inst$edges * netparams$inst$absdiff.age
  out$inst$absdiff_age <- absdiff_age

  # absdiff("sqrt.age") ---
  absdiff_sqrt.age <- out$inst$edges * netparams$inst$absdiff.sqrt.age
  out$inst$absdiff_sqrt.age <- absdiff_sqrt.age

  # nodefactor("risk.grp") ---
  nodefactor_risk.grp <- table(out$attr$risk.grp) * netparams$inst$nf.risk.grp
  out$inst$nodefactor_risk.grp <- unname(nodefactor_risk.grp)

  # nodefactor("deg.tot") ---
  nodefactor_deg.tot <- table(out$attr$deg.tot) * netparams$inst$nf.deg.tot
  out$inst$nodefactor_deg.tot <- unname(nodefactor_deg.tot)

  # nodefactor("diag.status") ---
  nodefactor_diag.status <- table(out$attr$diag.status) * netparams$inst$nf.diag.status
  out$inst$nodefactor_diag.status <- unname(nodefactor_diag.status)


  return(out)
}


#' Update mortality rates
#'
#' @description Replaces the default age- and sex-specific mortality rates used
#'              in \code{\link{build_netstats}} with user-specified values.
#'
#' @param netstats Output from \code{\link{build_netstats}}.
#' @param asmr_df A data frame with three columns (Race, Age, and DeathRate)
#'                specifying the desired mortality rates per time unit for each
#'                combination of race (Black, Hispanic, White) and age
#'                (1, 2, 3, ... 100).
#'
#' @details
#' \code{update_asmr} takes the output from \code{\link{build_netstats}} and
#' updates its mortality rates with the values specified in the input parameter
#' \code{asmr_df}. Mortality rates should be provided for each combination of
#' race (Black, Hispanic, White) and age (1, 2, 3, ... 100), for a total of 300
#' rates. All input mortality rates must be non-NA, numeric, and between 0 and
#' 1 (inclusive). The maximum mortality rate for each race must be 1,
#' representing deterministic departure from the network at that age.
#' \code{update_asmr} does not perform time unit conversions; the input
#' mortality rates should be pre-converted by the user to match
#' the \code{time.unit} used in \code{\link{build_netstats}}.
#'
#' @export
#'
#' @examples
#' epistats  <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
#'                             race = TRUE)
#' netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#' netstats  <- build_netstats(epistats, netparams)
#' asmr_df   <- data.frame(Race = rep(c("Black", "Hispanic", "White"),
#'                                    each = 100), Age = rep(1:100,3),
#'                         DeathRate = rep(c(runif(99), 1), 3))
#' netstats  <- update_asmr(netstats, asmr_df)

update_asmr <- function(netstats, asmr_df) {

  #Check that all input mortality rates are reasonable values
  if (sum(is.na(asmr_df$DeathRate)) > 0 | !is.numeric(asmr_df$DeathRate) |
      min(asmr_df$DeathRate < 0) | max(asmr_df$DeathRate > 1)) {
    stop("Ensure all mortality rates are non-NA, numeric, and between 0 and 1.")
  }

  asmr.B <- asmr_df[asmr_df$Race == "Black", 2:3]
  asmr.H <- asmr_df[asmr_df$Race == "Hispanic", 2:3]
  asmr.W <- asmr_df[asmr_df$Race == "White", 2:3]

  #Check for appropriate number of input mortality rates
  if (min(asmr.B$Age) != 1 | min(asmr.H$Age) != 1 | min(asmr.W$Age) != 1 |
      max(asmr.B$Age) != 100 | max(asmr.H$Age) != 100 | max(asmr.W$Age) != 100 |
      nrow(asmr.B) != 100 | nrow(asmr.H) != 100 | nrow(asmr.W) != 100) {
    stop("Provide mortality rates by race for 1-yr age groups (1 - 100).")
  }

  if (max(asmr.B$DeathRate) != 1 | max(asmr.H$DeathRate) != 1 |
      max(asmr.W$DeathRate) != 1) {
    stop("The mortality rate for at least one age group must be total (1).")
  }

  if (netstats$race == TRUE) {
    asmr <- data.frame(age = 1:100, asmr.B$DeathRate,
                       asmr.H$DeathRate, asmr.W$DeathRate)
  } else {
    asmr.O <-  rbind(asmr.B$DeathRate, asmr.H$DeathRate, asmr.W$DeathRate)
    asmr.O <- colMeans(asmr.O, na.rm = TRUE)
    asmr <- data.frame(age = 1:100, asmr.O, asmr.O, asmr.O)
  }

  message(paste0("The time unit used in build_netstats() was ",
                 netstats$time.unit,
                 ". Please ensure that that is consistent with your inputs."))

  netstats$demog$asmr <- asmr
  return(netstats)

}
