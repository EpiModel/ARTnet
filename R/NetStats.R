
#' Calculate Network Target Statistics
#'
#' @description Calculates the final target statistics for the network models by
#'              applying individual-level network statistics against the population
#'              size and structure, for use in the EpiModelHIV workflow.
#'
#' @param epistats Output from \code{\link{build_epistats}}.
#' @param netparams Output from \code{\link{build_netparams}}.
#' @param network.size Size of the starting network.
#' @param expect.mort Expected average mortality level to pass into
#'        \code{\link{dissolution_coefs}} function.
#' @param edges.avg Whether degree differences exist along race. TRUE
#'        or FALSE; default of FALSE.
#'
#' @details
#' \code{build_netstats} takes output from \code{\link{build_epistats}} and
#' \code{\link{build_netparams}} to build the relevant network statistics
#' that will be used in network estimation using package \link{EpiModel}.
#'
#' The parameter \code{edge.avg} allows a user set the network stated edges
#' to that estimated in \code{\link{build_netparams}} (divided by 2),
#' with \code{edges.avg = FALSE}, or, if sample proportions do not match
#' ARTnet population proportions, set to a weighted racial average
#' with \code{edges.avg = TRUE}.
#'
#' @export
#'
#' @examples
#' epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta")
#' netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#' netstats <- build_netstats(epistats, netparams)
#'
build_netstats <- function(epistats, netparams,
                           network.size = 10000,
                           expect.mort = 0.0001,
                           edges.avg = FALSE) {

  ## Data ##
  #NOTE: Not actually used
  d <- epistats$wide
  l <- epistats$long
  race.dist <- ARTnetData::race.dist

  ## Inputs ##
  geog.cat <- epistats$geog.cat
  geog.lvl <- epistats$geog.lvl
  race <- epistats$race
  age.grps <- epistats$age.grps
  age.limits <- epistats$age.limits


  # Demographic Initialization ----------------------------------------------

  out <- list()
  out$demog <- list()
  out$geog.lvl <- geog.lvl
  out$race <- race

  # Overall network size
  num <- out$demog$num <- network.size

  # Population size by race group
  # race.dist.3cat

  if (!is.null(geog.lvl)) {
  props <- race.dist[[geog.lvl]][which(race.dist[[geog.lvl]]$Geog == geog.cat), -c(1,2)]/100
  } else {
    props <- race.dist[["national"]][, -c(1,2)]/100
  }
  num.B <- out$demog$num.B <- round(num * props$Black)
  num.H <- out$demog$num.H <- round(num * props$Hispanic)
  num.W <- out$demog$num.W <- num - num.B - num.H

  ## Age-sex-specific mortality rates (B, H, W)
  #    in 5-year age decrments starting with age 15
  ages <- out$demog$ages <- age.limits[1]:age.limits[2]
  asmr.B <- c(0.00124, 0.00213, 0.00252, 0.00286, 0.00349,
              0.00422, 0.00578, 0.00870, 0.01366, 0.02052)
  asmr.H <- c(0.00062, 0.00114, 0.00127, 0.00132, 0.00154,
              0.00186, 0.00271, 0.00440, 0.00643, 0.00980)
  asmr.W <- c(0.00064, 0.00128, 0.00166, 0.00199, 0.00266,
              0.00272, 0.00382, 0.00591, 0.00889, 0.01266)

  if (race == TRUE) {
    # transformed to weekly rates
    trans.asmr.B <- 1 - (1 - asmr.B)^(1/52)
    trans.asmr.H <- 1 - (1 - asmr.H)^(1/52)
    trans.asmr.W <- 1 - (1 - asmr.W)^(1/52)

    # Null rate for 0-14, transformed rates, total rate for 65
    vec.asmr.B <- c(rep(0, 14), rep(trans.asmr.B, each = 5), 1)
    vec.asmr.H <- c(rep(0, 14), rep(trans.asmr.H, each = 5), 1)
    vec.asmr.W <- c(rep(0, 14), rep(trans.asmr.W, each = 5), 1)
    asmr <- data.frame(age = 1:65, vec.asmr.B, vec.asmr.H, vec.asmr.W)

    out$demog$asmr <- asmr
  } else {
    asmr.O <-  rbind(asmr.B, asmr.H, asmr.W)
    asmr.O <- colMeans(asmr.O)

    # transformed to weekly rates
    trans.asmr <- 1 - (1 - asmr.O)^(1/52)

    # Null rate for 0-14, transformed rates, total rate for 65
    vec.asmr <- c(rep(0,14), rep(trans.asmr, each = 5),1)
    asmr <- data.frame(age = 1:65, vec.asmr, vec.asmr, vec.asmr)
    out$demog$asmr <- asmr
  }

  # Nodal Attribute Initialization ------------------------------------------

  out$attr <- list()

  # age attributes
  attr_age <- runif(num, min = min(ages), max = max(ages) + (51/52))
  out$attr$age <- attr_age

  attr_sqrt.age <- sqrt(attr_age)
  out$attr$sqrt.age <- attr_sqrt.age

  age.breaks <- out$demog$age.breaks <- epistats$age.breaks
  attr_age.grp <- cut(attr_age, age.breaks, labels = FALSE,
                      right = FALSE, include.lowest = FALSE)
  out$attr$age.grp <- attr_age.grp

  # race attribute
  attr_race <- apportion_lr(num, 1:3, c(num.B/num, num.H/num, num.W/num),
                            shuffled = TRUE)
  out$attr$race <- attr_race

  # deg.casl attribute
  attr_deg.casl <- apportion_lr(num, 0:3, netparams$main$deg.casl.dist,
                                shuffled = TRUE)
  out$attr$deg.casl <- attr_deg.casl

  # deg main attribute
  attr_deg.main <- apportion_lr(num, 0:2, netparams$casl$deg.main.dist,
                                shuffled = TRUE)
  out$attr$deg.main <- attr_deg.main

  # deg tot 3 attribute
  attr_deg.tot <- apportion_lr(num, 0:3, netparams$inst$deg.tot.dist,
                               shuffled = TRUE)
  out$attr$deg.tot <- attr_deg.tot

  # risk group
  attr_risk.grp <- apportion_lr(num, 1:5, rep(0.2, 5), shuffled = TRUE)
  out$attr$risk.grp <- attr_risk.grp

  # role class
  attr_role.class <- apportion_lr(num, 0:2, netparams$all$role.type,
                                  shuffled = TRUE)
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
      exp.B <- ceiling(length(samp.B)*init.hiv.prev[1])
      samp.H <- which(attr_race == 2)
      exp.H <- ceiling(length(samp.H)*init.hiv.prev[2])
      samp.W <- which(attr_race == 3)
      exp.W <- ceiling(length(samp.W)*init.hiv.prev[3])

      attr_diag.status <- rep(0, network.size)

      attr_diag.status[sample(samp.B, exp.B)] <- 1
      attr_diag.status[sample(samp.H, exp.H)] <- 1
      attr_diag.status[sample(samp.W, exp.W)] <- 1

      out$attr$diag.status <- attr_diag.status

    } else {
      init.hiv.prev <- epistats$init.hiv.prev[1]
      samp.size <- ceiling(network.size*init.hiv.prev)
      attr_diag.status <- sample(1:network.size, samp.size)
      out$attr$diag.status <- rep(0, network.size)
      out$attr$diag.status[attr_diag.status] <- 1
    }
  }


  # Main Model -----------------------------------------------------------

  out$main <- list()

  ## edges
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$main$edges <- (netparams$main$md.main * num) / 2
    } else {
      out$main$edges <- sum(unname(table(out$attr$race)) * netparams$main$nf.race)/2
    }
    ## nodefactor("race")
    nodefactor_race <- table(out$attr$race) * netparams$main$nf.race
    out$main$nodefactor_race <- unname(nodefactor_race)

    ## nodematch("race")
    nodematch_race <- nodefactor_race/2 * netparams$main$nm.race
    out$main$nodematch_race <- unname(nodematch_race)

    ## nodematch("race", diff = FALSE)
    nodematch_race <- out$main$edges * netparams$main$nm.race_diffF
    out$main$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$main$edges <- (netparams$main$md.main * num) / 2
  }

  ## nodefactor("age.grp
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$main$nf.age.grp
  out$main$nodefactor_age.grp <- unname(nodefactor_age.grp)

  ## nodematch("age.grp")
  nodematch_age.grp <- nodefactor_age.grp/2 * netparams$main$nm.age.grp
  out$main$nodematch_age.grp <- unname(nodematch_age.grp)

  ## absdiff("age")
  absdiff_age <- out$main$edges * netparams$main$absdiff.age
  out$main$absdiff_age <- absdiff_age

  ## absdiff("sqrt.age")
  absdiff_sqrt.age <- out$main$edges * netparams$main$absdiff.sqrt.age
  out$main$absdiff_sqrt.age <- absdiff_sqrt.age

  ## nodefactor("deg.casl")
  out$main$nodefactor_deg.casl <-
    num * netparams$main$deg.casl.dist * netparams$main$nf.deg.casl

  ## concurrent
  out$main$concurrent <- num * netparams$main$concurrent

  ## nodefactor("diag.status")
  nodefactor_diag.status <-
    table(out$attr$diag.status) * netparams$main$nf.diag.status
  out$main$nodefactor_diag.status <- unname(nodefactor_diag.status)

  ## Dissolution
  out$main$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$main$durs.main.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$main$diss.byage <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("age.grp", diff = TRUE)),
                                           duration = netparams$main$durs.main.byage$mean.dur.adj,
                                           d.rate = expect.mort)



  # Casual Model ------------------------------------------------------------

  out$casl <- list()

  ## edges
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$casl$edges <- (netparams$casl$md.casl * num) / 2
    } else {
      out$casl$edges <- sum(unname(table(out$attr$race)) * netparams$casl$nf.race)/2
    }
    ## nodefactor("race")
    nodefactor_race <- table(out$attr$race) * netparams$casl$nf.race
    out$casl$nodefactor_race <- unname(nodefactor_race)

    ## nodematch("race")
    nodematch_race <- nodefactor_race/2 * netparams$casl$nm.race
    out$casl$nodematch_race <- unname(nodematch_race)

    ## nodematch("race", diff = FALSE)
    nodematch_race <- out$casl$edges * netparams$casl$nm.race_diffF
    out$casl$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$casl$edges <- (netparams$casl$md.casl * num) / 2
  }

  ## nodefactor("age.grp")
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$casl$nf.age.grp
  out$casl$nodefactor_age.grp <- unname(nodefactor_age.grp)

  ## nodematch("age.grp")
  nodematch_age.grp <- nodefactor_age.grp/2 * netparams$casl$nm.age.grp
  out$casl$nodematch_age.grp <- unname(nodematch_age.grp)

  ## absdiff("age")
  absdiff_age <- out$casl$edges * netparams$casl$absdiff.age
  out$casl$absdiff_age <- absdiff_age

  ## absdiff("sqrt.age")
  absdiff_sqrt.age <- out$casl$edges * netparams$casl$absdiff.sqrt.age
  out$casl$absdiff_sqrt.age <- absdiff_sqrt.age

  ## nodefactor("deg.main")
  out$casl$nodefactor_deg.main <- num * netparams$casl$deg.main.dist * netparams$casl$nf.deg.main

  ## concurrent
  out$casl$concurrent <- num * netparams$casl$concurrent

  ## nodefactor("diag.status")
  nodefactor_diag.status <- table(out$attr$diag.status) * netparams$casl$nf.diag.status
  out$casl$nodefactor_diag.status <- unname(nodefactor_diag.status)

  ## Dissolution
  out$casl$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$casl$durs.casl.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$casl$diss.byage <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("age.grp", diff = TRUE)),
                                           duration = netparams$casl$durs.casl.byage$mean.dur.adj,
                                           d.rate = expect.mort)


  # One-Time Model ----------------------------------------------------------

  out$inst <- list()

  ## edges
  if (race == TRUE) {
    if (edges.avg == FALSE) {
      out$inst$edges <- (netparams$inst$md.inst * num) / 2
    } else {
      out$inst$edges <- sum(unname(table(out$attr$race)) * netparams$inst$nf.race)/2
    }
    ## nodefactor("race")
    nodefactor_race <- table(out$attr$race) * netparams$inst$nf.race
    out$inst$nodefactor_race <- unname(nodefactor_race)

    ## nodematch("race")
    nodematch_race <- nodefactor_race/2 * netparams$inst$nm.race
    out$inst$nodematch_race <- unname(nodematch_race)

    ## nodematch("race", diff = FALSE)
    nodematch_race <- out$inst$edges * netparams$inst$nm.race_diffF
    out$inst$nodematch_race_diffF <- unname(nodematch_race)
  } else {
    out$inst$edges <- (netparams$inst$md.inst * num) / 2
  }

  ## nodefactor("age.grp")
  nodefactor_age.grp <- table(out$attr$age.grp) * netparams$inst$nf.age.grp
  out$inst$nodefactor_age.grp <- unname(nodefactor_age.grp)

  ## nodematch("age.grp")
  nodematch_age.grp <- nodefactor_age.grp/2 * netparams$inst$nm.age.grp
  out$inst$nodematch_age.grp <- unname(nodematch_age.grp)

  ## absdiff("age")
  absdiff_age <- out$inst$edges * netparams$inst$absdiff.age
  out$inst$absdiff_age <- absdiff_age

  ## absdiff("sqrt.age")
  absdiff_sqrt.age <- out$inst$edges * netparams$inst$absdiff.sqrt.age
  out$inst$absdiff_sqrt.age <- absdiff_sqrt.age

  ## nodefactor("risk.grp")
  nodefactor_risk.grp <- table(out$attr$risk.grp) * netparams$inst$nf.risk.grp
  out$inst$nodefactor_risk.grp <- unname(nodefactor_risk.grp)

  ## nodefactor("deg.tot")
  nodefactor_deg.tot <- table(out$attr$deg.tot) * netparams$inst$nf.deg.tot
  out$inst$nodefactor_deg.tot <- unname(nodefactor_deg.tot)

  ## nodefactor("diag.status")
  nodefactor_diag.status <- table(out$attr$diag.status) * netparams$inst$nf.diag.status
  out$inst$nodefactor_diag.status <- unname(nodefactor_diag.status)

  # Return Out

  return(out)
}
