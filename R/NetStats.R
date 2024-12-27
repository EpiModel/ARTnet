
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
#'        one. If `NULL`, then a default distribution based on 2020 NCHS data is used.
#' @param edges.avg If `TRUE`, calculates the overall edges target statistics as a weighted average
#'        of the statistics for edges by race/ethnicity group; if `FALSE`, takes the raw average.
#' @param race.prop A numerical vector of length 3, containing the proportion of the population with
#'        each of the three values for the nodal attribute "race" in order: Black,
#'        Hispanic, and White/Other).
#' @param young.prop The proportion of the population that should be below the age of sexual cessation.
#'        Default is NULL (meaning no re-weighting of the `age.pyramid` parameter is performed).
#'        This parameter is only used if the age of sexual cessation is less than the upper age bound.
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
#' data in `ARTnetData::race.dist` is not used. This may be the case `geog.lvl = "county"` or if
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
#' # Model with sexual cessation age < age limit, without age pyramid reweighting
#' epistats3 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams3 <- build_netparams(epistats3, smooth.main.dur = TRUE)
#' netstats3 <- build_netstats(epistats3, netparams3)
#'
#' # Model with sexual cessation age < age limit, with age pyramid reweighting
#' epistats4 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams4 <- build_netparams(epistats4, smooth.main.dur = TRUE)
#' netstats4 <- build_netstats(epistats4, netparams4, young.prop = 0.995)
#'
build_netstats <- function(epistats, netparams,
                           network.size = 10000,
                           expect.mort = 0.0001,
                           age.pyramid = NULL,
                           edges.avg = FALSE,
                           race.prop = NULL,
                           young.prop = NULL,
                           browser = FALSE) {
  # Ensures that ARTnetData is installed
  if (system.file(package = "ARTnetData") == "") stop(missing_data_msg)

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
  race.level <- epistats$race.level
  age.limits <- epistats$age.limits
  age.sexual.cessation <- epistats$age.sexual.cessation

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

  if (!is.null(race.prop)) {
    # Capitalize each race and join multiple with a `.`
    # c("white", "other") => "White.Other"
    flattened_race_level <- vapply(
      race.level,
      \(x) gsub("\\b(\\w)", "\\U\\1", paste0(x, collapse = "."), perl = TRUE),
      character(1)
    )
    props <- as.data.frame(t(race.prop))
    colnames(props) <- flattened_race_level
  } else {
    if (!is.null(geog.lvl) && geog.lvl != "county" && length(geog.cat) == 1) {
      props <- race.dist[[geog.lvl]][which(race.dist[[geog.lvl]]$Geog == geog.cat), -c(1, 2)] / 100
    } else {
      props <- race.dist[["national"]][, -c(1, 2)] / 100
    }
  }
  out$demog$props <- props
  race.num.vars <- list()
  total_remaining <- num
  for (i in seq_len(length(flattened_race_level) - 1)) {
    race_name <- flattened_race_level[i]
    race_num_var <- paste0("num.", race_name)
    race_num_value <- round(num * props[[race_name]])
    out$demog[[race_num_var]] <- race_num_value
    race.num.vars[[race_name]] <- race_num_value
    total_remaining <- total_remaining - race_num_value
  }

  # Assign the residual race group
  residual_race <- flattened_race_level[length(flattened_race_level)]
  residual_num_var <- paste0("num.", residual_race)
  out$demog[[residual_num_var]] <- total_remaining
  race.num.vars[[residual_race]] <- total_remaining

  for (race_name in names(race.num.vars)) {
    race_num_var <- paste0("num.", substr(race_name, 1, 1))
    assign(race_num_var, race.num.vars[[race_name]])
  }

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

    # Transformed rates, 85+ rate for ages 85 - 100
    vec.asmr.B <- c(trans.asmr.B, rep(tail(trans.asmr.B, n = 1), 15))
    vec.asmr.H <- c(trans.asmr.H, rep(tail(trans.asmr.H, n = 1), 15))
    vec.asmr.W <- c(trans.asmr.W, rep(tail(trans.asmr.W, n = 1), 15))

    asmr <- data.frame(age = 1:100, vec.asmr.B, vec.asmr.H, vec.asmr.W)
  } else {
    asmr.O <- props$Black * asmr.B + props$Hispanic * asmr.H +
      props$White.Other * asmr.W

    trans.asmr <- 1 - (1 - asmr.O)^(1 / (364 / time.unit))

    vec.asmr <- c(trans.asmr, rep(tail(trans.asmr, n = 1), 15))
    asmr <- data.frame(age = 1:100, vec.asmr, vec.asmr, vec.asmr)
  }

  # Setting deterministic mortality prob = 1 at upper age limit
  max.age <- age.limits[2]
  asmr[asmr$age >= max.age, -1] <- 1
  out$demog$asmr <- asmr


  # Nodal Attribute Initialization ------------------------------------------

  out$attr <- list()

  # age attributes
  # Currently uniform
  nAges <- age.limits[2] - age.limits[1]
  age.vals <- age.limits[1]:(age.limits[2] - 1)
  out$demog$ages <- age.vals
  if (!is.null(age.pyramid)) {
    if (length(age.pyramid) != nAges) {
      stop("Length of age.pyramid vector must be equal to length of unique age values: ", nAges)
    }
  } else {
    full.age.pyr <- c(0.01202, 0.01228, 0.01250, 0.01280, 0.01292, 0.01289,
                      0.01284, 0.01286, 0.01301, 0.01297, 0.01296, 0.01337,
                      0.01344, 0.01334, 0.01329, 0.01332, 0.01325, 0.01323,
                      0.01360, 0.01385, 0.01365, 0.01368, 0.01373, 0.01389,
                      0.01424, 0.01454, 0.01480, 0.01514, 0.01533, 0.01525,
                      0.01463, 0.01426, 0.01400, 0.01401, 0.01404, 0.01352,
                      0.01366, 0.01360, 0.01339, 0.01367, 0.01274, 0.01246,
                      0.01227, 0.01190, 0.01226, 0.01181, 0.01192, 0.01245,
                      0.01312, 0.01331, 0.01255, 0.01225, 0.01219, 0.01239,
                      0.01308, 0.01321, 0.01313, 0.01303, 0.01311, 0.01320,
                      0.01264, 0.01247, 0.01223, 0.01169, 0.01152, 0.01092,
                      0.01042, 0.00994, 0.00954, 0.00926, 0.00890, 0.00872,
                      0.00899, 0.00650, 0.00630, 0.00600, 0.00601, 0.00509,
                      0.00451, 0.00414, 0.00377, 0.00346, 0.00304, 0.00275,
                      0.00235, 0.00202, 0.00173, 0.00148, 0.00127, 0.00109,
                      0.00093, 0.00080, 0.00068, 0.00059, 0.00050, 0.00043,
                      0.00037, 0.00032, 0.00027, 0.00023)
    age.pyramid <- full.age.pyr[age.vals]
  }

  if (age.sexual.cessation < age.limits[2] && !is.null(young.prop)) {
    age.break <- age.sexual.cessation - (age.limits[1] - 1)
    age.pyramid <- reweight_age_pyr(age.pyramid, young.prop, age.break)
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
  out$attr$active.sex <- attr_active.sex

  # race attribute
  race_numbers <- sapply(flattened_race_level, function(race) {
    race_num_var <- paste0("num.", toupper(substr(race, 1, 1)))
    out$demog[[race_num_var]]
  })

  race_proportions <- race_numbers / num
  group_ids <- seq_along(flattened_race_level)
  attr_race <- apportion_lr(num, group_ids, race_proportions, shuffled = TRUE)
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
  attr_risk.grp <- apportion_lr(num, 1:nquants, rep(1 / nquants, nquants), shuffled = TRUE)
  out$attr$risk.grp <- attr_risk.grp

  # role class
  attr_role.class <- apportion_lr(num, 0:2, netparams$all$role.type, shuffled = TRUE)
  out$attr$role.class <- attr_role.class

  # diag status
  if (is.null(epistats$init.hiv.prev)) {
    if (race == TRUE) {
      xs <- data.frame(age = attr_age, race.cat.num = attr_race, geogYN = 1)
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
      attr_diag.status <- rep(0, network.size)

      for (i in seq_along(flattened_race_level)) {
        samp <- which(attr_race == i)
        exp <- ceiling(length(samp) * init.hiv.prev[i])
        attr_diag.status[sample(samp, exp)] <- 1
      }

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
#'              in [`build_netstats`] with user-specified values.
#'
#' @param netstats Output from [`build_netstats`].
#' @param asmr_df A data frame with three columns (Race, Age, and DeathRate)
#'                specifying the desired mortality rates per time unit for each
#'                combination of race (Black, Hispanic, White) and age
#'                (1, 2, 3, ... 100).
#'
#' @details
#' `update_asmr` takes the output from [`build_netstats`] and
#' updates its mortality rates with the values specified in the input parameter
#' `asmr_df`. Mortality rates should be provided for each combination of
#' race (Black, Hispanic, White) and age (1, 2, 3, ... 100), for a total of 300
#' rates. All input mortality rates must be non-NA, numeric, and between 0 and
#' 1 (inclusive). The maximum mortality rate for each race must be 1,
#' representing deterministic departure from the network at that age.
#' `update_asmr` does not perform time unit conversions; the input
#' mortality rates should be pre-converted by the user to match
#' the `time.unit` used in [`build_netstats`].
#'
#' @export
#'
#' @examples
#' epistats  <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
#'                             race = TRUE)
#' netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#' netstats  <- build_netstats(epistats, netparams)
#'
#' # Update mortality rates with random values
#' asmr_df   <- data.frame(Race = rep(c("Black", "Hispanic", "White"),
#'                                    each = 100), Age = rep(1:100,3),
#'                         DeathRate = rep(c(runif(99), 1), 3))
#' netstats  <- update_asmr(netstats, asmr_df)
#'
#' # Update mortality rates with 2020 NCHS data
#' asmr_df   <- read.csv(system.file("2020DeathRates.csv", package =  "ARTnet"))
#' netstats  <- update_asmr(netstats, asmr_df)
#'
update_asmr <- function(netstats, asmr_df) {
  # Check that all input mortality rates are reasonable values
  if (sum(is.na(asmr_df$DeathRate)) > 0 || !is.numeric(asmr_df$DeathRate) ||
      min(asmr_df$DeathRate) < 0 || max(asmr_df$DeathRate) > 1) {
    stop("Ensure all mortality rates are non-NA, numeric, and between 0 and 1.")
  }

  # Get unique races
  unique_races <- unique(asmr_df$Race)

  # Initialize a list to store mortality rates by race
  asmr_list <- list()

  # Validate mortality rates for each race
  for (race in unique_races) {
    race_data <- asmr_df[asmr_df$Race == race, 2:3]
    if (min(race_data$Age) != 1 || max(race_data$Age) != 100 || nrow(race_data) != 100) {
      stop(paste0("Provide mortality rates by race for 1-yr age groups (1 - 100) for race: ", race))
    }
    if (max(race_data$DeathRate) != 1) {
      stop(paste0("The mortality rate for at least one age group must be total (1) for race: ", race))
    }
    asmr_list[[race]] <- race_data$DeathRate
  }

  # If race-specific mortality rates should be included
  if (netstats$race == TRUE) {
    # Combine mortality rates for all races into a data frame
    asmr <- data.frame(
      age = 1:100,
      setNames(do.call(cbind, asmr_list), paste0("asmr.", unique_races))
    )
  } else {
    # Calculate combined mortality rates using demographic proportions
    props <- netstats$demog$props
    asmr_combined <- rep(0, 100)
    for (race in unique_races) {
      if (!race %in% names(props)) {
        stop(paste0("Demographic proportions for race: ", race, " are missing in netstats$demog$props."))
      }
      asmr_combined <- asmr_combined + props[[race]] * asmr_list[[race]]
    }
    # Create a combined mortality data frame
    asmr <- data.frame(age = 1:100, asmr_combined, asmr_combined, asmr_combined)
    colnames(asmr)[-1] <- c("asmr.O1", "asmr.O2", "asmr.O3")
  }

  # Display time unit consistency message
  message(paste0("The time unit used in build_netstats() was ",
                 netstats$time.unit,
                 ". Please ensure that that is consistent with your inputs."))

  # Assign the resulting asmr to netstats
  netstats$demog$asmr <- asmr
  return(netstats)
}

#' Re-weight age pyramid
#'
#' @description This function re-weights a given age pyramid (i.e., a vector
#' with one element per year in an age range that specifies what proportion of
#' the population should be  that age) so that the desired fraction of the
#' population falls below a specified age boundary.
#'
#' @param age.pyramid A vector of proportions.
#' @param young.prop The weight that should be applied to the portion of
#'                   `age.pyramid` that is below `age.break`
#' @param age.break The index at which the upper group should begin.
#'
#' @export
#'
#' @examples
#' unif.age.pyr <- rep(1/100, 100)
#' reweighted.age.pyr <- reweight_age_pyr(unif.age.pyr, 0.995, 65)
#'
reweight_age_pyr <- function(age.pyramid, young.prop, age.break) {

  #Check inputs
  if (age.break > length(age.pyramid) || young.prop > 1 || young.prop < 0 ||
      sum(is.na(age.pyramid)) > 0 || !is.numeric(age.pyramid) ||
      min(age.pyramid) < 0 || max(age.pyramid) > 1) {
    stop("Ensure that all input values are non-NA, numeric, and in range.")
  }

  #Re-weighting
  age1 <- age.pyramid[1:(age.break - 1)]
  age1 <- age1 / sum(age1) * young.prop

  age2 <- age.pyramid[age.break:length(age.pyramid)]
  age2 <- age2 / sum(age2) * (1 - young.prop)

  new.age.pyramid <- c(age1, age2)

  return(new.age.pyramid)
}

#' Reduces the Size of the Netstats Object by Removing the Networks it Contains
#'
#' Once `netest` has been built, the network object stored in `netstats` are
#' never used again and take up a lot of space. This function removes them to
#' reduce the stored size of `netstats`. Roughly 10MiB are saved for 100k nodes
#' networks.
#'
#' @param netstats the `netstats` object to be trimmed
#'
#' @return a trimmed `netstats`
#'
#' @export
#'
trim_netstats <- function(netstats) {
  netstats$main <- NULL
  netstats$casl <- NULL
  netstats$inst <- NULL
  return(netstats)
}
