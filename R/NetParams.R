
#' Build Network Parameters for Network Models
#'
#' @param epistats Output from \code{\link{build_epistats}}.
#' @param smooth.main.dur If \code{TRUE}, average main durations for oldest and second oldest age
#' groups. FALSE by default.
#'
#' @details
#' \code{build_netparams}is a helper function that constructs the neccessary network parameters
#' for use in building network models with \code{\link{build_netstats}}, building on models estimated
#' using \code{\link{build_epistats}}.
#'
#' The paramter \code{smooth.main.dur} is used when partnership
#' duration and mortality compete in the eldest age group; in such a case mean duration is averaged over
#' the oldest and second oldest age groups (as specified by \code{age.bks} in \code{\link{build_epistats}}).
#' Subsequently, this smoothing is only done if there are three or more age categories specified.
#'
#' @examples
#' epistats <- build_epistats(geog.lvl = "state", geog.cat = "GA", race = TRUE, age.lim = c(15, 65),
#' age.bks = c(20, 30, 40, 50, 60))
#' netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#'
#' @export

build_netparams <- function(epistats, smooth.main.dur = FALSE) {

  ## Data ##
  d <- ARTnet.wide
  l <- ARTnet.long

  #Append Data
  d$geog <- epistats$geog.d
  d$geogYN <- epistats$geogYN.d
  l$geog <- epistats$geog.l
  l$geogYN <- epistats$geogYN.l

  ## Inputs ##
  geog.lvl <- epistats$geog.lvl
  geog.cat <- epistats$geog.cat
  race <- epistats$race
  age.lim <- epistats$age.lim
  age.breaks <- age.bks <- epistats$age.bks


  # 0. Data Processing ------------------------------------------------------

  ## Age Processing ##

  l <- subset(l, age >= age.lim[1] & age <= age.lim[2])
  d <- subset(d, age >= age.lim[1] & age <= age.lim[2])

  l$comb.age <- l$age + l$p_age_imp
  l$diff.age <- abs(l$age - l$p_age_imp)

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

  # summary(d$deg.main)
  # summary(d$deg.casl)

  # recoding to truncate degree
  d$deg.casl <- ifelse(d$deg.casl > 3, 3, d$deg.casl)
  d$deg.main <- ifelse(d$deg.main > 2, 2, d$deg.main)

  d$deg.tot <- d$deg.main + d$deg.casl

  # md <- group_by(d, city2) %>%
  #   summarise(dm = mean(deg.main), dc = mean(deg.casl), dt = mean(deg.tot))
  # print(md, n = nrow(md))
  #
  # round(prop.table(table(d$deg.main, d$deg.casl)), 3)

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

  # head(data.frame(d$ai.part, d$count.mc.part, d$count.oo.part), 25)
  # summary(d$count.oo.part)

  # Truncated OO part
  d$count.oo.part.trunc <- ifelse(d$count.oo.part > 100, 100, d$count.oo.part)
  # summary(d$count.oo.part.trunc)
  # table(d$count.oo.part.trunc)


  if (race == TRUE) {
    ## Race/Ethnicity
    # table(d$race.cat)
    d$race.cat3 <- rep(NA, nrow(d))
    d$race.cat3[d$race.cat == "black"] <- 1
    d$race.cat3[d$race.cat == "hispanic"] <- 2
    d$race.cat3[d$race.cat %in% c("white", "other")] <- 3
    # table(d$race.cat, d$race.cat3)

    # table(l$race.cat, useNA = "always")
    # table(l$p_race.cat, useNA = "always")
    # table(l$race.cat, l$p_race.cat, useNA = "always")

    # l$race.cat3 <- rep(NA, nrow(l))
    l$race.cat3[l$race.cat == "black"] <- 1
    l$race.cat3[l$race.cat == "hispanic"] <- 2
    l$race.cat3[l$race.cat %in% c("white", "other")] <- 3
    # table(l$race.cat3, useNA = "always")

    # table(l$p_race.cat, useNA = "always")
    l$p_race.cat3 <- rep(NA, nrow(l))
    l$p_race.cat3[l$p_race.cat == "black"] <- 1
    l$p_race.cat3[l$p_race.cat == "hispanic"] <- 2
    l$p_race.cat3[l$p_race.cat %in% c("white", "other")] <- 3
    # table(l$p_race.cat3, useNA = "always")

    # redistribute NAs in proportion to non-missing partner races
    probs <- prop.table(table(l$race.cat3, l$p_race.cat3), 1)

    imp_black <- which(is.na(l$p_race.cat3) & l$race.cat3 == 1)
    l$p_race.cat3[imp_black] <- sample(1:3, length(imp_black), TRUE, probs[1, ])

    imp_hisp <- which(is.na(l$p_race.cat3) & l$race.cat3 == 2)
    l$p_race.cat3[imp_hisp] <- sample(1:3, length(imp_hisp), TRUE, probs[2, ])

    imp_white <- which(is.na(l$p_race.cat3) & l$race.cat3 == 3)
    l$p_race.cat3[imp_white] <- sample(1:3, length(imp_white), TRUE, probs[3, ])

    # table(l$race.cat3, l$p_race.cat3, useNA = "always")
  }

  ## HIV status

  l$p_hiv2 <- ifelse(l$p_hiv == 1, 1, 0)
  table(l$p_hiv, l$p_hiv2, useNA = "always")

  hiv.combo <- rep(NA, nrow(l))
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5
  # table(hiv.combo, useNA = "always")

  l$hiv.concord.pos <- ifelse(hiv.combo == 2, 1, 0)
  # table(l$hiv.concord.pos)


  ## Setup output list ##

  out <- list()


  # 1. Main Model -----------------------------------------------------------

  out$main <- list()
  lmain <- l[l$ptype == 1, ]


  ## edges ----
  if (geog.lvl != "all") {
    mod <- glm(deg.main ~ 1,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$md.main <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geog,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$md.main <- as.numeric(pred)
  }


  ## nodematch("age.grp") ----

  lmain$index.age.grp <- cut(lmain$age, age.breaks, labels = FALSE)
  lmain$part.age.grp <- cut(as.numeric(lmain$p_age_imp), age.breaks, labels = FALSE)
  # data.frame(lmain$age, lmain$index.age.grp, lmain$p_age_imp, lmain$part.age.grp)

  lmain$same.age.grp <- ifelse(lmain$index.age.grp == lmain$part.age.grp, 1, 0)

  if (geog.lvl != "all") {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = lmain, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geog + index.age.grp,
               data = lmain, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  lmain$ad <- abs(lmain$age - lmain$p_age_imp)
  lmain$ad.sr <- abs(sqrt(lmain$age) - sqrt(lmain$p_age_imp))

  if (geog.lvl != "all") {
    mod <- lm(ad ~ 1, data = lmain)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geog, data = lmain)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (geog.lvl != "all") {
    mod <- lm(ad.sr ~ 1, data = lmain)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geog, data = lmain)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE)

  if (geog.lvl != "all") {
    mod <- glm(deg.main ~ + age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geog + age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.age.grp <- as.numeric(pred)
  }


  ## nodematch("race", diff = TRUE) ----

  if (race == TRUE) {
    lmain$same.race <- ifelse(lmain$race.cat3 == lmain$p_race.cat3, 1, 0)

    if (geog.lvl != "all") {
      mod <- glm(same.race ~ as.factor(race.cat3),
                 data = lmain, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog + as.factor(race.cat3),
                 data = lmain, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (geog.lvl != "all") {
      mod <- glm(same.race ~ 1,
                 data = lmain, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog,
                 data = lmain, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nm.race_diffF <- as.numeric(pred)
    }


    ## nodefactor("race") ----

    if (geog.lvl != "all") {
      mod <- glm(deg.main ~ as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(deg.main ~ geog + as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$main$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("deg.casl") ----

  if (geog.lvl != "all") {
    mod <- glm(deg.main ~ deg.casl,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.casl = sort(unique(d$deg.casl)))
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.deg.casl <- as.numeric(pred)

    deg.casl.dist <- prop.table(table(d$deg.casl[d$geog == geog.cat]))
    out$main$deg.casl.dist <- as.numeric(deg.casl.dist)
  } else {
    mod <- glm(deg.main ~ geog + deg.casl,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.casl = sort(unique(d$deg.casl)))
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.deg.casl <- as.numeric(pred)

    deg.casl.dist <- prop.table(table(d$deg.casl[d$geog == geog.cat]))
    out$main$deg.casl.dist <- as.numeric(deg.casl.dist)
  }


  ## concurrent ----

  if (geog.lvl != "all") {
    mod <- glm(deg.main.conc ~ 1,
               data = d, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$concurrent <- as.numeric(pred)
  } else {
    mod <- glm(deg.main.conc ~ geog,
               data = d, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$concurrent <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (geog.lvl != "all") {
    mod <- glm(deg.main ~ hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(deg.main ~ geog + hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$main$nf.diag.status <- as.numeric(pred)
  }


  ## Durations ----

  # overall
  durs.main <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()

  # create city weights
  durs.main.geo <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    filter(geog == geog.cat) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()

  # city-specific weight based on ratio of medians
  wt <- durs.main.geo$median.dur/durs.main$median.dur

  # The weekly dissolution rate is function of the mean of the geometric distribution
  # which relates to the median as:
  durs.main$rates.main.adj <- 1 - (2^(-1/(wt*durs.main$median.dur)))

  # Mean duration associated with a geometric distribution that median:
  durs.main$mean.dur.adj <- 1/(1 - (2^(-1/(wt*durs.main$median.dur))))
  out$main$durs.main.homog <- durs.main

  # stratified by age

  # first, non-matched by age group
  durs.main.nonmatch <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 0) %>%
    # group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()
  durs.main.nonmatch$index.age.grp <- 0

  # then, matched within age-groups
  durs.main.matched <- lmain %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 1) %>%
    group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()
  durs.main.matched

  durs.main.all <- rbind(durs.main.nonmatch, durs.main.matched)

  durs.main.all$rates.main.adj <- 1 - (2^(-1/(wt*durs.main.all$median.dur)))
  durs.main.all$mean.dur.adj <- 1/(1 - (2^(-1/(wt*durs.main.all$median.dur))))

  durs.main.all <- durs.main.all[, c(3, 1, 2, 4, 5)]
  out$main$durs.main.byage <- durs.main.all

  if (smooth.main.dur == TRUE) {
    n2 <- nrow(durs.main.all)
    n1 <- n2-1
    if (n2 > 3) {
      out$main$durs.main.byage$mean.dur.adj[n2] <- mean(out$main$durs.main.byage$mean.dur.adj[n1:n2])
    }
  }


  # 2. Casual Model ---------------------------------------------------------

  out$casl <- list()
  lcasl <- l[l$ptype == 2, ]


  ## edges ----

  if (geog.lvl != "all") {
    mod <- glm(deg.casl ~ 1,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$md.casl <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geog,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$md.casl <- as.numeric(pred)
  }

  ## nodematch("age.grp") ----

  lcasl$index.age.grp <- cut(lcasl$age, age.breaks, labels = FALSE)
  lcasl$part.age.grp <- cut(as.numeric(lcasl$p_age_imp), age.breaks, labels = FALSE)
  # data.frame(lcasl$age, lcasl$index.age.grp, lcasl$p_age_imp, lcasl$part.age.grp)

  lcasl$same.age.grp <- ifelse(lcasl$index.age.grp == lcasl$part.age.grp, 1, 0)

  if (geog.lvl != "all") {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = lcasl, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geog + index.age.grp,
               data = lcasl, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  lcasl$ad <- abs(lcasl$age - lcasl$p_age_imp)
  lcasl$ad.sr <- abs(sqrt(lcasl$age) - sqrt(lcasl$p_age_imp))

  if (geog.lvl != "all") {
    mod <- lm(ad ~ 1, data = lcasl)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geog, data = lcasl)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (geog.lvl != "all") {
    mod <- lm(ad.sr ~ 1, data = lcasl)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geog, data = lcasl)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE)

  if (geog.lvl != "all") {
    mod <- glm(deg.casl ~ age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geog + age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.age.grp <- as.numeric(pred)
  }


  if (race == TRUE) {

    ## nodematch("race") ----

    # prop.table(table(lcasl$race.cat3, lcasl$p_race.cat3), 1)

    lcasl$same.race <- ifelse(lcasl$race.cat3 == lcasl$p_race.cat3, 1, 0)

    if (geog.lvl != "all") {
      mod <- glm(same.race ~ as.factor(race.cat3),
                 data = lcasl, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog + as.factor(race.cat3),
                 data = lcasl, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (geog.lvl == "all") {
      mod <- glm(same.race ~ 1,
                 data = lcasl, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog,
                 data = lcasl, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nm.race_diffF <- as.numeric(pred)
    }

    ## nodefactor("race", diff = TRUE) ----

    if (geog.lvl == "all") {
      mod <- glm(deg.casl ~ as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(deg.casl ~ geog + as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$casl$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("deg.main") ----

  if (geog.lvl == "all") {
    mod <- glm(deg.casl ~ deg.main,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.main = 0:2)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.deg.main <- as.numeric(pred)

    deg.main.dist <- prop.table(table(d$deg.main[d$geog == geog.cat]))
    out$casl$deg.main.dist <- as.numeric(deg.main.dist)
  } else {
    mod <- glm(deg.casl ~ geog + deg.main,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.main = 0:2)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.deg.main <- as.numeric(pred)

    deg.main.dist <- prop.table(table(d$deg.main[d$geog == geog.cat]))
    out$casl$deg.main.dist <- as.numeric(deg.main.dist)
  }


  ## concurrent ----

  if (geog.lvl == "all") {
    mod <- glm(deg.casl.conc ~ 1,
               data = d, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$concurrent <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl.conc ~ geog,
               data = d, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$concurrent <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (geog.lvl == "all") {
    mod <- glm(deg.casl ~ hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(deg.casl ~ geog + hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")

    out$casl$nf.diag.status <- as.numeric(pred)
  }


  ## Durations ----

  # overall
  durs.casl <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()

  # create city weights
  durs.casl.geo <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    filter(geog == geog.cat) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()

  # city-specific weight based on ratio of medians
  wt <- durs.casl.geo$median.dur/durs.casl$median.dur

  # The weekly dissolution rate is function of the mean of the geometric distribution
  # which relates to the median as:
  durs.casl$rates.casl.adj <- 1 - (2^(-1/(wt*durs.casl$median.dur)))

  # Mean duration associated with a geometric distribution that median:
  durs.casl$mean.dur.adj <- 1/(1 - (2^(-1/(wt*durs.casl$median.dur))))
  out$casl$durs.casl.homog <- durs.casl

  # stratified by age

  # first, non-matched by age group
  durs.casl.nonmatch <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(ongoing2 == 1) %>%
    filter(same.age.grp == 0) %>%
    # group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()
  durs.casl.nonmatch$index.age.grp <- 0

  # then, matched within age-groups
  durs.casl.matched <- lcasl %>%
    filter(RAI == 1 | IAI == 1) %>%
    filter(index.age.grp < 6) %>%
    filter(same.age.grp == 1) %>%
    filter(ongoing2 == 1) %>%
    group_by(index.age.grp) %>%
    summarise(mean.dur = mean(duration, na.rm = TRUE),
              median.dur = median(duration, na.rm = TRUE)) %>%
    as.data.frame()
  # durs.casl.matched

  durs.casl.all <- rbind(durs.casl.nonmatch, durs.casl.matched)
  # durs.casl.all

  durs.casl.all$rates.casl.adj <- 1 - (2^(-1/(wt*durs.casl.all$median.dur)))
  durs.casl.all$mean.dur.adj <- 1/(1 - (2^(-1/(wt*durs.casl.all$median.dur))))

  durs.casl.all <- durs.casl.all[, c(3, 1, 2, 4, 5)]
  out$casl$durs.casl.byage <- durs.casl.all


  # 3. One-off Model --------------------------------------------------------

  out$inst <- list()
  linst <- l[l$ptype == 3, ]

  ## edges ----

  head(d$count.oo.part, 25)
  # summary(d$count.oo.part)

  # weekly rate
  d$rate.oo.part <- d$count.oo.part/52
  # summary(d$rate.oo.part)

  if (geog.lvl != "all") {
    mod <- glm(count.oo.part ~ 1,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$md.inst <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geog,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$md.inst <- as.numeric(pred)
  }


  ## nodematch("age.grp") ----

  linst$index.age.grp <- cut(linst$age, age.breaks, labels = FALSE)
  linst$part.age.grp <- cut(as.numeric(linst$p_age_imp), age.breaks, labels = FALSE)
  # data.frame(linst$age, linst$index.age.grp, linst$p_age_imp, linst$part.age.grp)

  linst$same.age.grp <- ifelse(linst$index.age.grp == linst$part.age.grp, 1, 0)

  if (geog.lvl == "all") {
    mod <- glm(same.age.grp ~ index.age.grp,
               data = linst, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$nm.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(same.age.grp ~ geog + index.age.grp,
               data = linst, family = binomial())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, index.age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$nm.age.grp <- as.numeric(pred)
  }


  ## absdiff("age") ----

  linst$ad <- abs(linst$age - linst$p_age_imp)
  linst$ad.sr <- abs(sqrt(linst$age) - sqrt(linst$p_age_imp))

  if (geog.lvl == "all") {
    mod <- lm(ad ~ 1, data = linst)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.age <- as.numeric(pred)
  } else {
    mod <- lm(ad ~ geog, data = linst)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.age <- as.numeric(pred)
  }


  ## absdiff("sqrt.age") ----

  if (geog.lvl == "all") {
    mod <- lm(ad.sr ~ 1, data = linst)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.sqrt.age <- as.numeric(pred)
  } else {
    mod <- lm(ad.sr ~ geog, data = linst)
    # summary(mod)

    dat <- data.frame(geog = geog.cat)
    pred <- predict(mod, newdata = dat, type = "response")

    out$inst$absdiff.sqrt.age <- as.numeric(pred)
  }


  ## nodefactor("age.grp") ----

  d$age.grp <- cut(d$age, age.breaks, labels = FALSE)

  if (geog.lvl == "all") {
    mod <- glm(count.oo.part ~ age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.age.grp <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geog + age.grp + sqrt(age.grp),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, age.grp = 1:5)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.age.grp <- as.numeric(pred)
  }


  if (race == TRUE) {

    ## nodematch("race", diff = TRUE) ----

    prop.table(table(linst$race.cat3, linst$p_race.cat3), 1)

    linst$same.race <- ifelse(linst$race.cat3 == linst$p_race.cat3, 1, 0)

    if (geog.lvl != "all") {
      mod <- glm(same.race ~ as.factor(race.cat3),
                 data = linst, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog + as.factor(race.cat3),
                 data = linst, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race <- as.numeric(pred)
    }


    ## nodematch("race", diff = FALSE) ----

    if (geog.lvl == "all") {
      mod <- glm(same.race ~ 1,
                 data = linst, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race_diffF <- as.numeric(pred)
    } else {
      mod <- glm(same.race ~ geog,
                 data = linst, family = binomial())
      # summary(mod)

      dat <- data.frame(geog = geog.cat)
      pred <- predict(mod, newdata = dat, type = "response")

      out$inst$nm.race_diffF <- as.numeric(pred)
    }


    ## nodefactor("race") ----

    if (geog.lvl == "all") {
      mod <- glm(count.oo.part ~ as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")/52

      out$inst$nf.race <- as.numeric(pred)
    } else {
      mod <- glm(count.oo.part ~ geog + as.factor(race.cat3),
                 data = d, family = poisson())
      # summary(mod)

      dat <- data.frame(geog = geog.cat, race.cat3 = 1:3)
      pred <- predict(mod, newdata = dat, type = "response")/52

      out$inst$nf.race <- as.numeric(pred)
    }
  }

  ## nodefactor("risk.grp") ----

  # city-specific wts
  wt <- mean(d$rate.oo.part[d$geog == geog.cat], na.rm = TRUE)/mean(d$rate.oo.part, na.rm = TRUE)
  wt.rate <- d$rate.oo.part * wt

  nquants <- 5
  oo.quants <- rep(NA, nquants)
  sr <- sort(wt.rate)
  qsize <- floor(length(sr) / nquants)
  for (i in 1:nquants) {
    if (i == 1) {
      oo.quants[i] <- mean(sr[1:qsize])
    } else if (i > 1 & i < nquants) {
      oo.quants[i] <- mean(sr[(((i - 1)*qsize) + 1):(i*qsize)])
    } else if (i == nquants) {
      oo.quants[i] <- mean(sr[(((i - 1)*qsize) + 1):length(sr)])
    }
  }

  # Weekly acquisition rate
  # oo.quants

  # Yearly OO partner count
  # oo.quants * 52

  # New OO partner every X years
  # 1/(oo.quants * 52)

  # New OO partner every weeks
  # 1/oo.quants

  # Save it
  out$inst$nf.risk.grp <- oo.quants


  ## nodefactor("deg.tot") ----

  d$deg.tot3 <- ifelse(d$deg.tot >= 3, 3, d$deg.tot)
  # table(d$deg.tot3)

  deg.tot.dist <- prop.table(table(d$deg.tot3[d$geog == geog.cat]))
  out$inst$deg.tot.dist <- as.numeric(deg.tot.dist)

  if (geog.lvl != "all") {
    mod <- glm(count.oo.part ~ deg.tot3 + sqrt(deg.tot3),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.tot3 = 0:3)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.deg.tot <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geog + deg.tot3 + sqrt(deg.tot3),
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, deg.tot3 = 0:3)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.deg.tot <- as.numeric(pred)
  }


  ## nodefactor("diag.status") ----

  if (geog.lvl != "all") {
    mod <- glm(count.oo.part ~ hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.diag.status <- as.numeric(pred)
  } else {
    mod <- glm(count.oo.part ~ geog + hiv2,
               data = d, family = poisson())
    # summary(mod)

    dat <- data.frame(geog = geog.cat, hiv2 = 0:1)
    pred <- predict(mod, newdata = dat, type = "response")/52

    out$inst$nf.diag.status <- as.numeric(pred)
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
  # table(d$nRAIpart, useNA = "always")

  d <- l %>%
    filter(IAI == 1) %>%
    group_by(AMIS_ID) %>%
    count() %>%
    rename(nIAIpart = n) %>%
    right_join(d, by = "AMIS_ID") %>%
    as.data.frame()
  d$nIAIpart <- ifelse(is.na(d$nIAIpart), 0, d$nIAIpart)
  # table(d$nIAIpart, useNA = "always")

  # table(d$nRAIpart, d$nIAIpart, useNA = "always")

  # default NA for no AI
  roletype <- rep(NA, nrow(d))
  roletype[d$nRAIpart == 0 & d$nIAIpart > 0] <- 0
  roletype[d$nIAIpart == 0 & d$nRAIpart > 0] <- 1
  roletype[d$nIAIpart > 0 & d$nRAIpart > 0] <- 2
  # table(roletype, useNA = "always")

  out$all$role.type <- prop.table(table(roletype))



  # SAVE --------------------------------------------------------------------

  return(out)
}
