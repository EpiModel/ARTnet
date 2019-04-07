
##
## Epidemic parameters analysis for ARTnet Data
##

## Packages ##
library("tidyverse")
library("ARTnetData")

## Inputs ##
city_name <- "Atlanta"

## Data ##
d <- ARTnet.wide
l <- ARTnet.long

## Derivatives ##
coef_name <- paste0("city2", city_name)

names(l)


# Data Processing ---------------------------------------------------------

# Age
table(l$age, useNA = "always")
table(l$p_age_imp, useNA = "always")

l$comb.age <- l$age + l$p_age_imp
l$diff.age <- abs(l$age - l$p_age_imp)


# Race
table(d$race.cat)
d$race.cat3 <- rep(NA, nrow(d))
d$race.cat3[d$race.cat == "black"] <- 0
d$race.cat3[d$race.cat == "hispanic"] <- 1
d$race.cat3[d$race.cat %in% c("white", "other")] <- 2
table(d$race.cat, d$race.cat3)

table(l$race.cat, useNA = "always")
table(l$p_race.cat, useNA = "always")
table(l$race.cat, l$p_race.cat, useNA = "always")

l$race.cat3 <- rep(NA, nrow(l))
l$race.cat3[l$race.cat == "black"] <- 0
l$race.cat3[l$race.cat == "hispanic"] <- 1
l$race.cat3[l$race.cat %in% c("white", "other")] <- 2
table(l$race.cat3, useNA = "always")

table(l$p_race.cat, useNA = "always")
l$p_race.cat3 <- rep(NA, nrow(l))
l$p_race.cat3[l$p_race.cat == "black"] <- 0
l$p_race.cat3[l$p_race.cat == "hispanic"] <- 1
l$p_race.cat3[l$p_race.cat %in% c("white", "other")] <- 2
table(l$p_race.cat3, useNA = "always")

# redistribute NAs in proportion to non-missing partner races
probs <- prop.table(table(l$race.cat3, l$p_race.cat3), 1)

imp_black <- which(is.na(l$p_race.cat3) & l$race.cat3 == 0)
l$p_race.cat3[imp_black] <- sample(0:2, length(imp_black), TRUE, probs[1, ])

imp_hisp <- which(is.na(l$p_race.cat3) & l$race.cat3 == 1)
l$p_race.cat3[imp_hisp] <- sample(0:2, length(imp_hisp), TRUE, probs[2, ])

imp_white <- which(is.na(l$p_race.cat3) & l$race.cat3 == 2)
l$p_race.cat3[imp_white] <- sample(0:2, length(imp_white), TRUE, probs[3, ])

table(l$race.cat3, l$p_race.cat3, useNA = "always")

l$race.combo <- rep(NA, nrow(l))
l$race.combo[l$race.cat3 == 0 & l$p_race.cat3 == 0] <- 0
l$race.combo[l$race.cat3 == 0 & l$p_race.cat3 %in% 1:2] <- 1
l$race.combo[l$race.cat3 == 1 & l$p_race.cat3 %in% c(0, 2)] <- 2
l$race.combo[l$race.cat3 == 1 & l$p_race.cat3 == 1] <- 3
l$race.combo[l$race.cat3 == 2 & l$p_race.cat3 %in% 0:1] <- 4
l$race.combo[l$race.cat3 == 2 & l$p_race.cat3 == 2] <- 5

table(l$race.combo)
l <- select(l, -c(race.cat3, p_race.cat3))


# HIV
l$p_hiv2 <- ifelse(l$p_hiv == 1, 1, 0)
table(l$hiv2, l$p_hiv, useNA = "always")

hiv.combo <- rep(NA, nrow(l))
hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5
table(hiv.combo, useNA = "always")

l$hiv.combo <- hiv.combo
l$hiv.concord <- ifelse(hiv.combo %in% 0:1, 1, 0)
l$hiv.concord.pos <- ifelse(hiv.combo == 2, 1, 0)
table(l$hiv.concord)
table(l$hiv.concord.pos)

# PrEP
names(d)
table(d$PREP_REVISED, useNA = "always")
table(d$artnetPREP_CURRENT, useNA = "always")
table(d$PREP_REVISED, d$artnetPREP_CURRENT, useNA = "always")
d$prep <- ifelse(d$artnetPREP_CURRENT == 0 | is.na(d$artnetPREP_CURRENT), 0, 1)
table(d$prep, useNA = "always")

dlim <- select(d, c(AMIS_ID, survey.year, prep))
l <- left_join(l, dlim, by = "AMIS_ID")

# city
l$cityYN <- ifelse(l$city2 == city_name, 1, 0)
d$cityYN <- ifelse(d$city2 == city_name, 1, 0)


# Act Rates ---------------------------------------------------------------

# acts/per week/per partnership for main and casual partnerships

# Pull Data
la <- select(l, ptype, duration, comb.age, city = cityYN,
             race.combo, RAI, IAI, hiv.concord.pos, prep,
             acts = anal.acts.week, cp.acts = anal.acts.week.cp) %>%
  filter(ptype %in% 1:2) %>%
  filter(RAI == 1 | IAI == 1)
la <- select(la, -c(RAI, IAI))
head(la, 25)
nrow(la)


# Poisson Model
acts.mod <- glm(floor(acts*52) ~ duration + I(duration^2) + as.factor(race.combo) +
                as.factor(ptype) + duration*as.factor(ptype) + comb.age + I(comb.age^2) +
                hiv.concord.pos + city,
                family = poisson(), data = la)
summary(acts.mod)

x <- expand.grid(duration = 100,
                 ptype = 2,
                 race.combo = 0:5,
                 comb.age = 40,
                 hiv.concord.pos = 0,
                 city = 1)
pred <- predict(acts.mod, newdata = x, type = "response", se.fit = FALSE)
pred.acts <- cbind(x, pred = pred/52)
pred.acts


# Condom Use // Main Casual -----------------------------------------------

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(la$acts, la$cp.acts)
plot(la$acts, la$cp.acts, xlim = c(0, 10), ylim = c(0, 10))

summary(la$cp.acts)

la$prob.cond <- la$cp.acts / la$acts
head(la, 25)

table(la$acts, useNA = "always")

hist(la$prob.cond)
table(la$prob.cond)
summary(la$prob.cond)
summary(la$prob.cond[la$ptype == 1])
summary(la$prob.cond[la$ptype == 2])

la$any.cond <- ifelse(la$prob.cond > 0, 1, 0)
la$never.cond <- ifelse(la$prob.cond == 0, 1, 0)
table(la$never.cond)

cond.mc.mod <- glm(any.cond ~ duration + I(duration^2) + as.factor(race.combo) +
                     as.factor(ptype) + duration*as.factor(ptype) + comb.age + I(comb.age^2) +
                     hiv.concord.pos + prep + city,
                family = binomial(), data = la)
summary(cond.mc.mod)

x <- expand.grid(duration = 50,
                 ptype = 2,
                 race.combo = 0:5,
                 comb.age = 40,
                 hiv.concord.pos = 0,
                 prep = 0:1,
                 city = 1)
pred <- predict(cond.mc.mod, newdata = x, type = "response")
pred.cond <- cbind(x, pred)
pred.cond


# Condom Use // Inst ------------------------------------------------------

lb <- select(l, ptype, comb.age, city = cityYN,
         race.combo, hiv.concord.pos, prep,
         RAI, IAI, RECUAI, INSUAI) %>%
  filter(ptype == 3) %>%
  filter(RAI == 1 | IAI == 1)
head(lb, 40)

table(lb$RAI, lb$RECUAI, useNA = "always")
table(lb$IAI, lb$RAI)

lb$prob.cond <- rep(NA, nrow(lb))
lb$prob.cond[lb$RAI == 1 & lb$IAI == 0] <- lb$RECUAI[lb$RAI == 1 & lb$IAI == 0] /
                                              lb$RAI[lb$RAI == 1 & lb$IAI == 0]
lb$prob.cond[lb$RAI == 0 & lb$IAI == 1] <- lb$INSUAI[lb$RAI == 0 & lb$IAI == 1] /
                                              lb$IAI[lb$RAI == 0 & lb$IAI == 1]
lb$prob.cond[lb$RAI == 1 & lb$IAI == 1] <- (lb$RECUAI[lb$RAI == 1 & lb$IAI == 1] +
                                              lb$INSUAI[lb$RAI == 1 & lb$IAI == 1]) /
                                           (lb$RAI[lb$RAI == 1 & lb$IAI == 1] +
                                              lb$IAI[lb$RAI == 1 & lb$IAI == 1])
lb$prob.cond[which(lb$prob.cond == 0.5)] <- 0
lb$prob.cond[which(lb$prob.cond %in% c(88, 99, 44))] <- NA
table(lb$prob.cond)
lb <- select(lb, -c(RAI, IAI, RECUAI, INSUAI))
head(lb, 40)

cond.oo.mod <- glm(prob.cond ~ as.factor(race.combo) +
                     comb.age + I(comb.age^2) +
                     hiv.concord.pos + prep + city,
                   family = binomial(), data = lb)
summary(cond.oo.mod)

x <- expand.grid(race.combo = 0:5,
                 comb.age = 40,
                 hiv.concord.pos = 0,
                 prep = 0:1,
                 city = 1)
pred <- predict(cond.oo.mod, newdata = x, type = "response")
pred.cond <- cbind(x, pred)
pred.cond



# Init HIV Status ---------------------------------------------------------

d1 <- select(d, race.cat3, cityYN, age, hiv2)

hiv.mod <- glm(hiv2 ~ age + cityYN + as.factor(race.cat3) + cityYN*as.factor(race.cat3), 
               data = d1, family = binomial())
summary(hiv.mod)
x <- expand.grid(age = 15:65, race.cat3 = 0:2, cityYN = 0:1)
pred <- predict(hiv.mod, newdata = x)
pred <- cbind(x, est = plogis(pred))
pred

ggplot(pred, aes(age, est, color = as.factor(race.cat3), lty = as.factor(cityYN))) +
  geom_line() +
  scale_color_viridis_d() +
  theme_minimal()



# Save Out File -----------------------------------------------------------

out <- list()

out$acts.mod <- acts.mod
out$cond.mc.mod <- cond.mc.mod
out$cond.oo.mod <- cond.oo.mod
out$hiv.mod <- hiv.mod
fn <- paste("data/artnet.EpiStats", gsub(" ", "", city_name), "rda", sep = ".")
saveRDS(out, file = fn)
