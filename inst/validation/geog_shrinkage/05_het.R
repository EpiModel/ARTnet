suppressMessages({library(MASS); library(ARTnetData)})
res <- readRDS("base.rds"); co <- readRDS("coefs.rds")
options(width=250)
key <- c("joint_model","joint_concurrent_model","joint_duration_model","joint_nm_race_model")
y <- co[co$model %in% key, c("city","layer","model","b_cur","se_cur","se_quasi","b_flat","b_nest")]
y$z <- y$b_cur / ifelse(is.na(y$se_quasi), y$se_cur, y$se_quasi)
y[sapply(y,is.numeric)] <- lapply(y[sapply(y,is.numeric)], round, 3)
print(y[order(y$layer, y$model, y$city), ], row.names=FALSE)
# LRT of city heterogeneity (city2 factor, 23 df) on Detroit's formulas
geo <- ARTnet.wide[, c("AMIS_ID","city2")]; geo$city2 <- as.character(geo$city2)
np <- res$det$np
lr <- function(o, dat, kind) {
  f0 <- update(formula(o), . ~ . - geogYN); f1 <- update(f0, . ~ . + factor(city2))
  dat <- merge(dat[, setdiff(names(dat), "city2")], geo, by = "AMIS_ID")
  dat <- dat[complete.cases(dat[, intersect(all.vars(f1), names(dat))]), ]
  if (kind == "nb") { m0 <- glm.nb(f0, data = dat); m1 <- glm.nb(f1, data = dat) }
  else if (kind == "pois") { m0 <- glm(f0, data = dat, family = poisson); m1 <- glm(f1, data = dat, family = poisson) }
  else if (kind == "bin") { m0 <- glm(f0, data = dat, family = binomial); m1 <- glm(f1, data = dat, family = binomial) }
  else { m0 <- lm(f0, data = dat); m1 <- lm(f1, data = dat) }
  st <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))); df <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  c(chisq = round(st, 1), df = df, p = signif(pchisq(st, df, lower.tail = FALSE), 2))
}
ong <- function(L) { ld <- np[[L]]$joint_nm_age_model$data; ld[(ld$RAI == 1 | ld$IAI == 1) & ld$ongoing2 == 1 & !is.na(ld$duration.time) & ld$duration.time > 0, ] }
print(rbind(
  main_deg_pois = lr(np$main$joint_model, np$main$joint_model$data, "pois"),
  casl_deg_pois = lr(np$casl$joint_model, np$casl$joint_model$data, "pois"),
  casl_deg_nb   = lr(np$casl$joint_model, np$casl$joint_model$data, "nb"),
  inst_cnt_pois = lr(np$inst$joint_model, np$inst$joint_model$data, "pois"),
  inst_cnt_nb   = lr(np$inst$joint_model, np$inst$joint_model$data, "nb"),
  main_conc     = lr(np$main$joint_concurrent_model, np$main$joint_concurrent_model$data, "bin"),
  casl_conc     = lr(np$casl$joint_concurrent_model, np$casl$joint_concurrent_model$data, "bin"),
  main_dur      = lr(np$main$joint_duration_model, ong("main"), "lm"),
  casl_dur      = lr(np$casl$joint_duration_model, ong("casl"), "lm"),
  main_nmrace   = lr(np$main$joint_nm_race_model, np$main$joint_nm_race_model$data, "bin"),
  casl_nmrace   = lr(np$casl$joint_nm_race_model, np$casl$joint_nm_race_model$data, "bin"),
  inst_nmrace   = lr(np$inst$joint_nm_race_model, np$inst$joint_nm_race_model$data, "bin")))
