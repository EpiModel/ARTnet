suppressMessages({library(brms); library(ARTnetData)})
res <- readRDS("base.rds"); np <- res$det$np
geo <- ARTnet.wide[, c("AMIS_ID","city2","DIVCODE","REGCODE")]; geo$city2 <- as.character(geo$city2)
geo$DIVCODE <- paste0("D", geo$DIVCODE); geo$REGCODE <- paste0("R", geo$REGCODE)
labs <- c(atl="Atlanta", det="Detroit", la="Los Angeles", mia="Miami", sea="Seattle", sf="San Francisco")
fitone <- function(o, fam, tag) {
  dat <- merge(o$data[, setdiff(names(o$data), c("city2","DIVCODE","REGCODE"))], geo, by = "AMIS_ID")
  f0 <- update(formula(o), . ~ . - geogYN)
  dat <- dat[complete.cases(dat[, intersect(all.vars(f0), names(dat))]), ]
  f <- update(f0, . ~ . + (1 | REGCODE) + (1 | DIVCODE) + (1 | city2))
  pr <- c(prior(normal(0, 0.5), class = "sd"))
  m <- brm(f, data = dat, family = fam, prior = pr, chains = 4, iter = 1500, warmup = 750, cores = 4,
           seed = 1, refresh = 0, silent = 2, control = list(adapt_delta = 0.95))
  re <- ranef(m, summary = FALSE)
  draws_u <- function(rows) { u <- 0; for (g in c("REGCODE","DIVCODE","city2")) u <- u + re[[g]][, as.character(dat[[g]][rows]), 1]; u }
  out <- t(sapply(labs, function(cl) {
    inc <- dat$city2 == cl
    ui <- draws_u(which(inc)[1]); ur <- rowMeans(draws_u(which(!inc)))
    d <- ui - ur; c(mean = mean(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975)))
  }))
  sds <- posterior_summary(m, variable = c("sd_city2__Intercept","sd_DIVCODE__Intercept","sd_REGCODE__Intercept"))
  list(tag = tag, contrast = out, sds = sds, rhat = max(rhat(m), na.rm = TRUE))
}
r <- list(fitone(np$main$joint_model, poisson(), "main_deg"),
          fitone(np$casl$joint_model, negbinomial(), "casl_deg_nb"),
          fitone(np$inst$joint_model, negbinomial(), "inst_cnt_nb"))
saveRDS(r, "brms.rds")
for (x in r) { cat("==", x$tag, " max rhat", round(x$rhat, 3), "\n"); print(round(x$contrast, 3)); print(round(x$sds, 3)) }
