suppressMessages({library(glmmTMB); library(ARTnetData); library(nnet)})
res <- readRDS("base.rds")
geo <- ARTnet.wide[, c("AMIS_ID","city2","DIVCODE","REGCODE")]
geo$city2 <- as.character(geo$city2); geo$DIVCODE <- paste0("D", geo$DIVCODE); geo$REGCODE <- paste0("R", geo$REGCODE)

add_geo <- function(dat) {
  dat <- dat[, setdiff(names(dat), c("city2","DIVCODE","REGCODE")), drop = FALSE]
  out <- merge(dat, geo, by = "AMIS_ID", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(out) == nrow(dat), !anyNA(out$city2)); out
}
drop_geog <- function(f) update(f, . ~ . - geogYN)
re_flat   <- "(1 | city2)"
re_nest   <- "(1 | REGCODE) + (1 | DIVCODE) + (1 | city2)"
add_re <- function(f, re) update(f, as.formula(paste(". ~ . +", re)))

contrast <- function(fit, dat, city) {
  re <- ranef(fit)$cond
  u <- numeric(nrow(dat))
  for (g in names(re)) { v <- setNames(re[[g]][, 1], rownames(re[[g]])); u <- u + v[as.character(dat[[g]])] }
  inc <- dat$city2 == city
  c(contrast = mean(u[inc]) - mean(u[!inc]),
    sd = sapply(VarCorr(fit)$cond, function(x) attr(x, "stddev")))
}

one_model <- function(o, dat, city, fam_kind, tag) {
  f <- if (inherits(o, "lm") && !inherits(o, "glm")) formula(o) else formula(o)
  f0 <- drop_geog(f)
  vars <- all.vars(f)
  dat <- dat[complete.cases(dat[, intersect(vars, names(dat)), drop = FALSE]), , drop = FALSE]
  fam <- switch(fam_kind, count = nbinom2(), binary = binomial(), gauss = gaussian())
  cs <- summary(o)$coefficients
  row <- data.frame(model = tag, n = nrow(dat), n_city = sum(dat$city2 == city),
                    b_cur = cs["geogYN", 1], se_cur = cs["geogYN", 2])
  ctl <- glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
  fit_or_na <- function(ff) tryCatch(suppressWarnings(glmmTMB(ff, data = dat, family = fam, control = ctl,
                                     REML = fam_kind == "gauss")), error = function(e) NULL)
  if (fam_kind == "count") {
    fe <- fit_or_na(f); s <- summary(fe)$coefficients$cond
    row$b_nbfe <- s["geogYN", 1]; row$se_nbfe <- s["geogYN", 2]
    pq <- glm(f, data = dat, family = quasipoisson()); row$se_quasi <- summary(pq)$coefficients["geogYN", 2]
  } else { row$b_nbfe <- NA; row$se_nbfe <- NA; row$se_quasi <- NA }
  for (k in c("flat", "nest")) {
    fit <- fit_or_na(add_re(f0, if (k == "flat") re_flat else re_nest))
    if (is.null(fit)) { row[[paste0("b_", k)]] <- NA; next }
    cc <- contrast(fit, dat, city)
    row[[paste0("b_", k)]] <- cc[["contrast"]]
    row[[paste0("conv_", k)]] <- isTRUE(fit$sdr$pdHess)
    sds <- cc[-1]; names(sds) <- sub("\\..*$", "", sub("^sd\\.", "", names(sds)))
    row[[paste0("sdcity_", k)]] <- sds[["city2"]]
    if (k == "nest") { row$sddiv_nest <- sds[["DIVCODE"]]; row$sdreg_nest <- sds[["REGCODE"]] }
  }
  row
}

city_rows <- function(r) {
  city <- r$spec$label; np <- r$np; out <- list()
  for (L in c("main", "casl", "inst")) {
    for (m in c("joint_model", "joint_concurrent_model", "joint_nm_age_model", "joint_nm_race_model",
                "joint_absdiff_age_model", "joint_absdiff_sqrtage_model")) {
      o <- np[[L]][[m]]; if (is.null(o)) next
      kind <- if (o$family$family == "poisson") "count" else if (o$family$family == "binomial") "binary" else "gauss"
      out[[length(out) + 1]] <- cbind(layer = L, one_model(o, add_geo(o$data), city, kind, m))
    }
    # duration lm, refit data reconstructed from the long file
    o <- np[[L]]$joint_duration_model
    if (!is.null(o)) {
      ld <- np[[L]]$joint_nm_age_model$data
      ong <- ld[(ld$RAI == 1 | ld$IAI == 1) & ld$ongoing2 == 1 & !is.na(ld$duration.time) & ld$duration.time > 0, ]
      stopifnot(nrow(ong) == nobs(o), isTRUE(all.equal(coef(lm(formula(o), data = ong)), coef(o))))
      out[[length(out) + 1]] <- cbind(layer = L, one_model(o, add_geo(ong), city, "gauss", "joint_duration_model"))
    }
    # partner-race multinomial via baseline-category binary logits (ref = Black)
    o <- np[[L]]$joint_prace_model
    if (!is.null(o)) {
      ld <- np[[L]]$joint_nm_age_model$data
      ld$p_race_f <- factor(ld$p_race.cat.num, levels = 1:3)
      cf <- coef(o); se <- summary(o)$standard.errors
      for (k in 2:3) {
        sub <- add_geo(ld[ld$p_race_f %in% c(1, k), ]); sub$yk <- as.integer(sub$p_race_f == k)
        fk <- update(formula(o$terms), yk ~ .)
        bin <- glm(fk, data = sub, family = binomial())
        rr <- one_model(bin, sub, city, "binary", paste0("joint_prace_model_k", k))
        rr$b_cur <- cf[as.character(k), "geogYN"]; rr$se_cur <- se[as.character(k), "geogYN"]
        out[[length(out) + 1]] <- cbind(layer = L, rr)
      }
    }
  }
  cbind(city = r$spec$key, do.call(rbind, lapply(out, function(x) { x[setdiff(c("conv_flat","conv_nest","sdcity_flat","sdcity_nest","sddiv_nest","sdreg_nest","b_flat","b_nest"), names(x))] <- NA; x })))
}
t0 <- Sys.time()
tab <- parallel::mclapply(res, function(r) tryCatch(city_rows(r), error = function(e) conditionMessage(e)), mc.cores = 6)
print(Sys.time() - t0)
bad <- !sapply(tab, is.data.frame); if (any(bad)) print(tab[bad])
tab <- do.call(rbind, tab[!bad])
saveRDS(tab, "coefs.rds")
num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], round, 3)
print(tab[tab$city %in% c("det","mia"), ], row.names = FALSE)
