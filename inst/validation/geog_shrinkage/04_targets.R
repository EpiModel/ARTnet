suppressMessages({devtools::load_all("/Users/sjennes/git/ARTnet", quiet = TRUE); library(nnet)})
res <- readRDS("base.rds"); co <- readRDS("coefs.rds")
expect_mort <- 0.000478213; dur_cap <- 0.9/(2*expect_mort); seed <- 20261107L

registerS3method("predict", "shifted_multinom", function(object, newdata, type = "class", ...) {
  delta <- attr(object, "delta"); class(object) <- setdiff(class(object), "shifted_multinom")
  p <- predict(object, newdata = newdata, type = type, ...)
  if (type != "probs") stop("shifted_multinom: only type = 'probs'")
  p <- sweep(as.matrix(p), 2, exp(delta), "*"); p / rowSums(p)
})

# Dirichlet-multinomial shrinkage of a city's raw cross-layer degree distribution toward
# the national one; ICC from the one-way ANOVA estimator on each category indicator, over city2.
dm_shrink <- function(x, grp, city, levels) {
  n_g <- table(grp); C <- length(n_g); N <- length(x); n0 <- (N - sum(n_g^2) / N) / (C - 1)
  icc <- sapply(levels, function(k) {
    y <- as.numeric(x == k); m <- tapply(y, grp, mean)
    msb <- sum(n_g * (m - mean(y))^2) / (C - 1); msw <- sum(tapply(y, grp, function(v) sum((v - mean(v))^2))) / (N - C)
    max(0, (msb - msw) / (msb + (n0 - 1) * msw))
  })
  rho <- mean(icc); a0 <- if (rho > 0) (1 - rho) / rho else Inf
  nc <- sum(grp == city); w <- nc / (nc + a0)
  pc <- as.numeric(prop.table(table(factor(x[grp == city], levels))))
  pn <- as.numeric(prop.table(table(factor(x, levels))))
  list(p = w * pc + (1 - w) * pn, w = w, rho = rho, raw = pc, nat = pn)
}

shift_np <- function(np, cc, col) {
  for (L in c("main", "casl", "inst")) {
    rows <- cc[cc$layer == L, ]
    for (i in seq_len(nrow(rows))) {
      m <- rows$model[i]; b <- rows[[col]][i]
      if (grepl("^joint_prace_model_k", m)) next
      if (is.na(b)) stop("missing ", col, " for ", L, " ", m)
      np[[L]][[m]]$coefficients[["geogYN"]] <- b
    }
    pm <- np[[L]]$joint_prace_model
    if (!is.null(pm)) {
      k2 <- rows[rows$model == "joint_prace_model_k2", ]; k3 <- rows[rows$model == "joint_prace_model_k3", ]
      attr(pm, "delta") <- c(0, k2[[col]] - k2$b_cur, k3[[col]] - k3$b_cur)
      class(pm) <- c("shifted_multinom", class(pm)); np[[L]]$joint_prace_model <- pm
    }
  }
  np
}

summ <- function(ns, n, key, scen) {
  dur <- function(lay) {
    d <- ns[[lay]]$diss.byage$duration; nm <- ns[[lay]]$nodematch_age.grp
    w <- c(ns[[lay]]$edges - sum(nm), nm); c(nonmatch = d[1], wmean = sum(w * d) / sum(w), max = max(d))
  }
  dm <- dur("main"); dc <- dur("casl")
  data.frame(city = key, scen = scen,
    main = 1000 * ns$main$edges / n, casl = 1000 * ns$casl$edges / n, inst = 1000 * ns$inst$edges / n,
    conc_main = 1000 * ns$main$concurrent / n, conc_casl = 1000 * ns$casl$concurrent / n,
    dur_main = dm[["wmean"]], dur_main_nm = dm[["nonmatch"]], dur_main_max = dm[["max"]],
    dur_casl = dc[["wmean"]], dur_casl_nm = dc[["nonmatch"]],
    nmr_main_B = ns$main$nodematch_race[1], nmr_main_H = ns$main$nodematch_race[2], nmr_main_W = ns$main$nodematch_race[3],
    nmr_casl_B = ns$casl$nodematch_race[1], nmr_casl_H = ns$casl$nodematch_race[2], nmr_casl_W = ns$casl$nodematch_race[3],
    nmr_inst_B = ns$inst$nodematch_race[1], nmr_inst_H = ns$inst$nodematch_race[2], nmr_inst_W = ns$inst$nodematch_race[3],
    nmr_main_share = ns$main$nodematch_race_diffF / ns$main$edges, nmr_casl_share = ns$casl$nodematch_race_diffF / ns$casl$edges,
    nmr_inst_share = ns$inst$nodematch_race_diffF / ns$inst$edges)
}

out <- parallel::mclapply(res, function(r) {
  s <- r$spec; cc <- co[co$city == s$key, ]
  run <- function(np) suppressMessages(build_netstats(r$ep, np, expect.mort = expect_mort, network.size = s$n,
           race.prop = s$race.prop, method = "joint", seed = seed, duration.cap = dur_cap))
  d <- r$np$main$joint_model$data
  dmc <- dm_shrink(d$deg.casl, as.character(d$city2), s$label, 0:3)
  dmm <- dm_shrink(d$deg.main, as.character(d$city2), s$label, 0:2)
  np_dd <- function(np) { np$main$deg.casl.dist <- dmc$p; np$casl$deg.main.dist <- dmm$p; np }
  np_cur_chk <- shift_np(r$np, cc, "b_cur")
  rows <- rbind(summ(r$ns, s$n, s$key, "current"),
                summ(run(np_cur_chk), s$n, s$key, "current_rebuilt"),
                summ(run(np_dd(r$np)), s$n, s$key, "degdist_only"),
                summ(run(shift_np(r$np, cc, "b_flat")), s$n, s$key, "flat"),
                summ(run(shift_np(r$np, cc, "b_nest")), s$n, s$key, "nested"),
                summ(run(np_dd(shift_np(r$np, cc, "b_nest"))), s$n, s$key, "nested_dd"),
                summ(run(np_dd(shift_np(r$np, cc, "b_flat"))), s$n, s$key, "flat_dd"))
  list(rows = rows, dd = list(casl = dmc, main = dmm))
}, mc.cores = 6)
saveRDS(out, "targets.rds")
tab <- do.call(rbind, lapply(out, `[[`, "rows"))
options(width = 250); num <- sapply(tab, is.numeric); t2 <- tab; t2[num] <- lapply(t2[num], round, 3)
print(t2[, 1:12], row.names = FALSE)
for (o in out) cat(o$rows$city[1], "deg.casl w", round(o$dd$casl$w, 3), "rho", round(o$dd$casl$rho, 4), " deg.main w", round(o$dd$main$w, 3), "rho", round(o$dd$main$rho, 4), "\n")
