suppressMessages({devtools::load_all("/Users/sjennes/git/ARTnet", quiet = TRUE)})
res <- readRDS("base.rds"); em <- 0.000478213; cap <- 0.9/(2*em)
summ <- function(ns, n) { d <- function(l) { x <- ns[[l]]$diss.byage$duration; nm <- ns[[l]]$nodematch_age.grp; w <- c(ns[[l]]$edges - sum(nm), nm); sum(w*x)/sum(w) }
  c(main = 1000*ns$main$edges/n, casl = 1000*ns$casl$edges/n, inst = 1000*ns$inst$edges/n, dur_main = d("main"), dur_casl = d("casl"),
    nmr_main = ns$main$nodematch_race_diffF/ns$main$edges, nmr_casl = ns$casl$nodematch_race_diffF/ns$casl$edges, nmr_inst = ns$inst$nodematch_race_diffF/ns$inst$edges,
    riskgrp_over_2E = sum(ns$inst$nodefactor_risk.grp)/(2*ns$inst$edges)) }
out <- parallel::mclapply(res, function(r) { s <- r$spec
  np <- build_netparams(r$ep, smooth.main.dur = TRUE, method = "joint", duration.method = "joint_lm", seed = 20261107L, geog.pooling = "hierarchical")
  ns <- build_netstats(r$ep, np, expect.mort = em, network.size = s$n, race.prop = s$race.prop, method = "joint", seed = 20261107L, duration.cap = cap)
  rbind(current = summ(r$ns, s$n), pooled = summ(ns, s$n)) }, mc.cores = 6)
saveRDS(out, "six_pooled.rds")
for (k in names(out)) { x <- out[[k]]; cat("==", k, "\n"); print(round(rbind(x, pct = 100*(x[2,]/x[1,]-1)), 2)) }
