suppressMessages({devtools::load_all("/Users/sjennes/git/ARTnet", quiet = TRUE)})
res <- readRDS("base.rds"); em <- 0.000478213; cap <- 0.9/(2*em)
ok <- parallel::mclapply(res, function(r) { s <- r$spec
  np <- build_netparams(r$ep, smooth.main.dur = TRUE, method = "joint", duration.method = "joint_lm", seed = 20261107L)
  ns <- build_netstats(r$ep, np, expect.mort = em, network.size = s$n, race.prop = s$race.prop, method = "joint", seed = 20261107L, duration.cap = cap)
  strip <- function(x) { for (l in c("main","casl","inst")) x[[l]][grep("^joint_", names(x[[l]]))] <- NULL; x }
  c(netstats = isTRUE(all.equal(ns, r$ns, tolerance = 0)), netparams_nonmodel = isTRUE(all.equal(strip(np), strip(r$np), tolerance = 0)),
    coefs = all(sapply(c("main","casl","inst"), function(l) all(sapply(grep("^joint_", names(np[[l]]), value=TRUE), function(m) identical(coef(np[[l]][[m]]), coef(r$np[[l]][[m]])))))))
}, mc.cores = 6)
print(do.call(rbind, ok))
