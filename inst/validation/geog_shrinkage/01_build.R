suppressMessages({devtools::load_all("/Users/sjennes/git/ARTnet", quiet = TRUE)})
D <- "/Users/sjennes/git/DoxyPEP-Resist/data/input/"
cities <- c("atl","det","la","mia","sea","sf")
expect_mort <- 0.000478213; dur_cap <- 0.9/(2*expect_mort); seed <- 20261107L
spec <- lapply(cities, function(c) {
  t <- read.csv(paste0(D, c, "/city.csv"), colClasses = "character"); v <- setNames(t$value, t$param)
  list(key=c, label=v[["label"]], n=as.numeric(v[["msm.pop"]]), race.prop=as.numeric(v[paste0("race.prop_",1:3)]),
       geog.cat=v[["geog.cat"]], ihp=as.numeric(v[paste0("init.hiv.prev_",1:3)]))
}); names(spec) <- cities
res <- parallel::mclapply(spec, function(s) {
  t0 <- Sys.time()
  ep <- build_epistats(geog.lvl="city", geog.cat=s$geog.cat, init.hiv.prev=s$ihp, race=TRUE, time.unit=7)
  np <- build_netparams(ep, smooth.main.dur=TRUE, method="joint", duration.method="joint_lm", seed=seed)
  ns <- build_netstats(ep, np, expect.mort=expect_mort, network.size=s$n, race.prop=s$race.prop,
                       method="joint", seed=seed, duration.cap=dur_cap)
  list(spec=s, ep=ep, np=np, ns=ns, secs=as.numeric(Sys.time()-t0, units="secs"))
}, mc.cores = 6)
saveRDS(res, "base.rds")
for (r in res) cat(r$spec$key, round(r$secs), "main", round(1000*r$ns$main$edges/r$spec$n,1),
  "casl", round(1000*r$ns$casl$edges/r$spec$n,1), "inst", round(1000*r$ns$inst$edges/r$spec$n,1), "\n")
