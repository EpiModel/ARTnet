## ---- message=FALSE------------------------------------------------------
library("ARTnet")

## ------------------------------------------------------------------------
epistats2 <- build_epistats(geog.lvl = "state", geog.cat = "GA")

## ------------------------------------------------------------------------
epistats2 <- build_epistats(geog.lvl = "division", geog.cat = "5")
epistats3 <- build_epistats(geog.lvl = "region", geog.cat = "3")

## ------------------------------------------------------------------------
epistats4 <- build_epistats(geog.lvl = NULL)

## ------------------------------------------------------------------------
#Racial stratification
epistats4 <- build_epistats(geog.lvl = "state", geog.cat = "WA", race = TRUE)

#No racial stratification
epistats5 <- build_epistats(geog.lvl = "state", geog.cat = "GA", race = FALSE)

## ------------------------------------------------------------------------
#Set age limits
epistats6 <- build_epistats(geog.lvl = "state", geog.cat = "GA", 
                            race = TRUE, age.limits = c(15, 65),
                            age.breaks = c(20, 30, 40, 50, 60))

#Specify age categories: (0-20], (20, 30], . . . (60, 100]
netparams1 <- build_netparams(epistats = epistats6, smooth.main.dur = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  #1. Epistats: Specify geographic features, as well as race stratification and total age range
#  #under consideration
#  epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta", race = TRUE,
#                             age.limits = c(30, 50), age.breaks = c(35, 45))
#  
#  #2. Netparams: Specify age categories if needed, or let ARTnet determine age categories by number of
#  #categories desired
#  netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#  
#  #3. Netstats: Finalize network setup
#  netstats <- build_netstats(epistats, netparams, expect.mort = 0.0005,
#                             network.size = 1000, edges.avg =TRUE)
#  
#  #4. Initialize network using `netstats` object from previous step.
#  num <- netstats$demog$num
#  nw <- network::network.initialize(num, directed = FALSE)
#  
#  attr.names <- names(netstats$attr)
#  attr.values <- netstats$attr
#  nw <- network::set.vertex.attribute(nw, attr.names, attr.values)
#  nw_main <- nw_casl <- nw_inst <- nw
#  
#  # 5. Main Model
#  
#  #Formula:
#  model_main <- ~edges +
#    nodematch("age.grp", diff = TRUE) +
#    nodefactor("age.grp", base = 1) +
#    nodematch("race", diff = FALSE) + #race = TRUE; omit if FALSE
#    nodefactor("race", base = 1) +
#    nodefactor("deg.casl", base = 1) +
#    concurrent +
#    degrange(from = 3) +
#    nodematch("role.class", diff = TRUE, keep = 1:2)
#  
#  # Target Stats
#  netstats_main <- c(
#    edges = netstats$main$edges,
#    nodematch_age.grp = netstats$main$nodematch_age.grp,
#    nodefactor_age.grp = netstats$main$nodefactor_age.grp[-1],
#    nodematch_race = netstats$main$nodematch_race_diffF, #If race = FALSE, value will be NULL
#    nodefactor_race = netstats$main$nodefactor_race[-1],
#    nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
#    concurrent = netstats$main$concurrent,
#    degrange = 0,
#    nodematch_role.class = c(0, 0)
#  )
#  cbind(netstats_main)
#  netstats_main <- unname(netstats_main)
#  
#  # Fit model
#  fit_main <- netest(nw_main,
#                     formation = model_main,
#                     target.stats = netstats_main,
#                     coef.diss = netstats$main$diss.byage,
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500,
#                                                     SAN.maxit = 3,
#                                                     SAN.nsteps.times = 3),
#                     verbose = FALSE)
#  
#  # 6. Casual Model
#  
#  # Formula
#  model_casl <- ~edges +
#    nodematch("age.grp", diff = TRUE) +
#    nodefactor("age.grp", base = c(1,5)) +
#    nodefactor("deg.main", base = 3) +
#    concurrent +
#    degrange(from = 4) +
#    nodematch("role.class", diff = TRUE, keep = 1:2) +
#    #If race = TRUE:
#    nodematch("race", diff = FALSE) +
#    nodefactor("race", base = 1)
#  
#  # Target Stats
#  netstats_casl <- c(
#    edges = netstats$casl$edges,
#    nodematch_age.grp = netstats$casl$nodematch_age.grp,
#    nodefactor_age.grp = netstats$casl$nodefactor_age.grp[-c(1,5)],
#    nodefactor_deg.main = netstats$casl$nodefactor_deg.main[-3],
#    concurrent = netstats$casl$concurrent,
#    degrange = 0,
#    nodematch_role.class = c(0, 0),
#    #If race = TRUE:
#    nodematch_race = netstats$casl$nodematch_race_diffF,
#    nodefactor_race = netstats$casl$nodefactor_race[-1]
#  )
#  cbind(netstats_casl)
#  netstats_casl <- unname(netstats_casl)
#  
#  # Fit model
#  fit_casl <- netest(nw_casl,
#                     formation = model_casl,
#                     target.stats = netstats_casl,
#                     coef.diss = netstats$casl$diss.byage,
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500,
#                                                     SAN.maxit = 3,
#                                                     SAN.nsteps.times = 3),
#                     verbose = FALSE)
#  
#  # 7. One-Off Model
#  
#  # Formula
#  model_inst <- ~edges +
#    nodematch("age.grp", diff = FALSE) +
#    nodefactor("age.grp", base = 1) +
#    nodefactor("risk.grp", base = 5) +
#    nodefactor("deg.tot", base = 1) +
#    nodematch("role.class", diff = TRUE, keep = 1:2) +
#    #If race = TRUE
#    nodematch("race", diff = FALSE) +
#    nodefactor("race", base = 1)
#  
#  # Target Stats
#  netstats_inst <- c(
#    edges = netstats$inst$edges,
#    nodematch_age.grp = sum(netstats$inst$nodematch_age.grp),
#    nodefactor_age.grp = netstats$inst$nodefactor_age.grp[-1],
#    nodefactor_risk.grp = netstats$inst$nodefactor_risk.grp[-5],
#    nodefactor_deg.tot = netstats$inst$nodefactor_deg.tot[-1],
#    nodematch_role.class = c(0, 0),
#    #If race = TRUE
#    nodematch_race = netstats$inst$nodematch_race_diffF,
#    nodefactor_race = netstats$inst$nodefactor_race[-1]
#  )
#  cbind(netstats_inst)
#  netstats_inst <- unname(netstats_inst)
#  
#  # Fit model
#  fit_inst <- netest(nw_inst,
#                     formation = model_inst,
#                     target.stats = netstats_inst,
#                     coef.diss = dissolution_coefs(~offset(edges), 1),
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500,
#                                                     SAN.maxit = 3,
#                                                     SAN.nsteps.times = 3),
#                     verbose = FALSE)
#  
#  # 8. Save Data
#  
#  out <- list(fit_main, fit_casl, fit_inst)
#  
#  # saveRDS(out, file = "netest.rda")

