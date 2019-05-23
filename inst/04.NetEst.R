
##
## Network modeling for ARTnet Data
##

## Packages ##
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("ARTnet"))

epistats <- build_epistats(city_name = "Atlanta")
netparams <- build_netparams(epistats = epistats)
netstats <- build_netstats(epistats, netparams)


# 0. Initialize Network ---------------------------------------------------

num <- netstats$num
nw <- network::network.initialize(num, directed = FALSE)

attr.names <- names(netstats$attr)
attr.values <- netstats$attr
nw <- network::set.vertex.attribute(nw, attr.names, attr.values)
nw_main <- nw_casl <- nw_inst <- nw


# 1. Main Model -----------------------------------------------------------

# Formula
model_main <- ~edges +
               nodematch("age.grp", diff = TRUE) +
               nodefactor("age.grp", base = 1) +
               nodematch("race", diff = FALSE) +
               nodefactor("race", base = 1) +
               nodefactor("deg.casl", base = 1) +
               concurrent +
               nodefactor("diag.status", base = 1) +
               degrange(from = 3) +
               nodematch("role.class", diff = TRUE, keep = 1:2)

# Target Stats
netstats_main <- c(
  edges = netstats$main$edges,
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  nodefactor_age.grp = netstats$main$nodefactor_age.grp[-1],
  nodematch_race = netstats$main$nodematch_race_diffF,
  nodefactor_race = netstats$main$nodefactor_race[-1],
  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  concurrent = netstats$main$concurrent,
  nodefactor_diag.status = netstats$main$nodefactor_diag.status[-1],
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
cbind(netstats_main)
netstats_main <- unname(netstats_main)

# Fit model
fit_main <- netest(nw_main,
                   formation = model_main,
                   target.stats = netstats_main,
                   coef.diss = netstats$main$diss,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 5,
                                                   SAN.nsteps.times = 5),
                   verbose = FALSE)



# 2. Casual Model ---------------------------------------------------------

# Formula
model_casl <- ~edges +
               nodematch("age.grp", diff = TRUE) +
               nodefactor("age.grp", base = c(1,5)) +
               nodematch("race", diff = FALSE) +
               nodefactor("race", base = 1) +
               nodefactor("deg.main", base = 3) +
               concurrent +
               nodefactor("diag.status", base = 1) +
               degrange(from = 4) +
               nodematch("role.class", diff = TRUE, keep = 1:2)

# Target Stats
netstats_casl <- c(
  edges = netstats$casl$edges,
  nodematch_age.grp = netstats$casl$nodematch_age.grp,
  nodefactor_age.grp = netstats$casl$nodefactor_age.grp[-c(1,5)],
  nodematch_race = netstats$casl$nodematch_race_diffF,
  nodefactor_race = netstats$casl$nodefactor_race[-1],
  nodefactor_deg.main = netstats$casl$nodefactor_deg.main[-3],
  concurrent = netstats$casl$concurrent,
  nodefactor_diag.status = netstats$casl$nodefactor_diag.status[-1],
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
cbind(netstats_casl)
netstats_casl <- unname(netstats_casl)

# Fit model
fit_casl <- netest(nw_casl,
                   formation = model_casl,
                   target.stats = netstats_casl,
                   coef.diss = netstats$casl$diss,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 5,
                                                   SAN.nsteps.times = 5),
                   verbose = FALSE)


# 3. One-Off Model --------------------------------------------------------

# Formula
model_inst <- ~edges +
               nodematch("age.grp", diff = FALSE) +
               nodefactor("age.grp", base = 1) +
               nodematch("race", diff = FALSE) +
               nodefactor("race", base = 1) +
               nodefactor("risk.grp", base = 5) +
               nodefactor("deg.tot", base = 1) +
               nodefactor("diag.status", base = 1) +
               nodematch("role.class", diff = TRUE, keep = 1:2)

# Target Stats
netstats_inst <- c(
  edges = netstats$inst$edges,
  nodematch_age.grp = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race = netstats$inst$nodematch_race_diffF,
  nodefactor_race = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot = netstats$inst$nodefactor_deg.tot[-1],
  nodefactor_diag.status = netstats$inst$nodefactor_diag.status[-1],
  nodematch_role.class = c(0, 0)
)
cbind(netstats_inst)
netstats_inst <- unname(netstats_inst)

# Fit model
fit_inst <- netest(nw_inst,
                   formation = model_inst,
                   target.stats = netstats_inst,
                   coef.diss = dissolution_coefs(~offset(edges), 1),
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 5,
                                                   SAN.nsteps.times = 5),
                   verbose = FALSE)


# 4. Save Data ------------------------------------------------------------

out <- list(fit_main, fit_casl, fit_inst)

fns <- strsplit(fn, "[.]")[[1]]
fn.new <- paste(fns[1], "NetEst", fns[3], "rda", sep = ".")

saveRDS(out, file = fn.new)
