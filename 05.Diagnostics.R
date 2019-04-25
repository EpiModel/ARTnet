
##
## Network modeling diagnostics for ART-Net Data
## v1: 2018-08
##

## Packages ##
rm(list = ls())
suppressMessages(library("EpiModelHIV"))


## Inputs ##
city_name <- "Atlanta"


## Load Data ##
fn <- paste("data/artnet.NetEst", gsub(" ", "", city_name), "rda", sep = ".")
est <- readRDS(file = fn)

fn <- paste("data/artnet.NetStats", gsub(" ", "", city_name), "rda", sep = ".")
nstats <- readRDS(file = fn)


# Main --------------------------------------------------------------------

fit_main <- est[[1]]

model_main_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = 1:5) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = 1:3) +
  nodefactor("deg.casl", levels = 0:3) +
  degrange(from = 3) +
  concurrent +
  nodefactor("diag.status", levels = 0:1) +
  nodematch("role.class", diff = TRUE, levels = 1:3) +
  degree(0:3)
dx_main <- netdx(fit_main, nsims = 10, ncores = 6, nsteps = 500,
            nwstats.formula = model_main_dx, skip.dissolution = TRUE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx_main)
plot(dx_main)

nstats$main


# Casual ------------------------------------------------------------------

fit_casl <- est[[2]]

model_casl_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = 1:5) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = 1:3) +
  nodefactor("deg.main", levels = 0:2) +
  degrange(from = 4) +
  concurrent +
  nodefactor("diag.status", levels = 0:1) +
  nodematch("role.class", diff = TRUE, levels = 1:3) +
  degree(0:4)
dx_casl <- netdx(fit_casl, nsims = 10, ncores = 6, nsteps = 500,
            nwstats.formula = model_casl_dx, skip.dissolution = TRUE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx_casl)
plot(dx_casl)

nstats$casl


# One-Off -----------------------------------------------------------------

fit_inst <- est[[3]]

model_inst_dx <- ~edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = 1:5) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = 1:3) +
  nodefactor("risk.grp", levels = 1:5) +
  nodefactor("deg.tot", levels = 0:5) +
  nodefactor("diag.status", levels = 0:1) +
  nodematch("role.class", diff = TRUE, levels = 1:3) +
  degree(0:4)
dx_inst <- netdx(fit_inst, nsims = 10000, dynamic = FALSE,
            nwstats.formula = model_inst_dx,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))

print(dx_inst)

plot(dx_inst, sim.lines = TRUE, sim.lwd = 0.05)

nstats$inst
