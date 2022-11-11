
## testing ARTnet workflow on full-scale MSM model with max age 65

## should work on EpiModelHIV-p@main

# Setup  -----------------------------------------------------------------------
library("EpiModelHIV")
library("ARTnet")

networks_size <- 10000

# 0. Initialize Network --------------------------------------------------------
epistats <- build_epistats(
  geog.lvl = "city",
  geog.cat = "Atlanta",
  init.hiv.prev = c(0.33, 0.137, 0.084),
  time.unit = 7,
  age.limits = c(15, 100),
  age.sexual.cessation = 65
)

netparams <- build_netparams(
  epistats = epistats,
  smooth.main.dur = TRUE,
  cessation.dissolve.edges = TRUE
)

# Weighted uniform age distribution
age1 <- rep(1 / length(15:64) * 0.999, length(15:64))
age2 <- rep(1 / length(65:99) * 0.001, length(65:99))
age.pyramid <- c(age1, age2)

netstats <- build_netstats(
  epistats,
  netparams,
  expect.mort = 0.000478213,
  age.pyramid = age.pyramid,
  network.size = networks_size,
)

table(floor(netstats$attr$age))
sum(netstats$attr$age >= 65)
hist(netstats$attr$age)

netstats$main$nodefactor_age.grp
sum(netstats$main$nodefactor_age.grp) / 2
netstats$main$edges

num <- netstats$demog$num
nw <- EpiModel::network_initialize(num)

attr_names <- names(netstats$attr)
attr_values <- netstats$attr

nw_main <- EpiModel::set_vertex_attribute(nw, attr_names, attr_values)
nw_casl <- nw_main
nw_inst <- nw_main

# 1. Main Model ----------------------------------------------------------------

# Formula
model_main <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  degrange(from = 3) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_main <- c(
  edges                = netstats$main$edges,
  nodematch_age.grp    = netstats$main$nodematch_age.grp,
  nodefactor_age.grp   = netstats$main$nodefactor_age.grp[-1],
  nodematch_race       = netstats$main$nodematch_race_diffF,
  nodefactor_race      = netstats$main$nodefactor_race[-1],
  nodefactor_deg.casl  = netstats$main$nodefactor_deg.casl[-1],
  concurrent           = netstats$main$concurrent,
  degrange             = 0,
  nodematch_role.class = c(0, 0)
)
netstats_main <- unname(netstats_main)

# Fit model
fit_main <- netest(
  nw_main,
  formation = model_main,
  target.stats = netstats_main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 1e4,
    MCMC.interval = 5e3,
    parallel = 1
  ),
  verbose = FALSE
)
fit_main <- trim_netest(fit_main)

# 2. Casual Model ---------------------------------------------------------

# Formula
model_casl <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -5) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.main", levels = -3) +
  concurrent +
  degrange(from = 4) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_casl <- c(
  edges                = netstats$casl$edges,
  nodematch_age.grp    = netstats$casl$nodematch_age.grp,
  nodefactor_age.grp   = netstats$casl$nodefactor_age.grp[-5],
  nodematch_race       = netstats$casl$nodematch_race_diffF,
  nodefactor_race      = netstats$casl$nodefactor_race[-1],
  nodefactor_deg.main  = netstats$casl$nodefactor_deg.main[-3],
  concurrent           = netstats$casl$concurrent,
  degrange             = 0,
  nodematch_role.class = c(0, 0)
)
netstats_casl <- unname(netstats_casl)

# Fit model
fit_casl <- netest(
  nw_casl,
  formation = model_casl,
  target.stats = netstats_casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 1e4,
    MCMC.interval = 5e3,
    parallel = 1
  ),
  verbose = FALSE
)
fit_casl <- trim_netest(fit_casl)

# 3. One-Off Model -------------------------------------------------------------

# Formula
model_inst <- ~ edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("risk.grp", levels = -5) +
  nodefactor("deg.tot", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_inst <- c(
  edges                = netstats$inst$edges,
  nodematch_age.grp    = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp   = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race       = netstats$inst$nodematch_race_diffF,
  nodefactor_race      = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp  = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot   = netstats$inst$nodefactor_deg.tot[-1],
  nodematch_role.class = c(0, 0)
)
netstats_inst <- unname(netstats_inst)

# Fit model
fit_inst <- netest(
  nw_inst,
  formation = model_inst,
  target.stats = netstats_inst,
  coef.diss = dissolution_coefs(~ offset(edges), 1),
  set.control.ergm = control.ergm(
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 1e4,
    MCMC.interval = 5e3,
    parallel = 1
  ),
  verbose = FALSE
)
fit_inst <- trim_netest(fit_inst)


# 4. Run Diagnostics --------------------------------------------------------------------------

nsims <- 10
ncores <- 10
nsteps <- 1000

model_main_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("deg.casl", levels = TRUE) +
  degrange(from = 3) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:3)

dx_main <- netdx(
  fit_main,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_main_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_main_static <- netdx(
  fit_main,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_main_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

print(dx_main)
plot(dx_main, targ.col = "grey50")
plot(dx_main, type = "duration", targ.col = "grey50")
plot(dx_main, type = "dissolution", targ.col = "grey50")

print(dx_main_static)
plot(dx_main_static, sim.lines = TRUE, sim.lwd = 0.1, mean.lwd = 0.5)

# Casual
model_casl_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("deg.main", levels = TRUE) +
  degrange(from = 4) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)

dx_casl <- netdx(
  fit_casl,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_casl_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_casl_static <- netdx(
  fit_casl,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_casl_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

print(dx_casl)
plot(dx_casl, targ.col = "grey50")
plot(dx_casl, type = "duration", targ.col = "grey50")
plot(dx_casl, type = "dissolution", targ.col = "grey50")

print(dx_casl_static)
plot(dx_casl_static, sim.lines = TRUE, sim.lwd = 0.1, mean.lwd = 0.5)


# One Off
model_inst_dx <- ~edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = TRUE) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = TRUE) +
  nodefactor("risk.grp", levels = TRUE) +
  nodefactor("deg.tot", levels = TRUE) +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)

dx_inst <- netdx(
  fit_inst,
  nsims = 50000,
  dynamic = FALSE,
  nwstats.formula = model_inst_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

print(dx_inst)
plot(dx_inst, sim.lines = TRUE, sim.lwd = 0.1, mean.lwd = 0.5)

nw_main
get_vertex_attribute(nw_main, "active.sex")
