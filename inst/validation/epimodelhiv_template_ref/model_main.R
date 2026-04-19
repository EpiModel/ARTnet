## REFERENCE COPY (2026-04-19) of EpiModelHIV-Template/R/A-networks/model_main.R
## DO NOT EDIT — pins the ERGM specification consuming netstats$main.

## Define and fit the *main* network  model
##
## This script should not be run directly. But `sourced` by `1-estimation.R`

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
) |> unname()

# Fit model
fit_main <- EpiModel::netest(
  nw_main,
  formation = model_main,
  target.stats = netstats_main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> EpiModel::trim_netest()

# Keep only the necessary objects
rm(model_main, netstats_main)
