## REFERENCE COPY (2026-04-19) of EpiModelHIV-Template/R/A-networks/model_casl.R
## DO NOT EDIT — pins the ERGM specification consuming netstats$casl.

## Define and fit the *casual* network  model
##
## This script should not be run directly. But `sourced` by `1-estimation.R`

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
) |> unname()

# Fit model
fit_casl <- EpiModel::netest(
  nw_casl,
  formation = model_casl,
  target.stats = netstats_casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> trim_netest()

# Keep only the necessary objects
rm(model_casl, netstats_casl)
