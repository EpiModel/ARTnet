## REFERENCE COPY (2026-04-19) of EpiModelHIV-Template/R/A-networks/model_ooff.R
## DO NOT EDIT — pins the ERGM specification consuming netstats$inst.

## Define and fit the *one-off* network  model
##
## This script should not be run directly. But `sourced` by `1-estimation.R`

# Formula
model_ooff <- ~ edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("risk.grp", levels = -5) +
  nodefactor("deg.tot", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_ooff <- c(
  edges                = netstats$inst$edges,
  nodematch_age.grp    = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp   = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race       = netstats$inst$nodematch_race_diffF,
  nodefactor_race      = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp  = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot   = netstats$inst$nodefactor_deg.tot[-1],
  nodematch_role.class = c(0, 0)
) |> unname()

# Fit model
fit_ooff <- EpiModel::netest(
  nw_inst,
  formation = model_ooff,
  target.stats = netstats_ooff,
  coef.diss = dissolution_coefs(~ offset(edges), 1),
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> trim_netest()

# Keep only the necessary objects
rm(model_ooff, netstats_ooff)
