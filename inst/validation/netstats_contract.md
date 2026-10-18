# netstats Backward-Compatibility Contract

The `netstats` object returned by `build_netstats()` is consumed by the
ERGM estimation scripts in `EpiModelHIV-Template/R/A-networks/` (verbatim
copy pinned in `epimodelhiv_template_ref/`). These are the fields the
template reads — they must remain byte-identical under
`method = "existing"` (or whatever we name the legacy flag).

Snapshot taken 2026-04-19 against EpiModelHIV-Template@main.

## `initialize.R`
- `netstats$demog$num` — network size (scalar integer)
- `netstats$attr` — named list of vertex attributes (age, sqrt.age,
  age.grp, active.sex, race, deg.casl, deg.main, deg.tot, risk.grp,
  role.class, diag.status). The attribute vectors are constructed via
  `sample()` / `apportion_lr()` / `rbinom()` and depend on RNG state.
  Validation must `set.seed()` before comparison.

## `model_main.R`
- `netstats$main$edges`
- `netstats$main$nodematch_age.grp` (vector, one per age group)
- `netstats$main$nodefactor_age.grp` (vector, one per age group)
- `netstats$main$nodematch_race_diffF` (scalar)
- `netstats$main$nodefactor_race` (vector, one per race group)
- `netstats$main$nodefactor_deg.casl` (vector, one per deg.casl level)
- `netstats$main$concurrent` (scalar)
- `netstats$main$diss.byage` — `dissolution_coefs` S3 object

## `model_casl.R`
- `netstats$casl$edges`
- `netstats$casl$nodematch_age.grp`
- `netstats$casl$nodefactor_age.grp`
- `netstats$casl$nodematch_race_diffF`
- `netstats$casl$nodefactor_race`
- `netstats$casl$nodefactor_deg.main`
- `netstats$casl$concurrent`
- `netstats$casl$diss.byage`

## `model_ooff.R`
- `netstats$inst$edges`
- `netstats$inst$nodematch_age.grp`
- `netstats$inst$nodefactor_age.grp`
- `netstats$inst$nodematch_race_diffF`
- `netstats$inst$nodefactor_race`
- `netstats$inst$nodefactor_risk.grp`
- `netstats$inst$nodefactor_deg.tot`

## Not directly consumed but still part of the contract
Anything else currently in `netstats$*` — `nodematch_race`,
`absdiff_age`, `absdiff_sqrt.age`, etc. — is also part of the contract
by default because the package ships it publicly. The validation script
does a full-object diff rather than checking only the fields above.

## `netparams` contract
The validation also captures `netparams` whole (the input to
`build_netstats`). Joint models are *additive* outputs (new
`$joint_model` fields), so:
- existing `netparams$main$*`, `$casl$*`, `$inst$*`, `$all$*` fields must be
  byte-identical under `method = "existing"`;
- new fields (`$joint_model`) are ignored during comparison.
