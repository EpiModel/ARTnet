## ARTnet v2.9.2

Two fixes to the sexual-cessation configuration of `build_netstats()`, where
`age.limits[2]` is set above `age.sexual.cessation` so the population keeps
aging past the point where it stops forming partnerships.

* `active.sex` is now zeroed on the age group the breaks declare rather than on
  the largest age group realized in the sampled population. When the
  post-cessation group came out empty, the old behavior marked the oldest
  *sexually active* group inactive and zeroed its degree on all three layers,
  with no error and no warning. This affected both the sampling path and the
  `target_pop` data frame path, where a projected population whose oldest
  member sits below the cessation age would have hit the same problem.
* `young.prop` must now be strictly between 0 and 1 when a cessation band
  exists, and is rejected with a message pointing at a usable value. Either
  endpoint empties one side of the age split by construction, and an age group
  with no members carries no `ergm` term, so the failure used to surface much
  later as a `target.stats` length mismatch in `netest()`. A post-cessation
  group that comes out empty for other reasons (a `network.size` too small for
  the seeded share, or a `target_pop` with no members above the cessation age)
  now warns for the same reason.

Neither change affects the default configuration, where the population upper
limit and the cessation age coincide and there is no post-cessation band.

## ARTnet v2.9.0

* `ARTnetData` is now a suggested package and example `EpiStats` and `Netstats`
are provided. This allows experimenting with `EpiModelHIV` without being part
of the EpiModel organisation.

## ARTnet v2.8.0

* Modify the order of the `race.prop` argument to `build_netstats`. The new
order is: "Black, Hispanic, and White/Other" to be consistent with the rest of
ARTnet and EpiModelHIV. (Previously, "White/Other" was first).

## ARTnet v2.7.0

### Age Range Expansion
* ARTnet now allows for specifying an age range that spans outside the range of the age eligibility criteria of ARTnet itself. For example, it may be desired to have a population ranging from age 15 to 100 represented in the epidemic model but only to model sexual activity within 15 to 65 year olds (corresponding to ARTnet eligibility). This approach is motivated by cost-effectiveness analysis models that require tracking demographics through natural mortality rather than the age of sexual cessation. Examples are provided in the documentation of `build_epistats`, `build_netparams`, and `build_netstats`.
* There is related flexibility in passing in the age distribution structure in `build_netstats`. This allows passing in different starting age distributions necessary for the TERGM model fitting.

### Updated Mortality Rate Parameterization
* Default age-specific mortality rates from 2020 census are now included. Also provided is a general updating function, `update_asmr` to allow passing in other mortality rate data. See the help page for details and examples.


## ARTnet v2.6.0

* Primary update is handling of geographical stratification to allow for multiple jurisdictions to be used (see #40), contributed by @sgoodreau.
* Major cosmetic clean-up of scripts.


## ARTnet v2.5.6

* Adds a new `time.unit` parameter to `build_epistats` that allows for rescaling of time units from 1 (days) to 30 (months), instead of everything hard coded as weeks.
* Fixes error in casual duration estimates.


## ARTnet v2.5.0

* Further revision to ARTnet workflow to include age stratification as well as compatibility with revisions to EpiModelHIV.


## ARTnet v2.0.0

* Revised workflow and parameterization based on expanded geographic/racial parameterizations.


## ARTnet v1.0.0

* Initial package publication.
