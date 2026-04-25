# Parse the target_pop argument to build_netstats. Returns a list with
# $form in {"default", "list", "data.frame"} plus form-specific payload.
# Raises an informative error for the not-yet-implemented character form.
.parse_target_pop <- function(target_pop, race) {
  if (is.null(target_pop)) return(list(form = "default"))
  if (is.character(target_pop)) {
    stop("Built-in reference target populations (target_pop = '",
         target_pop, "') are not yet implemented in ARTnet. The planned ",
         "set is geography-specific general male population demographics ",
         "(NCHS age pyramid + ARTnetData::race.dist by geography); when ",
         "those bundles ship, this argument will accept their names. ",
         "For now, pass a list of marginal distributions or a data.frame. ",
         "See issue #64.", call. = FALSE)
  }
  if (is.data.frame(target_pop)) {
    req <- c("age", "deg.casl", "deg.main", "role.class", "risk.grp")
    if (isTRUE(race)) req <- c(req, "race")
    missing_cols <- setdiff(req, names(target_pop))
    if (length(missing_cols)) {
      stop("target_pop data.frame missing required columns: ",
           paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
    if (nrow(target_pop) < 1) {
      stop("target_pop data.frame must have at least one row.",
           call. = FALSE)
    }
    return(list(form = "data.frame", df = target_pop))
  }
  if (is.list(target_pop)) {
    known <- c("age.pyramid", "race.prop", "race.props",
               "deg.casl", "deg.main", "deg.tot",
               "role.class", "risk.grp")
    extra <- setdiff(names(target_pop), known)
    if (length(extra)) {
      stop("target_pop list has unknown elements: ",
           paste(extra, collapse = ", "),
           ". Allowed names: ", paste(known, collapse = ", "), ".",
           call. = FALSE)
    }
    # Normalize the alternate name spelling shown in #64's example.
    if (!is.null(target_pop$race.props) && is.null(target_pop$race.prop)) {
      target_pop$race.prop <- target_pop$race.props
    }
    return(list(form = "list", overrides = target_pop))
  }
  stop("target_pop must be NULL, a list, a data.frame, or a character ",
       "string. Got class: ", paste(class(target_pop), collapse = ", "), ".",
       call. = FALSE)
}


# Aggregate synth-population stratum-level durations under joint_lm. Per-ego
# predicted log(duration) from the joint_lm fit, marginalized over partner-race
# uncertainty using joint_nm_race_model when race = TRUE, then exponentiated
# and median-aggregated within each (same.age.grp x index.age.grp) stratum. The
# returned vector matches the row layout of `netparams$<layer>$durs.<layer>.byage`:
# nonmatch first (same.age.grp = 0), then matched-within-age-group 1..N, plus
# an optional deterministic "post-cessation" row when sex.cess.mod = TRUE.
#
# Used in `build_netstats(method = "joint")` to override the within-ARTnet
# stratum medians that build_netparams emits under duration.method = "joint_lm",
# so dissolution offsets reflect the synthetic target population's joint
# attribute distribution rather than ARTnet's. See issue #73.
#
# Returns NULL if joint_dur_model is NULL (caller falls back to netparams values).
.aggregate_synth_byage_durations <- function(joint_dur_model,
                                             joint_nm_race_model,
                                             synth_data,
                                             n_age_grps,
                                             sex_cess_extra_row = FALSE) {
  if (is.null(joint_dur_model)) return(NULL)

  if (!is.null(joint_nm_race_model)) {
    p_same_race <- predict(joint_nm_race_model, newdata = synth_data,
                           type = "response")
  } else {
    p_same_race <- rep(0, nrow(synth_data))
  }

  predict_stratum_median <- function(same_age, age_grp_select) {
    sel <- if (is.na(age_grp_select)) {
      rep(TRUE, nrow(synth_data))
    } else {
      synth_data$age.grp == age_grp_select
    }
    if (!any(sel)) return(NA_real_)
    sub <- synth_data[sel, , drop = FALSE]
    sub$same.age.grp <- as.integer(same_age)

    sub$same.race <- 0L
    pred_log_0 <- predict(joint_dur_model, newdata = sub)
    sub$same.race <- 1L
    pred_log_1 <- predict(joint_dur_model, newdata = sub)

    p <- p_same_race[sel]
    dur_marg <- p * exp(pred_log_1) + (1 - p) * exp(pred_log_0)
    median(dur_marg, na.rm = TRUE)
  }

  medians <- c(predict_stratum_median(0, NA),
               vapply(seq_len(n_age_grps),
                      function(k) predict_stratum_median(1, k),
                      numeric(1)))

  mean_dur_adj <- ifelse(is.na(medians) | medians <= 0,
                         NA_real_,
                         1 / (1 - 2^(-1 / medians)))

  if (isTRUE(sex_cess_extra_row)) mean_dur_adj <- c(mean_dur_adj, 1)
  mean_dur_adj
}


#' Calculate Network Target Statistics
#'
#' @description Calculates the final target statistics for the network models by applying
#'              individual-level network statistics against the population size and structure, for
#'              use in the EpiModelHIV workflow.
#'
#' @param epistats Output from [`build_epistats`].
#' @param netparams Output from [`build_netparams`].
#' @param network.size Size of the starting network.
#' @param expect.mort Expected average mortality level to pass into [`dissolution_coefs`] function.
#' @param age.pyramid Numerical vector of length equal to the length of the age range specified in
#'        [`build_epistats`], containing probability distribution of each year of age, summing to
#'        one. If `NULL`, then a default distribution based on 2020 NCHS data is used.
#' @param edges.avg If `TRUE`, calculates the overall edges target statistics as a weighted average
#'        of the statistics for edges by race/ethnicity group; if `FALSE`, takes the raw average.
#' @param race.prop A numerical vector of length 3, containing the proportion of the population with
#'        each of the three values for the nodal attribute "race" in order: Black,
#'        Hispanic, and White/Other).
#' @param young.prop The proportion of the population that should be below the age of sexual cessation.
#'        Default is NULL (meaning no re-weighting of the `age.pyramid` parameter is performed).
#'        This parameter is only used if the age of sexual cessation is less than the upper age bound.
#' @param target_pop Optional specification of the synthetic target population. Defaults to NULL,
#'        which uses the legacy patchwork of reference sources (NCHS age pyramid +
#'        `ARTnetData::race.dist` + ARTnet's own degree / role / risk-quintile distributions).
#'        Supports three forms:
#'        \itemize{
#'          \item A **named list** of marginal distribution overrides. Allowed names:
#'            `age.pyramid` (vector of length `nAges`), `race.prop` (length matching
#'            `race.level`), `deg.casl` (length 4), `deg.main` (length 3), `deg.tot` (length 4),
#'            `role.class` (length 3), `risk.grp` (length matching `nf.risk.grp`). Each
#'            element overrides the corresponding default source; absent names fall through
#'            to legacy defaults. Equivalent in form to passing the older `age.pyramid` /
#'            `race.prop` arguments, but extends the override surface to the per-attribute
#'            distributions previously sourced from `netparams`.
#'          \item A **data.frame** with one row per node. Required columns: `age`, `deg.casl`,
#'            `deg.main`, `role.class`, `risk.grp` (plus `race` when `epistats$race = TRUE`).
#'            Optional columns: `sqrt.age`, `age.grp`, `active.sex`, `deg.tot`, `diag.status`
#'            (derived from required columns if absent). When supplied, attribute sampling is
#'            bypassed entirely and `network.size` is set to `nrow(target_pop)`. This form is
#'            for users with a fully-specified joint synthetic population (e.g., post-stratified
#'            to NHBS or AMIS demographics).
#'          \item A **character string** naming a built-in reference population. Currently
#'            raises an informative error. The planned set is geography-specific general male
#'            population demographics (NCHS age pyramid + `ARTnetData::race.dist` by
#'            geography) — bundles like `"atlanta"` or `"us_msm_male"` packaged from data
#'            already in the package, no NHBS or other restricted data required.
#'        }
#' @param method Character. Either `"existing"` (default) or `"joint"`. `"existing"` reproduces
#'        the pre-refactor behavior byte-for-byte: target statistics for edges, nodefactor, and
#'        concurrent are computed layer-by-layer from the univariate marginal fits stored on
#'        `netparams`. `"joint"` uses the joint Poisson GLM fit at `netparams$<layer>$joint_model`
#'        (set by `build_netparams(..., method = "joint")`) to predict expected degree for each
#'        synthetic-population node and aggregates via g-computation: `edges = sum(pred) / 2`
#'        and `nodefactor_<attr>[level] = sum(pred[attr == level])`. `concurrent` uses a parallel
#'        joint binomial GLM on the `deg > 1` indicator. `nodematch_*` (age.grp and race) and
#'        `absdiff_*` (age and sqrt.age) are aggregated from per-ego dyad-level predictions:
#'        `sum(pred_deg * pred_dyad) / 2`, where `pred_dyad` is a per-ego expected value from a
#'        partnership-level GLM fit on long-format ARTnet data with ego attributes on the RHS
#'        (`netparams$<layer>$joint_nm_{age,race}_model` and
#'        `netparams$<layer>$joint_absdiff_{age,sqrtage}_model`). Under `"joint"`, edges and all
#'        nodefactor target stats are internally consistent by construction
#'        (`sum(nodefactor_<attr>) = 2 * edges`), so `edges.avg` has no effect. The
#'        `diss.byage` dissolution offset is computed from synth-aggregated stratum-level
#'        durations when `build_netparams(..., duration.method = "joint_lm")` was used:
#'        per-ego predicted log(duration) from `joint_duration_model`, marginalized over
#'        partner-race uncertainty via `joint_nm_race_model`, then median-aggregated within
#'        stratum. `nodefactor_risk.grp` and `diss.homog` still use the within-ARTnet
#'        univariate / aggregated values — those are not consumed by the standard
#'        EpiModelHIV-Template dissolution offset.
#' @param browser If `TRUE`, run `build_netparams` in interactive browser mode.
#'
#' @details
#' This function takes output from [`build_epistats`] and [`build_netparams`] to build the relevant
#' population-level network statistics that will be used in network estimation using the [EpiModel]
#' package.
#'
#' The parameter `edge.avg` allows a user specify how the network edges statistics as estimated from
#' ARTnet with [`build_netparams`] should be calculated. With `edges.avg = FALSE`, the edges count
#' is the overall mean degree (divided by 2) times the `network.size`. However, this may create model
#' fitting issues if sample proportions by race/ethnicity in ARTnet do not match the population
#' proportions (i.e., the overall average may be inconsistent with the weighted average due to
#' different race denominator sizes). The `edges.avg = TRUE` setting calculates the overall edges as
#' the sum of mean degrees in each race/ethnicity group (divided by 2) times the size of each group.
#'
#' The `race.prop` argument only needs to be specified if the built-in race/ethnicity distribution
#' data in `ARTnetData::race.dist` is not used. This may be the case `geog.lvl = "county"` or if
#' `geog.cat` has length >1; otherwise, the values can be obtained automatically. If `race.prop` is
#' not supplied in either of these cases, national US values will be used..
#'
#' @export
#'
#' @examples
#' # Standard model with default age stratification
#' epistats <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta")
#' netparams <- build_netparams(epistats, smooth.main.dur = TRUE)
#' netstats <- build_netstats(epistats, netparams)
#'
#' # Restricted age stratification
#' epistats2 <- build_epistats(geog.lvl = "state", geog.cat = "GA",
#'                             age.limits = c(20, 50),
#'                             age.breaks = c(20, 30, 40))
#' netparams2 <- build_netparams(epistats2, smooth.main.dur = TRUE)
#' netstats2 <- build_netstats(epistats2, netparams2)
#'
#' # Model with sexual cessation age < age limit, without age pyramid reweighting
#' epistats3 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams3 <- build_netparams(epistats3, smooth.main.dur = TRUE)
#' netstats3 <- build_netstats(epistats3, netparams3)
#'
#' # Model with sexual cessation age < age limit, with age pyramid reweighting
#' epistats4 <- build_epistats(geog.lvl = "city",
#'                             geog.cat = "Atlanta",
#'                             race = TRUE,
#'                             age.limits = c(15, 100),
#'                             age.breaks = c(25, 35, 45, 55),
#'                             age.sexual.cessation = 65)
#' netparams4 <- build_netparams(epistats4, smooth.main.dur = TRUE)
#' netstats4 <- build_netstats(epistats4, netparams4, young.prop = 0.995)
#'
build_netstats <- function(epistats, netparams,
                           network.size = 10000,
                           expect.mort = 0.0001,
                           age.pyramid = NULL,
                           edges.avg = FALSE,
                           race.prop = NULL,
                           young.prop = NULL,
                           method = c("existing", "joint"),
                           target_pop = NULL,
                           browser = FALSE) {
  method <- match.arg(method)

  # Ensures that ARTnetData is installed
  if (system.file(package = "ARTnetData") == "") stop(missing_data_msg)

  # target_pop API (#64): NULL preserves the legacy patchwork-default
  # behavior; a list overrides specific marginals; a data.frame supplies
  # the synthetic population directly. Character-form built-in references
  # are not yet implemented.
  .tp <- .parse_target_pop(target_pop, race = epistats$race)
  if (.tp$form == "list") {
    if (!is.null(.tp$overrides$age.pyramid)) age.pyramid <- .tp$overrides$age.pyramid
    if (!is.null(.tp$overrides$race.prop))   race.prop   <- .tp$overrides$race.prop
  }
  if (.tp$form == "data.frame") {
    network.size <- nrow(.tp$df)
  }

  if (method == "joint") {
    missing_joint <- vapply(c("main", "casl", "inst"),
      function(layer) is.null(netparams[[layer]]$joint_model),
      logical(1))
    if (any(missing_joint)) {
      stop("method = 'joint' requires netparams fit with method = 'joint'. ",
           "No joint_model on layer(s): ",
           paste(names(missing_joint)[missing_joint], collapse = ", "), ". ",
           "Re-run build_netparams(..., method = 'joint').")
    }
  }

  if (browser == TRUE) {
    browser()
  }

  ## Data ##
  race.dist <- ARTnetData::race.dist

  ## Inputs ##
  geog.cat <- epistats$geog.cat
  geog.lvl <- epistats$geog.lvl
  sex.cess.mod <- epistats$sex.cess.mod
  race <- epistats$race
  race.level <- epistats$race.level
  age.limits <- epistats$age.limits
  age.sexual.cessation <- epistats$age.sexual.cessation

  time.unit <- epistats$time.unit


  # Demographic Initialization ----------------------------------------------

  out <- list()
  out$demog <- list()
  out$geog.lvl <- geog.lvl
  out$race <- race
  out$time.unit <- time.unit

  # Overall network size
  num <- out$demog$num <- network.size

  # Capitalize each race and join multiple with a `.`
  # c("white", "other") => "White.Other"
  flattened_race_level <- vapply(
    race.level,
    \(x) gsub("\\b(\\w)", "\\U\\1", paste0(x, collapse = "."), perl = TRUE),
    character(1)
  )

  # Population size by race group
  if (!is.null(race.prop)) {
    props <- as.data.frame(t(race.prop))
    colnames(props) <- flattened_race_level
  } else {
    if (!is.null(geog.lvl) && geog.lvl != "county" && length(geog.cat) == 1) {
      props <- race.dist[[geog.lvl]][which(race.dist[[geog.lvl]]$Geog == geog.cat), -c(1, 2)] / 100
    } else {
      props <- race.dist[["national"]][, -c(1, 2)] / 100
    }
  }
  out$demog$props <- props
  race.num.vars <- list()
  total_remaining <- num
  for (i in seq_len(length(flattened_race_level) - 1)) {
    race_name <- flattened_race_level[i]
    race_num_var <- paste0("num.", race_name)
    race_num_value <- round(num * props[[race_name]])
    out$demog[[race_num_var]] <- race_num_value
    race.num.vars[[race_name]] <- race_num_value
    total_remaining <- total_remaining - race_num_value
  }

  # Assign the residual race group
  residual_race <- flattened_race_level[length(flattened_race_level)]
  residual_num_var <- paste0("num.", residual_race)
  out$demog[[residual_num_var]] <- total_remaining
  race.num.vars[[residual_race]] <- total_remaining

  for (race_name in names(race.num.vars)) {
    race_num_var <- paste0("num.", substr(race_name, 1, 1))
    assign(race_num_var, race.num.vars[[race_name]])
  }

  ## Age-sex-specific mortality rates (B, H, W)
  #  in 1-year age decrements starting with age 1
  #  from CDC NCHS Underlying Cause of Death database (for 2020)
  asmr.B <- c(0.00079, 0.00046, 0.00030, 0.00025, 0.00024, 0.00025, 0.00019,
              0.00019, 0.00021, 0.00020, 0.00026, 0.00026, 0.00038, 0.00056,
              0.00077, 0.00100, 0.00151, 0.00227, 0.00271, 0.00264, 0.00297,
              0.00302, 0.00315, 0.00319, 0.00322, 0.00319, 0.00336, 0.00337,
              0.00330, 0.00363, 0.00396, 0.00392, 0.00407, 0.00428, 0.00411,
              0.00453, 0.00485, 0.00486, 0.00533, 0.00513, 0.00575, 0.00580,
              0.00628, 0.00671, 0.00669, 0.00750, 0.00773, 0.00858, 0.00934,
              0.00947, 0.00999, 0.01141, 0.01216, 0.01360, 0.01432, 0.01517,
              0.01699, 0.01853, 0.02021, 0.02099, 0.02366, 0.02547, 0.02877,
              0.02979, 0.03104, 0.03467, 0.03653, 0.03941, 0.04114, 0.04320,
              0.04487, 0.04879, 0.05100, 0.05678, 0.05611, 0.06384, 0.06891,
              0.07399, 0.07682, 0.08209, 0.08938, 0.09737, 0.10400, 0.11336,
              0.16336)
  asmr.H <- c(0.00032, 0.00021, 0.00018, 0.00011, 0.00011, 0.00009, 0.00010,
              0.00009, 0.00009, 0.00012, 0.00013, 0.00015, 0.00016, 0.00025,
              0.00036, 0.00058, 0.00076, 0.00106, 0.00125, 0.00134, 0.00145,
              0.00156, 0.00164, 0.00166, 0.00164, 0.00159, 0.00176, 0.00172,
              0.00201, 0.00198, 0.00192, 0.00191, 0.00202, 0.00204, 0.00219,
              0.00223, 0.00251, 0.00246, 0.00272, 0.00272, 0.00298, 0.00307,
              0.00321, 0.00351, 0.00367, 0.00391, 0.00442, 0.00484, 0.00512,
              0.00521, 0.00616, 0.00649, 0.00714, 0.00790, 0.00863, 0.00938,
              0.00992, 0.01094, 0.01222, 0.01217, 0.01464, 0.01483, 0.01630,
              0.01731, 0.01850, 0.02054, 0.02269, 0.02321, 0.02515, 0.02734,
              0.02937, 0.03064, 0.03349, 0.03670, 0.03980, 0.04387, 0.04724,
              0.05151, 0.05591, 0.05902, 0.06345, 0.07317, 0.07849, 0.08617,
              0.13436)
  asmr.W <- c(0.00034, 0.00023, 0.00019, 0.00014, 0.00014, 0.00010, 0.00010,
              0.00009, 0.00009, 0.00012, 0.00014, 0.00015, 0.00022, 0.00028,
              0.00036, 0.00050, 0.00059, 0.00082, 0.00096, 0.00104, 0.00126,
              0.00128, 0.00134, 0.00144, 0.00153, 0.00163, 0.00172, 0.00186,
              0.00194, 0.00205, 0.00220, 0.00225, 0.00238, 0.00245, 0.00247,
              0.00264, 0.00274, 0.00280, 0.00306, 0.00312, 0.00324, 0.00329,
              0.00344, 0.00354, 0.00371, 0.00405, 0.00442, 0.00479, 0.00511,
              0.00547, 0.00599, 0.00653, 0.00706, 0.00768, 0.00827, 0.00922,
              0.00978, 0.01065, 0.01151, 0.01235, 0.01349, 0.01437, 0.01548,
              0.01664, 0.01730, 0.01879, 0.01986, 0.02140, 0.02263, 0.02419,
              0.02646, 0.02895, 0.03031, 0.03625, 0.03753, 0.04268, 0.04631,
              0.05235, 0.05724, 0.06251, 0.06934, 0.07589, 0.08669, 0.09582,
              0.16601)

  if (race == TRUE) {
    # transformed to rates by time unit
    trans.asmr.B <- 1 - (1 - asmr.B)^(1 / (364 / time.unit))
    trans.asmr.H <- 1 - (1 - asmr.H)^(1 / (364 / time.unit))
    trans.asmr.W <- 1 - (1 - asmr.W)^(1 / (364 / time.unit))

    # Transformed rates, 85+ rate for ages 85 - 100
    vec.asmr.B <- c(trans.asmr.B, rep(tail(trans.asmr.B, n = 1), 15))
    vec.asmr.H <- c(trans.asmr.H, rep(tail(trans.asmr.H, n = 1), 15))
    vec.asmr.W <- c(trans.asmr.W, rep(tail(trans.asmr.W, n = 1), 15))

    asmr <- data.frame(age = 1:100, vec.asmr.B, vec.asmr.H, vec.asmr.W)
  } else {
    asmr.O <- props$Black * asmr.B + props$Hispanic * asmr.H +
      props$White.Other * asmr.W

    trans.asmr <- 1 - (1 - asmr.O)^(1 / (364 / time.unit))

    vec.asmr <- c(trans.asmr, rep(tail(trans.asmr, n = 1), 15))
    asmr <- data.frame(age = 1:100, vec.asmr, vec.asmr, vec.asmr)
  }

  # Setting deterministic mortality prob = 1 at upper age limit
  max.age <- age.limits[2]
  asmr[asmr$age >= max.age, -1] <- 1
  out$demog$asmr <- asmr


  # Nodal Attribute Initialization ------------------------------------------

  out$attr <- list()

  # age attributes
  # Currently uniform
  nAges <- age.limits[2] - age.limits[1]
  age.vals <- age.limits[1]:(age.limits[2] - 1)
  out$demog$ages <- age.vals
  age.breaks <- out$demog$age.breaks <- epistats$age.breaks
  nquants <- length(netparams$inst$nf.risk.grp)

  # List-form distribution overrides for the sampling path. NULL means
  # "use the existing default source" (current behavior).
  .ov <- if (.tp$form == "list") .tp$overrides else list()
  .dist_deg.casl   <- if (!is.null(.ov$deg.casl))   .ov$deg.casl   else netparams$main$deg.casl.dist
  .dist_deg.main   <- if (!is.null(.ov$deg.main))   .ov$deg.main   else netparams$casl$deg.main.dist
  .dist_deg.tot    <- if (!is.null(.ov$deg.tot))    .ov$deg.tot    else netparams$inst$deg.tot.dist
  .dist_role.class <- if (!is.null(.ov$role.class)) .ov$role.class else netparams$all$role.type
  .dist_risk.grp   <- if (!is.null(.ov$risk.grp))   .ov$risk.grp   else rep(1 / nquants, nquants)

  if (.tp$form == "data.frame") {
    # ---- target_pop = data.frame: pull attributes directly from user df ----
    df <- .tp$df
    attr_age <- df$age
    attr_sqrt.age <- if (!is.null(df$sqrt.age)) df$sqrt.age else sqrt(attr_age)
    attr_age.grp <- if (!is.null(df$age.grp)) df$age.grp else
      cut(attr_age, age.breaks, labels = FALSE, right = FALSE, include.lowest = FALSE)
    attr_active.sex <- if (!is.null(df$active.sex)) {
      as.integer(df$active.sex)
    } else {
      as_active <- rep(1L, num)
      if (sex.cess.mod == TRUE) {
        as_active[attr_age.grp == max(attr_age.grp, na.rm = TRUE)] <- 0L
      }
      as_active
    }
    if (!is.null(df$race)) {
      attr_race <- df$race
    } else {
      # race column is required when race = TRUE (validated in
      # .parse_target_pop); otherwise sample from race.dist for parity
      # with the legacy code path which always populates out$attr$race.
      race_numbers <- vapply(flattened_race_level,
        function(r) out$demog[[paste0("num.", r)]], numeric(1))
      attr_race <- apportion_lr(num, seq_along(flattened_race_level),
                                race_numbers / num, shuffled = TRUE)
    }
    attr_deg.casl   <- df$deg.casl
    attr_deg.main   <- df$deg.main
    attr_deg.tot    <- if (!is.null(df$deg.tot)) df$deg.tot else
      pmin(attr_deg.main + attr_deg.casl, 3L)
    attr_risk.grp   <- df$risk.grp
    attr_role.class <- df$role.class
  } else {
    # ---- Existing sampling path (default + list-form marginal overrides) ----
    if (!is.null(age.pyramid)) {
      if (length(age.pyramid) != nAges) {
        stop("Length of age.pyramid vector must be equal to length of unique age values: ", nAges)
      }
    } else {
      full.age.pyr <- c(0.01202, 0.01228, 0.01250, 0.01280, 0.01292, 0.01289,
                        0.01284, 0.01286, 0.01301, 0.01297, 0.01296, 0.01337,
                        0.01344, 0.01334, 0.01329, 0.01332, 0.01325, 0.01323,
                        0.01360, 0.01385, 0.01365, 0.01368, 0.01373, 0.01389,
                        0.01424, 0.01454, 0.01480, 0.01514, 0.01533, 0.01525,
                        0.01463, 0.01426, 0.01400, 0.01401, 0.01404, 0.01352,
                        0.01366, 0.01360, 0.01339, 0.01367, 0.01274, 0.01246,
                        0.01227, 0.01190, 0.01226, 0.01181, 0.01192, 0.01245,
                        0.01312, 0.01331, 0.01255, 0.01225, 0.01219, 0.01239,
                        0.01308, 0.01321, 0.01313, 0.01303, 0.01311, 0.01320,
                        0.01264, 0.01247, 0.01223, 0.01169, 0.01152, 0.01092,
                        0.01042, 0.00994, 0.00954, 0.00926, 0.00890, 0.00872,
                        0.00899, 0.00650, 0.00630, 0.00600, 0.00601, 0.00509,
                        0.00451, 0.00414, 0.00377, 0.00346, 0.00304, 0.00275,
                        0.00235, 0.00202, 0.00173, 0.00148, 0.00127, 0.00109,
                        0.00093, 0.00080, 0.00068, 0.00059, 0.00050, 0.00043,
                        0.00037, 0.00032, 0.00027, 0.00023)
      age.pyramid <- full.age.pyr[age.vals]
    }

    if (age.sexual.cessation < age.limits[2] && !is.null(young.prop)) {
      age.break <- age.sexual.cessation - (age.limits[1] - 1)
      age.pyramid <- reweight_age_pyr(age.pyramid, young.prop, age.break)
    }

    attr_age <- sample(x = age.vals, size = num, prob = age.pyramid, replace = TRUE)
    age_noise <- runif(num)
    attr_age <- attr_age + age_noise

    attr_sqrt.age <- sqrt(attr_age)

    attr_age.grp <- cut(attr_age, age.breaks, labels = FALSE, right = FALSE, include.lowest = FALSE)

    # sexually active attribute
    attr_active.sex <- rep(1L, num)
    if (sex.cess.mod == TRUE) {
      attr_active.sex[attr_age.grp == max(attr_age.grp)] <- 0L
    }

    # race attribute
    race_numbers <- vapply(
      flattened_race_level,
      function(race) {
        race_num_var <- paste0("num.", race)
        out$demog[[race_num_var]]
        },
      numeric(1)
      )

    race_proportions <- race_numbers / num
    group_ids <- seq_along(flattened_race_level)
    attr_race <- apportion_lr(num, group_ids, race_proportions, shuffled = TRUE)

    # deg.casl attribute
    attr_deg.casl <- apportion_lr(num, 0:3, .dist_deg.casl, shuffled = TRUE)
    if (sex.cess.mod == TRUE) {
      attr_deg.casl[attr_active.sex == 0] <- 0
    }

    # deg main attribute
    attr_deg.main <- apportion_lr(num, 0:2, .dist_deg.main, shuffled = TRUE)
    if (sex.cess.mod == TRUE) {
      attr_deg.main[attr_active.sex == 0] <- 0
    }

    # deg tot 3 attribute
    attr_deg.tot <- apportion_lr(num, 0:3, .dist_deg.tot, shuffled = TRUE)
    if (sex.cess.mod == TRUE) {
      attr_deg.tot[attr_active.sex == 0] <- 0
    }

    # risk group
    attr_risk.grp <- apportion_lr(num, 1:nquants, .dist_risk.grp, shuffled = TRUE)

    # role class
    attr_role.class <- apportion_lr(num, 0:2, .dist_role.class, shuffled = TRUE)
  }

  # Common attr assignments (both paths) -----------------------------------
  out$attr$age <- attr_age
  out$attr$sqrt.age <- attr_sqrt.age
  out$attr$age.grp <- attr_age.grp
  out$attr$active.sex <- attr_active.sex
  out$attr$race <- attr_race
  out$attr$deg.casl <- attr_deg.casl
  out$attr$deg.main <- attr_deg.main
  out$attr$deg.tot <- attr_deg.tot
  out$attr$risk.grp <- attr_risk.grp
  out$attr$role.class <- attr_role.class

  # diag status
  if (.tp$form == "data.frame" && !is.null(.tp$df$diag.status)) {
    # User-supplied diag.status takes precedence over epistats-based draw.
    out$attr$diag.status <- as.integer(.tp$df$diag.status)
  } else if (is.null(epistats$init.hiv.prev)) {
    if (race == TRUE) {
      xs <- data.frame(age = attr_age, race.cat.num = attr_race, geogYN = 1)
      preds <- predict(epistats$hiv.mod, newdata = xs, type = "response")
      attr_diag.status <- rbinom(num, 1, preds)
      out$attr$diag.status <- attr_diag.status
    }  else {
      xs <- data.frame(age = attr_age, geogYN = 1)
      preds <- predict(epistats$hiv.mod, newdata = xs, type = "response")
      attr_diag.status <- rbinom(num, 1, preds)
      out$attr$diag.status <- attr_diag.status
    }
  } else {
    if (race == TRUE) {
      init.hiv.prev <- epistats$init.hiv.prev
      attr_diag.status <- rep(0, network.size)

      # Randomly applies HIV diagnosis based on race
      for (i in seq_along(flattened_race_level)) {
        samp <- which(attr_race == i)
        exp <- ceiling(length(samp) * init.hiv.prev[i])
        attr_diag.status[sample(samp, exp)] <- 1
      }

      out$attr$diag.status <- attr_diag.status

    } else {
      init.hiv.prev <- epistats$init.hiv.prev[1]
      samp.size <- ceiling(network.size * init.hiv.prev)
      attr_diag.status <- sample(1:network.size, samp.size)
      out$attr$diag.status <- rep(0, network.size)
      out$attr$diag.status[attr_diag.status] <- 1
    }
  }
  if (!is.integer(out$attr$diag.status)) {
    out$attr$diag.status <- as.integer(out$attr$diag.status)
  }


  # Joint g-computation predictions (method = "joint" only) -----------------
  # Predict per-synthetic-node expected degree from the joint Poisson GLMs
  # fit in build_netparams(..., method = "joint"). Aggregating these per-node
  # predictions gives target statistics that are internally consistent by
  # construction (sum(nodefactor_<attr>) == 2 * edges). Inactive-age nodes
  # (sex.cess.mod) are zeroed out so they contribute nothing to any layer's
  # edges or nodefactor counts.
  # Initialize synth-aggregated duration overrides; populated below under
  # method = "joint" + duration.method = "joint_lm". When NULL, the dissolution
  # offsets fall back to the within-ARTnet values from build_netparams.
  synth_dur_main_byage <- NULL
  synth_dur_casl_byage <- NULL

  if (method == "joint") {
    synth <- data.frame(
      age.grp      = out$attr$age.grp,
      race.cat.num = out$attr$race,
      deg.main     = out$attr$deg.main,
      deg.casl     = out$attr$deg.casl,
      hiv2         = out$attr$diag.status,
      geogYN       = 1L
    )
    synth$deg.tot3 <- pmin(out$attr$deg.tot, 3)
    # Alias for the joint duration model, which uses index.age.grp on the RHS
    # rather than age.grp (a leftover from the partnership-level fit on lmain
    # in build_netparams). Joint_nm_*_model uses age.grp; we keep both columns
    # available so predict() works for both families of models.
    synth$index.age.grp <- synth$age.grp

    pred_deg_main   <- predict(netparams$main$joint_model, newdata = synth, type = "response")
    pred_deg_casl   <- predict(netparams$casl$joint_model, newdata = synth, type = "response")
    pred_count_inst <- predict(netparams$inst$joint_model, newdata = synth, type = "response")
    pred_deg_inst   <- pred_count_inst / (364 / time.unit)

    # Concurrency probabilities from the joint binomial GLMs on the
    # deg.{main,casl}.conc indicators (fit alongside the Poisson in
    # build_netparams under method = "joint"). No concurrent target for
    # the one-off layer.
    pred_conc_main <- predict(netparams$main$joint_concurrent_model,
                              newdata = synth, type = "response")
    pred_conc_casl <- predict(netparams$casl$joint_concurrent_model,
                              newdata = synth, type = "response")

    # Dyad-level predictions (nodematch / absdiff). Each is a per-ego
    # expected value of a partnership property (P(same.attr), E[|age-p_age|])
    # from a partnership-level GLM with only ego attributes on the RHS.
    # Aggregation over edges: sum(pred_deg * pred_dyad) / 2 (each edge is
    # counted from both endpoints).
    .predict_dyad <- function(layer) {
      list(
        pred_nm_age_grp  = predict(netparams[[layer]]$joint_nm_age_model,
                                   newdata = synth, type = "response"),
        pred_nm_race     = if (race) predict(netparams[[layer]]$joint_nm_race_model,
                                             newdata = synth, type = "response") else NULL,
        pred_ad_age      = predict(netparams[[layer]]$joint_absdiff_age_model,
                                   newdata = synth, type = "response"),
        pred_ad_sqrtage  = predict(netparams[[layer]]$joint_absdiff_sqrtage_model,
                                   newdata = synth, type = "response")
      )
    }
    dyad_main <- .predict_dyad("main")
    dyad_casl <- .predict_dyad("casl")
    dyad_inst <- .predict_dyad("inst")

    if (sex.cess.mod == TRUE) {
      inactive <- out$attr$active.sex == 0L
      pred_deg_main[inactive]  <- 0
      pred_deg_casl[inactive]  <- 0
      pred_deg_inst[inactive]  <- 0
      pred_conc_main[inactive] <- 0
      pred_conc_casl[inactive] <- 0
      # dyad predictions are multiplied by pred_deg downstream, so
      # zeroing pred_deg above already suppresses their contribution.
    }

    # Synth-aggregated stratum durations (#73). Override the within-ARTnet
    # joint_lm aggregation that build_netparams emits with synth-population
    # aggregation. Only fires when duration.method = "joint_lm" was used in
    # build_netparams (otherwise joint_duration_model is NULL).
    n_age_grps_main <-
      nrow(netparams$main$durs.main.byage) - 1L - as.integer(sex.cess.mod)
    n_age_grps_casl <-
      nrow(netparams$casl$durs.casl.byage) - 1L - as.integer(sex.cess.mod)
    synth_dur_main_byage <- .aggregate_synth_byage_durations(
      joint_dur_model     = netparams$main$joint_duration_model,
      joint_nm_race_model = netparams$main$joint_nm_race_model,
      synth_data          = synth,
      n_age_grps          = n_age_grps_main,
      sex_cess_extra_row  = sex.cess.mod
    )
    synth_dur_casl_byage <- .aggregate_synth_byage_durations(
      joint_dur_model     = netparams$casl$joint_duration_model,
      joint_nm_race_model = netparams$casl$joint_nm_race_model,
      synth_data          = synth,
      n_age_grps          = n_age_grps_casl,
      sex_cess_extra_row  = sex.cess.mod
    )
  }


  # Main Model -----------------------------------------------------------

  out$main <- list()

  if (method == "existing") {

    # edges ---
    if (race == TRUE) {
      if (edges.avg == FALSE) {
        out$main$edges <- (netparams$main$md.main * num) / 2
      } else {
        out$main$edges <- sum(unname(table(out$attr$race)) * netparams$main$nf.race) / 2
      }
      # nodefactor("race") ---
      nodefactor_race <- table(out$attr$race) * netparams$main$nf.race
      out$main$nodefactor_race <- unname(nodefactor_race)

      # nodematch("race") ---
      nodematch_race <- nodefactor_race / 2 * netparams$main$nm.race
      out$main$nodematch_race <- unname(nodematch_race)

      # nodematch("race", diff = FALSE) ---
      nodematch_race <- out$main$edges * netparams$main$nm.race_diffF
      out$main$nodematch_race_diffF <- unname(nodematch_race)
    } else {
      out$main$edges <- (netparams$main$md.main * num) / 2
    }

    # nodefactor("age.grp") ---
    nodefactor_age.grp <- table(out$attr$age.grp) * netparams$main$nf.age.grp
    if (sex.cess.mod == TRUE) {
      nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
    }
    out$main$nodefactor_age.grp <- unname(nodefactor_age.grp)

    # nodematch("age.grp") ---
    nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$main$nm.age.grp
    out$main$nodematch_age.grp <- unname(nodematch_age.grp)

    # absdiff("age") ---
    absdiff_age <- out$main$edges * netparams$main$absdiff.age
    out$main$absdiff_age <- absdiff_age

    # absdiff("sqrt.age") ---
    absdiff_sqrt.age <- out$main$edges * netparams$main$absdiff.sqrt.age
    out$main$absdiff_sqrt.age <- absdiff_sqrt.age

    # nodefactor("deg.casl") ---
    out$main$nodefactor_deg.casl <- num * netparams$main$deg.casl.dist * netparams$main$nf.deg.casl

    # concurrent ---
    out$main$concurrent <- num * netparams$main$concurrent

    # nodefactor("diag.status") ---
    nodefactor_diag.status <- table(out$attr$diag.status) * netparams$main$nf.diag.status
    out$main$nodefactor_diag.status <- unname(nodefactor_diag.status)

  } else {  # method == "joint"

    # edges from summed per-node predictions --------------------------------
    out$main$edges <- sum(pred_deg_main) / 2

    # Per-ego expected edge contributions to dyad-level sums. edge_wt[i]
    # is what ego i contributes to edges × expected-value-per-partnership,
    # so that sum(edge_wt * dyad_pred) / 2 aggregates across all edges
    # with correct double-counting from both endpoints.
    edge_wt_main <- pred_deg_main

    # nodefactor("race") + nodematch("race") --------------------------------
    if (race == TRUE) {
      n_race <- length(netparams$main$nf.race)
      out$main$nodefactor_race <- vapply(seq_len(n_race),
        function(r) sum(pred_deg_main[out$attr$race == r]), numeric(1))
      # nodematch_race[r] = expected count of race-r-to-race-r edges
      race_match_contrib <- edge_wt_main * dyad_main$pred_nm_race
      out$main$nodematch_race <- vapply(seq_len(n_race),
        function(r) sum(race_match_contrib[out$attr$race == r]) / 2,
        numeric(1))
      out$main$nodematch_race_diffF <- sum(race_match_contrib) / 2
    }

    # nodefactor("age.grp") + nodematch("age.grp") --------------------------
    n_age <- length(netparams$main$nf.age.grp)
    out$main$nodefactor_age.grp <- vapply(seq_len(n_age),
      function(k) sum(pred_deg_main[out$attr$age.grp == k]), numeric(1))
    age_match_contrib <- edge_wt_main * dyad_main$pred_nm_age_grp
    out$main$nodematch_age.grp <- vapply(seq_len(n_age),
      function(k) sum(age_match_contrib[out$attr$age.grp == k]) / 2,
      numeric(1))

    # absdiff -- sum over edges of expected |age_i - age_j| --------------
    out$main$absdiff_age <- sum(edge_wt_main * dyad_main$pred_ad_age) / 2
    out$main$absdiff_sqrt.age <- sum(edge_wt_main * dyad_main$pred_ad_sqrtage) / 2

    # nodefactor("deg.casl") ------------------------------------------------
    out$main$nodefactor_deg.casl <- vapply(0:3,
      function(d) sum(pred_deg_main[out$attr$deg.casl == d]), numeric(1))

    # concurrent from the joint binomial GLM on deg.main.conc ----------------
    out$main$concurrent <- sum(pred_conc_main)

    # nodefactor("diag.status") ---------------------------------------------
    out$main$nodefactor_diag.status <- vapply(0:1,
      function(h) sum(pred_deg_main[out$attr$diag.status == h]), numeric(1))
  }

  # Dissolution -----------------------------------------------------------
  # diss.byage uses synth-aggregated durations under method = "joint" +
  # duration.method = "joint_lm" (#73). diss.homog still uses the within-
  # ARTnet aggregation from build_netparams; it is not consumed by
  # EpiModelHIV-Template's tergm offset, and its synth analog can be added
  # later without changing the byage interface that matters for production.
  out$main$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$main$durs.main.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$main$diss.byage <- dissolution_coefs(
    dissolution = ~offset(edges) + offset(nodematch("age.grp", diff = TRUE)),
    duration = if (!is.null(synth_dur_main_byage)) synth_dur_main_byage
               else netparams$main$durs.main.byage$mean.dur.adj,
    d.rate = expect.mort
  )



  # Casual Model ------------------------------------------------------------

  out$casl <- list()

  if (method == "existing") {

    # edges ---
    if (race == TRUE) {
      if (edges.avg == FALSE) {
        out$casl$edges <- (netparams$casl$md.casl * num) / 2
      } else {
        out$casl$edges <- sum(unname(table(out$attr$race)) * netparams$casl$nf.race) / 2
      }
      # nodefactor("race") ---
      nodefactor_race <- table(out$attr$race) * netparams$casl$nf.race
      out$casl$nodefactor_race <- unname(nodefactor_race)

      # nodematch("race") ---
      nodematch_race <- nodefactor_race / 2 * netparams$casl$nm.race
      out$casl$nodematch_race <- unname(nodematch_race)

      # nodematch("race", diff = FALSE) ---
      nodematch_race <- out$casl$edges * netparams$casl$nm.race_diffF
      out$casl$nodematch_race_diffF <- unname(nodematch_race)
    } else {
      out$casl$edges <- (netparams$casl$md.casl * num) / 2
    }

    # nodefactor("age.grp") ---
    nodefactor_age.grp <- table(out$attr$age.grp) * netparams$casl$nf.age.grp
    if (sex.cess.mod == TRUE) {
      nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
    }
    out$casl$nodefactor_age.grp <- unname(nodefactor_age.grp)

    # nodematch("age.grp") ---
    nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$casl$nm.age.grp
    out$casl$nodematch_age.grp <- unname(nodematch_age.grp)

    # absdiff("age") ---
    absdiff_age <- out$casl$edges * netparams$casl$absdiff.age
    out$casl$absdiff_age <- absdiff_age

    # absdiff("sqrt.age") ---
    absdiff_sqrt.age <- out$casl$edges * netparams$casl$absdiff.sqrt.age
    out$casl$absdiff_sqrt.age <- absdiff_sqrt.age

    # nodefactor("deg.main") ---
    out$casl$nodefactor_deg.main <- num * netparams$casl$deg.main.dist * netparams$casl$nf.deg.main

    # concurrent ---
    out$casl$concurrent <- num * netparams$casl$concurrent

    # nodefactor("diag.status") ---
    nodefactor_diag.status <- table(out$attr$diag.status) * netparams$casl$nf.diag.status
    out$casl$nodefactor_diag.status <- unname(nodefactor_diag.status)

  } else {  # method == "joint"

    out$casl$edges <- sum(pred_deg_casl) / 2
    edge_wt_casl <- pred_deg_casl

    if (race == TRUE) {
      n_race <- length(netparams$casl$nf.race)
      out$casl$nodefactor_race <- vapply(seq_len(n_race),
        function(r) sum(pred_deg_casl[out$attr$race == r]), numeric(1))
      race_match_contrib <- edge_wt_casl * dyad_casl$pred_nm_race
      out$casl$nodematch_race <- vapply(seq_len(n_race),
        function(r) sum(race_match_contrib[out$attr$race == r]) / 2,
        numeric(1))
      out$casl$nodematch_race_diffF <- sum(race_match_contrib) / 2
    }

    n_age <- length(netparams$casl$nf.age.grp)
    out$casl$nodefactor_age.grp <- vapply(seq_len(n_age),
      function(k) sum(pred_deg_casl[out$attr$age.grp == k]), numeric(1))
    age_match_contrib <- edge_wt_casl * dyad_casl$pred_nm_age_grp
    out$casl$nodematch_age.grp <- vapply(seq_len(n_age),
      function(k) sum(age_match_contrib[out$attr$age.grp == k]) / 2,
      numeric(1))

    out$casl$absdiff_age <- sum(edge_wt_casl * dyad_casl$pred_ad_age) / 2
    out$casl$absdiff_sqrt.age <- sum(edge_wt_casl * dyad_casl$pred_ad_sqrtage) / 2

    # nodefactor("deg.main") -- truncated at 2 in build_netparams
    out$casl$nodefactor_deg.main <- vapply(0:2,
      function(d) sum(pred_deg_casl[out$attr$deg.main == d]), numeric(1))

    # concurrent from the joint binomial GLM on deg.casl.conc ----------------
    out$casl$concurrent <- sum(pred_conc_casl)

    out$casl$nodefactor_diag.status <- vapply(0:1,
      function(h) sum(pred_deg_casl[out$attr$diag.status == h]), numeric(1))
  }

  # Dissolution (see note on diss.byage / diss.homog at the main layer block)
  out$casl$diss.homog <- dissolution_coefs(dissolution = ~offset(edges),
                                           duration = netparams$casl$durs.casl.homog$mean.dur.adj,
                                           d.rate = expect.mort)
  out$casl$diss.byage <- dissolution_coefs(
    dissolution = ~offset(edges) + offset(nodematch("age.grp", diff = TRUE)),
    duration = if (!is.null(synth_dur_casl_byage)) synth_dur_casl_byage
               else netparams$casl$durs.casl.byage$mean.dur.adj,
    d.rate = expect.mort
  )


  # One-Time Model ----------------------------------------------------------

  out$inst <- list()

  if (method == "existing") {

    # edges ---
    if (race == TRUE) {
      if (edges.avg == FALSE) {
        out$inst$edges <- (netparams$inst$md.inst * num) / 2
      } else {
        out$inst$edges <- sum(unname(table(out$attr$race)) * netparams$inst$nf.race) / 2
      }
      # nodefactor("race") ---
      nodefactor_race <- table(out$attr$race) * netparams$inst$nf.race
      out$inst$nodefactor_race <- unname(nodefactor_race)

      # nodematch("race") ---
      nodematch_race <- nodefactor_race / 2 * netparams$inst$nm.race
      out$inst$nodematch_race <- unname(nodematch_race)

      # nodematch("race", diff = FALSE) ---
      nodematch_race <- out$inst$edges * netparams$inst$nm.race_diffF
      out$inst$nodematch_race_diffF <- unname(nodematch_race)
    } else {
      out$inst$edges <- (netparams$inst$md.inst * num) / 2
    }

    # nodefactor("age.grp") ---
    nodefactor_age.grp <- table(out$attr$age.grp) * netparams$inst$nf.age.grp
    if (sex.cess.mod == TRUE) {
      nodefactor_age.grp[length(nodefactor_age.grp)] <- 0
    }
    out$inst$nodefactor_age.grp <- unname(nodefactor_age.grp)

    # nodematch("age.grp") ---
    nodematch_age.grp <- nodefactor_age.grp / 2 * netparams$inst$nm.age.grp
    out$inst$nodematch_age.grp <- unname(nodematch_age.grp)

    # absdiff("age") ---
    absdiff_age <- out$inst$edges * netparams$inst$absdiff.age
    out$inst$absdiff_age <- absdiff_age

    # absdiff("sqrt.age") ---
    absdiff_sqrt.age <- out$inst$edges * netparams$inst$absdiff.sqrt.age
    out$inst$absdiff_sqrt.age <- absdiff_sqrt.age

    # nodefactor("risk.grp") ---
    nodefactor_risk.grp <- table(out$attr$risk.grp) * netparams$inst$nf.risk.grp
    out$inst$nodefactor_risk.grp <- unname(nodefactor_risk.grp)

    # nodefactor("deg.tot") ---
    nodefactor_deg.tot <- table(out$attr$deg.tot) * netparams$inst$nf.deg.tot
    out$inst$nodefactor_deg.tot <- unname(nodefactor_deg.tot)

    # nodefactor("diag.status") ---
    nodefactor_diag.status <- table(out$attr$diag.status) * netparams$inst$nf.diag.status
    out$inst$nodefactor_diag.status <- unname(nodefactor_diag.status)

  } else {  # method == "joint"

    out$inst$edges <- sum(pred_deg_inst) / 2
    edge_wt_inst <- pred_deg_inst

    if (race == TRUE) {
      n_race <- length(netparams$inst$nf.race)
      out$inst$nodefactor_race <- vapply(seq_len(n_race),
        function(r) sum(pred_deg_inst[out$attr$race == r]), numeric(1))
      race_match_contrib <- edge_wt_inst * dyad_inst$pred_nm_race
      out$inst$nodematch_race <- vapply(seq_len(n_race),
        function(r) sum(race_match_contrib[out$attr$race == r]) / 2,
        numeric(1))
      out$inst$nodematch_race_diffF <- sum(race_match_contrib) / 2
    }

    n_age <- length(netparams$inst$nf.age.grp)
    out$inst$nodefactor_age.grp <- vapply(seq_len(n_age),
      function(k) sum(pred_deg_inst[out$attr$age.grp == k]), numeric(1))
    age_match_contrib <- edge_wt_inst * dyad_inst$pred_nm_age_grp
    out$inst$nodematch_age.grp <- vapply(seq_len(n_age),
      function(k) sum(age_match_contrib[out$attr$age.grp == k]) / 2,
      numeric(1))

    out$inst$absdiff_age <- sum(edge_wt_inst * dyad_inst$pred_ad_age) / 2
    out$inst$absdiff_sqrt.age <- sum(edge_wt_inst * dyad_inst$pred_ad_sqrtage) / 2

    # nodefactor("risk.grp") -- quantile-scheme preserved from univariate
    nodefactor_risk.grp <- table(out$attr$risk.grp) * netparams$inst$nf.risk.grp
    out$inst$nodefactor_risk.grp <- unname(nodefactor_risk.grp)

    # nodefactor("deg.tot") -- truncated at 3 in build_netparams
    out$inst$nodefactor_deg.tot <- vapply(0:3,
      function(d) sum(pred_deg_inst[out$attr$deg.tot == d]), numeric(1))

    out$inst$nodefactor_diag.status <- vapply(0:1,
      function(h) sum(pred_deg_inst[out$attr$diag.status == h]), numeric(1))
  }


  return(out)
}


#' Update mortality rates
#'
#' @description Replaces the default age- and sex-specific mortality rates used
#'              in [`build_netstats`] with user-specified values.
#'
#' @param netstats Output from [`build_netstats`].
#' @param asmr_df A data frame with three columns (Race, Age, and DeathRate)
#'                specifying the desired mortality rates per time unit for each
#'                combination of race (Black, Hispanic, White) and age
#'                (1, 2, 3, ... 100).
#'
#' @details
#' `update_asmr` takes the output from [`build_netstats`] and
#' updates its mortality rates with the values specified in the input parameter
#' `asmr_df`. Mortality rates should be provided for each combination of
#' race (Black, Hispanic, White) and age (1, 2, 3, ... 100), for a total of 300
#' rates. All input mortality rates must be non-NA, numeric, and between 0 and
#' 1 (inclusive). The maximum mortality rate for each race must be 1,
#' representing deterministic departure from the network at that age.
#' `update_asmr` does not perform time unit conversions; the input
#' mortality rates should be pre-converted by the user to match
#' the `time.unit` used in [`build_netstats`].
#'
#' @export
#'
#' @examples
#' epistats  <- build_epistats(geog.lvl = "city", geog.cat = "Atlanta",
#'                             race = TRUE)
#' netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)
#' netstats  <- build_netstats(epistats, netparams)
#'
#' # Update mortality rates with random values
#' asmr_df   <- data.frame(Race = rep(c("Black", "Hispanic", "White"),
#'                                    each = 100), Age = rep(1:100,3),
#'                         DeathRate = rep(c(runif(99), 1), 3))
#' netstats  <- update_asmr(netstats, asmr_df)
#'
#' # Update mortality rates with 2020 NCHS data
#' asmr_df   <- read.csv(system.file("2020DeathRates.csv", package =  "ARTnet"))
#' netstats  <- update_asmr(netstats, asmr_df)
#'
update_asmr <- function(netstats, asmr_df) {
  # Check that all input mortality rates are reasonable values
  if (sum(is.na(asmr_df$DeathRate)) > 0 || !is.numeric(asmr_df$DeathRate) ||
      min(asmr_df$DeathRate) < 0 || max(asmr_df$DeathRate) > 1) {
    stop("Ensure all mortality rates are non-NA, numeric, and between 0 and 1.")
  }

  # Get unique races
  unique_races <- unique(asmr_df$Race)

  # Initialize a list to store mortality rates by race
  asmr_list <- list()

  # Validate mortality rates for each race
  for (race in unique_races) {
    race_data <- asmr_df[asmr_df$Race == race, 2:3]
    if (min(race_data$Age) != 1 || max(race_data$Age) != 100 || nrow(race_data) != 100) {
      stop(paste0("Provide mortality rates by race for 1-yr age groups (1 - 100) for race: ", race))
    }
    if (max(race_data$DeathRate) != 1) {
      stop(paste0("The mortality rate for at least one age group must be total (1) for race: ", race))
    }
    asmr_list[[race]] <- race_data$DeathRate
  }

  # If race-specific mortality rates should be included
  if (netstats$race == TRUE) {
    # Combine mortality rates for all races into a data frame
    asmr <- data.frame(
      age = 1:100,
      do.call(cbind, asmr_list)
    )

    colnames(asmr)[-1] <- paste0("asmr.", names(asmr_list))

  } else {
    # Calculate combined mortality rates using demographic proportions
    props <- netstats$demog$props
    asmr_combined <- rep(0, 100)
    for (race in unique_races) {
      if (!race %in% names(props)) {
        stop(paste0("Demographic proportions for race: ", race, " are missing in netstats$demog$props."))
      }
      asmr_combined <- asmr_combined + props[[race]] * asmr_list[[race]]
    }
    # Create a combined mortality data frame
    asmr <- data.frame(age = 1:100, asmr_combined, asmr_combined, asmr_combined)
    colnames(asmr)[-1] <- c("asmr.O1", "asmr.O2", "asmr.O3")
  }

  # Display time unit consistency message
  message(paste0("The time unit used in build_netstats() was ",
                 netstats$time.unit,
                 ". Please ensure that that is consistent with your inputs."))

  # Assign the resulting asmr to netstats
  netstats$demog$asmr <- asmr
  return(netstats)
}

#' Re-weight age pyramid
#'
#' @description This function re-weights a given age pyramid (i.e., a vector
#' with one element per year in an age range that specifies what proportion of
#' the population should be  that age) so that the desired fraction of the
#' population falls below a specified age boundary.
#'
#' @param age.pyramid A vector of proportions.
#' @param young.prop The weight that should be applied to the portion of
#'                   `age.pyramid` that is below `age.break`
#' @param age.break The index at which the upper group should begin.
#'
#' @export
#'
#' @examples
#' unif.age.pyr <- rep(1/100, 100)
#' reweighted.age.pyr <- reweight_age_pyr(unif.age.pyr, 0.995, 65)
#'
reweight_age_pyr <- function(age.pyramid, young.prop, age.break) {

  #Check inputs
  if (age.break > length(age.pyramid) || young.prop > 1 || young.prop < 0 ||
      sum(is.na(age.pyramid)) > 0 || !is.numeric(age.pyramid) ||
      min(age.pyramid) < 0 || max(age.pyramid) > 1) {
    stop("Ensure that all input values are non-NA, numeric, and in range.")
  }

  #Re-weighting
  age1 <- age.pyramid[1:(age.break - 1)]
  age1 <- age1 / sum(age1) * young.prop

  age2 <- age.pyramid[age.break:length(age.pyramid)]
  age2 <- age2 / sum(age2) * (1 - young.prop)

  new.age.pyramid <- c(age1, age2)

  return(new.age.pyramid)
}

#' Reduces the Size of the Netstats Object by Removing the Networks it Contains
#'
#' Once `netest` has been built, the network object stored in `netstats` are
#' never used again and take up a lot of space. This function removes them to
#' reduce the stored size of `netstats`. Roughly 10MiB are saved for 100k nodes
#' networks.
#'
#' @param netstats the `netstats` object to be trimmed
#'
#' @return a trimmed `netstats`
#'
#' @export
#'
trim_netstats <- function(netstats) {
  netstats$main <- NULL
  netstats$casl <- NULL
  netstats$inst <- NULL
  return(netstats)
}