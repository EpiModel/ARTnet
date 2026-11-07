## Method comparison: marginal vs joint g-computation across scenarios.
## Phase 1.5 of the joint g-comp refactor (issue #65). Produces
## side-by-side target-stat comparison tables under matched parameter
## sets, identifies stats whose values shift materially under the
## method correction, and writes the result as a markdown report.
##
## Usage:
##   source(system.file("validation/method_comparison.R", package = "ARTnet"))
##   res <- compare_methods()
##   summarize_comparison(res)
##   render_comparison_report(res, file = "inst/validation/method_comparison.md")
##
## Requires ARTnetData. The full run takes ~30s on a recent laptop
## (4 scenarios x 2 methods x build_epistats + build_netparams + build_netstats).

# ---- Scenario definitions ----------------------------------------------------

# Each scenario: a list with $name and the args to build_{epistats,netparams,netstats}.
# Four city defaults selected because they are the cities we have applied
# ARTnet for EpiModelHIV modeling with. Each uses ARTnetData::race.dist to
# define the synthetic population's race composition (no fabricated mixes).
# Race composition (W/Other / Black / Hispanic), by city:
#   Atlanta:  43.9 / 51.5 /  4.6   -- far from ARTnet sample (~81% W/O)
#   Boston:   59.7 / 20.9 / 19.4   -- moderate distance
#   Chicago:  41.8 / 29.2 / 29.0   -- moderate distance, balanced
#   Seattle:  87.4 /  6.1 /  6.5   -- closest to ARTnet sample composition
.make_city_scenario <- function(city) {
  list(
    name = paste0(tolower(gsub(" ", "_", city)), "_default"),
    epistats = list(geog.lvl = "city", geog.cat = city,
                    init.hiv.prev = c(0.33, 0.137, 0.084),
                    race = TRUE, time.unit = 7),
    netparams = list(smooth.main.dur = TRUE),
    # expect.mort = 0.00025 (slightly lower than the EpiModelHIV-Template
    # default 0.000478213) so dissolution_coefs accommodates Seattle's
    # noisy empirical matched.5 main duration of ~1381 weeks. The
    # absolute coef values shift, but the marginal-vs-joint comparison
    # is consistent across scenarios at the same expect.mort.
    netstats  = list(expect.mort = 0.00025, network.size = 5000)
  )
}
COMPARISON_SCENARIOS <- lapply(
  c("Atlanta", "Boston", "Chicago", "Seattle"),
  .make_city_scenario
)

.COMPARISON_SEED <- 20260420L


# ---- Utilities --------------------------------------------------------------

.require_artnetdata <- function() {
  if (system.file(package = "ARTnetData") == "") {
    stop("ARTnetData not installed; method comparison cannot run.")
  }
}

# Walk a netstats object and produce a long-format data.frame of every
# target statistic with columns scenario, method, layer, stat, level, value.
# `level` is NA for scalars; integer index for vector-valued stats.
.extract_target_stats <- function(netstats, scenario, method) {
  rows <- list()
  add_scalar <- function(layer, stat, x) {
    if (is.null(x) || length(x) != 1) return(NULL)
    data.frame(scenario = scenario, method = method, layer = layer,
               stat = stat, level = NA_integer_, value = as.numeric(x),
               stringsAsFactors = FALSE)
  }
  add_vec <- function(layer, stat, v) {
    if (is.null(v) || length(v) < 1) return(NULL)
    data.frame(scenario = scenario, method = method, layer = layer,
               stat = stat, level = seq_along(v), value = unname(as.numeric(v)),
               stringsAsFactors = FALSE)
  }
  for (layer in c("main", "casl", "inst")) {
    L <- netstats[[layer]]
    if (is.null(L)) next
    rows <- c(rows, list(
      add_scalar(layer, "edges", L$edges),
      add_scalar(layer, "concurrent", L$concurrent),
      add_scalar(layer, "nodematch_race_diffF", L$nodematch_race_diffF),
      add_scalar(layer, "absdiff_age", L$absdiff_age),
      add_scalar(layer, "absdiff_sqrt.age", L$absdiff_sqrt.age),
      add_vec(layer, "nodefactor_race", L$nodefactor_race),
      add_vec(layer, "nodefactor_age.grp", L$nodefactor_age.grp),
      add_vec(layer, "nodefactor_diag.status", L$nodefactor_diag.status),
      add_vec(layer, "nodefactor_deg.casl", L$nodefactor_deg.casl),
      add_vec(layer, "nodefactor_deg.main", L$nodefactor_deg.main),
      add_vec(layer, "nodefactor_deg.tot",  L$nodefactor_deg.tot),
      add_vec(layer, "nodefactor_risk.grp", L$nodefactor_risk.grp),
      add_vec(layer, "nodematch_race", L$nodematch_race),
      add_vec(layer, "nodematch_age.grp", L$nodematch_age.grp)
    ))
    # Stratum-level dissolution durations (main / casl only — inst is
    # ~offset(edges) with a fixed 1).
    if (layer != "inst" && !is.null(L$diss.byage$duration) &&
        length(L$diss.byage$duration) > 1) {
      rows <- c(rows, list(
        add_vec(layer, "dissolution_duration", L$diss.byage$duration)
      ))
    }
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

# Build the netstats object for one scenario under one method. The seed is
# reset before each build_* call so the three stochastic steps (HIV model,
# partner-race imputation, attribute sampling) are reproducible and the
# same across methods. The feasibility check is off so that infeasible
# target sets (Seattle under method = "existing") come back and can be
# reported rather than stopping the run; fit_comparison_ergms() runs the
# check itself. `floor = TRUE` turns on the random-mixing floor (#9).
.build_scenario_netstats <- function(scenario, method, floor = FALSE) {
  dur_method <- if (method == "joint") "joint_lm" else "empirical"
  set.seed(.COMPARISON_SEED)
  epistats <- do.call(ARTnet::build_epistats, scenario$epistats)
  np_args <- c(list(epistats = epistats), scenario$netparams,
               list(method = method, duration.method = dur_method,
                    seed = .COMPARISON_SEED))
  netparams <- do.call(ARTnet::build_netparams, np_args)
  ns_args <- c(list(epistats = epistats, netparams = netparams),
               scenario$netstats,
               list(method = method, seed = .COMPARISON_SEED,
                    nodematch.race.floor = floor,
                    check.feasibility = FALSE))
  suppressMessages(do.call(ARTnet::build_netstats, ns_args))
}

# Run one scenario through both methods, return both extracted stat tables.
.run_one_scenario <- function(scenario) {
  list(
    existing = .extract_target_stats(.build_scenario_netstats(scenario, "existing"),
                                     scenario$name, "existing"),
    joint    = .extract_target_stats(.build_scenario_netstats(scenario, "joint"),
                                     scenario$name, "joint")
  )
}


# ---- Public entry points -----------------------------------------------------

#' Run the method comparison across scenarios.
#'
#' @param scenarios List of scenario specs (defaults to COMPARISON_SCENARIOS).
#' @return Long-format data.frame with columns scenario, layer, stat, level,
#'   existing, joint, abs_diff, pct_diff (set to NA when existing == 0).
#'   Each row is a single (scenario, layer, stat, level) cell with both
#'   methods' values side-by-side.
compare_methods <- function(scenarios = COMPARISON_SCENARIOS) {
  .require_artnetdata()
  message("Running ", length(scenarios), " scenarios x 2 methods (",
          "this is build_*x6, expect ~30s)...")
  pieces <- list()
  for (s in scenarios) {
    message("  ", s$name)
    res <- .run_one_scenario(s)
    pieces[[s$name]] <- res
  }
  # Wide on method
  long_existing <- do.call(rbind, lapply(pieces, `[[`, "existing"))
  long_joint    <- do.call(rbind, lapply(pieces, `[[`, "joint"))
  key <- c("scenario", "layer", "stat", "level")
  wide <- merge(long_existing[, c(key, "value")],
                long_joint[, c(key, "value")],
                by = key, all = TRUE,
                suffixes = c("_existing", "_joint"))
  names(wide)[names(wide) == "value_existing"] <- "existing"
  names(wide)[names(wide) == "value_joint"]    <- "joint"
  wide$abs_diff <- wide$joint - wide$existing
  wide$pct_diff <- ifelse(abs(wide$existing) > 1e-12,
                          100 * wide$abs_diff / wide$existing,
                          NA_real_)
  wide[order(wide$scenario, wide$layer, wide$stat,
             ifelse(is.na(wide$level), -1L, wide$level)), , drop = FALSE]
}


#' Print a high-level summary of the comparison.
summarize_comparison <- function(comparison, threshold_pct = 5) {
  cat(sprintf("\nTotal cells: %d (across %d scenarios)\n",
              nrow(comparison), length(unique(comparison$scenario))))
  ok <- !is.na(comparison$pct_diff)
  cat(sprintf("Cells with |%%diff| > %g%%: %d\n",
              threshold_pct,
              sum(abs(comparison$pct_diff) > threshold_pct, na.rm = TRUE)))
  for (s in unique(comparison$scenario)) {
    sub <- comparison[comparison$scenario == s, , drop = FALSE]
    nbig <- sum(abs(sub$pct_diff) > threshold_pct, na.rm = TRUE)
    cat(sprintf("\n=== %s: %d cells, %d materially shifted (>%g%%) ===\n",
                s, nrow(sub), nbig, threshold_pct))
    if (nbig == 0) {
      cat("  (no material shifts)\n")
      next
    }
    big <- sub[!is.na(sub$pct_diff) & abs(sub$pct_diff) > threshold_pct, ]
    big <- big[order(-abs(big$pct_diff)), , drop = FALSE]
    n_show <- min(10, nrow(big))
    cat(sprintf("  Top %d by |%%diff|:\n", n_show))
    for (i in seq_len(n_show)) {
      r <- big[i, ]
      level_str <- if (is.na(r$level)) "" else sprintf("[%d]", r$level)
      cat(sprintf("    %-7s %-22s%-5s existing=%9.2f  joint=%9.2f  (%+0.1f%%)\n",
                  r$layer, r$stat, level_str, r$existing, r$joint, r$pct_diff))
    }
  }
  invisible(comparison)
}


#' Render a markdown report of the comparison results.
render_comparison_report <- function(comparison,
                                     file = "inst/validation/method_comparison.md",
                                     threshold_pct = 5) {
  con <- file(file, "w")
  on.exit(close(con))
  out <- function(...) cat(..., "\n", sep = "", file = con)

  out("# Method comparison: marginal vs joint g-computation")
  out()
  out("Generated by `inst/validation/method_comparison.R` on ", as.character(Sys.Date()),
      ". Phase 1.5 of the joint g-comp refactor; closes part of issue #65.")
  out()
  out("ARTnet version: ", as.character(packageVersion("ARTnet")), ". Seed: ",
      .COMPARISON_SEED, ". Network size: 5000.")
  out()
  out("## Scenarios")
  out()
  out("Derived from the scenarios actually present in the comparison object, ",
      "so this table always matches the data below.")
  out()
  out("| Scenario | Description |")
  out("|---|---|")
  for (s in unique(comparison$scenario)) {
    city <- gsub("_default$", "", s)
    city <- paste0(toupper(substring(city, 1, 1)), substring(city, 2))
    out("| `", s, "` | City config (", city,
        ", race = TRUE, init.hiv.prev = c(0.33, 0.137, 0.084)) |")
  }
  out()
  out("## High-level summary")
  out()
  ok <- !is.na(comparison$pct_diff)
  total_big <- sum(abs(comparison$pct_diff) > threshold_pct, na.rm = TRUE)
  out("- Total target-stat cells across scenarios: ", nrow(comparison))
  out("- Cells where |joint vs existing %diff| > ", threshold_pct, "%: **",
      total_big, "**")
  out()
  for (s in unique(comparison$scenario)) {
    sub <- comparison[comparison$scenario == s, , drop = FALSE]
    nbig <- sum(abs(sub$pct_diff) > threshold_pct, na.rm = TRUE)
    out("- `", s, "`: ", nrow(sub), " cells, ", nbig, " materially shifted (>",
        threshold_pct, "%)")
  }
  out()
  out("## Per-scenario top shifts (by |% diff|)")
  out()
  for (s in unique(comparison$scenario)) {
    sub <- comparison[comparison$scenario == s, , drop = FALSE]
    big <- sub[!is.na(sub$pct_diff) & abs(sub$pct_diff) > threshold_pct, ]
    big <- big[order(-abs(big$pct_diff)), , drop = FALSE]
    out("### ", s)
    out()
    if (nrow(big) == 0) {
      out("_No cells exceed |", threshold_pct, "%| threshold._")
      out()
      next
    }
    out("| Layer | Stat | Level | Existing | Joint | %Δ |")
    out("|---|---|---:|---:|---:|---:|")
    for (i in seq_len(min(15, nrow(big)))) {
      r <- big[i, ]
      level_str <- if (is.na(r$level)) "—" else as.character(r$level)
      out(sprintf("| %s | %s | %s | %.2f | %.2f | %+.1f%% |",
                  r$layer, r$stat, level_str,
                  r$existing, r$joint, r$pct_diff))
    }
    if (nrow(big) > 15) {
      out()
      out("_...and ", nrow(big) - 15, " more cells over threshold._")
    }
    out()
  }
  out("## Edges / concurrent comparison (all scenarios)")
  out()
  out("| Scenario | Layer | Existing edges | Joint edges | %Δ | Existing concurrent | Joint concurrent | %Δ |")
  out("|---|---|---:|---:|---:|---:|---:|---:|")
  edges <- comparison[comparison$stat == "edges", ]
  conc  <- comparison[comparison$stat == "concurrent", ]
  for (s in unique(comparison$scenario)) {
    for (l in c("main", "casl", "inst")) {
      e <- edges[edges$scenario == s & edges$layer == l, ]
      c1 <- conc[conc$scenario == s & conc$layer == l, ]
      if (nrow(e) == 0) next
      conc_existing <- if (nrow(c1) > 0) sprintf("%.2f", c1$existing) else "—"
      conc_joint    <- if (nrow(c1) > 0) sprintf("%.2f", c1$joint) else "—"
      conc_pct      <- if (nrow(c1) > 0 && !is.na(c1$pct_diff))
                         sprintf("%+.1f%%", c1$pct_diff) else "—"
      out(sprintf("| %s | %s | %.2f | %.2f | %+.1f%% | %s | %s | %s |",
                  s, l, e$existing, e$joint, e$pct_diff,
                  conc_existing, conc_joint, conc_pct))
    }
  }
  out()
  out("## Reproducibility")
  out()
  out("```r")
  out("source(system.file(\"validation/method_comparison.R\", package = \"ARTnet\"))")
  out("res <- compare_methods()")
  out("summarize_comparison(res)")
  out("render_comparison_report(res)")
  out("```")
  invisible(file)
}


# ---- ERGM fits with MCMLE (issue #12) ---------------------------------------
#
# Fits the three EpiModelHIV-Template ERGMs (formulas and target mappings
# as pinned in inst/validation/epimodelhiv_template_ref/) for every
# scenario x arm with MCMLE, the production estimation method, after running
# check_race_feasibility() on each target set. Stochastic-Approximation
# settles on the nearest reachable network when targets are infeasible, so
# it could not have caught the Seattle problem (#9); MCMLE fails on it.
#
# Arms:
#   existing        legacy targets
#   existing_floor  legacy targets with nodematch.race.floor = TRUE
#   joint           joint g-computation with balanced partner-race mixing
#
# Infeasible target sets are recorded and not fit. Each fit records whether
# netest() returned, whether MCMLE converged, whether ergm warned about a
# linear dependence (nodefactor.race = edges + nodematch.race, the
# feasibility boundary seen from the other side), the fit time, and from a
# static netdx the largest |Z| and |% diff| over the formation targets.
#
# Usage:
#   fits <- fit_comparison_ergms(ncores = 4)   # ~15-20 min on a laptop
#   render_fit_report(fits)                    # writes method_comparison_fits.md
#
# Requires EpiModel (and through it ergm) in addition to ARTnetData.

.ERGM_CONTROL <- function() {
  ergm::control.ergm(
    main.method = "MCMLE",
    MCMLE.maxit = 500,
    SAN.maxit = 3,
    SAN.nsteps.times = 4,
    MCMC.samplesize = 1e4,
    MCMC.interval = 5e3,
    # 0, not 1: ergm starts a one-worker PSOCK cluster at parallel = 1, and
    # concurrent jobs under mclapply then collide on the cluster port.
    parallel = 0
  )
}

# Formation formula, target.stats vector and dissolution coefs for one
# layer, verbatim from the template (tests/workflows/workflow-standard.R).
.template_layer_spec <- function(netstats, layer) {
  ns <- netstats[[layer]]
  if (layer == "main") {
    list(
      formation = ~ edges +
        nodematch("age.grp", diff = TRUE) +
        nodefactor("age.grp", levels = -1) +
        nodematch("race", diff = FALSE) +
        nodefactor("race", levels = -1) +
        nodefactor("deg.casl", levels = -1) +
        concurrent +
        degrange(from = 3) +
        nodematch("role.class", diff = TRUE, levels = c(1, 2)),
      target.stats = unname(c(
        ns$edges, ns$nodematch_age.grp, ns$nodefactor_age.grp[-1],
        ns$nodematch_race_diffF, ns$nodefactor_race[-1],
        ns$nodefactor_deg.casl[-1], ns$concurrent, 0, 0, 0)),
      coef.diss = ns$diss.byage
    )
  } else if (layer == "casl") {
    list(
      formation = ~ edges +
        nodematch("age.grp", diff = TRUE) +
        nodefactor("age.grp", levels = -5) +
        nodematch("race", diff = FALSE) +
        nodefactor("race", levels = -1) +
        nodefactor("deg.main", levels = -3) +
        concurrent +
        degrange(from = 4) +
        nodematch("role.class", diff = TRUE, levels = c(1, 2)),
      target.stats = unname(c(
        ns$edges, ns$nodematch_age.grp, ns$nodefactor_age.grp[-5],
        ns$nodematch_race_diffF, ns$nodefactor_race[-1],
        ns$nodefactor_deg.main[-3], ns$concurrent, 0, 0, 0)),
      coef.diss = ns$diss.byage
    )
  } else {
    list(
      formation = ~ edges +
        nodematch("age.grp", diff = FALSE) +
        nodefactor("age.grp", levels = -1) +
        nodematch("race", diff = FALSE) +
        nodefactor("race", levels = -1) +
        nodefactor("risk.grp", levels = -5) +
        nodefactor("deg.tot", levels = -1) +
        nodematch("role.class", diff = TRUE, levels = c(1, 2)),
      target.stats = unname(c(
        ns$edges, sum(ns$nodematch_age.grp), ns$nodefactor_age.grp[-1],
        ns$nodematch_race_diffF, ns$nodefactor_race[-1],
        ns$nodefactor_risk.grp[-5], ns$nodefactor_deg.tot[-1], 0, 0)),
      coef.diss = EpiModel::dissolution_coefs(~ offset(edges), 1)
    )
  }
}

.empty_fit_row <- function() {
  data.frame(fit_ok = NA, converged = NA, lin_dep_warning = NA,
             fit_minutes = NA_real_, netdx_max_abs_z = NA_real_,
             netdx_max_abs_pct = NA_real_, error = NA_character_,
             stringsAsFactors = FALSE)
}

# Fit one layer with MCMLE and run a static netdx on the formation targets.
.fit_one_layer <- function(netstats, layer, nsims_dx = 1000) {
  spec <- .template_layer_spec(netstats, layer)
  nw <- EpiModel::network_initialize(netstats$demog$num)
  nw <- EpiModel::set_vertex_attribute(nw, names(netstats$attr), netstats$attr)
  warns <- character(0)
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(
      EpiModel::netest(nw, formation = spec$formation,
                       target.stats = spec$target.stats,
                       coef.diss = spec$coef.diss,
                       set.control.ergm = .ERGM_CONTROL(),
                       verbose = FALSE),
      error = function(e) {
        structure(list(message = conditionMessage(e)), class = "fit_error")
      }
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  row <- .empty_fit_row()
  row$fit_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  row$lin_dep_warning <- any(grepl("linear", warns, ignore.case = TRUE))
  if (inherits(fit, "fit_error")) {
    row$fit_ok <- FALSE
    row$converged <- FALSE
    row$error <- gsub("\\s+", " ", fit$message)
    return(row)
  }
  row$fit_ok <- TRUE
  row$converged <- !any(grepl("did not converge", warns, ignore.case = TRUE))
  dx <- tryCatch(
    withCallingHandlers(
      EpiModel::netdx(fit, dynamic = FALSE, nsims = nsims_dx,
                      nwstats.formula = "formation",
                      set.control.ergm = ergm::control.simulate.formula(MCMC.burnin = 1e5),
                      verbose = FALSE),
      message = function(m) invokeRestart("muffleMessage"),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  if (!is.null(dx)) {
    st <- dx$stats.table.formation
    z <- st[, grep("^Z", colnames(st))[1]]
    pct <- st[, grep("^Pct", colnames(st))[1]]
    ok <- is.finite(z) & is.finite(pct)
    row$netdx_max_abs_z <- round(max(abs(z[ok])), 2)
    row$netdx_max_abs_pct <- round(max(abs(pct[ok])), 2)
  }
  row
}

# Build netstats for one scenario x arm, check feasibility, fit what is
# feasible. Returns a list(rows = data.frame, mixmat = list or NULL).
.run_one_fit_job <- function(scenario, arm, nsims_dx = 1000) {
  method <- if (arm == "joint") "joint" else "existing"
  ns <- .build_scenario_netstats(scenario, method, floor = arm == "existing_floor")
  feas <- ARTnet::check_race_feasibility(ns)
  rows <- lapply(c("main", "casl", "inst"), function(l) {
    f <- feas[feas$layer == l, , drop = FALSE]
    worst <- which.max(f$slack)
    base <- data.frame(
      scenario = scenario$name, arm = arm, layer = l,
      feasible = all(f$feasible),
      worst_slack = round(max(f$slack), 1),
      worst_race = f$race[worst],
      edges = round(ns[[l]]$edges, 1),
      nodematch_race_diffF = round(ns[[l]]$nodematch_race_diffF, 1),
      stringsAsFactors = FALSE
    )
    if (!all(f$feasible)) {
      r <- .empty_fit_row()
      r$fit_ok <- FALSE
      r$converged <- FALSE
      r$error <- "skipped: infeasible race targets"
    } else {
      r <- .fit_one_layer(ns, l, nsims_dx)
    }
    cbind(base, r)
  })
  mixmat <- NULL
  if (method == "joint") {
    mixmat <- lapply(c(main = "main", casl = "casl", inst = "inst"), function(l) {
      list(reported = ns[[l]]$mixmat_race_reported, balanced = ns[[l]]$mixmat_race)
    })
  }
  list(rows = do.call(rbind, rows), mixmat = mixmat)
}


#' Fit the template ERGMs with MCMLE for every scenario and arm.
#'
#' @param scenarios List of scenario specs (defaults to COMPARISON_SCENARIOS).
#' @param arms Character vector of arms to run; see the header comment.
#' @param ncores Number of cores for parallel::mclapply over scenario x arm
#'   jobs. Each ergm fit itself runs on one core.
#' @param nsims_dx Number of static netdx simulations per fitted layer.
#' @param save_to Optional path for saveRDS() of the result.
#' @return A list with `fits` (one row per scenario x arm x layer), `mixmats`
#'   (the joint arm's reported and balanced race mixing matrices by scenario
#'   and layer), and run metadata.
fit_comparison_ergms <- function(scenarios = COMPARISON_SCENARIOS,
                                 arms = c("existing", "existing_floor", "joint"),
                                 ncores = 1,
                                 nsims_dx = 1000,
                                 save_to = NULL) {
  .require_artnetdata()
  if (!requireNamespace("EpiModel", quietly = TRUE)) {
    stop("EpiModel is required to fit the ERGMs.")
  }
  jobs <- expand.grid(scenario = seq_along(scenarios), arm = arms,
                      stringsAsFactors = FALSE)
  message("Fitting ", nrow(jobs), " scenario x arm jobs (3 layers each) with MCMLE ",
          "on ", ncores, " core(s)...")
  run_job <- function(i) {
    s <- scenarios[[jobs$scenario[i]]]
    arm <- jobs$arm[i]
    tryCatch(
      .run_one_fit_job(s, arm, nsims_dx),
      error = function(e) {
        rows <- data.frame(scenario = s$name, arm = arm,
                           layer = c("main", "casl", "inst"),
                           feasible = NA, worst_slack = NA_real_, worst_race = NA_integer_,
                           edges = NA_real_, nodematch_race_diffF = NA_real_,
                           stringsAsFactors = FALSE)
        r <- .empty_fit_row()[rep(1, 3), ]
        r$fit_ok <- FALSE
        r$error <- paste("job error:", gsub("\\s+", " ", conditionMessage(e)))
        list(rows = cbind(rows, r), mixmat = NULL)
      }
    )
  }
  res <- if (ncores > 1) {
    parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = ncores)
  } else {
    lapply(seq_len(nrow(jobs)), function(i) {
      message("  ", scenarios[[jobs$scenario[i]]]$name, " / ", jobs$arm[i])
      run_job(i)
    })
  }
  fits <- do.call(rbind, lapply(res, `[[`, "rows"))
  rownames(fits) <- NULL
  mixmats <- list()
  for (i in seq_along(res)) {
    if (!is.null(res[[i]]$mixmat)) {
      mixmats[[scenarios[[jobs$scenario[i]]]$name]] <- res[[i]]$mixmat
    }
  }
  out <- list(
    fits = fits,
    mixmats = mixmats,
    network.size = scenarios[[1]]$netstats$network.size,
    seed = .COMPARISON_SEED,
    artnet_version = as.character(packageVersion("ARTnet")),
    ergm_version = as.character(packageVersion("ergm")),
    date = as.character(Sys.Date())
  )
  if (!is.null(save_to)) saveRDS(out, save_to)
  out
}


.fmt_matrix_md <- function(M, digits = 1) {
  if (is.null(M)) return("_not available_")
  labs <- rownames(M)
  if (is.null(labs)) labs <- as.character(seq_len(nrow(M)))
  lines <- c(paste0("| | ", paste(labs, collapse = " | "), " |"),
             paste0("|---|", paste(rep("---:", ncol(M)), collapse = "|"), "|"))
  for (i in seq_len(nrow(M))) {
    lines <- c(lines, paste0("| ", labs[i], " | ",
                             paste(formatC(M[i, ], format = "f", digits = digits),
                                   collapse = " | "), " |"))
  }
  paste(lines, collapse = "\n")
}

#' Render a markdown report of the ERGM fit results.
render_fit_report <- function(fits_result,
                              file = "inst/validation/method_comparison_fits.md") {
  fits <- fits_result$fits
  con <- file(file, "w")
  on.exit(close(con))
  out <- function(...) cat(..., "\n", sep = "", file = con)

  out("# Method comparison: MCMLE fits of the template ERGMs by scenario city")
  out()
  out("Generated by `fit_comparison_ergms()` in `inst/validation/method_comparison.R` on ",
      fits_result$date, ". Addresses issue #12.")
  out()
  out("ARTnet version: ", fits_result$artnet_version, ". ergm version: ",
      fits_result$ergm_version, ". Seed: ", fits_result$seed,
      ". Network size: ", fits_result$network.size, ".")
  out()
  out("Every scenario x arm target set is first run through `check_race_feasibility()`. ",
      "Feasible sets are fit for all three layers with the template formulas and ",
      "`control.ergm(main.method = \"MCMLE\", MCMLE.maxit = 500, SAN.maxit = 3, ",
      "SAN.nsteps.times = 4, MCMC.samplesize = 1e4, MCMC.interval = 5e3)`, then checked ",
      "with a static `netdx` on the formation targets. Infeasible sets are skipped.")
  out()
  out("Arms: `existing` is the legacy per-attribute method; `existing_floor` is the same ",
      "with `nodematch.race.floor = TRUE`; `joint` is the joint g-computation with ",
      "partner-race mixing balanced to the target population.")
  out()
  out("## Summary")
  out()
  n_fit <- sum(fits$fit_ok, na.rm = TRUE)
  n_conv <- sum(fits$converged, na.rm = TRUE)
  n_skip <- sum(!is.na(fits$feasible) & !fits$feasible)
  out("- Target sets (scenario x arm x layer): ", nrow(fits))
  out("- Skipped as infeasible: ", n_skip)
  out("- Fits returned: ", n_fit, "; MCMLE converged: ", n_conv)
  by_arm <- split(fits, fits$arm)
  for (a in names(by_arm)) {
    f <- by_arm[[a]]
    out("- `", a, "`: ", sum(!is.na(f$feasible) & !f$feasible), " infeasible, ",
        sum(f$converged, na.rm = TRUE), " of ", nrow(f), " converged, max netdx |Z| ",
        if (all(is.na(f$netdx_max_abs_z))) "n/a" else sprintf("%.2f", max(f$netdx_max_abs_z, na.rm = TRUE)),
        ", max |% diff| ",
        if (all(is.na(f$netdx_max_abs_pct))) "n/a" else sprintf("%.2f%%", max(f$netdx_max_abs_pct, na.rm = TRUE)))
  }
  out()
  out("## Fits by scenario")
  out()
  tf <- function(x) ifelse(is.na(x), "n/a", ifelse(x, "yes", "no"))
  num <- function(x, d = 1) ifelse(is.na(x), "n/a", formatC(x, format = "f", digits = d))
  for (s in unique(fits$scenario)) {
    out("### ", s)
    out()
    out("| Arm | Layer | Feasible | Worst slack | Edges | Same-race target | Fit | Converged | ",
        "Lin. dep. warning | Minutes | netdx max abs Z | netdx max abs % | Note |")
    out("|---|---|---|---:|---:|---:|---|---|---|---:|---:|---:|---|")
    sub <- fits[fits$scenario == s, , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      out(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                  r$arm, r$layer, tf(r$feasible), num(r$worst_slack), num(r$edges),
                  num(r$nodematch_race_diffF), tf(r$fit_ok), tf(r$converged),
                  tf(r$lin_dep_warning), num(r$fit_minutes, 2), num(r$netdx_max_abs_z, 2),
                  num(r$netdx_max_abs_pct, 2),
                  if (is.na(r$error)) "" else r$error))
    }
    out()
    mm <- fits_result$mixmats[[s]]
    if (!is.null(mm)) {
      out("Joint arm race mixing matrices (rows = ego race). Reported: expected partnership ",
          "endpoints by ego race and partner race at the synthetic population, standardized ",
          "on the ego side only; row sums are `nodefactor_race`. Balanced: expected edge ",
          "counts after symmetric scaling with `nodefactor_race` held fixed; the diagonal is ",
          "`nodematch_race`.")
      out()
      for (l in c("main", "casl", "inst")) {
        out("**", l, ", reported endpoints**")
        out()
        out(.fmt_matrix_md(mm[[l]]$reported))
        out()
        out("**", l, ", balanced edges**")
        out()
        out(.fmt_matrix_md(mm[[l]]$balanced))
        out()
      }
    }
  }
  out("## Reproducibility")
  out()
  out("```r")
  out("source(system.file(\"validation/method_comparison.R\", package = \"ARTnet\"))")
  out("fits <- fit_comparison_ergms(ncores = 4)")
  out("render_fit_report(fits)")
  out("```")
  invisible(file)
}
