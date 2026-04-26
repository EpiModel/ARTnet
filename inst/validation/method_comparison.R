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

# Run one scenario through both methods, return both extracted stat tables.
.run_one_scenario <- function(scenario) {
  do_run <- function(np_method, dur_method) {
    set.seed(.COMPARISON_SEED)
    epistats <- do.call(ARTnet::build_epistats, scenario$epistats)
    set.seed(.COMPARISON_SEED)
    np_args <- c(list(epistats = epistats), scenario$netparams,
                 list(method = np_method, duration.method = dur_method))
    netparams <- do.call(ARTnet::build_netparams, np_args)
    set.seed(.COMPARISON_SEED)
    ns_args <- c(list(epistats = epistats, netparams = netparams),
                 scenario$netstats, list(method = np_method))
    do.call(ARTnet::build_netstats, ns_args)
  }
  list(
    existing = .extract_target_stats(do_run("existing", "empirical"),
                                     scenario$name, "existing"),
    joint    = .extract_target_stats(do_run("joint", "joint_lm"),
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
  out("| Scenario | Description |")
  out("|---|---|")
  out("| `atlanta_default` | Baseline EpiModelHIV-Template config (Atlanta, race = TRUE) |")
  out("| `national_no_geog` | No geographic stratification (sanity check) |")
  out("| `atlanta_nhbs_shifted` | Atlanta config with `race.prop = c(0.35, 0.25, 0.40)` (NHBS-MSM-like) |")
  out("| `atlanta_no_race` | `race = FALSE` path (sanity check) |")
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
