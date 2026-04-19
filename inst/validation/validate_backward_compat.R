## Backward-compatibility validation harness for the joint g-computation
## refactor (issues #61-#65). See inst/validation/README.md for the workflow.
##
## Two entry points:
##   capture_snapshot()         - run on pre-refactor `main` to save reference
##   compare_to_snapshot(...)   - run on refactor branch to diff against reference
##
## Both functions iterate over the parameter sets in PARAM_SETS (edit below to
## add coverage). Each set defines the args to build_epistats / build_netparams
## / build_netstats. Keep them small — this is a regression harness, not a
## simulation.

# ---- Configuration -----------------------------------------------------------

# Where snapshots are stored, relative to the package install (or the repo
# root when using devtools::load_all()).
.snapshot_dir <- function() {
  candidates <- c(
    # Dev mode: ARTnet repo root / inst / validation / snapshots
    file.path(getwd(), "inst", "validation", "snapshots"),
    # Installed package
    system.file("validation", "snapshots", package = "ARTnet")
  )
  hit <- candidates[nzchar(candidates) & dir.exists(dirname(candidates))]
  if (length(hit) == 0) {
    stop("Cannot locate inst/validation/. Run from the ARTnet repo root or ",
         "install the package.")
  }
  dir <- hit[1]
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

# Parameter sets to cover. Add more as edge cases surface.
# Each entry: a list with $name (snapshot key), $epistats (args to build_epistats),
# $netparams (args to build_netparams), $netstats (args to build_netstats).
PARAM_SETS <- list(
  list(
    name = "atlanta_default",
    epistats = list(
      geog.lvl = "city",
      geog.cat = "Atlanta",
      init.hiv.prev = c(0.33, 0.137, 0.084),
      race = TRUE,
      time.unit = 7
    ),
    netparams = list(smooth.main.dur = TRUE),
    netstats = list(expect.mort = 0.000478213, network.size = 5000)
  ),
  list(
    name = "national_no_geog",
    epistats = list(race = TRUE, time.unit = 7),
    netparams = list(smooth.main.dur = TRUE),
    netstats = list(expect.mort = 0.000478213, network.size = 5000)
  ),
  list(
    name = "atlanta_no_race",
    epistats = list(
      geog.lvl = "city",
      geog.cat = "Atlanta",
      init.hiv.prev = c(0.33, 0.137, 0.084),
      race = FALSE,
      time.unit = 7
    ),
    netparams = list(smooth.main.dur = TRUE),
    netstats = list(expect.mort = 0.000478213, network.size = 5000)
  )
)

# Fixed seed so the stochastic bits of build_netstats (sample/rbinom/runif)
# are reproducible across runs.
.VALIDATION_SEED <- 20260419L


# ---- Utilities ---------------------------------------------------------------

.require_artnetdata <- function() {
  if (system.file(package = "ARTnetData") == "") {
    stop("ARTnetData not installed; validation cannot run. ",
         "See https://github.com/EpiModel/ARTnet#artnetdata-dependency")
  }
}

# Strip fields that are new/additive under the refactor so they don't cause
# spurious diffs against a pre-refactor snapshot. Extend this list as new
# fields are added (e.g., $joint_model).
.strip_additive <- function(netparams) {
  for (layer in c("main", "casl", "inst", "all")) {
    if (is.null(netparams[[layer]])) next
    netparams[[layer]]$joint_model <- NULL
  }
  netparams
}

# Run one parameter set end-to-end. `netparams_extra` and `netstats_extra`
# let the caller pass e.g. method = "existing" post-refactor without
# affecting the pre-capture. The convenience `method` arg threads a single
# value to both functions (the common case).
.run_one <- function(set, netparams_extra = list(), netstats_extra = list()) {
  set.seed(.VALIDATION_SEED)
  epistats <- do.call(ARTnet::build_epistats, set$epistats)

  netparams_args <- c(list(epistats = epistats), set$netparams, netparams_extra)
  netparams <- do.call(ARTnet::build_netparams, netparams_args)

  netstats_args <- c(list(epistats = epistats, netparams = netparams),
                     set$netstats, netstats_extra)
  netstats <- do.call(ARTnet::build_netstats, netstats_args)

  list(netparams = netparams, netstats = netstats)
}


# ---- Public entry points -----------------------------------------------------

#' Capture golden-reference snapshots from the current code
#'
#' Call once on the pre-refactor `main` branch. Saves one `.rds` per entry in
#' `PARAM_SETS` under `inst/validation/snapshots/`.
#'
#' @param overwrite If TRUE, overwrite existing snapshot files.
capture_snapshot <- function(overwrite = FALSE) {
  .require_artnetdata()
  dir <- .snapshot_dir()

  for (set in PARAM_SETS) {
    path <- file.path(dir, paste0(set$name, ".rds"))
    if (file.exists(path) && !overwrite) {
      message("SKIP (exists): ", path, "  -- pass overwrite = TRUE to replace")
      next
    }
    message("CAPTURE: ", set$name, " -> ", path)
    result <- .run_one(set)
    saveRDS(result, path)
  }
  invisible(TRUE)
}


#' Compare current code output against the captured snapshots
#'
#' Call on the refactor branch. Reports per-parameter-set diffs between the
#' current code's output and the snapshot saved by `capture_snapshot()`. The
#' joint g-comp refactor should pass this with zero diffs when the legacy
#' code path is selected (e.g. `method = "existing"`).
#'
#' @param method If non-NULL, a single value passed as `method = <val>` to
#'   **both** `build_netparams()` and `build_netstats()`. This is the common
#'   case (symmetric method) — pass `method = "existing"` to compare legacy
#'   behavior against snapshots. For asymmetric testing, use the low-level
#'   `netparams_extra` / `netstats_extra` args.
#' @param netparams_extra,netstats_extra Named lists of extra args forwarded
#'   to the respective functions.
#' @param tolerance Numeric tolerance for `all.equal()`. Default 0 (exact).
#' @return Invisibly: TRUE iff all sets match; FALSE otherwise.
compare_to_snapshot <- function(method = NULL,
                                netparams_extra = list(),
                                netstats_extra = list(),
                                tolerance = 0) {
  .require_artnetdata()
  dir <- .snapshot_dir()
  if (!is.null(method)) {
    netparams_extra$method <- method
    netstats_extra$method  <- method
  }

  overall_ok <- TRUE
  for (set in PARAM_SETS) {
    path <- file.path(dir, paste0(set$name, ".rds"))
    if (!file.exists(path)) {
      warning("No snapshot for ", set$name, " at ", path,
              " -- did you forget to run capture_snapshot()?")
      overall_ok <- FALSE
      next
    }
    message("COMPARE: ", set$name)
    ref <- readRDS(path)
    cur <- .run_one(set, netparams_extra = netparams_extra,
                    netstats_extra = netstats_extra)

    np_ref <- .strip_additive(ref$netparams)
    np_cur <- .strip_additive(cur$netparams)
    np_diff <- all.equal(np_ref, np_cur, tolerance = tolerance)
    ns_diff <- all.equal(ref$netstats, cur$netstats, tolerance = tolerance)

    np_ok <- isTRUE(np_diff)
    ns_ok <- isTRUE(ns_diff)
    if (np_ok && ns_ok) {
      message("  OK  (netparams + netstats identical)")
    } else {
      overall_ok <- FALSE
      if (!np_ok) {
        message("  FAIL netparams:")
        message(paste("   ", np_diff, collapse = "\n"))
      }
      if (!ns_ok) {
        message("  FAIL netstats:")
        message(paste("   ", ns_diff, collapse = "\n"))
      }
    }
  }

  if (overall_ok) {
    message("\n==============================")
    message("ALL MATCH (", length(PARAM_SETS), " parameter sets)")
    message("==============================")
  } else {
    message("\n==============================")
    message("REGRESSION DETECTED -- see diffs above")
    message("==============================")
  }
  invisible(overall_ok)
}


#' Show which snapshots currently exist on disk.
list_snapshots <- function() {
  dir <- .snapshot_dir()
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    message("No snapshots in ", dir)
    return(invisible(character(0)))
  }
  info <- file.info(files)
  out <- data.frame(
    name = basename(files),
    size_kb = round(info$size / 1024, 1),
    mtime = info$mtime,
    row.names = NULL
  )
  print(out)
  invisible(files)
}
