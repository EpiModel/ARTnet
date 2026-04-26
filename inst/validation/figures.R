## Render the three figures referenced from method_refactor_report.md.
## Run from the ARTnet repo root via:
##   source("inst/validation/figures.R")
##   render_all_figures()
##
## Outputs: inst/validation/figures/{fig1_joint_vs_existing.png,
##                                    fig2_dissolution_durations.png,
##                                    fig3_race_composition.png}
##
## PNGs are .gitignored — regenerate locally as needed.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

.fig_dir <- function() {
  d <- file.path(getwd(), "inst", "validation", "figures")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}


# ---- Figure 1: joint vs existing scatter (all 363 cells) ---------------------

render_fig1 <- function(comparison) {
  d <- comparison
  d <- d[is.finite(d$existing) & is.finite(d$joint) &
         d$existing > 0 & d$joint > 0, , drop = FALSE]
  d$layer <- factor(d$layer, levels = c("main", "casl", "inst"),
                    labels = c("Main", "Casual", "One-time"))

  p <- ggplot(d, aes(x = existing, y = joint, colour = layer)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60",
                linetype = "dashed", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.6) +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    scale_colour_manual(values = c(Main = "#1f77b4", Casual = "#ff7f0e",
                                   `One-time` = "#2ca02c"), name = "Layer") +
    facet_wrap(~ scenario, ncol = 2) +
    labs(
      x = "Target stat under method = \"existing\"",
      y = "Target stat under method = \"joint\"",
      title = "Joint vs existing target statistics",
      subtitle = paste0("Each point is one (scenario, layer, stat, level) cell. ",
                        "Dashed line is y = x (no correction).")
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title.position = "plot")

  ggsave(file.path(.fig_dir(), "fig1_joint_vs_existing.png"),
         plot = p, width = 8.0, height = 7.5, dpi = 300,
         bg = "white")
  invisible(p)
}


# ---- Figure 2: per-stratum dissolution duration, main layer -----------------

render_fig2 <- function() {
  set.seed(20260420L)
  ep <- ARTnet::build_epistats(
    geog.lvl = "city", geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE, time.unit = 7
  )
  set.seed(20260420L)
  np_e <- ARTnet::build_netparams(ep, smooth.main.dur = TRUE,
                                  duration.method = "empirical")
  set.seed(20260420L)
  np_j <- ARTnet::build_netparams(ep, smooth.main.dur = TRUE,
                                  duration.method = "joint_lm")
  set.seed(20260420L)
  np_e_ns <- ARTnet::build_netparams(ep, smooth.main.dur = FALSE,
                                     duration.method = "empirical")

  strata <- c("nonmatch", paste0("matched.", 1:5))
  rows <- rbind(
    data.frame(stratum = strata,
               mean_dur_adj = np_e_ns$main$durs.main.byage$mean.dur.adj,
               method = "empirical (no smoothing)"),
    data.frame(stratum = strata,
               mean_dur_adj = np_e$main$durs.main.byage$mean.dur.adj,
               method = "empirical (default = smoothed)"),
    data.frame(stratum = strata,
               mean_dur_adj = np_j$main$durs.main.byage$mean.dur.adj,
               method = "joint_lm (smoothed)")
  )
  rows$stratum <- factor(rows$stratum, levels = strata)
  rows$method <- factor(rows$method,
                        levels = c("empirical (no smoothing)",
                                   "empirical (default = smoothed)",
                                   "joint_lm (smoothed)"))

  p <- ggplot(rows, aes(x = stratum, y = mean_dur_adj,
                        colour = method, group = method, shape = method)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("empirical (no smoothing)" = "#7f7f7f",
                                   "empirical (default = smoothed)" = "#1f77b4",
                                   "joint_lm (smoothed)" = "#d62728")) +
    scale_shape_manual(values = c(1, 16, 17)) +
    labs(
      x = "Age-match stratum (main layer)",
      y = "mean.dur.adj (weeks; passed to dissolution_coefs)",
      title = "Stratum-level dissolution duration: empirical vs joint_lm",
      subtitle = "Atlanta, race = TRUE. matched.5 is where the methods diverge most."
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.title.position = "plot")

  ggsave(file.path(.fig_dir(), "fig2_dissolution_durations.png"),
         plot = p, width = 8.0, height = 5.0, dpi = 300,
         bg = "white")
  invisible(p)
}


# ---- Figure 3: ARTnet vs Atlanta race composition ---------------------------

render_fig3 <- function() {
  d <- data.frame(
    race   = factor(rep(c("Black", "Hispanic", "White / Other"), 2),
                    levels = c("Black", "Hispanic", "White / Other")),
    source = factor(rep(c("ARTnet sample", "Atlanta MSM target"), each = 3),
                    levels = c("ARTnet sample", "Atlanta MSM target")),
    prop   = c(0.055, 0.138, 0.807,
               0.515, 0.046, 0.439)
  )

  p <- ggplot(d, aes(x = race, y = prop, fill = source)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.55) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * prop)),
              position = position_dodge(width = 0.65),
              vjust = -0.4, size = 3.4) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values = c("ARTnet sample"      = "#1f77b4",
                                 "Atlanta MSM target" = "#d62728"),
                      name = NULL) +
    labs(
      x = NULL, y = "Population share",
      title = "Race composition: ARTnet sample vs Atlanta MSM target population",
      subtitle = paste0("ARTnet sample is overwhelmingly White / Other; Atlanta's actual MSM ",
                        "population is roughly half Black.")
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title.position = "plot")

  ggsave(file.path(.fig_dir(), "fig3_race_composition.png"),
         plot = p, width = 7.5, height = 4.3, dpi = 300,
         bg = "white")
  invisible(p)
}


# ---- Convenience: render all three -----------------------------------------

render_all_figures <- function() {
  source(file.path(getwd(), "inst", "validation", "method_comparison.R"))
  message("Running compare_methods() across 4 scenarios... (~30s)")
  res <- compare_methods()
  message("Rendering figure 1 (joint vs existing scatter)")
  render_fig1(res)
  message("Rendering figure 2 (dissolution durations)")
  render_fig2()
  message("Rendering figure 3 (race composition)")
  render_fig3()
  message("Done. PNGs at ", .fig_dir())
  invisible(NULL)
}
