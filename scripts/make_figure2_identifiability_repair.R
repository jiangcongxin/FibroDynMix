#!/usr/bin/env Rscript

# Figure 2: controlled diagnosis and stabilization of state-recovery drift.
#
# This figure deliberately distinguishes three evidence levels:
#   (1) shared-initialization fitting trajectories,
#   (2) simulator-oracle decompositions, and
#   (3) a held-out, initializer-anchored stability repair.
# It does not claim that an oracle prior is available for real data.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "scripts/make_figure2_identifiability_repair.R"
}
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- Sys.getenv(
  "FDM_FIGURE_OUT",
  unset = file.path(ROOT, "figures", "figure2_identifiability_repair")
)
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure2_identifiability_repair")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
QC <- file.path(OUT, "qc")
for (directory in c(EXPORTS, SOURCE, QC)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

source(file.path(ROOT, "scripts", "figure_style_fdm2.R"))

input_path <- function(...) file.path(ROOT, ...)
read_required_tsv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required Figure 2 input is missing: %s", path), call. = FALSE)
  }
  fdm_read_tsv(path)
}
require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(sprintf("%s is missing required columns: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

# Return a transparent mean and 95% t interval. The raw paired observations are
# also exported, so the plotting summaries are not the only available evidence.
mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- if (n > 1L) stats::sd(x) else NA_real_
  half_width <- if (n > 1L) stats::qt(0.975, df = n - 1L) * sd_x / sqrt(n) else NA_real_
  c(n = n, mean = mean_x, sd = sd_x, lower = mean_x - half_width, upper = mean_x + half_width)
}

summarise_ci <- function(data, by, value) {
  require_columns(data, c(by, value), "summary input")
  split_key <- interaction(
    lapply(data[by], as.character),
    drop = TRUE,
    lex.order = TRUE
  )
  pieces <- lapply(split(data, split_key), function(group) {
    stat <- mean_ci(group[[value]])
    out <- group[1L, by, drop = FALSE]
    out$n <- unname(stat[["n"]])
    out$mean <- unname(stat[["mean"]])
    out$sd <- unname(stat[["sd"]])
    out$lower <- unname(stat[["lower"]])
    out$upper <- unname(stat[["upper"]])
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

extract_endpoint_loss <- function(metrics, end_iteration = 20L) {
  key <- c("scenario", "replicate", "method")
  require_columns(metrics, c(key, "iteration", "z_rmse"), "trajectory metrics")
  initial <- metrics[metrics$iteration == 0L, c(key, "z_rmse"), drop = FALSE]
  final <- metrics[metrics$iteration == end_iteration, c(key, "z_rmse"), drop = FALSE]
  if (anyDuplicated(initial[key]) || anyDuplicated(final[key])) {
    stop("Expected exactly one initial and one final trajectory row per matched run.", call. = FALSE)
  }
  names(initial)[names(initial) == "z_rmse"] <- "z_rmse_initial"
  names(final)[names(final) == "z_rmse"] <- "z_rmse_final"
  out <- merge(initial, final, by = key, all = FALSE, sort = FALSE)
  out$rmse_loss <- out$z_rmse_final - out$z_rmse_initial
  out
}

scenario_label <- c(
  continuous = "Continuous",
  discrete = "Discrete",
  batch_confounding = "Batch-confounded",
  rare_transition = "Rare transition"
)
method_label <- c(
  fibrodynmix_nb = "NB",
  fibrodynmix_nb_study = "NB + study",
  fibrodynmix_nb_study_donor = "NB + study + donor"
)

palette <- fdm_editorial_palette()
ink <- unname(palette[["ink"]])
signal <- unname(palette[["signal"]])
signal_deep <- unname(palette[["signal_deep"]])
selected <- unname(palette[["selected"]])
selected_deep <- unname(palette[["selected_deep"]])
neutral <- unname(palette[["neutral"]])
neutral_dark <- unname(palette[["neutral_dark"]])
grid_col <- unname(palette[["grid"]])

# All typography stays black; colours identify fitted conditions or evidence
# classes only. Panel letters are embedded at the start of each title, so their
# publication-scale placement cannot overlap a neighbouring panel.
theme_set(fdm_theme(base_size = 7.4))
black_text_theme <- theme(
  text = element_text(colour = ink),
  plot.title = element_text(colour = ink),
  plot.subtitle = element_text(colour = ink),
  plot.caption = element_text(colour = ink),
  axis.title = element_text(colour = ink),
  axis.text = element_text(colour = ink),
  legend.title = element_text(colour = ink),
  legend.text = element_text(colour = ink),
  strip.text = element_text(colour = ink)
)
tag_panel <- function(plot, tag) {
  plot + theme(
    plot.title = element_text(margin = margin(b = 1)),
    plot.subtitle = element_text(margin = margin(b = 5)),
    plot.margin = margin(t = 12, r = 10, b = 8, l = 10)
  )
}

# ---- Inputs ----------------------------------------------------------------

current_metrics_path <- input_path("analysis", "shared_initialization_trajectory_current_v1", "trajectory_metrics.tsv")
current_objective_path <- input_path("analysis", "shared_initialization_trajectory_current_v1", "objective_components.tsv")
constraint_metrics_path <- input_path("analysis", "shared_initialization_trajectory_constraint_v1", "trajectory_metrics.tsv")
oracle_metrics_path <- input_path("analysis", "oracle_identifiability_diagnostic_v1", "diagnostic_metrics.tsv")
holdout_none_path <- input_path("analysis", "shared_initialization_trajectory_holdout_none_v1", "trajectory_metrics.tsv")
holdout_anchor_path <- input_path("analysis", "shared_initialization_trajectory_holdout_initializer_cell_sd01_v1", "trajectory_metrics.tsv")

current_metrics <- read_required_tsv(current_metrics_path)
current_objective <- read_required_tsv(current_objective_path)
constraint_metrics <- read_required_tsv(constraint_metrics_path)
oracle_metrics <- read_required_tsv(oracle_metrics_path)
holdout_none_metrics <- read_required_tsv(holdout_none_path)
holdout_anchor_metrics <- read_required_tsv(holdout_anchor_path)

require_columns(current_objective,
  c("scenario", "replicate", "method", "iteration", "complete_objective_average"),
  "current objective components"
)

# ---- Panel A: shared-initialization trajectories --------------------------

trajectory_key <- c("scenario", "replicate", "method", "iteration")
traj <- merge(
  current_metrics[, c(trajectory_key, "z_rmse"), drop = FALSE],
  current_objective[, c(trajectory_key, "complete_objective_average"), drop = FALSE],
  by = trajectory_key,
  all = FALSE,
  sort = FALSE
)
traj_initial <- traj[traj$iteration == 0L,
  c("scenario", "replicate", "method", "z_rmse", "complete_objective_average"), drop = FALSE]
names(traj_initial)[names(traj_initial) == "z_rmse"] <- "z_rmse_initial"
names(traj_initial)[names(traj_initial) == "complete_objective_average"] <- "objective_initial"
traj <- merge(traj, traj_initial, by = c("scenario", "replicate", "method"), all.x = TRUE, sort = FALSE)
traj$objective_reduction_pct <- 100 * (traj$objective_initial - traj$complete_objective_average) / traj$objective_initial
traj$rmse_increase_pct <- 100 * (traj$z_rmse - traj$z_rmse_initial) / traj$z_rmse_initial
traj_summary <- summarise_ci(traj, c("scenario", "method", "iteration"), "objective_reduction_pct")
names(traj_summary)[names(traj_summary) %in% c("n", "mean", "sd", "lower", "upper")] <-
  paste0("objective_", names(traj_summary)[names(traj_summary) %in% c("n", "mean", "sd", "lower", "upper")])
rmse_summary <- summarise_ci(traj, c("scenario", "method", "iteration"), "rmse_increase_pct")
names(rmse_summary)[names(rmse_summary) %in% c("n", "mean", "sd", "lower", "upper")] <-
  paste0("rmse_", names(rmse_summary)[names(rmse_summary) %in% c("n", "mean", "sd", "lower", "upper")])
traj_plot <- merge(traj_summary, rmse_summary, by = c("scenario", "method", "iteration"), all = TRUE, sort = FALSE)
traj_plot$scenario_label <- factor(scenario_label[traj_plot$scenario], levels = unname(scenario_label))
traj_plot$method_label <- factor(method_label[traj_plot$method], levels = unname(method_label[c("fibrodynmix_nb", "fibrodynmix_nb_study")]))
traj_plot <- traj_plot[order(traj_plot$scenario_label, traj_plot$method_label, traj_plot$iteration), , drop = FALSE]
fdm_write_tsv(traj_plot, file.path(SOURCE, "fig2a_shared_initialization_phase_paths.tsv"))

endpoint_points <- traj_plot[traj_plot$iteration %in% c(0L, 20L), , drop = FALSE]
panel_a <- ggplot(
  traj_plot,
  aes(x = objective_mean, y = rmse_mean, colour = method_label, linetype = method_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = neutral_dark) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = neutral_dark) +
  geom_path(linewidth = 0.72, alpha = 0.95, arrow = arrow(length = grid::unit(0.055, "in"), type = "closed", ends = "last")) +
  geom_point(data = endpoint_points[endpoint_points$iteration == 0L, , drop = FALSE],
    shape = 21, fill = "white", size = 1.8, stroke = 0.5
  ) +
  geom_point(data = endpoint_points[endpoint_points$iteration == 20L, , drop = FALSE],
    shape = 16, size = 1.95
  ) +
  facet_wrap(~ scenario_label, nrow = 1) +
  scale_colour_manual(values = c("NB" = neutral_dark, "NB + study" = selected_deep), name = NULL) +
  scale_linetype_manual(values = c("NB" = "solid", "NB + study" = "dashed"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "A. Shared-start trajectories: count fit versus state recovery",
    subtitle = "Updates 0 to 20; open = common initializer and filled = final update\nMean of 5 matched simulations per scenario and method",
    x = "Reduction in complete penalized objective (%)",
    y = "Increase in state-weight RMSE (%)"
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.margin = margin(b = 1),
    panel.grid.major = element_line(linewidth = 0.22, colour = grid_col),
    strip.text = element_text(size = 7.1, face = "bold")
  ) +
  black_text_theme
panel_a <- tag_panel(panel_a, "A")

# ---- Panel B: oracle decomposition ----------------------------------------

oracle_order <- c(
  "initializer",
  "oracle_z_mle_all_truth",
  "oracle_z_map_all_truth",
  "oracle_z_posterior_mean_all_truth",
  "fixed_beta_estimated_nuisance",
  "joint_current_nb"
)
oracle_labels <- c(
  initializer = "Marker initializer",
  oracle_z_mle_all_truth = "Conditional z MLE\n(true count parameters)",
  oracle_z_map_all_truth = "Exact logistic-normal MAP [oracle]",
  oracle_z_posterior_mean_all_truth = "Posterior mean [oracle]",
  fixed_beta_estimated_nuisance = "True beta; estimated nuisance",
  joint_current_nb = "Joint current NB fit"
)
oracle_class <- c(
  initializer = "Initializer",
  oracle_z_mle_all_truth = "Conditional likelihood only",
  oracle_z_map_all_truth = "Oracle prior",
  oracle_z_posterior_mean_all_truth = "Oracle prior",
  fixed_beta_estimated_nuisance = "Incomplete nuisance model",
  joint_current_nb = "Joint current model"
)
require_columns(oracle_metrics, c("scenario", "replicate", "method", "rmse"), "oracle diagnostic metrics")
oracle_metrics <- oracle_metrics[oracle_metrics$method %in% oracle_order & oracle_metrics$status == "ok", , drop = FALSE]
oracle_summary <- summarise_ci(oracle_metrics, c("scenario", "method"), "rmse")
oracle_summary$method_label <- factor(oracle_labels[oracle_summary$method], levels = rev(unname(oracle_labels[oracle_order])))
oracle_summary$evidence_class <- unname(oracle_class[oracle_summary$method])
oracle_summary$scenario_label <- factor(
  scenario_label[oracle_summary$scenario],
  levels = unname(scenario_label[c("continuous", "batch_confounding")])
)
oracle_summary$y <- as.numeric(oracle_summary$method_label)
oracle_summary <- oracle_summary[order(oracle_summary$scenario_label, oracle_summary$y), , drop = FALSE]
oracle_export <- oracle_summary
oracle_export$method_label <- gsub("\n", "; ", as.character(oracle_export$method_label), fixed = TRUE)
fdm_write_tsv(oracle_export, file.path(SOURCE, "fig2b_oracle_decomposition.tsv"))

oracle_colours <- c(
  "Initializer" = neutral_dark,
  "Conditional likelihood only" = signal,
  "Oracle prior" = selected_deep,
  "Incomplete nuisance model" = signal_deep,
  "Joint current model" = ink
)
panel_b <- ggplot(oracle_summary, aes(x = mean, y = y, colour = evidence_class)) +
  geom_segment(aes(x = lower, xend = upper, yend = y), linewidth = 0.72, colour = ink) +
  geom_point(size = 2.15) +
  facet_wrap(~ scenario_label, nrow = 1) +
  scale_colour_manual(values = oracle_colours, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(rev(unname(oracle_labels[oracle_order]))),
    labels = rev(unname(oracle_labels[oracle_order])),
    expand = expansion(add = c(0.45, 0.45))
  ) +
  scale_x_continuous(limits = c(0, 0.285), breaks = seq(0, 0.25, by = 0.05), expand = c(0, 0)) +
  labs(
    title = "B. Oracle decomposition: conditional likelihood versus simulator prior",
    subtitle = "State-weight RMSE; mean and 95% t interval (n = 5)\n[oracle] uses the known simulated donor logistic-normal prior and is not a real-data result",
    x = "State-weight RMSE",
    y = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.22, colour = grid_col),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 7.1, face = "bold")
  ) +
  black_text_theme
panel_b <- tag_panel(panel_b, "B")

# ---- Panel C: sum-to-zero constraint negative control ---------------------

current_end <- extract_endpoint_loss(current_metrics)
constraint_end <- extract_endpoint_loss(constraint_metrics)
constraint_pair <- merge(
  current_end[, c("scenario", "replicate", "method", "rmse_loss"), drop = FALSE],
  constraint_end[, c("scenario", "replicate", "method", "rmse_loss"), drop = FALSE],
  by = c("scenario", "replicate", "method"),
  suffixes = c("_unconstrained", "_sum_to_zero"),
  all = FALSE,
  sort = FALSE
)
constraint_pair$constraint_effect <- constraint_pair$rmse_loss_sum_to_zero - constraint_pair$rmse_loss_unconstrained
constraint_pair$constraint_effect_x1000 <- 1000 * constraint_pair$constraint_effect
constraint_summary <- summarise_ci(constraint_pair, c("scenario", "method"), "constraint_effect_x1000")
constraint_summary$condition <- paste(
  scenario_label[constraint_summary$scenario],
  method_label[constraint_summary$method],
  sep = " - "
)
constraint_condition_order <- c(
  outer(unname(scenario_label), unname(method_label[c("fibrodynmix_nb", "fibrodynmix_nb_study")]), paste, sep = " - ")
)
constraint_summary$condition <- factor(constraint_summary$condition, levels = rev(constraint_condition_order))
constraint_summary$y <- as.numeric(constraint_summary$condition)
constraint_summary <- constraint_summary[order(constraint_summary$y), , drop = FALSE]
fdm_write_tsv(constraint_pair, file.path(SOURCE, "fig2c_sum_to_zero_paired_observations.tsv"))
fdm_write_tsv(constraint_summary, file.path(SOURCE, "fig2c_sum_to_zero_constraint_effect.tsv"))

constraint_limit <- max(0.12, max(abs(c(constraint_summary$lower, constraint_summary$upper)), na.rm = TRUE) * 1.18)
panel_c <- ggplot(constraint_summary, aes(x = mean, y = y)) +
  geom_vline(xintercept = 0, linewidth = 0.36, colour = ink) +
  geom_segment(aes(x = lower, xend = upper, yend = y), linewidth = 0.72, colour = ink) +
  geom_point(size = 2.05, colour = signal) +
  scale_y_continuous(
    breaks = seq_along(levels(constraint_summary$condition)),
    labels = levels(constraint_summary$condition),
    expand = expansion(add = c(0.42, 0.42))
  ) +
  scale_x_continuous(
    limits = c(-constraint_limit, constraint_limit),
    labels = function(x) sprintf("%+.2f", x),
    expand = c(0, 0)
  ) +
  labs(
    title = "C. Sum-to-zero constraint\nalone is insufficient",
    subtitle = "At update 20: constrained minus unconstrained\nMatched starts and seeds; mean and 95% t interval",
    x = "Change in RMSE loss\n(x10^-3)",
    y = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.22, colour = grid_col),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.2)
  ) +
  black_text_theme
panel_c <- tag_panel(panel_c, "C")

# ---- Panel D: held-out initializer-cell anchor ----------------------------

holdout_none_end <- extract_endpoint_loss(holdout_none_metrics)
holdout_anchor_end <- extract_endpoint_loss(holdout_anchor_metrics)
anchor_pair <- merge(
  holdout_none_end[, c("scenario", "replicate", "method", "rmse_loss"), drop = FALSE],
  holdout_anchor_end[, c("scenario", "replicate", "method", "rmse_loss"), drop = FALSE],
  by = c("scenario", "replicate", "method"),
  suffixes = c("_none", "_initializer_anchor"),
  all = FALSE,
  sort = FALSE
)
anchor_long <- rbind(
  data.frame(anchor_pair[, c("scenario", "replicate", "method"), drop = FALSE], mode = "No anchor", rmse_loss = anchor_pair$rmse_loss_none),
  data.frame(anchor_pair[, c("scenario", "replicate", "method"), drop = FALSE], mode = "Initializer-cell anchor", rmse_loss = anchor_pair$rmse_loss_initializer_anchor)
)
anchor_summary <- summarise_ci(anchor_long, c("scenario", "method", "mode"), "rmse_loss")
anchor_summary$condition <- paste(
  scenario_label[anchor_summary$scenario],
  method_label[anchor_summary$method],
  sep = " - "
)
anchor_condition_order <- c(
  outer(
    unname(scenario_label[c("continuous", "batch_confounding")]),
    unname(method_label[c("fibrodynmix_nb", "fibrodynmix_nb_study", "fibrodynmix_nb_study_donor")]),
    paste,
    sep = " - "
  )
)
anchor_summary$condition <- factor(anchor_summary$condition, levels = rev(anchor_condition_order))
anchor_summary$mode <- factor(anchor_summary$mode, levels = c("No anchor", "Initializer-cell anchor"))
anchor_summary$y <- as.numeric(anchor_summary$condition)
anchor_summary <- anchor_summary[order(anchor_summary$y, anchor_summary$mode), , drop = FALSE]

anchor_none_summary <- anchor_summary[anchor_summary$mode == "No anchor", c("condition", "y", "mean"), drop = FALSE]
names(anchor_none_summary)[names(anchor_none_summary) == "mean"] <- "mean_no_anchor"
anchor_cell_summary <- anchor_summary[anchor_summary$mode == "Initializer-cell anchor", c("condition", "y", "mean"), drop = FALSE]
names(anchor_cell_summary)[names(anchor_cell_summary) == "mean"] <- "mean_initializer_cell_anchor"
anchor_segments <- merge(
  anchor_none_summary,
  anchor_cell_summary,
  by = c("condition", "y"),
  all = TRUE,
  sort = FALSE
)
fdm_write_tsv(anchor_pair, file.path(SOURCE, "fig2d_heldout_anchor_paired_observations.tsv"))
fdm_write_tsv(anchor_summary, file.path(SOURCE, "fig2d_heldout_anchor_repair.tsv"))

panel_d <- ggplot(anchor_summary, aes(x = mean, y = y, colour = mode, shape = mode)) +
  geom_vline(xintercept = 0, linewidth = 0.36, colour = ink) +
  geom_segment(
    data = anchor_segments,
    aes(x = mean_no_anchor, xend = mean_initializer_cell_anchor, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.72,
    colour = neutral
  ) +
  geom_segment(aes(x = lower, xend = upper, yend = y), linewidth = 0.72) +
  geom_point(size = 2.2, fill = "white", stroke = 0.72) +
  scale_colour_manual(values = c("No anchor" = signal, "Initializer-cell anchor" = selected_deep), labels = c("No anchor", "Cell anchor"), name = NULL) +
  scale_shape_manual(values = c("No anchor" = 21, "Initializer-cell anchor" = 24), labels = c("No anchor", "Cell anchor"), name = NULL) +
  scale_y_continuous(
    breaks = seq_along(levels(anchor_summary$condition)),
    labels = levels(anchor_summary$condition),
    expand = expansion(add = c(0.42, 0.42))
  ) +
  scale_x_continuous(limits = c(-0.016, 0.125), breaks = seq(0, 0.12, by = 0.03), expand = c(0, 0)) +
  labs(
    title = "D. Cell anchor limits\nRMSE drift",
    subtitle = "20 updates; n = 5; 95% t interval\nAnchor SD = 0.10 (reference logits)",
    x = "Increase in state-weight RMSE\nfrom shared initializer",
    y = NULL
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.x = element_line(linewidth = 0.22, colour = grid_col),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.2)
  ) +
  black_text_theme
panel_d <- tag_panel(panel_d, "D")

# ---- Assemble and export ----------------------------------------------------

figure2 <- panel_a / panel_b / (panel_c | panel_d) +
  plot_layout(heights = c(1.06, 1.10, 0.98), widths = c(1, 1)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

fdm_save_plot(figure2, EXPORTS, STEM, width = fdm_main_width(), height = 8.8, dpi = 600)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = "Figure 2",
    primary_claim = "Shared-start trajectories show that a lower fitted count objective can coincide with worse state recovery; oracle and held-out anchor analyses distinguish conditional likelihood, a simulator-specific prior, and cell-coordinate stabilization.",
    claim_boundary = "The exact logistic-normal prior is an explicit simulator oracle. The held-out initializer-cell anchor demonstrates optimization stability relative to the shared marker initializer, not external biological validation, a universal prior choice, or improved recovery beyond initialization.",
    current_trajectory_replicates = 5L,
    oracle_replicates = 5L,
    heldout_replicates = 5L,
    anchor_reference_logit_sd = 0.10,
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C", "D"),
    source_data = c(
      "source_data/fig2a_shared_initialization_phase_paths.tsv",
      "source_data/fig2b_oracle_decomposition.tsv",
      "source_data/fig2c_sum_to_zero_constraint_effect.tsv; source_data/fig2c_sum_to_zero_paired_observations.tsv",
      "source_data/fig2d_heldout_anchor_repair.tsv; source_data/fig2d_heldout_anchor_paired_observations.tsv"
    ),
    claim = c(
      "From a common marker-guided start, the complete penalized objective declines while known-truth state-weight RMSE increases.",
      "In five-replicate means, conditional likelihood optimization does not improve RMSE relative to the marker initializer; the exact simulator logistic-normal prior substantially reduces mean RMSE.",
      "Removing the alpha-beta intercept/program gauge with a sum-to-zero constraint alone has a negligible effect on later RMSE drift.",
      "On a held-out seed set, a cell-level anchor to the marker initializer keeps later RMSE close to its starting value."
    ),
    evidence_boundary = c(
      "Controlled simulation trajectory; no causal mechanism is inferred from this panel alone.",
      "Oracle comparison: exact prior is known only because data were simulated.",
      "Matched negative control; it tests one exact gauge constraint, not all forms of factorization ambiguity.",
      "Optimization-stability result on simulated data; not a claim of external biological validation."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    source_file = c(
      current_metrics_path,
      current_objective_path,
      constraint_metrics_path,
      oracle_metrics_path,
      holdout_none_path,
      holdout_anchor_path
    ),
    role = c(
      "Panel A shared-start state-recovery trajectories",
      "Panel A complete penalized objective trajectories",
      "Panel C sum-to-zero matched negative control",
      "Panel B conditional-likelihood and oracle decomposition",
      "Panel D held-out trajectories without the cell-specific marker-logit anchor",
      "Panel D held-out initializer-cell anchor trajectories"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "analysis_input_manifest.tsv")
)

writeLines(c(
  "# Figure 2. Shared-start trajectories and cell-coordinate stabilization",
  "",
  "**A.** In five matched simulations per scenario, FibroDynMix NB and NB plus study-effect fits begin from the same marker-guided initializer. Across 20 outer updates, the complete penalized count objective decreases while known-truth state-weight RMSE increases. Paths are plotted in objective-reduction and RMSE-increase coordinates; open and filled points mark updates 0 and 20, respectively.",
  "",
  "**B.** Oracle decomposition of two generative scenarios. In the five-replicate mean, conditional per-cell NB maximum likelihood, even with all count parameters fixed to their simulated truth, does not improve RMSE relative to the marker initializer. MAP and posterior-mean estimates under the exact donor-specific logistic-normal prior substantially reduce mean RMSE. These two oracle conditions use information unavailable in real data and are labeled accordingly.",
  "",
  "**C.** A sum-to-zero constraint removes the exact alpha-beta intercept/program gauge, but has a negligible matched-seed effect on the RMSE loss accumulated after 20 updates.",
  "",
  "**D.** On a held-out set of matched simulation seeds, a cell-level reference-logit anchor centered on the common marker initializer (SD 0.10) keeps the later state-recovery error near its starting value across NB, NB plus study, and NB plus study-plus-donor fits. This panel demonstrates optimization stability in simulation; it is not an external biological-validation result or an oracle comparison.",
  "",
  "All intervals are 95% t intervals across five simulations. State-weight RMSE is evaluated against the known simulated mixture."
), file.path(OUT, "main_figure_legends.md"))

writeLines(c(
  "# Reproducibility note",
  "",
  "Run `Rscript scripts/make_figure2_identifiability_repair.R` from the repository root.",
  "",
  "The script only reads the six analysis TSV inputs listed in `analysis_input_manifest.tsv` and writes this output directory. Raw matched observations and plotted summaries are preserved under `source_data/`."
), file.path(OUT, "README.md"))

message(sprintf("Figure 2 coordinate-stabilization package written to %s", OUT))
