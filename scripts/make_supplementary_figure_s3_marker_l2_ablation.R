#!/usr/bin/env Rscript

# Supplementary Figure S3: program-level marker penalties versus the
# cell-specific state coordinate.  The figure deliberately separates fixed
# penalty-strength sensitivity from a likelihood-normalized control.  It never
# compares penalized objectives across different penalty definitions.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "scripts/make_supplementary_figure_s3_marker_l2_ablation.R"
}
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
OUT <- Sys.getenv(
  "FDM_SUPPLEMENTARY_FIGURE_OUT",
  unset = file.path(ROOT, "figures", "figure_s3_marker_l2_ablation")
)
STEM <- Sys.getenv("FDM_SUPPLEMENTARY_FIGURE_STEM", unset = "Figure_S3")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
QC <- file.path(OUT, "qc")
for (directory in c(EXPORTS, SOURCE, QC)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

source(file.path(ROOT, "scripts", "figure_style_fdm2.R"))

read_required_tsv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required Figure S3 input is missing: %s", path), call. = FALSE)
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

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(c(n = 0, mean = NA_real_, sd = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  mean_x <- mean(x)
  sd_x <- if (n > 1L) stats::sd(x) else NA_real_
  half_width <- if (n > 1L) stats::qt(0.975, df = n - 1L) * sd_x / sqrt(n) else NA_real_
  c(n = n, mean = mean_x, sd = sd_x, lower = mean_x - half_width, upper = mean_x + half_width)
}

summarise_ci <- function(data, by, value) {
  require_columns(data, c(by, value), "summary input")
  key <- interaction(lapply(data[by], as.character), drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, key), function(group) {
    statistic <- mean_ci(group[[value]])
    out <- group[1L, by, drop = FALSE]
    out$n <- unname(statistic[["n"]])
    out$mean <- unname(statistic[["mean"]])
    out$sd <- unname(statistic[["sd"]])
    out$lower <- unname(statistic[["lower"]])
    out$upper <- unname(statistic[["upper"]])
    out
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

scenario_label <- c(
  continuous = "Continuous",
  batch_confounding = "Batch-confounded"
)
method_label <- c(
  fibrodynmix_nb = "NB",
  fibrodynmix_nb_study_donor = "NB + study + donor"
)
lambda_spec <- data.frame(
  directory = c("lambda_000", "lambda_005", "lambda_160", "lambda_1600"),
  marker_l2 = c(0, 0.05, 1.6, 16),
  marker_l2_label = c("0", "0.05", "1.6", "16"),
  stringsAsFactors = FALSE
)
fixed_root <- file.path(ROOT, "analysis", "marker_l2_fixed_scaling_ablation_v1")
normalized_root <- file.path(ROOT, "analysis", "marker_l2_likelihood_normalized_v1")

read_fixed_run <- function(directory, marker_l2, marker_l2_label) {
  path <- file.path(fixed_root, directory, "trajectory_metrics.tsv")
  metrics <- read_required_tsv(path)
  require_columns(
    metrics,
    c("scenario", "replicate", "method", "iteration", "z_rmse", "beta_rmse_direct"),
    basename(path)
  )
  metrics$marker_l2 <- marker_l2
  metrics$marker_l2_label <- marker_l2_label
  metrics$input_path <- file.path("analysis", "marker_l2_fixed_scaling_ablation_v1", directory, "trajectory_metrics.tsv")
  metrics
}

fixed_metrics <- do.call(rbind, lapply(seq_len(nrow(lambda_spec)), function(index) {
  read_fixed_run(
    lambda_spec$directory[[index]],
    lambda_spec$marker_l2[[index]],
    lambda_spec$marker_l2_label[[index]]
  )
}))
row.names(fixed_metrics) <- NULL

fixed_end <- fixed_metrics[fixed_metrics$iteration == 20L, , drop = FALSE]
fixed_start <- fixed_metrics[fixed_metrics$iteration == 0L & fixed_metrics$marker_l2 == 0.05, , drop = FALSE]
expected_conditions <- expand.grid(
  scenario = names(scenario_label),
  method = names(method_label),
  replicate = seq_len(5L),
  stringsAsFactors = FALSE
)
if (nrow(fixed_end) != 80L || anyDuplicated(fixed_end[c("scenario", "method", "replicate", "marker_l2")])) {
  stop("Figure S3 requires 80 unique fixed-scale endpoints (4 strengths x 4 conditions x 5 replicates).", call. = FALSE)
}
if (nrow(fixed_start) != 20L || anyDuplicated(fixed_start[c("scenario", "method", "replicate")])) {
  stop("Figure S3 requires 20 common-initializer rows from the fixed-scale analysis.", call. = FALSE)
}
if (!identical(
  sort(interaction(fixed_start$scenario, fixed_start$method, fixed_start$replicate)),
  sort(interaction(expected_conditions$scenario, expected_conditions$method, expected_conditions$replicate))
)) {
  stop("The fixed-scale initializer rows do not match the planned scenarios, methods, and replicates.", call. = FALSE)
}

normalized_metrics_path <- file.path(normalized_root, "trajectory_metrics.tsv")
normalized_components_path <- file.path(normalized_root, "objective_components.tsv")
normalized_metrics <- read_required_tsv(normalized_metrics_path)
normalized_components <- read_required_tsv(normalized_components_path)
require_columns(
  normalized_metrics,
  c("scenario", "replicate", "method", "iteration", "z_rmse"),
  basename(normalized_metrics_path)
)
require_columns(
  normalized_components,
  c("scenario", "replicate", "method", "iteration", "marker_l2_effective"),
  basename(normalized_components_path)
)
normalized_end <- normalized_metrics[normalized_metrics$iteration == 20L, , drop = FALSE]
if (nrow(normalized_end) != 20L || anyDuplicated(normalized_end[c("scenario", "method", "replicate")])) {
  stop("Figure S3 requires 20 unique likelihood-normalized endpoints.", call. = FALSE)
}

add_labels <- function(data) {
  data$scenario_label <- factor(scenario_label[data$scenario], levels = unname(scenario_label))
  data$method_label <- factor(method_label[data$method], levels = unname(method_label))
  data
}
fixed_end <- add_labels(fixed_end)
fixed_start <- add_labels(fixed_start)
normalized_end <- add_labels(normalized_end)
normalized_components <- add_labels(normalized_components)
fixed_end$marker_l2_label <- factor(fixed_end$marker_l2_label, levels = lambda_spec$marker_l2_label)

fixed_z_summary <- summarise_ci(
  fixed_end,
  c("scenario", "method", "scenario_label", "method_label", "marker_l2", "marker_l2_label"),
  "z_rmse"
)
fixed_beta_summary <- summarise_ci(
  fixed_end,
  c("scenario", "method", "scenario_label", "method_label", "marker_l2", "marker_l2_label"),
  "beta_rmse_direct"
)
initializer_summary <- summarise_ci(
  fixed_start,
  c("scenario", "method", "scenario_label", "method_label"),
  "z_rmse"
)

normalized_final <- normalized_end[, c("scenario", "method", "replicate", "scenario_label", "method_label", "z_rmse"), drop = FALSE]
normalized_final$treatment <- "Likelihood-normalized\n1% target"
fixed_reference <- fixed_end[fixed_end$marker_l2 == 0.05,
  c("scenario", "method", "replicate", "scenario_label", "method_label", "z_rmse"), drop = FALSE]
fixed_reference$treatment <- "Fixed\n0.05"
initializer_reference <- fixed_start[, c("scenario", "method", "replicate", "scenario_label", "method_label", "z_rmse"), drop = FALSE]
initializer_reference$treatment <- "Common marker\ninitializer"
normalized_control <- rbind(initializer_reference, fixed_reference, normalized_final)
normalized_control$treatment <- factor(
  normalized_control$treatment,
  levels = c("Common marker\ninitializer", "Fixed\n0.05", "Likelihood-normalized\n1% target")
)
normalized_control_summary <- summarise_ci(
  normalized_control,
  c("scenario", "method", "scenario_label", "method_label", "treatment"),
  "z_rmse"
)

normalized_effective <- normalized_components[
  normalized_components$iteration %in% c(0L, 20L),
  c("scenario", "method", "replicate", "iteration", "marker_l2_effective"),
  drop = FALSE
]

fdm_write_tsv(
  fixed_end[, c("scenario", "method", "replicate", "marker_l2", "marker_l2_label", "z_rmse", "beta_rmse_direct", "z_rmse_delta_from_initializer")],
  file.path(SOURCE, "figs3_fixed_scale_endpoints.tsv")
)
fdm_write_tsv(fixed_z_summary, file.path(SOURCE, "figs3a_fixed_scale_state_weight_summary.tsv"))
fdm_write_tsv(fixed_beta_summary, file.path(SOURCE, "figs3b_fixed_scale_program_summary.tsv"))
fdm_write_tsv(initializer_summary, file.path(SOURCE, "figs3_common_initializer_summary.tsv"))
fdm_write_tsv(normalized_control, file.path(SOURCE, "figs3c_likelihood_normalized_control_endpoints.tsv"))
fdm_write_tsv(normalized_control_summary, file.path(SOURCE, "figs3c_likelihood_normalized_control_summary.tsv"))
fdm_write_tsv(normalized_effective, file.path(SOURCE, "figs3c_effective_marker_l2_at_endpoints.tsv"))

input_manifest <- data.frame(
  analysis_input = c(
    file.path("analysis", "marker_l2_fixed_scaling_ablation_v1", lambda_spec$directory, "trajectory_metrics.tsv"),
    file.path("analysis", "marker_l2_likelihood_normalized_v1", "trajectory_metrics.tsv"),
    file.path("analysis", "marker_l2_likelihood_normalized_v1", "objective_components.tsv")
  ),
  purpose = c(
    rep("Fixed marker_l2 endpoint sensitivity", nrow(lambda_spec)),
    "Likelihood-normalized endpoint control",
    "Effective likelihood-normalized marker_l2 record"
  ),
  stringsAsFactors = FALSE
)
input_manifest$exists <- file.exists(file.path(ROOT, input_manifest$analysis_input))
input_manifest$bytes <- ifelse(
  input_manifest$exists,
  file.info(file.path(ROOT, input_manifest$analysis_input))$size,
  NA_real_
)
fdm_write_tsv(input_manifest, file.path(OUT, "analysis_input_manifest.tsv"))

palette <- fdm_editorial_palette()
ink <- unname(palette[["ink"]])
neutral_dark <- unname(palette[["neutral_dark"]])
signal <- unname(palette[["signal"]])
signal_deep <- unname(palette[["signal_deep"]])
selected <- unname(palette[["selected"]])
selected_deep <- unname(palette[["selected_deep"]])
grid_col <- unname(palette[["grid"]])
fixed_colours <- c("0" = neutral_dark, "0.05" = selected, "1.6" = signal, "16" = signal_deep)
control_colours <- c(
  "Common marker\ninitializer" = neutral_dark,
  "Fixed\n0.05" = selected,
  "Likelihood-normalized\n1% target" = signal
)

theme_set(fdm_theme(base_size = 7.4))
black_text_theme <- theme(
  text = element_text(colour = "#000000"),
  plot.title = element_text(colour = "#000000"),
  plot.subtitle = element_text(colour = "#000000"),
  axis.title = element_text(colour = "#000000"),
  axis.text = element_text(colour = "#000000"),
  legend.title = element_text(colour = "#000000"),
  legend.text = element_text(colour = "#000000"),
  strip.text = element_text(colour = "#000000")
)
panel_theme <- theme(
  panel.grid.major.y = element_line(linewidth = 0.22, colour = grid_col),
  strip.text = element_text(size = 7.2, face = "bold"),
  strip.background = element_rect(fill = "white", colour = ink, linewidth = 0.28),
  plot.title = element_text(margin = margin(b = 1)),
  plot.subtitle = element_text(margin = margin(b = 5)),
  plot.margin = margin(t = 12, r = 10, b = 8, l = 10)
)

panel_a <- ggplot(fixed_end, aes(x = marker_l2_label, y = z_rmse, colour = marker_l2_label)) +
  geom_hline(
    data = initializer_summary,
    aes(yintercept = mean),
    inherit.aes = FALSE,
    colour = ink,
    linetype = "dashed",
    linewidth = 0.34
  ) +
  geom_point(
    position = position_jitter(width = 0.065, height = 0),
    size = 1.25,
    alpha = 0.74
  ) +
  geom_errorbar(
    data = fixed_z_summary,
    aes(x = marker_l2_label, ymin = lower, ymax = upper, colour = marker_l2_label),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.52
  ) +
  geom_point(
    data = fixed_z_summary,
    aes(x = marker_l2_label, y = mean, colour = marker_l2_label),
    inherit.aes = FALSE,
    size = 2.05
  ) +
  facet_grid(scenario_label ~ method_label) +
  scale_colour_manual(values = fixed_colours, name = expression(marker_l2), guide = "none") +
  labs(
    title = "A. Fixed program-level marker penalties do not stabilize state weights",
    subtitle = "Final state-weight RMSE after 20 updates; points = 5 matched replicates; dashed = common marker initializer",
    x = expression(marker_l2),
    y = "State-weight RMSE"
  ) +
  panel_theme +
  black_text_theme

panel_b <- ggplot(fixed_end, aes(x = marker_l2_label, y = beta_rmse_direct, colour = marker_l2_label)) +
  geom_point(
    position = position_jitter(width = 0.065, height = 0),
    size = 1.25,
    alpha = 0.74
  ) +
  geom_errorbar(
    data = fixed_beta_summary,
    aes(x = marker_l2_label, ymin = lower, ymax = upper, colour = marker_l2_label),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.52
  ) +
  geom_point(
    data = fixed_beta_summary,
    aes(x = marker_l2_label, y = mean, colour = marker_l2_label),
    inherit.aes = FALSE,
    size = 2.05
  ) +
  facet_grid(scenario_label ~ method_label) +
  scale_colour_manual(values = fixed_colours, name = expression(marker_l2), guide = "none") +
  scale_y_log10() +
  labs(
    title = "B. Program-coefficient recovery and cell-coordinate recovery separate",
    subtitle = "Direct state-program coefficient RMSE after 20 updates; logarithmic y-axis; points = 5 matched replicates",
    x = expression(marker_l2),
    y = "Direct program-coefficient RMSE (log scale)"
  ) +
  panel_theme +
  black_text_theme

panel_c <- ggplot(normalized_control, aes(x = treatment, y = z_rmse, colour = treatment)) +
  geom_point(
    position = position_jitter(width = 0.065, height = 0),
    size = 1.25,
    alpha = 0.74
  ) +
  geom_errorbar(
    data = normalized_control_summary,
    aes(x = treatment, ymin = lower, ymax = upper, colour = treatment),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.52
  ) +
  geom_point(
    data = normalized_control_summary,
    aes(x = treatment, y = mean, colour = treatment),
    inherit.aes = FALSE,
    size = 2.05
  ) +
  facet_grid(scenario_label ~ method_label) +
  scale_colour_manual(values = control_colours, name = NULL, guide = "none") +
  labs(
    title = "C. A likelihood-normalized program penalty does not replace a cell anchor",
    subtitle = "The dynamic coefficient targets a 1% marker-penalty/NB-likelihood ratio before each program update; points = 5 matched replicates",
    x = NULL,
    y = "State-weight RMSE"
  ) +
  panel_theme +
  black_text_theme

combined <- panel_a / panel_b / panel_c + plot_layout(heights = c(1, 1, 1))
fdm_save_plot(combined, EXPORTS, STEM, width = fdm_main_width(), height = 10.8, dpi = 600)
fdm_write_export_qc(EXPORTS, QC, STEM)

legend_text <- paste(
  "# Figure S3. Program-level marker penalties do not stabilize cell-state coordinates in the tested simulations",
  "",
  "**(A)** Fixed-strength sensitivity analysis for `marker_l2 = 0`, `0.05`, `1.6`, or `16`. State-weight RMSE was evaluated after 20 updates in continuous and batch-confounded simulations, with NB and NB + study + donor fits. Every fit used `z_anchor = \"none\"` (no cell-specific marker-logit anchor) and began from the matched marker-guided initializer. Points are five matched replicate fits, error bars are 95% *t* intervals, and dashed lines mark the corresponding common initializer. **(B)** Direct state-program coefficient RMSE from the same, sum-to-zero-canonicalized fits. Stronger program-level targeting reduced coefficient error without stabilizing cell-state weights in these simulations. **(C)** A dynamic likelihood-normalized control. Before each program update, the effective `marker_l2` was selected to target a 1% marker-penalty/negative-binomial-likelihood ratio; it still did not return final state-weight RMSE to the common marker initializer. This control changes the penalty weight between updates and is not interpreted by cross-iteration objective values. These are small known-truth simulations, not a general penalty-tuning recommendation or a real-data validation result.",
  sep = "\n"
)
writeLines(legend_text, file.path(OUT, "main_figure_legend.md"), useBytes = TRUE)

panel_manifest <- data.frame(
  panel = c("A", "B", "C"),
  source_data = c(
    "source_data/figs3_fixed_scale_endpoints.tsv; source_data/figs3a_fixed_scale_state_weight_summary.tsv",
    "source_data/figs3_fixed_scale_endpoints.tsv; source_data/figs3b_fixed_scale_program_summary.tsv",
    "source_data/figs3c_likelihood_normalized_control_endpoints.tsv; source_data/figs3c_likelihood_normalized_control_summary.tsv; source_data/figs3c_effective_marker_l2_at_endpoints.tsv"
  ),
  results_anchor = c(
    "Fixed program-level marker-penalty scaling did not stabilize z in the tested simulations.",
    "Program-coefficient recovery and cell-coordinate recovery separated.",
    "Likelihood-normalized program penalty did not replace a cell anchor."
  ),
  claim_boundary = c(
    "Known-truth simulation without the cell-specific marker-logit anchor; no cross-penalty objective comparison.",
    "Direct beta comparison uses simulator truth after the same sum-to-zero canonicalization.",
    "Dynamic scale control; not a learned prior, general tuning rule, or real-data validation."
  ),
  stringsAsFactors = FALSE
)
fdm_write_tsv(panel_manifest, file.path(OUT, "panel_source_data_manifest.tsv"))

figure_manifest <- data.frame(
  figure = "Figure S3",
  title = "Program-level marker penalties do not stabilize the cell-state coordinate",
  panels = "A-C",
  exports = paste0("exports/", STEM, c(".pdf", ".svg", ".png", ".tiff"), collapse = "; "),
  all_nondata_text_black = TRUE,
  source_data_manifest = "panel_source_data_manifest.tsv",
  stringsAsFactors = FALSE
)
fdm_write_tsv(figure_manifest, file.path(OUT, "figure_manifest.tsv"))

message("Supplementary Figure S3 written to: ", normalizePath(OUT))
