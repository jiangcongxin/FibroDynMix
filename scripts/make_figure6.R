#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure6.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "figures", "figure6")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
PALETTES <- file.path(OUT, "palettes")
QC <- file.path(OUT, "qc")

dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(PALETTES, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "figure_style_fdm.R"))

required <- file.path(
  ROOT,
  c(
    "analysis/bootstrap_uncertainty/sample_composition_uncertainty.tsv",
    "analysis/bootstrap_uncertainty/cell_state_uncertainty.tsv",
    "analysis/vi_posterior/vi_cell_state_intervals.tsv",
    "analysis/vi_posterior/vi_cell_entropy_intervals.tsv",
    "analysis/vi_benchmark/vi_benchmark_summary.tsv",
    "analysis/selected_nb_final_summaries/selected_nb_vi_benchmark_summary.tsv"
  )
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop(sprintf("Missing uncertainty inputs for Figure 6: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

state_levels <- c("resident", "inflammatory", "myofibroblast", "ECM-remodeling", "antigen-presenting", "IFN-stress")
state_labels <- c(
  resident = "Resident",
  inflammatory = "Inflamm.",
  myofibroblast = "MyoFB",
  "ECM-remodeling" = "ECM",
  "antigen-presenting" = "AP",
  "IFN-stress" = "IFN"
)
state_palette <- fdm_state_palette()[state_levels]
fdm_write_tsv(data.frame(name = names(state_palette), colour = unname(state_palette)), file.path(PALETTES, "figure6_state_palette.tsv"))

bootstrap_comp <- fdm_read_tsv(file.path(ROOT, "analysis", "bootstrap_uncertainty", "sample_composition_uncertainty.tsv"))
bootstrap_cell <- fdm_read_tsv(file.path(ROOT, "analysis", "bootstrap_uncertainty", "cell_state_uncertainty.tsv"))
vi_state <- fdm_read_tsv(file.path(ROOT, "analysis", "vi_posterior", "vi_cell_state_intervals.tsv"))
vi_entropy <- fdm_read_tsv(file.path(ROOT, "analysis", "vi_posterior", "vi_cell_entropy_intervals.tsv"))
vi_benchmark <- fdm_read_tsv(file.path(ROOT, "analysis", "vi_benchmark", "vi_benchmark_summary.tsv"))
selected_vi <- fdm_read_tsv(file.path(ROOT, "analysis", "selected_nb_final_summaries", "selected_nb_vi_benchmark_summary.tsv"))

bootstrap_comp$state <- factor(bootstrap_comp$state, levels = state_levels)
bootstrap_comp$state_label <- factor(state_labels[as.character(bootstrap_comp$state)], levels = rev(state_labels[state_levels]))
bootstrap_comp$sample_id <- factor(bootstrap_comp$sample_id, levels = unique(bootstrap_comp$sample_id))

bootstrap_cell$interval_width <- bootstrap_cell$upper - bootstrap_cell$lower
bootstrap_cell_summary <- aggregate(interval_width ~ state, data = bootstrap_cell, FUN = mean)
bootstrap_cell_summary$state <- factor(bootstrap_cell_summary$state, levels = state_levels)

vi_state$interval_width <- vi_state$upper - vi_state$lower
vi_state$state <- factor(vi_state$state, levels = state_levels)
vi_state$state_label <- factor(state_labels[as.character(vi_state$state)], levels = state_labels[state_levels])

vi_entropy$entropy_width <- vi_entropy$upper - vi_entropy$lower
vi_cell_width <- aggregate(interval_width ~ cell_id, data = vi_state, FUN = mean)
vi_uncertainty <- merge(
  vi_cell_width,
  vi_entropy[, c("cell_id", "mean", "entropy_width")],
  by = "cell_id",
  all.x = TRUE
)
colnames(vi_uncertainty)[colnames(vi_uncertainty) == "mean"] <- "mean_entropy"

coverage <- selected_vi[
  selected_vi$method == "fibrodynmix_vi" &
    !is.na(selected_vi$vi_interval_coverage_mean),
  ,
  drop = FALSE
]
coverage_raw <- data.frame(
  scenario = coverage$scenario,
  interval_type = "Raw VI",
  coverage = coverage$vi_interval_coverage_mean,
  width = coverage$vi_mean_interval_width_mean,
  selected_n_outer = coverage$selected_n_outer,
  stringsAsFactors = FALSE
)
coverage_cal <- data.frame(
  scenario = coverage$scenario,
  interval_type = "Calibrated VI",
  coverage = coverage$vi_calibrated_interval_coverage_mean,
  width = coverage$vi_calibrated_mean_interval_width_mean,
  selected_n_outer = coverage$selected_n_outer,
  stringsAsFactors = FALSE
)
coverage_plot <- rbind(coverage_raw, coverage_cal)
coverage_plot$interval_type <- factor(coverage_plot$interval_type, levels = c("Raw VI", "Calibrated VI"))
coverage_plot$scenario <- factor(
  coverage_plot$scenario,
  levels = c("continuous", "batch_confounding", "rare_transition"),
  labels = c("Continuous", "Batch", "Rare")
)

fdm_write_tsv(bootstrap_comp, file.path(SOURCE, "fig6_bootstrap_sample_composition_intervals.tsv"))
fdm_write_tsv(bootstrap_cell_summary, file.path(SOURCE, "fig6_bootstrap_cell_state_interval_width.tsv"))
fdm_write_tsv(vi_state, file.path(SOURCE, "fig6_vi_cell_state_intervals.tsv"))
fdm_write_tsv(vi_uncertainty, file.path(SOURCE, "fig6_vi_uncertainty_vs_entropy.tsv"))
fdm_write_tsv(coverage_plot, file.path(SOURCE, "fig6_vi_interval_coverage_width.tsv"))

theme_set(fdm_theme(base_size = 9))

panel_a <- ggplot(bootstrap_comp, aes(x = mean, y = state_label, xmin = lower, xmax = upper, colour = state)) +
  geom_segment(aes(x = lower, xend = upper, yend = state_label), linewidth = 0.35, alpha = 0.85) +
  geom_point(size = 1.5) +
  facet_wrap(~ sample_id, nrow = 2) +
  scale_colour_manual(values = state_palette, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, NA)) +
  labs(
    title = "Bootstrap composition intervals",
    subtitle = "Sample-level state composition, mean and 95% interval",
    x = "State fraction",
    y = NULL
  )

panel_b <- ggplot(vi_state, aes(x = state_label, y = interval_width, fill = state)) +
  geom_violin(width = 0.78, linewidth = 0.25, scale = "width", alpha = 0.82) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, alpha = 0.92) +
  scale_fill_manual(values = state_palette, guide = "none") +
  labs(
    title = "VI state-weight interval width",
    subtitle = "Cell-level logistic-normal posterior intervals",
    x = NULL,
    y = "95% interval width"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

panel_c <- ggplot(coverage_plot, aes(x = width, y = coverage, colour = interval_type, shape = scenario)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.3, colour = "#6B7280") +
  geom_point(size = 2.5, stroke = 0.6) +
  geom_line(aes(group = scenario), linewidth = 0.35, colour = "#B0B0B0") +
  scale_colour_manual(values = c("Raw VI" = "#386CB0", "Calibrated VI" = "#984EA3"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1)) +
  labs(
    title = "Interval coverage calibration",
    subtitle = "Dashed line marks nominal 95% coverage",
    x = "Mean interval width",
    y = "Empirical coverage",
    shape = "Scenario"
  ) +
  theme(legend.position = "right")

panel_d <- ggplot(vi_uncertainty, aes(x = mean_entropy, y = interval_width)) +
  geom_point(size = 1.8, alpha = 0.72, colour = "#287C6F") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45, colour = "#555555") +
  labs(
    title = "Uncertainty tracks mixed-state entropy",
    subtitle = "Cells with higher entropy have wider intervals",
    x = "VI mean entropy",
    y = "Mean interval width"
  )

figure6 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(widths = c(1.12, 0.88), heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A") &
  fdm_panel_tag_theme()

fdm_save_plot(figure6, EXPORTS, "figure6", width = 8.2, height = 6.4, dpi = 600)
fdm_write_export_qc(EXPORTS, QC, "figure6")

fdm_write_tsv(
  data.frame(
    figure = "Figure 6",
    primary_claim = "FibroDynMix reports uncertainty through bootstrap state-composition intervals, cell-level VI state-weight intervals, simulation-calibrated VI coverage diagnostics, and entropy-linked uncertainty summaries.",
    claim_boundary = "Uncertainty summaries combine bootstrap and lightweight logistic-normal VI diagnostics; they are not a full all-parameter Bayesian posterior or real-data calibration guarantee.",
    bootstrap_source = "analysis/bootstrap_uncertainty",
    vi_source = "analysis/vi_posterior; analysis/selected_nb_final_summaries",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C", "D"),
    source_data = c(
      "source_data/fig6_bootstrap_sample_composition_intervals.tsv",
      "source_data/fig6_vi_cell_state_intervals.tsv",
      "source_data/fig6_vi_interval_coverage_width.tsv",
      "source_data/fig6_vi_uncertainty_vs_entropy.tsv"
    ),
    claim = c(
      "Bootstrap resampling provides sample-level state-composition intervals.",
      "The lightweight VI posterior returns cell-level state-weight interval widths by state.",
      "Simulation-calibrated VI intervals improve empirical coverage at the cost of wider intervals.",
      "Cells with greater mixed-state entropy tend to show wider state-weight uncertainty."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 6. Uncertainty quantification and VI calibration",
    "",
    "**A.** Bootstrap sample-level state-composition means and 95% intervals.",
    "**B.** Cell-level VI state-weight interval widths by state.",
    "**C.** Raw and simulation-calibrated VI interval coverage versus mean interval width; dashed line marks nominal 95% coverage.",
    "**D.** Relationship between cell-level VI mean entropy and mean state-weight interval width.",
    "",
    "Bootstrap intervals and lightweight logistic-normal VI intervals quantify uncertainty around state-mixture inference. These diagnostics do not constitute a full all-parameter Bayesian posterior or real-data calibration guarantee."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Figure 6 uncertainty figure written to: ", OUT)
