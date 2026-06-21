#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure3_redesign.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE_IN <- file.path(ROOT, "figures", "figure3", "source_data")
OUT <- Sys.getenv("FDM_FIGURE_OUT", unset = file.path(ROOT, "figures", "figure3_redesign"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure3_redesign")
FIGURE_LABEL <- Sys.getenv("FDM_FIGURE_LABEL", unset = "Figure 3 redesign")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
QC <- file.path(OUT, "qc")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "figure_style_fdm.R"))

scenario_levels <- c("continuous", "discrete", "batch_confounding", "rare_transition")
scenario_labels <- c(
  continuous = "Continuous",
  discrete = "Discrete",
  batch_confounding = "Batch",
  rare_transition = "Rare"
)
method_levels <- c("marker_scoring", "fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_vi")
method_labels <- c(
  marker_scoring = "Marker",
  fibrodynmix_initializer = "FDM init.",
  fibrodynmix_nb = "FDM NB",
  fibrodynmix_vi = "FDM VI"
)
method_palette <- fdm_method_palette()[method_levels]

benchmark <- fdm_read_tsv(file.path(SOURCE_IN, "fig3_benchmark_metrics.tsv"))
benchmark <- benchmark[
  benchmark$scenario %in% scenario_levels &
    benchmark$method %in% method_levels &
    !is.na(benchmark$rmse),
  ,
  drop = FALSE
]

entropy <- fdm_read_tsv(file.path(SOURCE_IN, "fig3_panel_d_entropy_calibration.tsv"))
entropy <- entropy[entropy$scenario %in% scenario_levels & entropy$method %in% method_levels, , drop = FALSE]
marker <- fdm_read_tsv(file.path(SOURCE_IN, "fig3_panel_f_marker_recovery_summary.tsv"))
marker <- marker[marker$scenario %in% scenario_levels & marker$method %in% method_levels, , drop = FALSE]
transition <- fdm_read_tsv(file.path(SOURCE_IN, "fig3_panel_e_transition_scores.tsv"))
transition <- transition[transition$method %in% method_levels, , drop = FALSE]

benchmark$scenario <- factor(benchmark$scenario, levels = scenario_levels)
benchmark$method <- factor(benchmark$method, levels = rev(method_levels))
benchmark$scenario_label <- factor(scenario_labels[as.character(benchmark$scenario)], levels = scenario_labels[scenario_levels])
benchmark$method_label <- factor(method_labels[as.character(benchmark$method)], levels = rev(method_labels[method_levels]))

summary_core <- aggregate(
  cbind(rmse, dominant_accuracy, downstream_balanced_accuracy, downstream_macro_f1) ~ scenario + scenario_label + method + method_label,
  data = benchmark,
  FUN = mean,
  na.rm = TRUE
)
summary_core$rmse_label <- sprintf("%.2f", summary_core$rmse)

entropy_plot <- entropy[entropy$scenario %in% c("continuous", "rare_transition") &
  entropy$method %in% c("marker_scoring", "fibrodynmix_nb", "fibrodynmix_vi"), , drop = FALSE]
set.seed(930)
if (nrow(entropy_plot) > 850L) {
  entropy_plot <- entropy_plot[sample(seq_len(nrow(entropy_plot)), 850L), , drop = FALSE]
}
entropy_plot$scenario_label <- factor(scenario_labels[entropy_plot$scenario], levels = scenario_labels[c("continuous", "rare_transition")])
entropy_plot$method <- factor(entropy_plot$method, levels = method_levels)

marker$scenario <- factor(marker$scenario, levels = scenario_levels)
marker$scenario_label <- factor(scenario_labels[as.character(marker$scenario)], levels = scenario_labels[scenario_levels])
marker$method <- factor(marker$method, levels = method_levels)
marker_summary <- aggregate(marker_auprc ~ scenario_label + method, data = marker, FUN = mean, na.rm = TRUE)

transition_plot <- transition[transition$method %in% c("marker_scoring", "fibrodynmix_nb", "fibrodynmix_vi"), , drop = FALSE]
transition_plot$method <- factor(transition_plot$method, levels = c("marker_scoring", "fibrodynmix_nb", "fibrodynmix_vi"))
transition_plot$transition_class <- factor(
  ifelse(as.character(transition_plot$is_transition) == "TRUE", "Rare transition", "Other cells"),
  levels = c("Other cells", "Rare transition")
)

fdm_write_tsv(summary_core, file.path(SOURCE, "fig3_redesign_core_summary.tsv"))
fdm_write_tsv(entropy_plot, file.path(SOURCE, "fig3_redesign_entropy_sample.tsv"))
fdm_write_tsv(marker_summary, file.path(SOURCE, "fig3_redesign_marker_summary.tsv"))
fdm_write_tsv(transition_plot, file.path(SOURCE, "fig3_redesign_transition_scores.tsv"))

theme_set(fdm_theme(base_size = 8))

panel_a <- ggplot(summary_core, aes(x = scenario_label, y = method_label, fill = rmse)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = rmse_label), size = 2.25, fontface = "bold", colour = "#1F2937") +
  scale_fill_gradient(low = "#2C7A7B", high = "#F8E7DF", name = "RMSE\nlower is better") +
  labs(
    title = "State-weight recovery heatmap",
    subtitle = "Mean RMSE across known-truth simulations",
    x = NULL,
    y = NULL
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

panel_b <- ggplot(summary_core, aes(x = rmse, y = downstream_balanced_accuracy, colour = method, shape = scenario_label)) +
  geom_point(aes(size = dominant_accuracy), alpha = 0.9, stroke = 0.35) +
  scale_colour_manual(values = method_palette, labels = method_labels, name = NULL) +
  scale_shape_discrete(name = "Scenario") +
  scale_size_continuous(range = c(1.7, 4.2), name = "Dominant\naccuracy") +
  labs(
    title = "Recovery vs downstream utility",
    subtitle = "Each point is a scenario-method mean",
    x = "State-weight RMSE",
    y = "Balanced accuracy"
  ) +
  theme(legend.position = "right")

panel_c <- ggplot(entropy_plot, aes(x = true_entropy, y = pred_entropy, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.25, colour = "#6B7280") +
  geom_point(size = 0.55, alpha = 0.28) +
  facet_wrap(~ scenario_label, nrow = 1) +
  scale_colour_manual(values = method_palette, labels = method_labels, name = NULL) +
  labs(
    title = "Entropy calibration scatter",
    subtitle = "Sampled cells; dashed line is equality",
    x = "True entropy",
    y = "Predicted entropy"
  ) +
  theme(legend.position = "bottom")

panel_d <- ggplot(marker_summary, aes(x = scenario_label, y = marker_auprc, colour = method, group = method)) +
  geom_line(linewidth = 0.35, alpha = 0.75) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = method_palette, labels = method_labels, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Marker-program recovery",
    x = NULL,
    y = "AUPRC"
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")

panel_e <- ggplot(transition_plot, aes(x = method, y = entropy_score, fill = transition_class)) +
  geom_violin(width = 0.78, linewidth = 0.25, scale = "width", alpha = 0.82) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, alpha = 0.85) +
  scale_x_discrete(labels = method_labels[c("marker_scoring", "fibrodynmix_nb", "fibrodynmix_vi")]) +
  scale_fill_manual(values = c("Other cells" = "#BDBDBD", "Rare transition" = fdm_condition_palette()[["disease"]]), name = NULL) +
  labs(
    title = "Rare-transition entropy signal",
    x = NULL,
    y = "Entropy score"
  ) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "top")

figure3 <- (panel_a | panel_b) / panel_c / (panel_d | panel_e) +
  plot_layout(heights = c(1.08, 1.0, 0.95), widths = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A") &
  fdm_panel_tag_theme()

fdm_save_plot(figure3, EXPORTS, STEM, width = 8.4, height = 7.2)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = FIGURE_LABEL,
    primary_claim = "Known-truth simulations summarize core state-weight recovery, downstream utility, entropy calibration, marker-program recovery, and rare-transition signal while preserving baseline-competitive boundaries.",
    claim_boundary = "Simulation benchmark only; not patient-tissue causal validation.",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C", "D", "E"),
    source_data = c(
      "source_data/fig3_redesign_core_summary.tsv",
      "source_data/fig3_redesign_core_summary.tsv",
      "source_data/fig3_redesign_entropy_sample.tsv",
      "source_data/fig3_redesign_marker_summary.tsv",
      "source_data/fig3_redesign_transition_scores.tsv"
    ),
    claim = c(
      "Mean state-weight RMSE is summarized across methods and known-truth simulation scenarios.",
      "State recovery is compared against downstream balanced accuracy and dominant-state accuracy.",
      "Entropy calibration is shown in selected continuous-mixture and rare-transition scenarios.",
      "Marker-program recovery is quantified by AUPRC across scenarios.",
      "Rare-transition cells show elevated entropy signal relative to other cells in the rare-transition simulation."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 3. Known-truth simulation benchmark",
    "",
    "**A.** Heatmap of mean state-weight RMSE across simulation scenarios and methods.",
    "**B.** Scenario-method means comparing state-weight RMSE with downstream balanced accuracy; point size encodes dominant-state accuracy.",
    "**C.** True-versus-predicted entropy scatter for sampled cells in continuous and rare-transition simulations.",
    "**D.** Marker-program recovery summarized as AUPRC across scenarios.",
    "**E.** Rare-transition entropy-score distributions for selected methods.",
    "",
    "These are simulation benchmarks with known latent truth. Marker scoring and other baselines remain competitive in some settings, so the claim is not uniform superiority."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Figure 3 redesign written to: ", file.path(EXPORTS, paste0(STEM, ".png")))
