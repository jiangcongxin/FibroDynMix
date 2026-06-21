#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure4_redesign.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE_IN <- file.path(ROOT, "figures", "figure4", "source_data")
OUT <- Sys.getenv("FDM_FIGURE_OUT", unset = file.path(ROOT, "figures", "figure4_redesign"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure4_redesign")
FIGURE_LABEL <- Sys.getenv("FDM_FIGURE_LABEL", unset = "Figure 4 redesign")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
QC <- file.path(OUT, "qc")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "figure_style_fdm.R"))

tradeoff <- fdm_read_tsv(file.path(SOURCE_IN, "fig4_study_effect_tradeoff.tsv"))
selected_lines <- readLines(file.path(SOURCE_IN, "fig4_selected_penalty.txt"), warn = FALSE)
selected_study <- as.numeric(sub("recommended_study_l2: ", "", selected_lines[grep("recommended_study_l2", selected_lines)]))
selected_marker <- as.numeric(sub("recommended_marker_l2: ", "", selected_lines[grep("recommended_marker_l2", selected_lines)]))

tradeoff$study_l2_factor <- factor(as.character(tradeoff$study_l2), levels = c("0.05", "0.1", "0.5", "1", "5"))
tradeoff$marker_l2_factor <- factor(as.character(tradeoff$marker_l2), levels = c("0.05", "0.1"))
tradeoff$is_selected <- tradeoff$study_l2 == selected_study & tradeoff$marker_l2 == selected_marker
selected_row <- tradeoff[tradeoff$is_selected, , drop = FALSE]

fdm_write_tsv(tradeoff, file.path(SOURCE, "fig4_redesign_tradeoff.tsv"))

marker_palette <- c("0.05" = fdm_method_palette()[["fibrodynmix_nb"]], "0.1" = fdm_method_palette()[["fibrodynmix_initializer"]])
theme_set(fdm_theme(base_size = 8))

panel_a <- ggplot(tradeoff, aes(x = study_l2_factor, y = marker_l2_factor, fill = selection_score)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", selection_score)), size = 2.5, fontface = "bold", colour = "#1F2937") +
  geom_point(
    data = selected_row,
    aes(x = study_l2_factor, y = marker_l2_factor),
    shape = 21,
    size = 6.0,
    fill = NA,
    colour = "#111111",
    stroke = 0.85
  ) +
  scale_fill_gradient(low = "#2C7A7B", high = "#F6E5D8", name = "Selection\nscore") +
  labs(
    tag = "A",
    title = "Penalty-selection heatmap",
    subtitle = "Lower score preferred; ring marks selected grid point",
    x = "Study-effect penalty",
    y = "Marker penalty"
  )

panel_b <- ggplot(tradeoff, aes(x = study_l2, y = rmse_mean, colour = marker_l2_factor, group = marker_l2_factor)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.9) +
  geom_point(data = selected_row, shape = 21, size = 3.8, fill = "white", stroke = 0.8) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.5, 1, 5)) +
  scale_colour_manual(values = marker_palette, name = "Marker penalty") +
  labs(
    tag = "B",
    title = "State-recovery calibration",
    subtitle = "Stronger study shrinkage lowers RMSE",
    x = "Study-effect penalty",
    y = "Mean RMSE"
  ) +
  theme(legend.position = "top")

panel_c <- ggplot(tradeoff, aes(x = study_l2, y = study_effect_l2_norm_mean, colour = marker_l2_factor, group = marker_l2_factor)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.9) +
  geom_point(data = selected_row, shape = 21, size = 3.8, fill = "white", stroke = 0.8) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.5, 1, 5)) +
  scale_colour_manual(values = marker_palette, name = "Marker penalty") +
  labs(
    tag = "C",
    title = "Study-effect magnitude is regularized",
    x = "Study-effect penalty",
    y = "Study-effect L2 norm"
  ) +
  theme(legend.position = "none")

figure4 <- panel_a | (panel_b / panel_c) +
  plot_layout(widths = c(1.12, 0.88)) &
  fdm_panel_tag_theme()

fdm_save_plot(figure4, EXPORTS, STEM, width = 7.8, height = 5.0)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = FIGURE_LABEL,
    primary_claim = "A two-dimensional penalty grid selects a study-effect and marker-orientation setting that balances state recovery and study-effect shrinkage.",
    claim_boundary = "Simulation calibration only; not universal hyperparameter proof.",
    selected_study_l2 = selected_study,
    selected_marker_l2 = selected_marker,
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C"),
    source_data = c(
      "source_data/fig4_redesign_tradeoff.tsv",
      "source_data/fig4_redesign_tradeoff.tsv",
      "source_data/fig4_redesign_tradeoff.tsv"
    ),
    claim = c(
      "The two-dimensional penalty grid summarizes the selection score across study-effect and marker-orientation penalties.",
      "State-weight RMSE improves as study-effect shrinkage strengthens in the selected marker-penalty branch.",
      "Study-effect L2 norm decreases as ridge penalty increases."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 4. Study-effect penalty calibration",
    "",
    "**A.** Two-dimensional penalty-selection heatmap. Lower scores are preferred and the ring marks the selected penalty pair.",
    "**B.** Mean state-weight RMSE across study-effect penalties and marker-orientation penalties.",
    "**C.** Mean fitted study-effect L2 norm across the same penalty grid.",
    "",
    "This figure is a simulation calibration of penalty selection under batch confounding, not a universal hyperparameter proof."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Figure 4 redesign written to: ", file.path(EXPORTS, paste0(STEM, ".png")))
