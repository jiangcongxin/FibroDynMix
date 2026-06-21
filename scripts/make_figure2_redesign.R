#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure2_redesign.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE <- file.path(ROOT, "figures", "figure2", "source_data")
OUT <- Sys.getenv("FDM_FIGURE_OUT", unset = file.path(ROOT, "figures", "figure2_redesign"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure2_redesign")
FIGURE_LABEL <- Sys.getenv("FDM_FIGURE_LABEL", unset = "Figure 2 redesign")
EXPORTS <- file.path(OUT, "exports")
QC <- file.path(OUT, "qc")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "figure_style_fdm.R"))

required <- file.path(
  SOURCE,
  c(
    "fig2_condition_cell_counts.tsv",
    "fig2_nb_objective.tsv",
    "fig2_state_composition.tsv",
    "fig2_state_composition_delta.tsv",
    "fig2_cell_state_weight_heatmap.tsv",
    "fig2_cell_order_condition_annotation.tsv",
    "fig2_cell_fpi.tsv",
    "fig2_transition_flow.tsv"
  )
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop(sprintf("Missing Figure 2 source data files: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

read_tsv <- fdm_read_tsv
cell_counts <- read_tsv(file.path(SOURCE, "fig2_condition_cell_counts.tsv"))
objective <- read_tsv(file.path(SOURCE, "fig2_nb_objective.tsv"))
composition <- read_tsv(file.path(SOURCE, "fig2_state_composition.tsv"))
composition_delta <- read_tsv(file.path(SOURCE, "fig2_state_composition_delta.tsv"))
weights_long <- read_tsv(file.path(SOURCE, "fig2_cell_state_weight_heatmap.tsv"))
condition_blocks <- read_tsv(file.path(SOURCE, "fig2_cell_order_condition_annotation.tsv"))
fpi <- read_tsv(file.path(SOURCE, "fig2_cell_fpi.tsv"))
flow <- read_tsv(file.path(SOURCE, "fig2_transition_flow.tsv"))

state_levels <- c("resident", "inflammatory", "myofibroblast", "ECM-remodeling", "antigen-presenting", "IFN-stress")
state_short <- c(
  resident = "Resident",
  inflammatory = "Inflamm.",
  myofibroblast = "MyoFB",
  "ECM-remodeling" = "ECM",
  "antigen-presenting" = "AP",
  "IFN-stress" = "IFN"
)
state_palette <- fdm_state_palette()[state_levels]
condition_palette <- fdm_condition_palette()[c("normal", "disease")]

cell_counts$condition <- factor(cell_counts$condition, levels = c("normal", "disease"))
composition$condition <- factor(composition$condition, levels = c("normal", "disease"))
composition$state <- factor(composition$state, levels = state_levels)
fpi$condition <- factor(fpi$condition, levels = c("normal", "disease"))
condition_blocks$condition <- factor(condition_blocks$condition, levels = c("normal", "disease"))
flow$source_state <- factor(flow$source_state, levels = state_levels)
flow$target_state <- factor(flow$target_state, levels = state_levels)

delta_order <- composition_delta[order(composition_delta$delta_disease_minus_normal), ]
delta_order$state <- factor(delta_order$state, levels = delta_order$state)
delta_order$direction <- ifelse(delta_order$delta_disease_minus_normal >= 0, "disease", "normal")
delta_order$label_x <- delta_order$delta_disease_minus_normal +
  ifelse(delta_order$delta_disease_minus_normal >= 0, 0.035, -0.035)
delta_order$hjust <- ifelse(delta_order$delta_disease_minus_normal >= 0, 0, 1)

weights_long$state <- factor(weights_long$state, levels = rev(state_levels))
weights_long$state_label <- factor(state_short[as.character(weights_long$state)], levels = rev(state_short[state_levels]))

objective_improvement <- objective$nb_objective[objective$step == "Initial"] -
  objective$nb_objective[objective$step == "Best"]
objective_pct <- objective_improvement / objective$nb_objective[objective$step == "Initial"] * 100
n_cells <- sum(cell_counts$n_cells)
n_genes <- 1217

theme_set(fdm_theme(base_size = 8))

panel_a_data <- data.frame(
  xmin = c(0.02, 0.35, 0.68),
  xmax = c(0.29, 0.62, 0.98),
  ymin = c(0.30, 0.30, 0.30),
  ymax = c(0.78, 0.78, 0.78),
  fill = c(condition_palette[["normal"]], condition_palette[["disease"]], "#F3F4F6"),
  label = c(
    sprintf("Normal\nn=%s", cell_counts$n_cells[cell_counts$condition == "normal"]),
    sprintf("Disease\nn=%s", cell_counts$n_cells[cell_counts$condition == "disease"]),
    sprintf("NB fit\nobjective -%.1f%%\n%s genes", objective_pct, format(n_genes, big.mark = ","))
  ),
  label_colour = c("white", "white", "#4B5563")
)

panel_a <- ggplot(panel_a_data) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "#2F2F2F",
    linewidth = 0.25
  ) +
  geom_segment(
    data = data.frame(x = c(0.30, 0.63), xend = c(0.34, 0.67), y = 0.54, yend = 0.54),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = grid::unit(0.06, "in"), type = "closed"),
    linewidth = 0.25,
    colour = "#4B5563"
  ) +
  geom_text(aes(x = (xmin + xmax) / 2, y = 0.54, label = label, colour = label_colour), size = 2.45, lineheight = 0.88, fontface = "bold") +
  scale_fill_identity() +
  scale_colour_identity() +
  annotate("text", x = 0.02, y = 0.92, label = "Public raw-count run", hjust = 0, fontface = "bold", size = 3.0) +
  annotate("text", x = 0.02, y = 0.15, label = sprintf("Balanced subset: %s cells", n_cells), hjust = 0, size = 2.2, colour = "#555555") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_size = 8) +
  theme(plot.margin = margin(8, 4, 5, 4))

panel_b <- ggplot(delta_order, aes(x = delta_disease_minus_normal, y = state)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EEF3F9", alpha = 0.75) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FBEDEE", alpha = 0.75) +
  geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = delta_disease_minus_normal, yend = state), linewidth = 0.65, colour = "#555555") +
  geom_point(aes(fill = direction), shape = 21, size = 2.6, colour = "#333333", stroke = 0.25) +
  geom_text(aes(x = label_x, label = sprintf("%+.2f", delta_disease_minus_normal), hjust = hjust), size = 2.1, colour = "#333333") +
  scale_fill_manual(values = condition_palette, guide = "none") +
  scale_x_continuous(limits = c(-0.48, 0.58), breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  labs(
    title = "Disease subset gains ECM and myofibroblast programs",
    subtitle = "Disease-minus-normal fitted state weight",
    x = "Delta state weight",
    y = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.18, colour = "#E5E7EB"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

panel_c <- ggplot(weights_long, aes(x = cell_order, y = state_label, fill = state_weight)) +
  geom_rect(
    data = condition_blocks,
    aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = length(state_levels) + 0.62, ymax = length(state_levels) + 0.82, fill = NULL),
    inherit.aes = FALSE,
    fill = NA,
    colour = NA
  ) +
  geom_rect(
    data = condition_blocks[condition_blocks$condition == "normal", , drop = FALSE],
    aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = length(state_levels) + 0.64, ymax = length(state_levels) + 0.82),
    inherit.aes = FALSE,
    fill = condition_palette[["normal"]],
    colour = NA
  ) +
  geom_rect(
    data = condition_blocks[condition_blocks$condition == "disease", , drop = FALSE],
    aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = length(state_levels) + 0.64, ymax = length(state_levels) + 0.82),
    inherit.aes = FALSE,
    fill = condition_palette[["disease"]],
    colour = NA
  ) +
  geom_text(
    data = condition_blocks,
    aes(x = (xmin + xmax) / 2, y = length(state_levels) + 0.92, label = label),
    inherit.aes = FALSE,
    size = 2.1,
    colour = "#333333"
  ) +
  geom_tile(width = 1, height = 0.92) +
  scale_fill_gradient(low = "#F4F4EF", high = "#244A86", name = "Weight") +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(ylim = c(0.5, length(state_levels) + 1.05), clip = "off") +
  labs(
    title = "Cell-level mixture landscape",
    subtitle = "Cells ordered by condition and dominant state; rows are fitted latent state weights",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(8, 5, 5, 5)
  )

panel_fpi_supplement_candidate <- ggplot(fpi, aes(x = condition, y = fpi, fill = condition)) +
  geom_violin(width = 0.72, linewidth = 0.25, scale = "width", alpha = 0.88) +
  geom_jitter(aes(colour = condition), width = 0.12, height = 0, size = 0.45, alpha = 0.22, show.legend = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, alpha = 0.92) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 1.8, fill = "white", colour = "#222222") +
  scale_fill_manual(values = condition_palette, guide = "none") +
  scale_colour_manual(values = condition_palette, guide = "none") +
  labs(
    title = "FPI is lower in the disease subset",
    subtitle = "Raw cell-level values with median marker",
    x = NULL,
    y = "FPI"
  )

flow_plot <- flow[flow$flow >= 0.009, , drop = FALSE]
flow_plot$source_y <- match(as.character(flow_plot$source_state), state_levels)
flow_plot$target_y <- match(as.character(flow_plot$target_state), state_levels)
flow_plot$target_state_chr <- as.character(flow_plot$target_state)
flow_plot <- flow_plot[order(flow_plot$flow), ]
state_axis <- data.frame(
  state = state_levels,
  label = state_short[state_levels],
  y = seq_along(state_levels),
  stringsAsFactors = FALSE
)

panel_d <- ggplot() +
  geom_curve(
    data = flow_plot,
    aes(x = 0.20, xend = 0.80, y = source_y, yend = target_y, linewidth = flow, colour = target_state_chr),
    curvature = 0.22,
    alpha = 0.45,
    lineend = "round"
  ) +
  geom_point(data = state_axis, aes(x = 0.13, y = y), size = 2.2, shape = 21, fill = "white", colour = "#555555", stroke = 0.3) +
  geom_point(data = state_axis, aes(x = 0.87, y = y, fill = state), size = 2.4, shape = 21, colour = "#333333", stroke = 0.25) +
  geom_text(data = state_axis, aes(x = 0.09, y = y, label = label), hjust = 1, size = 2.1, colour = "#333333") +
  geom_text(data = state_axis, aes(x = 0.91, y = y, label = label), hjust = 0, size = 2.1, colour = "#333333") +
  annotate("text", x = 0.13, y = 6.65, label = "normal", size = 2.2, fontface = "bold", colour = condition_palette[["normal"]]) +
  annotate("text", x = 0.87, y = 6.65, label = "disease", size = 2.2, fontface = "bold", colour = condition_palette[["disease"]]) +
  scale_colour_manual(values = state_palette, guide = "none") +
  scale_fill_manual(values = state_palette, guide = "none") +
  scale_linewidth(range = c(0.35, 4.4), guide = "none") +
  coord_cartesian(xlim = c(-0.02, 1.08), ylim = c(-0.75, 7.05), clip = "on") +
  labs(
    title = "Cross-sectional state-flow summary",
    subtitle = "Top nonzero flows only; descriptive, not lineage tracing",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, face = "bold", colour = "#111111"),
    plot.subtitle = element_text(size = 7.3, colour = "#555555"),
    plot.margin = margin(8, 8, 14, 16)
  )

figure2_redesign <- (
  (panel_a | panel_b) /
    panel_c /
    panel_d
) +
  plot_layout(heights = c(0.90, 1.08, 1.02), widths = c(0.75, 1.45)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = "#111111"))

fdm_save_plot(figure2_redesign, EXPORTS, STEM, width = 8.8, height = 7.35)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = FIGURE_LABEL,
    primary_claim = "FibroDynMix runs on public raw-count fibroblast data and returns interpretable condition-level composition shifts, cell-level mixture weights, and cross-sectional state-flow summaries.",
    claim_boundary = "Public pooled-count method demonstration; not a disease-mechanism result, clinical biomarker, or lineage-tracing claim.",
    n_cells = n_cells,
    n_genes = n_genes,
    objective_improvement = objective_improvement,
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C", "D"),
    source_data = c(
      "source_data/fig2_condition_cell_counts.tsv; source_data/fig2_nb_objective.tsv",
      "source_data/fig2_state_composition_delta.tsv",
      "source_data/fig2_cell_state_weight_heatmap.tsv; source_data/fig2_cell_order_condition_annotation.tsv",
      "source_data/fig2_transition_flow.tsv; source_data/fig2_transition_summary.tsv"
    ),
    claim = c(
      "A balanced normal/disease public raw-count subset is fitted under the NB model.",
      "Disease-minus-normal composition shifts summarize fitted state-weight changes.",
      "Cell-level latent state weights form a continuous condition-ordered mixture landscape.",
      "Cross-sectional state flow summarizes normal-to-disease composition redistribution without lineage interpretation."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 2. Public raw-count method demonstration",
    "",
    "**A.** Balanced public normal and disease fibroblast raw-count subset and NB objective improvement after fitting.",
    "**B.** Disease-minus-normal fitted latent state-weight differences, ordered by effect size.",
    "**C.** Cell-level latent state-mixture heatmap ordered by condition and dominant state.",
    "**D.** Cross-sectional state-flow summary showing top nonzero normal-to-disease state flows.",
    "",
    "This figure demonstrates real raw-count execution and interpretable FibroDynMix readouts. State flow is cross-sectional and descriptive, not lineage tracing."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Redesigned Figure 2 written to: ", file.path(EXPORTS, paste0(STEM, ".png")))
