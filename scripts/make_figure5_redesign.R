#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure5_redesign.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
SOURCE_IN <- file.path(ROOT, "figures", "figure6", "source_data")
OUT <- Sys.getenv("FDM_FIGURE_OUT", unset = file.path(ROOT, "figures", "figure5_redesign"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure5_redesign")
FIGURE_LABEL <- Sys.getenv("FDM_FIGURE_LABEL", unset = "Figure 5 redesign")
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
  SOURCE_IN,
  c(
    "fig6_transition_cost_long.tsv",
    "fig6_transition_flow_long.tsv",
    "fig6_cell_fpi.tsv",
    "fig6_condition_state_composition.tsv",
    "fig6_transition_flow_summary.tsv"
  )
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop(sprintf("Missing FPI inputs for Figure 5: %s", paste(missing, collapse = ", ")), call. = FALSE)
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
condition_labels <- c(normal = "Normal", disease = "Disease")
condition_palette_soft <- c(normal = "#5D8DB6", disease = "#C66B68")
condition_palette_display <- setNames(condition_palette_soft, condition_labels[names(condition_palette_soft)])
roma_fallback <- c("#4A4E87", "#F7F4EF", "#9A3F47")

fdm_write_tsv(
  data.frame(name = names(condition_palette_soft), colour = unname(condition_palette_soft)),
  file.path(PALETTES, "condition_palette_soft.tsv")
)

cost_long <- fdm_read_tsv(file.path(SOURCE_IN, "fig6_transition_cost_long.tsv"))
flow_long <- fdm_read_tsv(file.path(SOURCE_IN, "fig6_transition_flow_long.tsv"))
cell_fpi <- fdm_read_tsv(file.path(SOURCE_IN, "fig6_cell_fpi.tsv"))
composition <- fdm_read_tsv(file.path(SOURCE_IN, "fig6_condition_state_composition.tsv"))
flow_summary <- fdm_read_tsv(file.path(SOURCE_IN, "fig6_transition_flow_summary.tsv"))

cost_mat <- xtabs(cost ~ source_state + target_state, data = cost_long)
cost_mat <- cost_mat[state_levels, state_levels]
cost_for_cluster <- cost_mat
diag(cost_for_cluster) <- median(cost_mat[row(cost_mat) != col(cost_mat)], na.rm = TRUE)
state_order <- rownames(cost_for_cluster)[stats::hclust(stats::as.dist(cost_for_cluster), method = "average")$order]
state_order <- intersect(state_order, state_levels)

affinity_long <- cost_long
non_diag <- affinity_long$source_state != affinity_long$target_state
center <- median(affinity_long$cost[non_diag], na.rm = TRUE)
scale_cost <- max(abs(affinity_long$cost[non_diag] - center), na.rm = TRUE)
if (!is.finite(scale_cost) || scale_cost == 0) {
  scale_cost <- 1
}
affinity_long$relative_affinity <- (center - affinity_long$cost) / scale_cost
affinity_long$source_idx <- match(affinity_long$source_state, state_order)
affinity_long$target_idx <- match(affinity_long$target_state, state_order)
affinity_long$keep <- affinity_long$source_idx >= affinity_long$target_idx
heatmap_data <- affinity_long[affinity_long$keep, , drop = FALSE]
heatmap_data$source_state <- factor(heatmap_data$source_state, levels = rev(state_order))
heatmap_data$target_state <- factor(heatmap_data$target_state, levels = state_order)
heatmap_data$diag_label <- ifelse(as.character(heatmap_data$source_state) == as.character(heatmap_data$target_state),
  state_labels[as.character(heatmap_data$target_state)],
  ""
)

flow_plot <- flow_long[flow_long$source_state != flow_long$target_state & flow_long$flow > 1e-8, , drop = FALSE]
flow_plot <- flow_plot[order(flow_plot$flow, decreasing = TRUE), , drop = FALSE]
flow_plot <- head(flow_plot, 12L)
flow_plot$x <- 1
flow_plot$xend <- 2
flow_plot$y <- match(flow_plot$source_state, rev(state_order))
flow_plot$yend <- match(flow_plot$target_state, rev(state_order))
flow_plot$source_label <- state_labels[flow_plot$source_state]
flow_plot$target_label <- state_labels[flow_plot$target_state]

composition$state <- factor(composition$state, levels = state_order)
composition$condition <- factor(composition$condition, levels = c("normal", "disease"))

cell_fpi$condition <- factor(cell_fpi$disease, levels = c("normal", "disease"))
cell_fpi$condition_label <- factor(condition_labels[as.character(cell_fpi$condition)], levels = condition_labels)
metric_long <- rbind(
  data.frame(cell_id = cell_fpi$cell_id, condition = cell_fpi$condition, condition_label = cell_fpi$condition_label, metric = "FPI", value = cell_fpi$fpi),
  data.frame(cell_id = cell_fpi$cell_id, condition = cell_fpi$condition, condition_label = cell_fpi$condition_label, metric = "Entropy", value = cell_fpi$entropy),
  data.frame(cell_id = cell_fpi$cell_id, condition = cell_fpi$condition, condition_label = cell_fpi$condition_label, metric = "Transition\npotential", value = cell_fpi$transition_potential)
)
metric_long$metric <- factor(metric_long$metric, levels = c("FPI", "Entropy", "Transition\npotential"))

centroids <- aggregate(cbind(entropy, transition_potential, fpi) ~ condition, data = cell_fpi, FUN = mean)
centroids$condition_label <- factor(condition_labels[as.character(centroids$condition)], levels = condition_labels)
direction <- data.frame(
  x = centroids$entropy[centroids$condition == "normal"],
  y = centroids$transition_potential[centroids$condition == "normal"],
  xend = centroids$entropy[centroids$condition == "disease"],
  yend = centroids$transition_potential[centroids$condition == "disease"]
)

fdm_write_tsv(heatmap_data, file.path(SOURCE, "fig5_state_relative_affinity_lower_triangle.tsv"))
fdm_write_tsv(flow_plot, file.path(SOURCE, "fig5_transition_flow_curves.tsv"))
fdm_write_tsv(metric_long, file.path(SOURCE, "fig5_fpi_component_raincloud.tsv"))
fdm_write_tsv(cell_fpi, file.path(SOURCE, "fig5_fpi_decomposition_cells.tsv"))
fdm_write_tsv(centroids, file.path(SOURCE, "fig5_fpi_decomposition_centroids.tsv"))
fdm_write_tsv(flow_summary, file.path(SOURCE, "fig5_transition_flow_summary.tsv"))

theme_set(fdm_theme(base_size = 8.5))

panel_a <- ggplot(heatmap_data, aes(x = target_state, y = source_state, fill = relative_affinity)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = diag_label), size = 2.35, fontface = "bold", colour = "#222222") +
  scale_x_discrete(labels = state_labels) +
  scale_y_discrete(labels = state_labels) +
  labs(
    tag = "A",
    title = "State landscape ordered by transcriptomic proximity",
    subtitle = "Lower triangle removes matrix redundancy; diagonal carries state labels",
    x = NULL,
    y = NULL,
    fill = "Relative\naffinity"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.border = element_blank(),
    legend.position = "right"
  )
if (requireNamespace("scico", quietly = TRUE)) {
  panel_a <- panel_a + scico::scale_fill_scico(palette = "roma", limits = c(-1, 1), oob = scales::squish)
} else {
  panel_a <- panel_a + scale_fill_gradient2(
    low = roma_fallback[1],
    mid = roma_fallback[2],
    high = roma_fallback[3],
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  )
}

node_data <- data.frame(
  x = rep(c(1, 2), each = length(state_order)),
  y = rep(seq_along(rev(state_order)), times = 2),
  state = rep(rev(state_order), times = 2),
  side = rep(c("Normal mix", "Disease mix"), each = length(state_order)),
  stringsAsFactors = FALSE
)
node_data <- merge(
  node_data,
  composition[, c("state", "condition", "composition")],
  by = "state",
  all.x = TRUE
)
node_data <- node_data[
  (node_data$x == 1 & node_data$condition == "normal") |
    (node_data$x == 2 & node_data$condition == "disease"),
  ,
  drop = FALSE
]
node_data$side <- factor(node_data$side, levels = c("Normal mix", "Disease mix"))

panel_b <- ggplot() +
  geom_curve(
    data = flow_plot,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = flow),
    curvature = 0.28,
    colour = "#3E6C88",
    alpha = 0.58,
    lineend = "round"
  ) +
  geom_point(
    data = node_data,
    aes(x = x, y = y, size = composition, fill = side),
    shape = 21,
    colour = "white",
    stroke = 0.35
  ) +
  geom_text(
    data = node_data[node_data$x == 1, , drop = FALSE],
    aes(x = x - 0.05, y = y, label = state_labels[state]),
    hjust = 1,
    size = 2.65,
    colour = "#222222"
  ) +
  geom_text(
    data = node_data[node_data$x == 2, , drop = FALSE],
    aes(x = x + 0.05, y = y, label = state_labels[state]),
    hjust = 0,
    size = 2.65,
    colour = "#222222"
  ) +
  scale_fill_manual(values = c("Normal mix" = condition_palette_soft[["normal"]], "Disease mix" = condition_palette_soft[["disease"]]), name = NULL) +
  scale_linewidth_continuous(range = c(0.35, 3.0), guide = "none") +
  scale_size_continuous(range = c(2.1, 7.2), guide = "none") +
  coord_cartesian(xlim = c(0.54, 2.46), ylim = c(0.45, length(state_order) + 0.55), clip = "off") +
  labs(
    tag = "B",
    title = "Dominant cross-condition transition flow",
    subtitle = sprintf("Expected cost %.2f; OT entropy %.2f", flow_summary$expected_cost[[1]], flow_summary$entropy[[1]]),
    x = NULL,
    y = NULL
  ) +
  fdm_blank_theme(base_size = 8.5) +
  theme(legend.position = "top")

panel_c <- ggplot(metric_long, aes(x = condition_label, y = value, fill = condition_label, colour = condition_label)) +
  geom_violin(width = 0.78, alpha = 0.42, linewidth = 0.28, trim = FALSE) +
  ggbeeswarm::geom_quasirandom(width = 0.18, size = 0.65, alpha = 0.48, stroke = 0) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.78, linewidth = 0.28, colour = "#222222") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = condition_palette_display, name = NULL) +
  scale_colour_manual(values = condition_palette_display, guide = "none") +
  labs(
    tag = "C",
    title = "FPI separates entropy from transition potential",
    subtitle = "Violin, box, and individual cells show distributional structure",
    x = NULL,
    y = "Cell-level score"
  ) +
  theme(legend.position = "none")

panel_d <- ggplot(cell_fpi, aes(x = entropy, y = transition_potential)) +
  geom_abline(intercept = seq(0.35, 1.25, by = 0.15), slope = -1, linewidth = 0.2, colour = "#D5D5D5") +
  geom_point(aes(fill = condition_label), shape = 21, size = 1.75, alpha = 0.72, colour = "white", stroke = 0.12) +
  geom_point(
    data = centroids,
    aes(x = entropy, y = transition_potential, fill = condition_label),
    shape = 21,
    size = 4.2,
    colour = "#111111",
    stroke = 0.42
  ) +
  geom_segment(
    data = direction,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.07, "in"), type = "closed"),
    linewidth = 0.55,
    colour = "#111111"
  ) +
  annotate("text", x = direction$xend - 0.13, y = direction$yend + 0.07, label = "Disease shift:\ntransition potential loss", hjust = 0, vjust = 0.5, size = 2.55, colour = "#333333") +
  scale_fill_manual(values = condition_palette_display, name = NULL) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "D",
    title = "Low disease FPI is driven by transition-potential collapse",
    subtitle = "Grey diagonals mark constant FPI = entropy + transition potential",
    x = "Entropy",
    y = "Transition potential"
  ) +
  theme(legend.position = "top")

figure5 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(widths = c(1.02, 0.98), heights = c(0.92, 1.08)) &
  fdm_panel_tag_theme()

fdm_save_plot(figure5, EXPORTS, STEM, width = 8.6, height = 6.9, dpi = 600)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = FIGURE_LABEL,
    primary_claim = "FPI decomposes cell-level fibroblast plasticity into mixed-state entropy and transition potential; disease cells show lower FPI primarily through reduced transition potential.",
    claim_boundary = "Transition flow and FPI are cross-sectional model-derived summaries, not lineage tracing or causal dynamics.",
    transition_source = "figures/figure6/source_data",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = c("A", "B", "C", "D"),
    source_data = c(
      "source_data/fig5_state_relative_affinity_lower_triangle.tsv",
      "source_data/fig5_transition_flow_curves.tsv",
      "source_data/fig5_fpi_component_raincloud.tsv",
      "source_data/fig5_fpi_decomposition_cells.tsv"
    ),
    claim = c(
      "State relationships are shown as a clustered lower-triangle relative-affinity heatmap to reduce redundant matrix encoding.",
      "The strongest cross-condition transition-flow edges summarize model-derived state redistribution.",
      "FPI, entropy, and transition-potential distributions separate total plasticity from its components.",
      "Disease cells have lower transition potential at comparable entropy, explaining lower FPI."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 5. FPI decomposition reveals transition-potential loss in disease fibroblasts",
    "",
    "**A.** Clustered lower-triangle state-affinity heatmap derived from the transition cost matrix; diagonal labels identify fibroblast states.",
    "**B.** Curved two-column transition-flow summary from normal to disease state composition; node size reflects condition-level state composition.",
    "**C.** Raincloud distributions of cell-level FPI, entropy, and transition potential by condition.",
    "**D.** FPI decomposition scatter. Grey diagonals denote constant FPI, large points show condition centroids, and the arrow marks the normal-to-disease centroid shift.",
    "",
    "FPI is interpreted as a bounded cross-sectional readout: FPI = normalized state entropy + transition potential. The observed disease shift is driven mainly by reduced transition potential rather than reduced entropy."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Figure 5 FPI redesign written to: ", file.path(EXPORTS, paste0(STEM, ".png")))
