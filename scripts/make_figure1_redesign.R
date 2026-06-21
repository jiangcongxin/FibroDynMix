#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure1_redesign.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- Sys.getenv("FDM_FIGURE_OUT", unset = file.path(ROOT, "figures", "figure1_redesign"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", unset = "figure1_redesign")
FIGURE_LABEL <- Sys.getenv("FDM_FIGURE_LABEL", unset = "Figure 1 redesign")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
QC <- file.path(OUT, "qc")

dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "figure_style_fdm.R"))

blue <- "#4E79A7"
red <- "#E15759"
teal <- "#59A14F"
purple <- "#B07AA1"
orange <- "#F28E2B"
ink <- "#1F2937"
muted <- "#6B7280"
line <- "#374151"
soft_blue <- "#EEF3F9"
soft_red <- "#FBEDEE"
soft_green <- "#ECF5EA"
soft_gold <- "#FFF4D8"
soft_gray <- "#F4F5F7"

box_data <- data.frame(
  id = c("counts", "likelihood", "simplex", "outputs", "study", "markers"),
  label = c(
    "Raw UMI counts",
    "NB likelihood\nlibrary-size offset",
    "Latent simplex\nstate weights",
    "Composition\nUncertainty\nTransfer\nFlow / FPI",
    "Study / donor\nadjustment",
    "Weak marker\norientation"
  ),
  x = c(0.13, 0.36, 0.59, 0.84, 0.36, 0.59),
  y = c(0.55, 0.55, 0.55, 0.55, 0.29, 0.29),
  w = c(0.17, 0.18, 0.18, 0.19, 0.18, 0.18),
  h = c(0.23, 0.23, 0.23, 0.28, 0.16, 0.16),
  fill = c(soft_blue, soft_green, soft_gold, soft_gray, "#F0E7FA", "#FDE7C7"),
  stringsAsFactors = FALSE
)
box_data$label_x <- box_data$x
box_data$label_y <- box_data$y
box_data$label_y[box_data$id == "simplex"] <- 0.472
box_data$label_x[box_data$id == "outputs"] <- 0.815

edge_data <- data.frame(
  x = c(0.215, 0.45, 0.68, 0.36, 0.59, 0.59),
  xend = c(0.27, 0.50, 0.745, 0.36, 0.59, 0.745),
  y = c(0.55, 0.55, 0.55, 0.37, 0.37, 0.35),
  yend = c(0.55, 0.55, 0.55, 0.435, 0.435, 0.47),
  type = c("main", "main", "main", "support", "support", "support"),
  stringsAsFactors = FALSE
)

matrix_tiles <- expand.grid(i = seq_len(6), j = seq_len(8))
matrix_tiles$x <- 0.065 + matrix_tiles$j * 0.009
matrix_tiles$y <- 0.485 + matrix_tiles$i * 0.014
matrix_tiles$value <- c(
  0.88, 0.23, 0.11, 0.62, 0.05, 0.33, 0.42, 0.18,
  0.14, 0.71, 0.09, 0.25, 0.54, 0.06, 0.29, 0.82,
  0.35, 0.18, 0.72, 0.10, 0.21, 0.58, 0.11, 0.44,
  0.06, 0.41, 0.17, 0.81, 0.15, 0.24, 0.66, 0.12,
  0.51, 0.09, 0.31, 0.22, 0.78, 0.14, 0.19, 0.57,
  0.17, 0.68, 0.13, 0.46, 0.08, 0.73, 0.28, 0.16
)

simplex_bars <- data.frame(
  cell = rep(seq_len(5), each = 4),
  state = rep(c("resident", "inflamm.", "myoFB", "ECM"), times = 5),
  fraction = c(
    0.70, 0.12, 0.10, 0.08,
    0.28, 0.46, 0.14, 0.12,
    0.20, 0.14, 0.48, 0.18,
    0.10, 0.08, 0.24, 0.58,
    0.35, 0.10, 0.18, 0.37
  ),
  stringsAsFactors = FALSE
)
simplex_bars$xmin <- rep(seq(0.535, 0.635, length.out = 5), each = 4)
simplex_bars$xmax <- simplex_bars$xmin + 0.012
simplex_bars$ymin <- ave(simplex_bars$fraction, simplex_bars$cell, FUN = function(v) c(0, head(cumsum(v), -1)))
simplex_bars$ymax <- ave(simplex_bars$fraction, simplex_bars$cell, FUN = cumsum)
simplex_bars$ymin <- 0.545 + simplex_bars$ymin * 0.105
simplex_bars$ymax <- 0.545 + simplex_bars$ymax * 0.105
simplex_bars$fill <- c(
  resident = blue,
  "inflamm." = red,
  myoFB = teal,
  ECM = purple
)[simplex_bars$state]

output_points <- data.frame(
  x = c(0.895, 0.865, 0.925, 0.925),
  y = c(0.455, 0.630, 0.520, 0.675),
  fill = c(blue, red, teal, purple),
  stringsAsFactors = FALSE
)

visible_text <- c(
  "FibroDynMix: raw-count latent mixtures for fibroblast plasticity",
  "Counts are modeled directly; cell states are continuous; outputs stay cross-sectional.",
  box_data$label,
  "Cross-sectional summaries, not lineage tracing"
)
visible_word_count <- length(unlist(strsplit(gsub("[\n/;:.-]", " ", paste(visible_text, collapse = " ")), "\\s+")))

figure <- ggplot() +
  annotate("text", x = 0.05, y = 0.91, label = "FibroDynMix: raw-count latent mixtures for fibroblast plasticity", hjust = 0, size = 4.3, fontface = "bold", colour = ink) +
  annotate("text", x = 0.05, y = 0.845, label = "Counts are modeled directly; cell states are continuous; outputs stay cross-sectional.", hjust = 0, size = 2.7, colour = muted) +
  geom_rect(
    data = box_data,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill),
    colour = "#2F2F2F",
    linewidth = 0.28
  ) +
  geom_segment(
    data = edge_data,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = type),
    linewidth = 0.38,
    colour = line,
    arrow = arrow(length = unit(0.07, "in"), type = "closed")
  ) +
  geom_text(data = box_data, aes(x = label_x, y = label_y, label = label), size = 2.65, lineheight = 0.9, fontface = "bold", colour = ink) +
  geom_tile(
    data = matrix_tiles,
    aes(x = x, y = y, fill = value),
    width = 0.0075,
    height = 0.011,
    colour = "white",
    linewidth = 0.05
  ) +
  geom_rect(
    data = simplex_bars,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "white",
    linewidth = 0.12
  ) +
  geom_curve(
    aes(x = 0.87, xend = 0.93, y = 0.45, yend = 0.67),
    curvature = -0.35,
    linewidth = 1.2,
    colour = "#BFC7D5",
    alpha = 0.8,
    arrow = arrow(length = unit(0.06, "in"), type = "closed")
  ) +
  geom_point(data = output_points, aes(x = x, y = y, fill = fill), shape = 21, size = 3.0, colour = "white", stroke = 0.35) +
  annotate("rect", xmin = 0.72, xmax = 0.96, ymin = 0.12, ymax = 0.20, fill = soft_red, colour = "#F0C9CA", linewidth = 0.22) +
  annotate("text", x = 0.84, y = 0.16, label = "Cross-sectional summaries, not lineage tracing", size = 2.45, fontface = "bold", colour = "#7A2E2E") +
  annotate("text", x = 0.13, y = 0.19, label = "data", size = 2.1, colour = muted, fontface = "bold") +
  annotate("text", x = 0.36, y = 0.19, label = "count model", size = 2.1, colour = muted, fontface = "bold") +
  annotate("text", x = 0.59, y = 0.19, label = "state inference", size = 2.1, colour = muted, fontface = "bold") +
  annotate("text", x = 0.84, y = 0.25, label = "bounded readouts", size = 2.1, colour = muted, fontface = "bold") +
  scale_fill_identity() +
  scale_linetype_manual(values = c(main = "solid", support = "22"), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 0.96), clip = "off") +
  theme_void(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(12, 12, 12, 12)
  )

fdm_save_plot(figure, EXPORTS, STEM, width = 8.8, height = 4.2, dpi = 600)
fdm_write_export_qc(EXPORTS, QC, STEM)

fdm_write_tsv(
  data.frame(
    figure = FIGURE_LABEL,
    visible_word_count = visible_word_count,
    primary_claim = "Full-width graphical abstract for FibroDynMix as raw-count latent mixture inference with bounded cross-sectional readouts.",
    claim_boundary = "Graphical summary only; no causal transition, clinical biomarker, or complete disease-atlas claim.",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv")
)

fdm_write_tsv(
  box_data[, c("id", "label", "x", "y", "w", "h")],
  file.path(SOURCE, "fig1_redesign_components.tsv")
)

fdm_write_tsv(
  data.frame(
    panel = "A",
    source_data = "source_data/fig1_redesign_components.tsv",
    claim = "FibroDynMix links raw UMI counts to a negative-binomial likelihood, latent simplex state weights, study/donor adjustment, marker orientation, and bounded cross-sectional readouts.",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "panel_source_data_manifest.tsv")
)

writeLines(
  c(
    "# Figure 1. FibroDynMix graphical summary",
    "",
    "**A.** Full-width graphical summary of the FibroDynMix workflow. Raw UMI counts enter a negative-binomial likelihood with library-size offset, producing continuous latent simplex state weights. Study/donor adjustment and weak marker orientation support state inference, while composition, uncertainty, transfer, flow, and FPI are reported as bounded cross-sectional readouts.",
    "",
    "The figure summarizes model architecture only. It does not claim lineage tracing, clinical biomarker readiness, or a completed disease atlas."
  ),
  file.path(OUT, "main_figure_legends.md")
)

message("Figure 1 redesign written to: ", file.path(EXPORTS, paste0(STEM, ".png")))
