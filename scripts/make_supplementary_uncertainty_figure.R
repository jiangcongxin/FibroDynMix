#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_supplementary_uncertainty_figure.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "figures", "supplementary_uncertainty")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
PALETTES <- file.path(OUT, "palettes")
QC <- file.path(OUT, "qc")

dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)
dir.create(PALETTES, recursive = TRUE, showWarnings = FALSE)
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "R", "matrix_utils.R"))
source(file.path(ROOT, "R", "simulate_fibrodynmix.R"))
source(file.path(ROOT, "R", "benchmark_metrics.R"))
source(file.path(ROOT, "R", "baseline_marker_scoring.R"))
source(file.path(ROOT, "R", "fibrodynmix_initializer.R"))
source(file.path(ROOT, "R", "nb_likelihood.R"))
source(file.path(ROOT, "R", "fit_nb_model.R"))
source(file.path(ROOT, "R", "bootstrap_uncertainty.R"))

theme_set(theme_classic(base_size = 8, base_family = "Helvetica"))
theme_update(
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 7, colour = "black"),
  plot.title = element_text(size = 9, face = "bold", hjust = 0),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  strip.background = element_blank(),
  strip.text = element_text(size = 8, face = "bold")
)

state_palette <- c(
  resident = "#3B6EA8",
  inflammatory = "#B84A4A",
  myofibroblast = "#4F8B5B",
  "ECM-remodeling" = "#9A6B32",
  "antigen-presenting" = "#7A5EA8",
  "IFN-stress" = "#5A9AA8"
)
write.table(
  data.frame(name = names(state_palette), colour = unname(state_palette)),
  file.path(PALETTES, "supplementary_uncertainty_state_palette.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sim <- simulate_fibrodynmix(
  n_studies = 2,
  donors_per_study = 2,
  cells_per_donor = 10,
  n_genes = 100,
  marker_genes_per_state = 5,
  scenario = "batch_confounding",
  seed = 520
)

boot <- bootstrap_fibrodynmix(
  counts = sim$counts,
  marker_index = sim$parameters$marker_index,
  library_size = sim$cell_metadata$library_size,
  cell_metadata = sim$cell_metadata,
  sample_col = "donor_id",
  method = "nb_study",
  n_boot = 5,
  seed = 521,
  fit_args = list(
    n_outer = 1,
    initializer_args = list(n_iter = 2),
    study_l2 = 5,
    marker_l2 = 0.05,
    maxit_beta = 10,
    maxit_z = 8
  )
)

sample_summary <- boot$sample_summary
cell_summary <- boot$cell_summary
marker_summary <- boot$marker_summary
cell_draws <- boot$cell_draws

entropy_by_cell <- aggregate(entropy ~ cell_id + replicate, data = unique(cell_draws[, c("cell_id", "replicate", "entropy")]), FUN = mean)
entropy_summary <- summarize_draws(
  draws = transform(entropy_by_cell, state = "entropy"),
  group_cols = c("cell_id", "state"),
  value_col = "entropy",
  probs = c(0.025, 0.975)
)
entropy_summary$interval_width <- entropy_summary$upper - entropy_summary$lower
entropy_summary$mean_entropy <- entropy_summary$mean

z_width <- cell_summary
z_width$interval_width <- z_width$upper - z_width$lower
cell_uncertainty <- aggregate(interval_width ~ cell_id, data = z_width, FUN = mean)
cell_uncertainty <- merge(
  cell_uncertainty,
  entropy_summary[, c("cell_id", "mean_entropy")],
  by = "cell_id",
  all.x = TRUE
)

top_markers <- marker_summary[order(marker_summary$state, -marker_summary$mean), ]
top_markers <- do.call(rbind, lapply(split(top_markers, top_markers$state), function(df) head(df, 6)))

write.table(sample_summary, file.path(SOURCE, "fig5_sample_composition_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cell_summary, file.path(SOURCE, "fig5_cell_state_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(entropy_summary, file.path(SOURCE, "fig5_entropy_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cell_uncertainty, file.path(SOURCE, "fig5_uncertainty_vs_entropy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_markers, file.path(SOURCE, "fig5_marker_program_stability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

panel_a_data <- sample_summary[sample_summary$state %in% c("resident", "myofibroblast", "inflammatory"), ]
panel_a <- ggplot(panel_a_data, aes(x = sample_id, y = mean, ymin = lower, ymax = upper, colour = state)) +
  geom_pointrange(position = position_dodge(width = 0.6), linewidth = 0.35, size = 1.2) +
  scale_colour_manual(values = state_palette) +
  coord_cartesian(ylim = c(0, max(panel_a_data$upper, na.rm = TRUE) * 1.1)) +
  labs(title = "Sample-level state composition intervals", x = NULL, y = "Bootstrap mean and 95% interval", colour = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

panel_b <- ggplot(entropy_summary, aes(x = mean_entropy, y = interval_width)) +
  geom_point(size = 1.4, alpha = 0.75, colour = "#2F5E9E") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45, colour = "#555555") +
  labs(title = "Entropy uncertainty by cell", x = "Mean entropy", y = "95% interval width") +
  theme(legend.position = "none")

panel_c <- ggplot(top_markers, aes(x = reorder(gene, mean), y = mean, ymin = lower, ymax = upper, colour = state)) +
  geom_pointrange(linewidth = 0.3, size = 0.9) +
  facet_wrap(~ state, scales = "free_y", ncol = 2) +
  coord_flip() +
  scale_colour_manual(values = state_palette) +
  labs(title = "Marker-program bootstrap stability", x = NULL, y = "|beta| mean and 95% interval", colour = NULL) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6))

panel_d <- ggplot(cell_uncertainty, aes(x = mean_entropy, y = interval_width)) +
  geom_point(size = 1.4, alpha = 0.75, colour = "#287C6F") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45, colour = "#555555") +
  labs(title = "State-weight uncertainty across mixed states", x = "Mean entropy", y = "Mean z interval width") +
  theme(legend.position = "none")

supplementary_uncertainty <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold"))

save_plot <- function(plot, stem, width = 7.2, height = 6.0, dpi = 450) {
  ggsave(file.path(EXPORTS, paste0(stem, ".pdf")), plot, width = width, height = height, device = grDevices::pdf)
  ggsave(file.path(EXPORTS, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  ggsave(file.path(EXPORTS, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, device = ragg::agg_png)
  ggsave(file.path(EXPORTS, paste0(stem, ".tiff")), plot, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
}

save_plot(supplementary_uncertainty, "supplementary_uncertainty")

write.table(
  data.frame(
    figure = "Supplementary uncertainty figure",
    claim_boundary = "Bootstrap uncertainty summaries from simulated data; not full Bayesian posterior credible intervals.",
    primary_claim = "Cell bootstrap provides uncertainty summaries for state composition, cell entropy, and marker-program stability.",
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "figure_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

panel_manifest <- data.frame(
  panel = c("A", "B", "C", "D"),
  source_data = c(
    "source_data/fig5_sample_composition_uncertainty.tsv",
    "source_data/fig5_entropy_uncertainty.tsv",
    "source_data/fig5_marker_program_stability.tsv",
    "source_data/fig5_uncertainty_vs_entropy.tsv"
  ),
  claim = c(
    "Sample-level fibroblast state composition can be summarized with bootstrap intervals.",
    "Cells with higher entropy show quantifiable uncertainty in entropy estimates.",
    "State-gene program stability is summarized by bootstrap intervals on |beta|.",
    "State-weight interval width is evaluated against mixed-state entropy."
  ),
  stringsAsFactors = FALSE
)
write.table(panel_manifest, file.path(OUT, "panel_source_data_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

legend_text <- c(
  "# Supplementary Figure. Bootstrap uncertainty summaries for FibroDynMix state inference",
  "",
  "**A.** Sample-level bootstrap intervals for resident, inflammatory, and myofibroblast composition.",
  "**B.** Cell-level entropy interval width as a function of mean entropy.",
  "**C.** Bootstrap stability of learned state-gene programs, shown for the top six genes per state by mean absolute beta.",
  "**D.** Mean state-weight interval width evaluated against mixed-state entropy.",
  "",
  "Intervals are cell-bootstrap uncertainty summaries from simulated data. They are not full Bayesian posterior credible intervals."
)
writeLines(legend_text, file.path(OUT, "main_figure_legends.md"))

if (requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(file.path(EXPORTS, "supplementary_uncertainty.png"))
  magick::image_write(img, file.path(EXPORTS, "contact_sheet.png"))
}

qc_files <- file.path(EXPORTS, c("supplementary_uncertainty.pdf", "supplementary_uncertainty.svg", "supplementary_uncertainty.png", "supplementary_uncertainty.tiff", "contact_sheet.png"))
qc <- data.frame(
  file = basename(qc_files),
  exists = file.exists(qc_files),
  bytes = ifelse(file.exists(qc_files), file.info(qc_files)$size, NA_real_),
  stringsAsFactors = FALSE
)
if (requireNamespace("magick", quietly = TRUE)) {
  dims <- lapply(qc_files, function(path) {
    if (!file.exists(path) || !grepl("\\.(png|tiff)$", path)) {
      return(c(width = NA_real_, height = NA_real_))
    }
    info <- magick::image_info(magick::image_read(path))
    c(width = info$width[1], height = info$height[1])
  })
  dims <- do.call(rbind, dims)
  qc$width <- dims[, "width"]
  qc$height <- dims[, "height"]
}
write.table(qc, file.path(QC, "export_image_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Supplementary uncertainty figure written to: ", OUT)
