#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure3_incremental_utility.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- Sys.getenv("FDM_FIGURE_OUT", file.path(ROOT, "paper_rewriting_output", "figure_rebuild_iscience_marker_anchored_v20", "main", "Figure_3"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", "Figure_3")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)

replicate_path <- file.path(ROOT, "analysis", "incremental_utility_benchmark", "simulation_replicate_means.tsv")
coordinate_path <- file.path(ROOT, "analysis", "incremental_utility_benchmark", "state_coordinate_metrics.tsv")
split_path <- file.path(ROOT, "analysis", "incremental_utility_benchmark", "marker_split_variability.tsv")
for (path in c(replicate_path, coordinate_path, split_path)) if (!file.exists(path)) stop("Missing input: ", path)
simulation <- read.delim(replicate_path, check.names = FALSE)
coordinates <- read.delim(coordinate_path, check.names = FALSE)
splits <- read.delim(split_path, check.names = FALSE)
if (nrow(simulation) != 20L) stop("Figure 3 requires 20 independent simulation summaries.")

write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
write_tsv(simulation, file.path(SOURCE, "fig3_simulation_replicate_means.tsv"))
write_tsv(coordinates, file.path(SOURCE, "fig3_state_coordinate_metrics.tsv"))
write_tsv(splits, file.path(SOURCE, "fig3_marker_split_variability.tsv"))

ink <- "#202124"; blue <- "#2F6690"; teal <- "#3A7D78"; orange <- "#C46A32"
light_blue <- "#D9E6F2"; light_teal <- "#D8EAE7"
panel_label <- function(label) { u <- par("usr"); text(u[1], u[4], label, adj = c(-0.55, 1.2), xpd = NA, font = 2, cex = 1.2, col = ink) }

draw <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  layout(matrix(1:4, 2, byrow = TRUE))
  par(family = "sans", fg = ink, col.axis = ink, col.lab = ink, col.main = ink,
      mar = c(5.4, 5.2, 3.2, 1.2), mgp = c(3.05, 0.8, 0), tcl = -0.25)

  ll <- list(simulation$null_heldout_gene_loglik_per_observation,
             simulation$anchored_no_source_heldout_gene_loglik_per_observation,
             simulation$anchored_source_heldout_gene_loglik_per_observation)
  boxplot(ll, names = c("Intercept-only", "Anchored NB", "Anchored NB\n+ study"),
          col = c("#E5E5E5", light_blue, light_teal), border = ink, outline = FALSE,
          ylab = "Held-out log likelihood / observation", main = "Held-out count prediction")
  set.seed(3201); for (j in seq_along(ll)) points(jitter(rep(j, length(ll[[j]])), 0.08), ll[[j]], pch = 16, cex = 0.55, col = adjustcolor(ink, 0.55))
  panel_label("A")

  gains <- list(simulation$source_gain_over_null, simulation$source_gain_over_no_source)
  boxplot(gains, names = c("vs intercept", "vs no-source"), col = c(light_blue, light_teal),
          border = ink, outline = FALSE, ylab = "Source-aware gain / observation",
          main = "Incremental count-model utility")
  abline(h = 0, lty = 2, col = "#6B6B6B")
  set.seed(3202); for (j in seq_along(gains)) points(jitter(rep(j, length(gains[[j]])), 0.08), gains[[j]], pch = 16, cex = 0.55, col = adjustcolor(blue, 0.65))
  panel_label("B")

  coord_mean <- aggregate(cbind(state_weight_rmse, dominant_accuracy, drift_from_marker) ~ method, coordinates, mean)
  method_order <- c("marker_scoring", "anchored_no_source", "anchored_source_aware")
  coord_mean <- coord_mean[match(method_order, coord_mean$method), ]
  bp <- barplot(coord_mean$state_weight_rmse, names.arg = c("Marker score", "Anchored NB", "Anchored NB\n+ study"),
          col = c("#E5E5E5", light_blue, light_teal), border = ink,
          ylab = "State-weight RMSE", ylim = c(0, max(coord_mean$state_weight_rmse) * 1.22),
          main = "Coordinate recovery is retained, not improved")
  text(bp, coord_mean$state_weight_rmse,
       labels = sprintf("%.4f", coord_mean$state_weight_rmse), pos = 3, cex = 0.8)
  panel_label("C")

  yr <- range(c(simulation$marker_downsample_mean_abs_delta, simulation$anchored_projection_downsample_mean_abs_delta))
  plot(c(1, 2), yr, type = "n", xaxt = "n", xlab = "", ylab = "Mean absolute change after 50% thinning",
       main = "No downsampling-stability advantage")
  axis(1, at = c(1, 2), labels = c("Marker\nscore", "Anchored\nprojection"), cex.axis = 0.88)
  for (i in seq_len(nrow(simulation))) lines(c(1, 2), c(simulation$marker_downsample_mean_abs_delta[i], simulation$anchored_projection_downsample_mean_abs_delta[i]), col = adjustcolor("#7A7A7A", 0.45), lwd = 0.8)
  points(rep(1, nrow(simulation)), simulation$marker_downsample_mean_abs_delta, pch = 16, col = blue, cex = 0.65)
  points(rep(2, nrow(simulation)), simulation$anchored_projection_downsample_mean_abs_delta, pch = 16, col = orange, cex = 0.65)
  panel_label("D")
}

ragg::agg_png(file.path(EXPORTS, paste0(STEM, ".png")), 4800, 3600, units = "px", res = 400); draw(); dev.off()
ragg::agg_tiff(file.path(EXPORTS, paste0(STEM, ".tiff")), 4800, 3600, units = "px", res = 400, compression = "lzw"); draw(); dev.off()
pdf(file.path(EXPORTS, paste0(STEM, ".pdf")), 12, 9, useDingbats = FALSE); draw(); dev.off()
svglite::svglite(file.path(EXPORTS, paste0(STEM, ".svg")), 12, 9); draw(); dev.off()

writeLines(c(
  "# Figure 3. Anchored count modeling adds held-out prediction around a retained marker coordinate", "",
  "**(A)** Held-out state-associated-gene log likelihood for the intercept-only model, anchored NB without source terms, and source-aware anchored NB across 20 independent simulations; each point is the mean of four random marker splits. **(B)** Within-simulation gains of source-aware anchored NB over the intercept-only and no-source models; both contrasts were favorable in 19 of 20 simulations. **(C)** Mean state-weight RMSE for direct marker scoring and the two anchored models. Anchored fitting retained rather than improved the marker-defined coordinate. **(D)** Mean absolute coordinate change after 50% binomial thinning for direct marker scoring and anchored fixed-program projection; anchored projection was more stable in 0 of 20 simulations. Complementary markers were excluded from coordinate construction and used only for held-out count evaluation."
), file.path(OUT, "main_figure_legend.md"))
write_tsv(data.frame(panel = c("A", "B", "C", "D"), source_data = c(
  "source_data/fig3_simulation_replicate_means.tsv", "source_data/fig3_simulation_replicate_means.tsv",
  "source_data/fig3_state_coordinate_metrics.tsv", "source_data/fig3_simulation_replicate_means.tsv"),
  description = c("Held-out count likelihood", "Paired source-aware gains", "Coordinate-recovery comparison", "Paired thinning stability")),
  file.path(OUT, "panel_source_data_manifest.tsv"))
write_tsv(data.frame(figure = "Figure 3", title = "Anchored count modeling adds held-out prediction around a retained marker coordinate", panels = "A-D", exports = paste0("exports/", STEM, ".pdf; exports/", STEM, ".svg; exports/", STEM, ".png; exports/", STEM, ".tiff"), source_data_complete = TRUE, source_manifest = "panel_source_data_manifest.tsv"), file.path(OUT, "figure_manifest.tsv"))
