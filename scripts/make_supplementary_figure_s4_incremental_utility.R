#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/make_supplementary_figure_s4_incremental_utility.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

OUT <- Sys.getenv(
  "FDM_SUPPLEMENTARY_FIGURE_OUT",
  file.path(ROOT, "paper_rewriting_output", "figure_rebuild_iscience_identifiability_v18", "supplementary", "Figure_S4")
)
STEM <- Sys.getenv("FDM_SUPPLEMENTARY_FIGURE_STEM", "Figure_S4")
EXPORTS <- file.path(OUT, "exports")
SOURCE <- file.path(OUT, "source_data")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)

simulation_path <- file.path(ROOT, "analysis", "incremental_utility_benchmark", "simulation_replicate_means.tsv")
realdata_path <- file.path(ROOT, "analysis", "realdata_heldout_count_validation", "realdata_holdout_group_means.tsv")
for (path in c(simulation_path, realdata_path)) {
  if (!file.exists(path)) stop("Required Figure S4 input is missing: ", path, call. = FALSE)
}

simulation <- read.delim(simulation_path, check.names = FALSE, stringsAsFactors = FALSE)
realdata <- read.delim(realdata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(simulation) != 20L) stop("Figure S4 requires 20 independent simulation summaries.", call. = FALSE)
if (nrow(realdata) != 7L) stop("Figure S4 requires seven held-out public-data groups.", call. = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_tsv(simulation, file.path(SOURCE, "figs4_simulation_replicate_means.tsv"))
write_tsv(realdata, file.path(SOURCE, "figs4_realdata_holdout_group_means.tsv"))

ink <- "#202124"
blue <- "#2F6690"
teal <- "#3A7D78"
orange <- "#C46A32"
light_blue <- "#D9E6F2"
light_teal <- "#D8EAE7"
light_orange <- "#F2DED0"

panel_label <- function(label) {
  usr <- par("usr")
  text(usr[1], usr[4], label, adj = c(-0.6, 1.2), xpd = NA, font = 2, cex = 1.2, col = ink)
}

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))
  par(family = "sans", fg = ink, col.axis = ink, col.lab = ink, col.main = ink,
      mar = c(5.6, 5.2, 3.2, 1.2), mgp = c(3.1, 0.8, 0), tcl = -0.25)

  ll <- list(
    simulation$null_heldout_gene_loglik_per_observation,
    simulation$anchored_no_source_heldout_gene_loglik_per_observation,
    simulation$anchored_source_heldout_gene_loglik_per_observation
  )
  boxplot(ll, names = c("Intercept-only", "Anchored NB", "Anchored NB\n+ study"),
          col = c("#E5E5E5", light_blue, light_teal), border = ink, outline = FALSE,
          ylab = "Held-out log likelihood / observation", main = "Held-out count prediction")
  set.seed(2401)
  for (j in seq_along(ll)) points(jitter(rep(j, length(ll[[j]])), amount = 0.08), ll[[j]], pch = 16, cex = 0.55, col = adjustcolor(ink, 0.55))
  panel_label("A")

  gains <- list(simulation$source_gain_over_null, simulation$source_gain_over_no_source)
  boxplot(gains, names = c("vs intercept", "vs no-source"), col = c(light_blue, light_teal),
          border = ink, outline = FALSE, ylab = "Source-aware gain / observation",
          main = "Incremental count-model utility")
  abline(h = 0, lty = 2, col = "#6B6B6B")
  set.seed(2402)
  for (j in seq_along(gains)) points(jitter(rep(j, length(gains[[j]])), amount = 0.08), gains[[j]], pch = 16, cex = 0.55, col = adjustcolor(blue, 0.65))
  panel_label("B")

  y_range <- range(c(simulation$marker_downsample_mean_abs_delta, simulation$anchored_projection_downsample_mean_abs_delta))
  plot(c(1, 2), y_range, type = "n", xaxt = "n", xlab = "", ylab = "Mean absolute change after 50% thinning",
       main = "Downsampling stability")
  axis(1, at = c(1, 2), labels = c("Marker score", "Anchored projection"))
  for (i in seq_len(nrow(simulation))) {
    lines(c(1, 2), c(simulation$marker_downsample_mean_abs_delta[i], simulation$anchored_projection_downsample_mean_abs_delta[i]),
          col = adjustcolor("#7A7A7A", 0.45), lwd = 0.8)
  }
  points(rep(1, nrow(simulation)), simulation$marker_downsample_mean_abs_delta, pch = 16, col = blue, cex = 0.65)
  points(rep(2, nrow(simulation)), simulation$anchored_projection_downsample_mean_abs_delta, pch = 16, col = orange, cex = 0.65)
  panel_label("C")

  realdata <- realdata[order(realdata$dataset, realdata$holdout_group), , drop = FALSE]
  point_col <- ifelse(realdata$dataset == "GSE167339", blue, orange)
  plot(seq_len(nrow(realdata)), realdata$anchored_gain, type = "n", xaxt = "n", xlab = "Held-out group",
       ylab = "Anchored NB gain / observation", main = "Unseen-group count prediction",
       ylim = range(c(realdata$anchored_gain, 0)) + c(-0.02, 0.02))
  abline(h = 0, lty = 2, col = "#6B6B6B")
  axis(1, at = seq_len(nrow(realdata)), labels = realdata$holdout_group, las = 2, cex.axis = 0.82)
  points(seq_len(nrow(realdata)), realdata$anchored_gain, pch = 21, bg = point_col, col = ink, cex = 1.25)
  legend("topright", legend = c("GSE167339 donor", "GSE246215 cancer dataset"),
         pch = 21, pt.bg = c(blue, orange), col = ink, bty = "n", cex = 0.78)
  panel_label("D")
}

ragg::agg_png(file.path(EXPORTS, paste0(STEM, ".png")), width = 4800, height = 3600, units = "px", res = 400)
draw_figure()
dev.off()

ragg::agg_tiff(file.path(EXPORTS, paste0(STEM, ".tiff")), width = 4800, height = 3600, units = "px", res = 400, compression = "lzw")
draw_figure()
dev.off()

pdf(file.path(EXPORTS, paste0(STEM, ".pdf")), width = 12, height = 9, useDingbats = FALSE)
draw_figure()
dev.off()

svglite::svglite(file.path(EXPORTS, paste0(STEM, ".svg")), width = 12, height = 9)
draw_figure()
dev.off()

legend <- c(
  "# Figure S4. Anchored count fitting adds held-out count prediction but not improved coordinate recovery",
  "",
  "**(A)** Held-out state-associated-gene log likelihood for the intercept-only model, anchored NB without source terms, and source-aware anchored NB across 20 independent simulations; each point is the mean of four random marker splits. **(B)** Within-simulation gains of source-aware anchored NB over the intercept-only and no-source models. Gains were positive in 19 of 20 simulations for both contrasts. **(C)** Mean absolute coordinate change after 50% binomial thinning for direct marker scoring and anchored fixed-program projection. Anchored projection was more stable in 0 of 20 simulations. **(D)** Anchored NB gain over the intercept-only model for unseen GSE167339 donors and GSE246215 cancer datasets, averaged over four marker splits. Six of seven held-out groups were positive; HCC was negative. Complementary markers were excluded from coordinate construction and used for held-out count evaluation. Panels A--C are simulations; panel D is a numerical public-data generalization check, not biological ground-truth validation."
)
writeLines(legend, file.path(OUT, "main_figure_legend.md"))

write_tsv(data.frame(
  panel = c("A", "B", "C", "D"),
  source_data = c(
    "source_data/figs4_simulation_replicate_means.tsv",
    "source_data/figs4_simulation_replicate_means.tsv",
    "source_data/figs4_simulation_replicate_means.tsv",
    "source_data/figs4_realdata_holdout_group_means.tsv"
  ),
  description = c(
    "Held-out state-associated-gene log likelihood across simulation replicates",
    "Paired source-aware predictive gains across simulation replicates",
    "Paired coordinate change after 50% count thinning",
    "Held-out donor and cancer-dataset predictive gains"
  ), stringsAsFactors = FALSE
), file.path(OUT, "panel_source_data_manifest.tsv"))

write_tsv(data.frame(
  figure = "Figure S4",
  title = "Anchored count fitting adds held-out count prediction but not improved coordinate recovery",
  panels = "A-D",
  exports = paste0("exports/", STEM, ".pdf; exports/", STEM, ".svg; exports/", STEM, ".png; exports/", STEM, ".tiff"),
  source_data_complete = TRUE,
  source_manifest = "panel_source_data_manifest.tsv",
  stringsAsFactors = FALSE
), file.path(OUT, "figure_manifest.tsv"))

message("Supplementary Figure S4 written to: ", normalizePath(OUT))
