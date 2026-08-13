#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/make_figure4_realdata_holdout.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- Sys.getenv("FDM_FIGURE_OUT", file.path(ROOT, "paper_rewriting_output", "figure_rebuild_iscience_marker_anchored_v20", "main", "Figure_4"))
STEM <- Sys.getenv("FDM_FIGURE_STEM", "Figure_4")
EXPORTS <- file.path(OUT, "exports"); SOURCE <- file.path(OUT, "source_data")
dir.create(EXPORTS, recursive = TRUE, showWarnings = FALSE); dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)

metrics_path <- file.path(ROOT, "analysis", "realdata_heldout_count_validation", "realdata_heldout_metrics.tsv")
means_path <- file.path(ROOT, "analysis", "realdata_heldout_count_validation", "realdata_holdout_group_means.tsv")
metrics <- read.delim(metrics_path, check.names = FALSE); means <- read.delim(means_path, check.names = FALSE)
if (nrow(means) != 7L) stop("Figure 4 requires seven held-out public-data groups.")
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
write_tsv(metrics, file.path(SOURCE, "fig4_realdata_marker_split_metrics.tsv"))
write_tsv(means, file.path(SOURCE, "fig4_realdata_holdout_group_means.tsv"))

ink <- "#202124"; blue <- "#2F6690"; orange <- "#C46A32"; teal <- "#3A7D78"
panel_label <- function(label) { u <- par("usr"); text(u[1], u[4], label, adj = c(-0.55, 1.2), xpd = NA, font = 2, cex = 1.2, col = ink) }
draw <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  layout(matrix(1:4, 2, byrow = TRUE)); par(family = "sans", fg = ink, col.axis = ink, col.lab = ink, col.main = ink, mar = c(5.7, 5.2, 3.2, 1.2), mgp = c(3.05, 0.8, 0), tcl = -0.25)
  ord <- order(means$dataset, means$holdout_group); x <- means[ord, ]; cols <- ifelse(x$dataset == "GSE167339", blue, orange)
  plot(seq_len(nrow(x)), x$anchored_gain, type = "n", xaxt = "n", xlab = "Held-out donor or cancer dataset", ylab = "Anchored NB gain / observation", main = "Unseen-group count prediction", ylim = range(c(0, x$anchored_gain)) + c(-0.02, 0.02))
  abline(h = 0, lty = 2, col = "#6B6B6B"); axis(1, seq_len(nrow(x)), x$holdout_group, las = 2, cex.axis = 0.82)
  points(seq_len(nrow(x)), x$anchored_gain, pch = 21, bg = cols, col = ink, cex = 1.3)
  legend("topright", c("GSE167339 donor", "GSE246215 cancer dataset"), pch = 21, pt.bg = c(blue, orange), col = ink, bty = "n", cex = 0.78); panel_label("A")

  split_gain <- split(metrics$anchored_gain, metrics$dataset)
  boxplot(split_gain, col = c("#D9E6F2", "#F2DED0"), border = ink, outline = FALSE, ylab = "Anchored gain / observation", main = "Marker-split replication")
  abline(h = 0, lty = 2, col = "#6B6B6B"); set.seed(4402)
  for (j in seq_along(split_gain)) points(jitter(rep(j, length(split_gain[[j]])), 0.08), split_gain[[j]], pch = 16, cex = 0.55, col = adjustcolor(if (j == 1) blue else orange, 0.65)); panel_label("B")

  plot(means$train_drift_from_marker, means$anchored_gain, pch = 21, bg = ifelse(means$dataset == "GSE167339", blue, orange), col = ink, cex = 1.25,
       xlab = "Training-coordinate drift from marker score", ylab = "Held-out gain / observation", main = "Prediction gain around the anchored coordinate")
  abline(h = 0, lty = 2, col = "#6B6B6B"); text(means$train_drift_from_marker, means$anchored_gain, labels = means$holdout_group, pos = ifelse(means$holdout_group == "NSCLC", 1, 3), cex = 0.68); panel_label("C")

  dsum <- aggregate(cbind(anchored_gain, train_drift_from_marker) ~ dataset, means, mean)
  favourable <- c(sum(means$dataset == "GSE167339" & means$anchored_gain > 0), sum(means$dataset == "GSE246215" & means$anchored_gain > 0))
  total <- c(sum(means$dataset == "GSE167339"), sum(means$dataset == "GSE246215"))
  bp <- barplot(dsum$anchored_gain, names.arg = dsum$dataset, col = c("#D9E6F2", "#F2DED0"), border = ink, ylab = "Mean held-out gain / observation", main = "Dataset-level generalization summary", ylim = c(0, max(dsum$anchored_gain) * 1.18))
  text(bp, dsum$anchored_gain, labels = paste0(favourable, "/", total, " groups positive"), pos = 3, cex = 0.8); panel_label("D")
}
ragg::agg_png(file.path(EXPORTS, paste0(STEM, ".png")), 4800, 3600, units = "px", res = 400); draw(); dev.off()
ragg::agg_tiff(file.path(EXPORTS, paste0(STEM, ".tiff")), 4800, 3600, units = "px", res = 400, compression = "lzw"); draw(); dev.off()
pdf(file.path(EXPORTS, paste0(STEM, ".pdf")), 12, 9, useDingbats = FALSE); draw(); dev.off()
svglite::svglite(file.path(EXPORTS, paste0(STEM, ".svg")), 12, 9); draw(); dev.off()
writeLines(c("# Figure 4. Anchored count programs generalize to held-out donors and datasets", "", "**(A)** Anchored NB gain over an intercept-only model for three unseen GSE167339 donors and four unseen GSE246215 cancer datasets, averaged over four marker splits. Six of seven held-out groups were positive; HCC was negative. **(B)** Marker-split-level gains within each public dataset. **(C)** Group-level held-out gain plotted against training-coordinate drift from the marker score. **(D)** Dataset-level mean gains and favorable-group counts. In every split, half of each curated marker set defined the coordinate and the complementary markers were used only for held-out count evaluation. These panels test numerical generalization of count programs to unseen groups, not recovery of a biological ground-truth coordinate."), file.path(OUT, "main_figure_legend.md"))
write_tsv(data.frame(panel = c("A", "B", "C", "D"), source_data = c("source_data/fig4_realdata_holdout_group_means.tsv", "source_data/fig4_realdata_marker_split_metrics.tsv", "source_data/fig4_realdata_holdout_group_means.tsv", "source_data/fig4_realdata_holdout_group_means.tsv"), description = c("Held-out group gains", "Four marker-split replicates per group", "Gain versus anchored coordinate drift", "Dataset-level summary")), file.path(OUT, "panel_source_data_manifest.tsv"))
write_tsv(data.frame(figure = "Figure 4", title = "Anchored count programs generalize to held-out donors and datasets", panels = "A-D", exports = paste0("exports/", STEM, ".pdf; exports/", STEM, ".svg; exports/", STEM, ".png; exports/", STEM, ".tiff"), source_data_complete = TRUE, source_manifest = "panel_source_data_manifest.tsv"), file.path(OUT, "figure_manifest.tsv"))
