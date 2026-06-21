#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_vi_posterior.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "vi_posterior")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "simulate_fibrodynmix.R",
  "baseline_marker_scoring.R",
  "benchmark_metrics.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "bootstrap_uncertainty.R",
  "vi_posterior.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

sim <- simulate_fibrodynmix(
  n_studies = 2,
  donors_per_study = 2,
  cells_per_donor = 8,
  n_genes = 90,
  marker_genes_per_state = 4,
  scenario = "continuous",
  seed = 810
)
metadata <- sim$cell_metadata
rownames(metadata) <- colnames(sim$counts)

vi <- fit_fibrodynmix_vi(
  counts = sim$counts,
  marker_index = sim$parameters$marker_index,
  library_size = metadata$library_size,
  cell_metadata = metadata,
  sample_col = "donor_id",
  nb_args = list(
    study_id = metadata$study_id,
    donor_id = metadata$donor_id,
    fit_study_effect = TRUE,
    fit_donor_effect = TRUE,
    n_outer = 2,
    initializer_args = list(n_iter = 2),
    maxit_beta = 15,
    maxit_z = 20
  ),
  n_draws = 40,
  n_elbo_draws = 8,
  n_vi_iter = 3,
  posterior_scale = 0.25,
  seed = 811,
  keep_draws = FALSE
)

cell_z <- vi$cell_summary$z
cell_entropy <- vi$cell_summary$entropy
sample_summary <- vi$sample_summary
metrics <- evaluate_state_weights(sim$z, vi$z_mean)
scalar_metrics <- metrics[c("rmse", "mean_absolute_error", "dominant_accuracy", "mean_entropy_true", "mean_entropy_pred")]

write_tsv(vi$elbo_trace, file.path(OUT, "vi_elbo_trace.tsv"))
write_tsv(cell_z, file.path(OUT, "vi_cell_state_intervals.tsv"))
write_tsv(cell_entropy, file.path(OUT, "vi_cell_entropy_intervals.tsv"))
write_tsv(sample_summary, file.path(OUT, "vi_sample_composition_intervals.tsv"))
write_tsv(
  data.frame(
    metric = names(scalar_metrics),
    value = unlist(scalar_metrics, use.names = FALSE),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "vi_state_recovery_metrics.tsv")
)
write_tsv(
  data.frame(
    state = names(metrics$per_state_rmse),
    per_state_rmse = as.numeric(metrics$per_state_rmse),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "vi_per_state_rmse.tsv")
)

manifest <- data.frame(
  analysis = "vi_posterior",
  primary_claim = "A logistic-normal variational posterior layer returns raw-count NB posterior state summaries with credible intervals and an ELBO-like trace.",
  claim_boundary = "Lightweight mean-field posterior over cell state logits around the NB mode; not a full amortized VI or posterior over all hierarchical parameters.",
  n_cells = ncol(sim$counts),
  n_genes = nrow(sim$counts),
  n_states = ncol(sim$z),
  n_draws = 40,
  best_elbo = vi$best_elbo,
  state_weight_rmse = metrics$rmse,
  dominant_accuracy = metrics$dominant_accuracy,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "vi_manifest.tsv"))

message("VI posterior analysis written to: ", normalizePath(OUT))
