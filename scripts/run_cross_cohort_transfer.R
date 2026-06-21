#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_cross_cohort_transfer.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "cross_cohort_transfer")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "simulate_fibrodynmix.R",
  "baseline_marker_scoring.R",
  "benchmark_metrics.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "cross_cohort_transfer.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

benchmark <- run_cross_cohort_transfer_benchmark(
  n_replicates = 2,
  seed = 710,
  simulation_args = list(
    n_studies = 3,
    donors_per_study = 2,
    cells_per_donor = 6,
    n_genes = 90,
    marker_genes_per_state = 4,
    scenario = "batch_confounding"
  ),
  train_fit_args = list(
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 20,
    maxit_z = 20
  ),
  transfer_args = list(maxit_z = 50)
)

summary <- aggregate(
  cbind(transfer_rmse, transfer_mean_absolute_error, transfer_dominant_accuracy, transfer_nb_objective, transfer_z_convergence_rate) ~ holdout_study,
  data = benchmark,
  FUN = mean
)
summary$n_replicates <- as.integer(table(benchmark$holdout_study)[summary$holdout_study])

write.table(benchmark, file.path(OUT, "transfer_benchmark.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(OUT, "transfer_benchmark_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

manifest <- data.frame(
  analysis = "cross_cohort_transfer",
  claim_boundary = "Leave-study-out simulation transfer benchmark; not completed human disease-atlas or cross-atlas generalization.",
  n_replicates = 2,
  n_holdout_rows = nrow(benchmark),
  mean_transfer_rmse = mean(benchmark$transfer_rmse),
  mean_transfer_dominant_accuracy = mean(benchmark$transfer_dominant_accuracy),
  mean_transfer_z_convergence_rate = mean(benchmark$transfer_z_convergence_rate),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(OUT, "transfer_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Cross-cohort transfer analysis written to: ", normalizePath(OUT))
