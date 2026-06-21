#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_vi_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "vi_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "topic_nmf_baseline.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "bootstrap_uncertainty.R",
  "vi_posterior.R",
  "simulation_benchmark.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

benchmark <- run_simulation_benchmark(
  scenarios = c("continuous", "batch_confounding", "rare_transition"),
  n_replicates = 2,
  seed = 910,
  methods = c("fibrodynmix_nb", "fibrodynmix_vi"),
  simulation_args = list(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 8,
    n_genes = 90,
    marker_genes_per_state = 4
  ),
  nb_args = list(
    n_outer = 2,
    initializer_args = list(n_iter = 2),
    maxit_beta = 15,
    maxit_z = 20
  ),
  vi_args = list(
    n_draws = 30,
    n_elbo_draws = 6,
    n_vi_iter = 2,
    posterior_scale = 0.25,
    seed = 911,
    keep_draws = FALSE
  )
)

summary <- summarize_benchmark_results(
  benchmark,
  metrics = c(
    "rmse",
    "dominant_accuracy",
    "mean_entropy_pred",
    "vi_interval_coverage",
    "vi_mean_interval_width",
    "vi_calibrated_interval_coverage",
    "vi_calibrated_mean_interval_width"
  )
)
optimizer <- summarize_optimizer_diagnostics(benchmark)

write_tsv(benchmark, file.path(OUT, "vi_benchmark_metrics.tsv"))
write_tsv(summary, file.path(OUT, "vi_benchmark_summary.tsv"))
write_tsv(optimizer, file.path(OUT, "vi_optimizer_diagnostics.tsv"))

vi_rows <- benchmark[benchmark$method == "fibrodynmix_vi", , drop = FALSE]
manifest <- data.frame(
  analysis = "vi_benchmark",
  primary_claim = "FibroDynMix VI benchmark quantifies posterior interval coverage and state recovery under simulated latent truth.",
  claim_boundary = "Simulation calibration of the lightweight logistic-normal posterior; not full real-human atlas validation or fully hierarchical VI.",
  n_rows = nrow(benchmark),
  n_vi_rows = nrow(vi_rows),
  mean_vi_interval_coverage = mean(vi_rows$vi_interval_coverage, na.rm = TRUE),
  mean_vi_interval_width = mean(vi_rows$vi_mean_interval_width, na.rm = TRUE),
  mean_vi_calibrated_interval_coverage = mean(vi_rows$vi_calibrated_interval_coverage, na.rm = TRUE),
  mean_vi_calibrated_interval_width = mean(vi_rows$vi_calibrated_mean_interval_width, na.rm = TRUE),
  mean_vi_rmse = mean(vi_rows$rmse, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "vi_benchmark_manifest.tsv"))

message("VI benchmark analysis written to: ", normalizePath(OUT))
