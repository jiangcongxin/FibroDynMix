#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_extended_method_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "extended_method_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "matrix_utils.R",
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "topic_nmf_baseline.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "bootstrap_uncertainty.R",
  "vi_posterior.R",
  "scvi_baseline.R",
  "simulation_benchmark.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

topic_backend <- if (requireNamespace("NMF", quietly = TRUE)) "nmf" else "multiplicative_update"
scvi_available <- isTRUE(tryCatch(
  requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("scvi") &&
    reticulate::py_module_available("anndata"),
  error = function(e) FALSE
))
benchmark_methods <- c(
  "marker_scoring",
  "topic_nmf",
  "fibrodynmix_initializer",
  "fibrodynmix_nb",
  "fibrodynmix_nb_study",
  "fibrodynmix_vi"
)
if (scvi_available) {
  benchmark_methods <- append(benchmark_methods, "scvi_latent", after = match("topic_nmf", benchmark_methods))
}

benchmark <- run_simulation_benchmark(
  scenarios = c("continuous", "discrete", "batch_confounding", "rare_transition"),
  n_replicates = 2,
  seed = 720,
  methods = benchmark_methods,
  simulation_args = list(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 8,
    n_genes = 90,
    marker_genes_per_state = 4
  ),
  initializer_args = list(n_iter = 2),
  topic_nmf_args = list(
    backend = topic_backend,
    n_iter = 40
  ),
  scvi_args = list(
    max_epochs = 12
  ),
  nb_args = list(
    n_outer = 2,
    initializer_args = list(n_iter = 2),
    maxit_beta = 12,
    maxit_z = 12
  ),
  vi_args = list(
    n_draws = 18,
    n_elbo_draws = 4,
    n_vi_iter = 2,
    seed = 721,
    keep_draws = FALSE
  ),
  downstream_task = TRUE,
  downstream_label_col = "disease",
  downstream_group_col = "donor_id",
  downstream_folds = 3
)

summary <- summarize_benchmark_results(
  benchmark,
  metrics = c(
    "rmse",
    "mean_absolute_error",
    "dominant_accuracy",
    "mean_entropy_pred",
    "downstream_balanced_accuracy",
    "downstream_macro_f1",
    "downstream_macro_auroc",
    "vi_interval_coverage",
    "vi_calibrated_interval_coverage"
  )
)
optimizer <- summarize_optimizer_diagnostics(benchmark)

write_tsv(benchmark, file.path(OUT, "extended_method_benchmark_metrics.tsv"))
write_tsv(summary, file.path(OUT, "extended_method_benchmark_summary.tsv"))
write_tsv(optimizer, file.path(OUT, "extended_method_optimizer_diagnostics.tsv"))

method_rank <- do.call(rbind, lapply(split(benchmark, benchmark$scenario), function(df) {
  aggregate(
    cbind(rmse, dominant_accuracy, downstream_balanced_accuracy, downstream_macro_f1) ~ method,
    data = df,
    FUN = mean
  )
}))
method_rank$scenario <- rep(names(split(benchmark, benchmark$scenario)), vapply(split(benchmark, benchmark$scenario), function(x) length(unique(x$method)), integer(1)))
method_rank <- method_rank[, c("scenario", "method", "rmse", "dominant_accuracy", "downstream_balanced_accuracy", "downstream_macro_f1"), drop = FALSE]
method_rank <- method_rank[order(method_rank$scenario, method_rank$rmse), , drop = FALSE]
write_tsv(method_rank, file.path(OUT, "extended_method_rankings.tsv"))

fibro_rows <- benchmark[benchmark$method %in% c("fibrodynmix_nb", "fibrodynmix_nb_study", "fibrodynmix_vi"), , drop = FALSE]
baseline_rows <- benchmark[benchmark$method %in% c("marker_scoring", "topic_nmf", "scvi_latent"), , drop = FALSE]
manifest <- data.frame(
  analysis = "extended_method_benchmark",
  primary_claim = "FibroDynMix is benchmarked against marker scoring, NMF/topic, optional scVI latent, and downstream disease-label prediction across simulated fibroblast state scenarios.",
  claim_boundary = "Simulation benchmark with known latent truth and simulated donor-level disease labels; not a replacement for independent biological validation. Topic baseline uses the NMF package when available and an internal KL-NMF fallback otherwise. scVI runs only when reticulate, scvi-tools, and anndata are available.",
  n_rows = nrow(benchmark),
  n_scenarios = length(unique(benchmark$scenario)),
  n_methods = length(unique(benchmark$method)),
  topic_backend = topic_backend,
  scvi_available = scvi_available,
  scvi_status = if (scvi_available) "run" else "not_available_skipped",
  mean_fibrodynmix_rmse = mean(fibro_rows$rmse, na.rm = TRUE),
  mean_baseline_rmse = mean(baseline_rows$rmse, na.rm = TRUE),
  mean_fibrodynmix_dominant_accuracy = mean(fibro_rows$dominant_accuracy, na.rm = TRUE),
  mean_baseline_dominant_accuracy = mean(baseline_rows$dominant_accuracy, na.rm = TRUE),
  mean_fibrodynmix_downstream_balanced_accuracy = mean(fibro_rows$downstream_balanced_accuracy, na.rm = TRUE),
  mean_baseline_downstream_balanced_accuracy = mean(baseline_rows$downstream_balanced_accuracy, na.rm = TRUE),
  mean_fibrodynmix_downstream_macro_f1 = mean(fibro_rows$downstream_macro_f1, na.rm = TRUE),
  mean_baseline_downstream_macro_f1 = mean(baseline_rows$downstream_macro_f1, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "extended_method_benchmark_manifest.tsv"))

message("Extended method benchmark written to: ", normalizePath(OUT))
