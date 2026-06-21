#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_selected_nb_final_summaries.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "selected_nb_final_summaries"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

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
  "simulation_benchmark.R",
  "cross_cohort_transfer.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))
source(file.path(ROOT, "scripts", "utils_selected_nb_defaults.R"))

defaults <- read_selected_nb_defaults(file.path(ROOT, "analysis", "validation_aware_nb_selection", "validation_aware_nb_selected_summary.tsv"))
defaults_used <- selected_nb_defaults_used(defaults)
write_tsv(defaults_used, file.path(OUT, "selected_nb_defaults_used.tsv"))

scenario_list <- c("continuous", "batch_confounding", "rare_transition")
vi_metrics <- list()
vi_summaries <- list()
vi_optimizers <- list()
for (scenario in scenario_list) {
  n_outer <- selected_nb_default(defaults, scenario, "nb", fallback = if (scenario == "batch_confounding") 10L else 2L)
  benchmark <- run_simulation_benchmark(
    scenarios = scenario,
    n_replicates = as.integer(get_arg("n-replicates", "2")),
    seed = as.integer(get_arg("seed", "202615")) + match(scenario, scenario_list) * 100L,
    methods = c("fibrodynmix_nb", "fibrodynmix_vi"),
    simulation_args = list(
      n_studies = 2,
      donors_per_study = 2,
      cells_per_donor = 8,
      n_genes = 90,
      marker_genes_per_state = 4
    ),
    nb_args = list(
      n_outer = n_outer,
      initializer_args = list(n_iter = 2),
      maxit_beta = 15,
      maxit_z = 20
    ),
    vi_args = list(
      n_draws = 30,
      n_elbo_draws = 6,
      n_vi_iter = 2,
      posterior_scale = 0.25,
      seed = as.integer(get_arg("seed", "202615")) + match(scenario, scenario_list) * 100L + 1L,
      keep_draws = FALSE
    )
  )
  benchmark$selected_n_outer <- n_outer
  vi_metrics[[scenario]] <- benchmark
  summary <- summarize_benchmark_results(
    benchmark,
    metrics = c(
      "rmse",
      "dominant_accuracy",
      "downstream_balanced_accuracy",
      "downstream_macro_f1",
      "mean_entropy_pred",
      "vi_interval_coverage",
      "vi_mean_interval_width",
      "vi_calibrated_interval_coverage",
      "vi_calibrated_mean_interval_width"
    )
  )
  summary$selected_n_outer <- n_outer
  vi_summaries[[scenario]] <- summary
  optimizer <- summarize_optimizer_diagnostics(benchmark)
  optimizer$selected_n_outer <- n_outer
  vi_optimizers[[scenario]] <- optimizer
}
vi_metrics <- do.call(rbind, vi_metrics)
vi_summary <- do.call(rbind, vi_summaries)
vi_optimizer <- do.call(rbind, vi_optimizers)
write_tsv(vi_metrics, file.path(OUT, "selected_nb_vi_benchmark_metrics.tsv"))
write_tsv(vi_summary, file.path(OUT, "selected_nb_vi_benchmark_summary.tsv"))
write_tsv(vi_optimizer, file.path(OUT, "selected_nb_vi_optimizer_diagnostics.tsv"))

transfer_n_outer <- selected_nb_default(defaults, "batch_confounding", "nb", fallback = 10L)
transfer <- run_cross_cohort_transfer_benchmark(
  n_replicates = as.integer(get_arg("transfer-replicates", "2")),
  seed = as.integer(get_arg("seed", "202615")) + 900L,
  simulation_args = list(
    n_studies = 3,
    donors_per_study = 2,
    cells_per_donor = 6,
    n_genes = 90,
    marker_genes_per_state = 4,
    scenario = "batch_confounding"
  ),
  train_fit_args = list(
    n_outer = transfer_n_outer,
    initializer_args = list(n_iter = 2),
    maxit_beta = 20,
    maxit_z = 20
  ),
  transfer_args = list(maxit_z = 50)
)
transfer$selected_n_outer <- transfer_n_outer
transfer_summary <- aggregate(
  cbind(transfer_rmse, transfer_mean_absolute_error, transfer_dominant_accuracy, transfer_nb_objective, transfer_z_convergence_rate) ~ holdout_study,
  data = transfer,
  FUN = mean
)
transfer_summary$n_replicates <- as.integer(table(transfer$holdout_study)[transfer_summary$holdout_study])
transfer_summary$selected_n_outer <- transfer_n_outer
write_tsv(transfer, file.path(OUT, "selected_nb_transfer_benchmark.tsv"))
write_tsv(transfer_summary, file.path(OUT, "selected_nb_transfer_benchmark_summary.tsv"))

downstream_n_outer <- selected_nb_default(defaults, "batch_confounding", "nb_study", fallback = 2L)
downstream_out <- file.path(OUT, "gse246215_downstream_selected_nb")
cmd <- file.path(R.home("bin"), "Rscript")
cmd_args <- c(
  file.path(ROOT, "scripts", "run_gse246215_downstream_benchmark.R"),
  paste0("--out=", downstream_out),
  paste0("--n-outer=", downstream_n_outer),
  "--initializer-iter=3",
  "--maxit-beta=12",
  "--maxit-z=12"
)
status <- system2(cmd, cmd_args)
if (!identical(status, 0L)) {
  stop("Selected-NB GSE246215 downstream rerun failed.", call. = FALSE)
}

downstream_manifest <- read_tsv(file.path(downstream_out, "gse246215_downstream_manifest.tsv"))
downstream_metrics <- read_tsv(file.path(downstream_out, "gse246215_downstream_classification_metrics.tsv"))
downstream_metrics$selected_n_outer <- downstream_n_outer
write_tsv(downstream_metrics, file.path(OUT, "selected_nb_gse246215_downstream_classification_metrics.tsv"))

vi_rows <- vi_metrics[vi_metrics$method == "fibrodynmix_vi", , drop = FALSE]
manifest <- data.frame(
  analysis = "selected_nb_final_summaries",
  primary_claim = "Final VI, transfer, and downstream summaries are anchored to validation-aware selected NB outer-iteration defaults rather than fixed smoke settings.",
  claim_boundary = "Simulation-calibrated selected defaults. GSE246215 downstream rerun remains a representation-utility case study from one GEO accession, not diagnostic validation.",
  vi_scenarios = paste(scenario_list, collapse = ";"),
  vi_selected_n_outer = paste(unique(paste(vi_summary$scenario, vi_summary$selected_n_outer, sep = "=")), collapse = ";"),
  mean_vi_rmse = mean(vi_rows$rmse, na.rm = TRUE),
  mean_vi_downstream_balanced_accuracy = mean(vi_rows$downstream_balanced_accuracy, na.rm = TRUE),
  transfer_selected_n_outer = transfer_n_outer,
  mean_transfer_rmse = mean(transfer$transfer_rmse, na.rm = TRUE),
  mean_transfer_z_convergence_rate = mean(transfer$transfer_z_convergence_rate, na.rm = TRUE),
  downstream_selected_n_outer = downstream_n_outer,
  downstream_best_method = downstream_manifest$best_method_by_balanced_accuracy[[1]],
  downstream_best_balanced_accuracy = downstream_manifest$best_balanced_accuracy[[1]],
  downstream_fibrodynmix_balanced_accuracy = downstream_manifest$fibrodynmix_balanced_accuracy[[1]],
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "selected_nb_final_summaries_manifest.tsv"))

message("Selected-NB final summaries written to: ", normalizePath(OUT))
