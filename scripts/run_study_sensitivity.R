#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_study_sensitivity.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "study_effect_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "R", "simulate_fibrodynmix.R"))
source(file.path(ROOT, "R", "benchmark_metrics.R"))
source(file.path(ROOT, "R", "baseline_marker_scoring.R"))
source(file.path(ROOT, "R", "topic_nmf_baseline.R"))
source(file.path(ROOT, "R", "fibrodynmix_initializer.R"))
source(file.path(ROOT, "R", "nb_likelihood.R"))
source(file.path(ROOT, "R", "fit_nb_model.R"))
source(file.path(ROOT, "R", "simulation_benchmark.R"))
source(file.path(ROOT, "R", "study_effect_sensitivity.R"))

results <- run_study_effect_sensitivity(
  study_l2_grid = c(0.05, 0.1, 0.5, 1, 5),
  marker_l2_grid = c(0.05, 0.1),
  n_replicates = 2,
  seed = 410,
  simulation_args = list(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 12,
    n_genes = 120,
    marker_genes_per_state = 6
  ),
  nb_args = list(
    n_outer = 2,
    initializer_args = list(n_iter = 3),
    maxit_beta = 15,
    maxit_z = 12
  )
)

write.table(
  results,
  file.path(OUT, "study_effect_sensitivity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- aggregate(
  cbind(rmse, dominant_accuracy, nb_best_objective, study_effect_l2_norm, study_effect_mean_abs) ~
    method + study_l2 + marker_l2,
  data = results,
  FUN = function(x) mean(x, na.rm = TRUE)
)
write.table(
  summary,
  file.path(OUT, "study_effect_sensitivity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

selected <- select_study_effect_penalty(results)
write.table(
  selected$tradeoff_table,
  file.path(OUT, "study_effect_penalty_tradeoff.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(
  c(
    paste0("recommended_study_l2: ", selected$recommended_study_l2),
    paste0("recommended_marker_l2: ", selected$recommended_marker_l2),
    paste0("selection_reason: ", selected$selection_reason)
  ),
  file.path(OUT, "study_effect_penalty_selection.txt")
)

message("Study-effect sensitivity written to: ", OUT)
