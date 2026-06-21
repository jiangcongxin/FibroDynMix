#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_bootstrap_uncertainty.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "bootstrap_uncertainty")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "R", "simulate_fibrodynmix.R"))
source(file.path(ROOT, "R", "benchmark_metrics.R"))
source(file.path(ROOT, "R", "baseline_marker_scoring.R"))
source(file.path(ROOT, "R", "fibrodynmix_initializer.R"))
source(file.path(ROOT, "R", "nb_likelihood.R"))
source(file.path(ROOT, "R", "fit_nb_model.R"))
source(file.path(ROOT, "R", "bootstrap_uncertainty.R"))

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

write.table(boot$sample_summary, file.path(OUT, "sample_composition_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(boot$cell_summary, file.path(OUT, "cell_state_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(boot$marker_summary, file.path(OUT, "marker_program_uncertainty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(boot$sample_draws, file.path(OUT, "sample_composition_draws.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Bootstrap uncertainty written to: ", OUT)
