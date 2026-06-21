#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

parse_int_vector <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || anyNA(out) || any(out < 1L)) {
    stop("Integer vector arguments must contain positive comma-separated integers.", call. = FALSE)
  }
  sort(unique(out))
}

parse_char_vector <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  if (length(out) == 0L) {
    stop("Character vector arguments cannot be empty.", call. = FALSE)
  }
  out
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_validation_aware_nb_selection.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "validation_aware_nb_selection"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "matrix_utils.R",
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "cross_cohort_transfer.R",
  "simulation_benchmark.R",
  "nb_model_selection.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

scenarios <- parse_char_vector(
  get_arg("scenarios"),
  c("continuous", "discrete", "batch_confounding", "rare_transition")
)
n_replicates <- as.integer(get_arg("n-replicates", "2"))
candidate_n_outer <- parse_int_vector(get_arg("n-outer-grid"), c(2L, 5L, 10L, 20L))
seed <- as.integer(get_arg("seed", "202614"))
n_studies <- as.integer(get_arg("n-studies", "2"))
donors_per_study <- as.integer(get_arg("donors-per-study", "2"))
cells_per_donor <- as.integer(get_arg("cells-per-donor", "8"))
n_genes <- as.integer(get_arg("n-genes", "90"))
marker_genes_per_state <- as.integer(get_arg("marker-genes-per-state", "4"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "12"))
maxit_z <- as.integer(get_arg("maxit-z", "12"))
holdout_fraction <- as.numeric(get_arg("holdout-fraction", "0.25"))
variants <- parse_char_vector(get_arg("variants"), c("nb", "nb_study"))

selection_rows <- list()
component_rows <- list()
split_rows <- list()
row_idx <- 1L
component_idx <- 1L
split_idx <- 1L

for (scenario in scenarios) {
  for (replicate in seq_len(n_replicates)) {
    replicate_seed <- seed + scenario_seed_offset(scenario) + replicate - 1L
    sim <- simulate_fibrodynmix(
      scenario = scenario,
      seed = replicate_seed,
      n_studies = n_studies,
      donors_per_study = donors_per_study,
      cells_per_donor = cells_per_donor,
      n_genes = n_genes,
      marker_genes_per_state = marker_genes_per_state
    )

    for (variant in variants) {
      if (!variant %in% c("nb", "nb_study")) {
        stop("Supported variants are `nb` and `nb_study`.", call. = FALSE)
      }
      nb_args <- list(
        initializer_args = list(n_iter = initializer_iter),
        maxit_beta = maxit_beta,
        maxit_z = maxit_z
      )
      if (variant == "nb_study") {
        nb_args$study_id <- sim$cell_metadata$study_id
        nb_args$fit_study_effect <- TRUE
        nb_args$study_l2 <- as.numeric(get_arg("study-l2", "5"))
      }
      selector <- select_fibrodynmix_nb_model(
        counts = sim$counts,
        marker_index = sim$parameters$marker_index,
        library_size = sim$cell_metadata$library_size,
        candidate_n_outer = candidate_n_outer,
        labels = sim$cell_metadata$disease,
        groups = sim$cell_metadata$donor_id,
        holdout_fraction = holdout_fraction,
        seed = replicate_seed + if (variant == "nb_study") 10000L else 0L,
        nb_args = nb_args,
        transfer_args = list(maxit_z = maxit_z),
        truth_z = sim$z
      )
      candidates <- selector$candidate_scores
      candidates$scenario <- scenario
      candidates$replicate <- replicate
      candidates$seed <- replicate_seed
      candidates$variant <- variant
      candidates$selected_n_outer <- selector$selected_n_outer
      selection_rows[[row_idx]] <- candidates
      row_idx <- row_idx + 1L

      components <- selector$selection_components
      components$scenario <- scenario
      components$replicate <- replicate
      components$seed <- replicate_seed
      components$variant <- variant
      component_rows[[component_idx]] <- components
      component_idx <- component_idx + 1L

      split <- selector$split
      split$scenario <- scenario
      split$replicate <- replicate
      split$seed <- replicate_seed
      split$variant <- variant
      split_rows[[split_idx]] <- split
      split_idx <- split_idx + 1L
    }
  }
}

selection <- do.call(rbind, selection_rows)
components <- do.call(rbind, component_rows)
splits <- do.call(rbind, split_rows)

write_tsv(selection, file.path(OUT, "validation_aware_nb_selection_candidates.tsv"))
write_tsv(components, file.path(OUT, "validation_aware_nb_marker_gradient_components.tsv"))
write_tsv(splits, file.path(OUT, "validation_aware_nb_splits.tsv"))

selected <- selection[selection$n_outer == selection$selected_n_outer, , drop = FALSE]
summary_keys <- unique(selected[, c("scenario", "variant", "selected_n_outer"), drop = FALSE])
summary_rows <- vector("list", nrow(summary_keys))
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}
for (i in seq_len(nrow(summary_keys))) {
  keep <- selected$scenario == summary_keys$scenario[i] &
    selected$variant == summary_keys$variant[i] &
    selected$selected_n_outer == summary_keys$selected_n_outer[i]
  df <- selected[keep, , drop = FALSE]
  summary_rows[[i]] <- data.frame(
    scenario = summary_keys$scenario[i],
    variant = summary_keys$variant[i],
    selected_n_outer = summary_keys$selected_n_outer[i],
    n_selected_rows = nrow(df),
    selection_score_mean = safe_mean(df$selection_score),
    heldout_nb_objective_mean = safe_mean(df$heldout_nb_objective),
    z_stability_delta_vs_previous_mean = safe_mean(df$z_stability_delta_vs_previous),
    marker_gradient_mean_spearman_mean = safe_mean(df$marker_gradient_mean_spearman),
    downstream_balanced_accuracy_mean = safe_mean(df$downstream_balanced_accuracy),
    truth_rmse_mean = safe_mean(df$truth_rmse),
    stringsAsFactors = FALSE
  )
}
selected_summary <- do.call(rbind, summary_rows)
write_tsv(selected_summary, file.path(OUT, "validation_aware_nb_selected_summary.tsv"))

manifest <- data.frame(
  analysis = "validation_aware_nb_selection",
  primary_claim = "FibroDynMix NB outer-iteration budgets are selected by held-out likelihood, z stability, marker-gradient preservation, and downstream validation rather than training objective alone.",
  claim_boundary = "Simulation-calibrated selector. Truth metrics are audit columns only and are not used in the validation-aware selection score.",
  scenarios = paste(scenarios, collapse = ";"),
  variants = paste(variants, collapse = ";"),
  n_replicates = n_replicates,
  candidate_n_outer = paste(candidate_n_outer, collapse = ";"),
  n_rows = nrow(selection),
  selected_n_outer_values = paste(sort(unique(selected$selected_n_outer)), collapse = ";"),
  holdout_fraction = holdout_fraction,
  cells_per_donor = cells_per_donor,
  n_genes = n_genes,
  seed = seed,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "validation_aware_nb_selection_manifest.tsv"))

message("Validation-aware NB selection written to: ", normalizePath(OUT))
