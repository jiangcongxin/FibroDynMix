#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_marker_stress_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "marker_stress_benchmark")
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
  "vi_posterior.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

rbind_fill <- function(rows) {
  all_cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_cols, colnames(row))
    for (col in missing) {
      row[[col]] <- NA
    }
    row[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

perturb_marker_index <- function(marker_index,
                                 gene_names,
                                 stress_mode,
                                 seed = NULL,
                                 corruption_fraction = 0.5,
                                 retention_fraction = 0.35,
                                 shared_n = 3L,
                                 swap_fraction = 0.5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  marker_index <- lapply(marker_index, function(x) {
    if (is.numeric(x)) gene_names[x] else as.character(x)
  })
  names(marker_index) <- names(marker_index)
  true_markers <- unique(unlist(marker_index, use.names = FALSE))
  nonmarkers <- setdiff(gene_names, true_markers)
  out <- marker_index

  if (stress_mode == "corrupted_prior") {
    out <- lapply(marker_index, function(markers) {
      n_replace <- max(1L, floor(length(markers) * corruption_fraction))
      keep <- setdiff(markers, sample(markers, n_replace))
      decoys <- sample(nonmarkers, n_replace)
      unique(c(keep, decoys))
    })
  } else if (stress_mode == "missing_markers") {
    out <- lapply(marker_index, function(markers) {
      n_keep <- max(1L, floor(length(markers) * retention_fraction))
      sample(markers, n_keep)
    })
  } else if (stress_mode == "shared_markers") {
    shared <- sample(nonmarkers, min(shared_n, length(nonmarkers)))
    out <- lapply(marker_index, function(markers) unique(c(markers, shared)))
  } else if (stress_mode == "swapped_markers") {
    states <- names(marker_index)
    out <- marker_index
    for (i in seq_along(states)) {
      donor_state <- states[ifelse(i == length(states), 1L, i + 1L)]
      markers <- marker_index[[states[i]]]
      n_swap <- max(1L, floor(length(markers) * swap_fraction))
      remove <- sample(markers, n_swap)
      add <- sample(marker_index[[donor_state]], n_swap)
      out[[states[i]]] <- unique(c(setdiff(markers, remove), add))
    }
  } else if (stress_mode == "clean_prior") {
    out <- marker_index
  } else {
    stop(sprintf("Unknown stress_mode: %s", stress_mode), call. = FALSE)
  }
  out
}

marker_prior_summary <- function(true_marker_index, perturbed_marker_index) {
  states <- names(true_marker_index)
  do.call(rbind, lapply(states, function(state) {
    truth <- as.character(true_marker_index[[state]])
    prior <- as.character(perturbed_marker_index[[state]])
    data.frame(
      state = state,
      n_true_markers = length(truth),
      n_prior_markers = length(prior),
      n_true_retained = length(intersect(truth, prior)),
      precision = length(intersect(truth, prior)) / max(1L, length(prior)),
      recall = length(intersect(truth, prior)) / max(1L, length(truth)),
      stringsAsFactors = FALSE
    )
  }))
}

inject_hidden_program <- function(sim,
                                  hidden_genes_per_state = 8L,
                                  effect_scale = 0.0007,
                                  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  counts <- sim$counts
  gene_names <- rownames(counts)
  marker_genes <- unique(unlist(lapply(sim$parameters$marker_index, function(idx) gene_names[idx]), use.names = FALSE))
  candidate_genes <- setdiff(gene_names, marker_genes)
  state_names <- colnames(sim$z)
  hidden_index <- vector("list", length(state_names))
  names(hidden_index) <- state_names
  cursor_candidates <- sample(candidate_genes)
  cursor <- 1L
  for (state in state_names) {
    idx <- cursor:(cursor + hidden_genes_per_state - 1L)
    hidden_index[[state]] <- cursor_candidates[idx]
    cursor <- cursor + hidden_genes_per_state
    lambda <- outer(
      rep(1, length(hidden_index[[state]])),
      sim$cell_metadata$library_size * sim$z[, state] * effect_scale
    )
    added <- matrix(
      stats::rpois(length(lambda), lambda = as.numeric(lambda)),
      nrow = length(hidden_index[[state]]),
      dimnames = list(hidden_index[[state]], colnames(counts))
    )
    counts[hidden_index[[state]], ] <- counts[hidden_index[[state]], , drop = FALSE] + added
  }
  sim$counts <- counts
  sim$parameters$hidden_program_index <- hidden_index
  sim$parameters$hidden_program_effect_scale <- effect_scale
  sim
}

metric_row <- function(stress_mode, replicate, seed, method, metrics, extra = list()) {
  base <- data.frame(
    stress_mode = stress_mode,
    replicate = replicate,
    seed = seed,
    method = method,
    rmse = metrics$rmse,
    mean_absolute_error = metrics$mean_absolute_error,
    dominant_accuracy = metrics$dominant_accuracy,
    mean_entropy_true = metrics$mean_entropy_true,
    mean_entropy_pred = metrics$mean_entropy_pred,
    stringsAsFactors = FALSE
  )
  if (length(extra) > 0L) {
    for (nm in names(extra)) {
      base[[nm]] <- extra[[nm]]
    }
  }
  base
}

stress_modes <- c(
  "clean_prior",
  "missing_markers",
  "corrupted_prior",
  "shared_markers",
  "swapped_markers",
  "hidden_program_corrupted_prior",
  "hidden_program_missing_markers"
)
methods <- c("marker_scoring", "topic_nmf", "fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_vi")
n_replicates <- 2L
base_seed <- 880
topic_backend <- if (requireNamespace("NMF", quietly = TRUE)) "nmf" else "multiplicative_update"

rows <- list()
prior_rows <- list()
row_idx <- 1L
prior_idx <- 1L

for (stress_mode in stress_modes) {
  for (replicate in seq_len(n_replicates)) {
    replicate_seed <- base_seed + match(stress_mode, stress_modes) * 1000L + replicate
    sim <- simulate_fibrodynmix(
      n_studies = 3,
      donors_per_study = 2,
      cells_per_donor = 8,
      n_genes = 120,
      marker_genes_per_state = 6,
      scenario = "batch_confounding",
      tau_high = 0.9,
      tau_low = 0.12,
      study_effect_sd = 0.2,
      donor_effect_sd = 0.1,
      seed = replicate_seed
    )
    if (grepl("^hidden_program", stress_mode)) {
      sim <- inject_hidden_program(
        sim,
        hidden_genes_per_state = 8L,
        effect_scale = 0.0009,
        seed = replicate_seed + 31L
      )
    }
    true_marker_index <- lapply(sim$parameters$marker_index, function(idx) rownames(sim$counts)[idx])
    marker_stress_mode <- switch(
      stress_mode,
      hidden_program_corrupted_prior = "corrupted_prior",
      hidden_program_missing_markers = "missing_markers",
      stress_mode
    )
    perturbed_marker_index <- perturb_marker_index(
      marker_index = true_marker_index,
      gene_names = rownames(sim$counts),
      stress_mode = marker_stress_mode,
      seed = replicate_seed + 17L
    )
    prior_summary <- marker_prior_summary(true_marker_index, perturbed_marker_index)
    prior_summary$stress_mode <- stress_mode
    prior_summary$replicate <- replicate
    prior_summary$seed <- replicate_seed
    prior_rows[[prior_idx]] <- prior_summary
    prior_idx <- prior_idx + 1L

    for (method in methods) {
      if (method == "marker_scoring") {
        fit <- score_marker_baseline(
          counts = sim$counts,
          marker_index = perturbed_marker_index,
          library_size = sim$cell_metadata$library_size,
          temperature = 1
        )
        metrics <- evaluate_state_weights(sim$z, fit$z_pred)
        rows[[row_idx]] <- metric_row(stress_mode, replicate, replicate_seed, method, metrics)
      } else if (method == "topic_nmf") {
        fit <- fit_topic_nmf_baseline(
          counts = sim$counts,
          marker_index = perturbed_marker_index,
          n_topics = length(perturbed_marker_index),
          backend = topic_backend,
          n_iter = 40,
          seed = replicate_seed
        )
        metrics <- evaluate_state_weights(sim$z, fit$z_pred)
        rows[[row_idx]] <- metric_row(stress_mode, replicate, replicate_seed, method, metrics)
      } else if (method == "fibrodynmix_initializer") {
        fit <- fit_fibrodynmix_initializer(
          counts = sim$counts,
          marker_index = perturbed_marker_index,
          library_size = sim$cell_metadata$library_size,
          n_iter = 3
        )
        metrics <- evaluate_state_weights(sim$z, fit$z_hat)
        rows[[row_idx]] <- metric_row(stress_mode, replicate, replicate_seed, method, metrics)
      } else if (method == "fibrodynmix_nb") {
        fit <- fit_fibrodynmix_nb(
          counts = sim$counts,
          marker_index = perturbed_marker_index,
          library_size = sim$cell_metadata$library_size,
          n_outer = 2,
          initializer_args = list(n_iter = 3),
          maxit_beta = 12,
          maxit_z = 12
        )
        metrics <- evaluate_state_weights(sim$z, fit$z_hat)
        rows[[row_idx]] <- metric_row(
          stress_mode, replicate, replicate_seed, method, metrics,
          extra = list(
            nb_initial_objective = fit$nb_objective_trace[1L],
            nb_best_objective = fit$best_objective,
            nb_objective_improvement = fit$nb_objective_trace[1L] - fit$best_objective,
            nb_stop_reason = fit$stop_reason
          )
        )
      } else if (method == "fibrodynmix_vi") {
        fit <- fit_fibrodynmix_vi(
          counts = sim$counts,
          marker_index = perturbed_marker_index,
          library_size = sim$cell_metadata$library_size,
          nb_args = list(
            n_outer = 2,
            initializer_args = list(n_iter = 3),
            maxit_beta = 12,
            maxit_z = 12
          ),
          n_draws = 18,
          n_elbo_draws = 4,
          n_vi_iter = 2,
          seed = replicate_seed + 101L,
          keep_draws = FALSE
        )
        metrics <- evaluate_state_weights(sim$z, fit$z_mean)
        posterior <- evaluate_posterior_intervals(sim$z, fit$cell_summary$z)
        calibration <- calibrate_posterior_interval_scale(sim$z, fit$cell_summary$z)
        rows[[row_idx]] <- metric_row(
          stress_mode, replicate, replicate_seed, method, metrics,
          extra = list(
            nb_initial_objective = fit$nb_fit$nb_objective_trace[1L],
            nb_best_objective = fit$nb_fit$best_objective,
            nb_objective_improvement = fit$nb_fit$nb_objective_trace[1L] - fit$nb_fit$best_objective,
            nb_stop_reason = fit$nb_fit$stop_reason,
            vi_interval_coverage = posterior$coverage,
            vi_calibrated_interval_coverage = calibration$calibrated_coverage,
            vi_calibrated_interval_scale = calibration$scale
          )
        )
      }
      row_idx <- row_idx + 1L
    }
  }
}

benchmark <- rbind_fill(rows)
prior_audit <- do.call(rbind, prior_rows)
write_tsv(benchmark, file.path(OUT, "marker_stress_benchmark_metrics.tsv"))
write_tsv(prior_audit, file.path(OUT, "marker_stress_prior_audit.tsv"))

summary <- aggregate(
  cbind(rmse, mean_absolute_error, dominant_accuracy, mean_entropy_pred) ~ stress_mode + method,
  data = benchmark,
  FUN = mean
)
colnames(summary)[-(1:2)] <- paste0(colnames(summary)[-(1:2)], "_mean")
write_tsv(summary, file.path(OUT, "marker_stress_benchmark_summary.tsv"))

rankings <- do.call(rbind, lapply(split(summary, summary$stress_mode), function(df) {
  df <- df[order(df$rmse_mean), , drop = FALSE]
  df$rmse_rank <- seq_len(nrow(df))
  df
}))
rankings <- rankings[, c("stress_mode", "rmse_rank", "method", "rmse_mean", "dominant_accuracy_mean", "mean_entropy_pred_mean"), drop = FALSE]
write_tsv(rankings, file.path(OUT, "marker_stress_rankings.tsv"))

baseline <- summary[summary$method == "marker_scoring", c("stress_mode", "rmse_mean", "dominant_accuracy_mean"), drop = FALSE]
colnames(baseline) <- c("stress_mode", "marker_rmse_mean", "marker_dominant_accuracy_mean")
contrast <- merge(summary, baseline, by = "stress_mode", sort = FALSE)
contrast$rmse_delta_vs_marker <- contrast$rmse_mean - contrast$marker_rmse_mean
contrast$dominant_accuracy_delta_vs_marker <- contrast$dominant_accuracy_mean - contrast$marker_dominant_accuracy_mean
write_tsv(contrast, file.path(OUT, "marker_stress_contrast_vs_marker.tsv"))

fibro_methods <- c("fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_vi")
fibro_rows <- contrast[contrast$method %in% fibro_methods, , drop = FALSE]
manifest <- data.frame(
  analysis = "marker_stress_benchmark",
  primary_claim = "FibroDynMix was stress-tested against marker scoring and NMF/topic baselines under deliberately imperfect marker priors.",
  claim_boundary = "Stress benchmark uses simulated truth and bounded replicates. It tests prior failure modes, not universal superiority across all simulation designs.",
  n_rows = nrow(benchmark),
  n_stress_modes = length(unique(benchmark$stress_mode)),
  n_methods = length(unique(benchmark$method)),
  n_replicates = n_replicates,
  topic_backend = topic_backend,
  n_modes_where_any_fibrodynmix_beats_marker_rmse = length(unique(fibro_rows$stress_mode[fibro_rows$rmse_delta_vs_marker < 0])),
  n_modes_where_any_fibrodynmix_beats_marker_accuracy = length(unique(fibro_rows$stress_mode[fibro_rows$dominant_accuracy_delta_vs_marker > 0])),
  mean_topic_nmf_rmse = mean(summary$rmse_mean[summary$method == "topic_nmf"]),
  mean_marker_rmse = mean(summary$rmse_mean[summary$method == "marker_scoring"]),
  mean_fibrodynmix_rmse = mean(summary$rmse_mean[summary$method %in% fibro_methods]),
  seed = base_seed,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "marker_stress_benchmark_manifest.tsv"))

message("Marker stress benchmark written to: ", normalizePath(OUT))
