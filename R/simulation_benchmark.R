#' Run simulation benchmarks across scenarios
#'
#' Runs repeated FibroDynMix simulations and evaluates the marker-scoring
#' baseline across one or more scenarios.
#'
#' @param scenarios Character vector of simulation scenarios.
#' @param n_replicates Number of replicates per scenario.
#' @param seed Optional base seed. Replicate seeds are deterministically derived
#'   from this value.
#' @param methods Character vector of methods to run. Implemented methods are
#'   `marker_scoring`, `seurat_cluster`, `topic_nmf`, `fibrodynmix_initializer`,
#'   `fibrodynmix_nb`, `fibrodynmix_nb_study`, `fibrodynmix_nb_donor`,
#'   `fibrodynmix_nb_study_donor`, and `fibrodynmix_vi`.
#' @param simulation_args Named list of arguments passed to
#'   `simulate_fibrodynmix()`.
#' @param temperature Softmax temperature for the marker-scoring baseline.
#' @param initializer_args Named list of arguments passed to
#'   `fit_fibrodynmix_initializer()`.
#' @param topic_nmf_args Named list of arguments passed to
#'   `fit_topic_nmf_baseline()`.
#' @param scvi_args Named list of arguments passed to the optional
#'   `fit_scvi_latent_baseline()`.
#' @param seurat_cluster_args Named list of arguments passed to the Seurat
#'   cluster baseline.
#' @param nb_args Named list of arguments passed to `fit_fibrodynmix_nb()`.
#' @param vi_args Named list of arguments passed to `fit_fibrodynmix_vi()`.
#' @param downstream_task Whether to evaluate a downstream label-prediction
#'   benchmark from each inferred state-feature matrix.
#' @param downstream_label_col Column in simulated cell metadata used as the
#'   downstream prediction label.
#' @param downstream_group_col Optional cell metadata column used for grouped
#'   cross-validation to avoid donor leakage.
#' @param downstream_folds Requested number of downstream CV folds.
#' @param keep_results Whether to retain full simulation and baseline objects.
#'
#' @return A data frame with one row per scenario, replicate, and method. If
#'   `keep_results = TRUE`, full run objects are attached as the `results`
#'   attribute.
#' @export
run_simulation_benchmark <- function(scenarios = c(
                                       "continuous",
                                       "discrete",
                                       "batch_confounding",
                                       "rare_transition"
                                     ),
                                     n_replicates = 3,
                                     seed = 1,
                                     methods = "marker_scoring",
                                     simulation_args = list(),
                                     temperature = 1,
                                     initializer_args = list(),
                                     topic_nmf_args = list(),
                                     scvi_args = list(),
                                     seurat_cluster_args = list(),
                                     nb_args = list(),
                                     vi_args = list(),
                                     downstream_task = TRUE,
                                     downstream_label_col = "disease",
                                     downstream_group_col = "donor_id",
                                     downstream_folds = 3,
                                     keep_results = FALSE) {
  if (length(scenarios) == 0L) {
    stop("`scenarios` must contain at least one scenario.", call. = FALSE)
  }
  assert_positive_integer(n_replicates, "n_replicates")
  methods <- match.arg(
    methods,
    choices = c("marker_scoring", "seurat_cluster", "topic_nmf", "scvi_latent", "fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_nb_study", "fibrodynmix_nb_donor", "fibrodynmix_nb_study_donor", "fibrodynmix_vi"),
    several.ok = TRUE
  )
  if (!isTRUE(downstream_task) && !identical(downstream_task, FALSE)) {
    stop("`downstream_task` must be TRUE or FALSE.", call. = FALSE)
  }
  assert_positive_integer(downstream_folds, "downstream_folds")

  rows <- list()
  full_results <- list()
  row_index <- 1L

  for (scenario in scenarios) {
    scenario <- match.arg(
      scenario,
      choices = c("continuous", "discrete", "batch_confounding", "rare_transition")
    )

    for (replicate_id in seq_len(n_replicates)) {
      replicate_seed <- if (is.null(seed)) {
        NULL
      } else {
        as.integer(seed + scenario_seed_offset(scenario) + replicate_id - 1L)
      }

      sim_args <- modifyList(
        list(scenario = scenario, seed = replicate_seed),
        simulation_args
      )
      sim_args$scenario <- scenario
      if (!is.null(seed)) {
        sim_args$seed <- replicate_seed
      }
      sim <- do.call(simulate_fibrodynmix, sim_args)

      for (method in methods) {
        if (method == "marker_scoring") {
          baseline <- score_marker_baseline(
            counts = sim$counts,
            marker_index = sim$parameters$marker_index,
            library_size = sim$cell_metadata$library_size,
            temperature = temperature
          )
          metrics <- evaluate_state_weights(sim$z, baseline$z_pred)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = baseline$z_pred,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = NULL,
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              baseline = baseline,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "seurat_cluster") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              seed = replicate_seed
            ),
            seurat_cluster_args
          )
          baseline <- do.call(fit_seurat_cluster_baseline, fit_args)
          metrics <- evaluate_state_weights(sim$z, baseline$z_pred)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = baseline$z_pred,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = NULL,
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              baseline = baseline,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "topic_nmf") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              n_topics = length(sim$parameters$marker_index),
              n_iter = 60,
              seed = replicate_seed
            ),
            topic_nmf_args
          )
          fit <- do.call(fit_topic_nmf_baseline, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_pred)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_pred,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = NULL,
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "scvi_latent") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              n_latent = length(sim$parameters$marker_index),
              max_epochs = 40,
              seed = replicate_seed
            ),
            scvi_args
          )
          fit <- do.call(fit_scvi_latent_baseline, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_pred)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_pred,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = NULL,
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_initializer") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size
            ),
            initializer_args
          )
          fit <- do.call(fit_fibrodynmix_initializer, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_hat)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_hat,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = NULL,
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_nb") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size
            ),
            nb_args
          )
          fit <- do.call(fit_fibrodynmix_nb, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_hat)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_hat,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = nb_optimizer_diagnostics(fit),
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_nb_study") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              study_id = sim$cell_metadata$study_id,
              fit_study_effect = TRUE
            ),
            nb_args
          )
          fit <- do.call(fit_fibrodynmix_nb, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_hat)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_hat,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = nb_optimizer_diagnostics(fit),
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_nb_donor") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              donor_id = sim$cell_metadata$donor_id,
              fit_donor_effect = TRUE
            ),
            nb_args
          )
          fit <- do.call(fit_fibrodynmix_nb, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_hat)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_hat,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = nb_optimizer_diagnostics(fit),
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_nb_study_donor") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              study_id = sim$cell_metadata$study_id,
              donor_id = sim$cell_metadata$donor_id,
              fit_study_effect = TRUE,
              fit_donor_effect = TRUE
            ),
            nb_args
          )
          fit <- do.call(fit_fibrodynmix_nb, fit_args)
          metrics <- evaluate_state_weights(sim$z, fit$z_hat)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = fit$z_hat,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = nb_optimizer_diagnostics(fit),
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = fit,
              metrics = metrics,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
        if (method == "fibrodynmix_vi") {
          fit_args <- modifyList(
            list(
              counts = sim$counts,
              marker_index = sim$parameters$marker_index,
              library_size = sim$cell_metadata$library_size,
              nb_args = nb_args,
              n_draws = 30,
              n_elbo_draws = 6,
              n_vi_iter = 2,
              keep_draws = FALSE
            ),
            vi_args
          )
          vi_fit <- do.call(fit_fibrodynmix_vi, fit_args)
          metrics <- evaluate_state_weights(sim$z, vi_fit$z_mean)
          downstream_metrics <- benchmark_downstream_diagnostics(
            sim = sim,
            z_pred = vi_fit$z_mean,
            enabled = downstream_task,
            label_col = downstream_label_col,
            group_col = downstream_group_col,
            n_folds = downstream_folds,
            seed = replicate_seed
          )
          posterior_metrics <- evaluate_posterior_intervals(sim$z, vi_fit$cell_summary$z)
          posterior_calibration <- calibrate_posterior_interval_scale(sim$z, vi_fit$cell_summary$z)

          rows[[row_index]] <- benchmark_metric_row(
            scenario = scenario,
            replicate_id = replicate_id,
            seed = replicate_seed,
            method = method,
            metrics = metrics,
            diagnostics = nb_optimizer_diagnostics(vi_fit$nb_fit),
            posterior_diagnostics = vi_posterior_diagnostics(vi_fit, posterior_metrics, posterior_calibration),
            downstream_diagnostics = downstream_metrics
          )
          if (keep_results) {
            full_results[[row_index]] <- list(
              simulation = sim,
              fit = vi_fit,
              metrics = metrics,
              posterior_metrics = posterior_metrics,
              posterior_calibration = posterior_calibration,
              downstream_metrics = downstream_metrics
            )
          }
          row_index <- row_index + 1L
        }
      }
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (keep_results) {
    attr(result, "results") <- full_results
  }
  result
}

#' Summarize benchmark results
#'
#' Aggregates replicate-level benchmark rows by scenario and method.
#'
#' @param benchmark_results Data frame returned by `run_simulation_benchmark()`.
#' @param metrics Character vector of metric columns to summarize.
#'
#' @return A data frame with mean, standard deviation, and replicate count for
#'   each metric.
#' @export
summarize_benchmark_results <- function(benchmark_results,
                                        metrics = c(
                                          "rmse",
                                          "mean_absolute_error",
                                          "dominant_accuracy",
                                          "mean_entropy_pred"
                                        )) {
  required <- c("scenario", "method", metrics)
  missing <- setdiff(required, colnames(benchmark_results))
  if (length(missing) > 0L) {
    stop(
      sprintf("`benchmark_results` is missing columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  groups <- unique(benchmark_results[, c("scenario", "method"), drop = FALSE])
  summary_rows <- vector("list", nrow(groups))

  for (i in seq_len(nrow(groups))) {
    keep <- benchmark_results$scenario == groups$scenario[i] &
      benchmark_results$method == groups$method[i]
    subset <- benchmark_results[keep, , drop = FALSE]

    row <- data.frame(
      scenario = groups$scenario[i],
      method = groups$method[i],
      n_replicates = nrow(subset),
      stringsAsFactors = FALSE
    )

    for (metric in metrics) {
      values <- subset[[metric]]
      row[[paste0(metric, "_mean")]] <- mean(values, na.rm = TRUE)
      row[[paste0(metric, "_sd")]] <- stats::sd(values, na.rm = TRUE)
    }

    summary_rows[[i]] <- row
  }

  result <- do.call(rbind, summary_rows)
  rownames(result) <- NULL
  result
}

#' Summarize FibroDynMix NB optimizer diagnostics
#'
#' Aggregates NB optimizer diagnostics emitted by `run_simulation_benchmark()`.
#' Non-NB methods are retained with `NA` diagnostics so benchmark tables remain
#' rectangular.
#'
#' @param benchmark_results Data frame returned by `run_simulation_benchmark()`.
#'
#' @return A data frame summarizing objective improvement, rollback frequency,
#'   early stopping frequency, and best iteration by scenario and method.
#' @export
summarize_optimizer_diagnostics <- function(benchmark_results) {
  required <- c(
    "scenario",
    "method",
    "nb_initial_objective",
    "nb_final_objective",
    "nb_best_objective",
    "nb_best_iteration",
    "nb_executed_iterations",
    "nb_stop_reason",
    "nb_any_rollback",
    "nb_rollback_count",
    "nb_study_effect_l2_norm",
    "nb_donor_effect_l2_norm"
  )
  missing <- setdiff(required, colnames(benchmark_results))
  if (length(missing) > 0L) {
    stop(
      sprintf("`benchmark_results` is missing columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  groups <- unique(benchmark_results[, c("scenario", "method"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    keep <- benchmark_results$scenario == groups$scenario[i] &
      benchmark_results$method == groups$method[i]
    subset <- benchmark_results[keep, , drop = FALSE]
    improvement <- subset$nb_initial_objective - subset$nb_best_objective
    n_nb_runs <- sum(!is.na(subset$nb_best_objective))
    rows[[i]] <- data.frame(
      scenario = groups$scenario[i],
      method = groups$method[i],
      n_replicates = nrow(subset),
      n_nb_runs = n_nb_runs,
      objective_improvement_mean = safe_mean(improvement),
      objective_improvement_sd = safe_sd(improvement),
      rollback_rate = safe_mean(subset$nb_any_rollback),
      rollback_count_mean = safe_mean(subset$nb_rollback_count),
      study_effect_l2_norm_mean = safe_mean(subset$nb_study_effect_l2_norm),
      donor_effect_l2_norm_mean = safe_mean(subset$nb_donor_effect_l2_norm),
      early_stopping_rate = safe_mean(subset$nb_stop_reason == "early_stopping"),
      median_best_iteration = safe_median(subset$nb_best_iteration),
      median_executed_iterations = safe_median(subset$nb_executed_iterations),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::median(x)
}

benchmark_metric_row <- function(scenario, replicate_id, seed, method, metrics, diagnostics = NULL, posterior_diagnostics = NULL, downstream_diagnostics = NULL) {
  if (is.null(diagnostics)) {
    diagnostics <- empty_nb_diagnostics()
  }
  if (is.null(posterior_diagnostics)) {
    posterior_diagnostics <- empty_vi_diagnostics()
  }
  if (is.null(downstream_diagnostics)) {
    downstream_diagnostics <- empty_downstream_diagnostics()
  }
  data.frame(
    scenario = scenario,
    replicate = replicate_id,
    seed = if (is.null(seed)) NA_integer_ else seed,
    method = method,
    rmse = metrics$rmse,
    mean_absolute_error = metrics$mean_absolute_error,
    dominant_accuracy = metrics$dominant_accuracy,
    mean_entropy_true = metrics$mean_entropy_true,
    mean_entropy_pred = metrics$mean_entropy_pred,
    nb_initial_objective = diagnostics$nb_initial_objective,
    nb_final_objective = diagnostics$nb_final_objective,
    nb_best_objective = diagnostics$nb_best_objective,
    nb_objective_improvement = diagnostics$nb_initial_objective - diagnostics$nb_best_objective,
    nb_best_iteration = diagnostics$nb_best_iteration,
    nb_executed_iterations = diagnostics$nb_executed_iterations,
    nb_stop_reason = diagnostics$nb_stop_reason,
    nb_any_rollback = diagnostics$nb_any_rollback,
    nb_rollback_count = diagnostics$nb_rollback_count,
    nb_study_effect_l2_norm = diagnostics$nb_study_effect_l2_norm,
    nb_donor_effect_l2_norm = diagnostics$nb_donor_effect_l2_norm,
    vi_best_elbo = posterior_diagnostics$vi_best_elbo,
    vi_interval_coverage = posterior_diagnostics$vi_interval_coverage,
    vi_mean_interval_width = posterior_diagnostics$vi_mean_interval_width,
    vi_median_interval_width = posterior_diagnostics$vi_median_interval_width,
    vi_posterior_mae = posterior_diagnostics$vi_posterior_mae,
    vi_n_interval_rows = posterior_diagnostics$vi_n_interval_rows,
    vi_calibrated_interval_scale = posterior_diagnostics$vi_calibrated_interval_scale,
    vi_calibrated_interval_coverage = posterior_diagnostics$vi_calibrated_interval_coverage,
    vi_calibrated_mean_interval_width = posterior_diagnostics$vi_calibrated_mean_interval_width,
    downstream_status = downstream_diagnostics$downstream_status,
    downstream_label_col = downstream_diagnostics$downstream_label_col,
    downstream_group_col = downstream_diagnostics$downstream_group_col,
    downstream_n_observations = downstream_diagnostics$downstream_n_observations,
    downstream_n_evaluated = downstream_diagnostics$downstream_n_evaluated,
    downstream_n_classes = downstream_diagnostics$downstream_n_classes,
    downstream_n_folds = downstream_diagnostics$downstream_n_folds,
    downstream_accuracy = downstream_diagnostics$downstream_accuracy,
    downstream_balanced_accuracy = downstream_diagnostics$downstream_balanced_accuracy,
    downstream_macro_f1 = downstream_diagnostics$downstream_macro_f1,
    downstream_macro_auroc = downstream_diagnostics$downstream_macro_auroc,
    stringsAsFactors = FALSE
  )
}

empty_nb_diagnostics <- function() {
  list(
    nb_initial_objective = NA_real_,
    nb_final_objective = NA_real_,
    nb_best_objective = NA_real_,
    nb_best_iteration = NA_integer_,
    nb_executed_iterations = NA_integer_,
    nb_stop_reason = NA_character_,
    nb_any_rollback = NA,
    nb_rollback_count = NA_integer_,
    nb_study_effect_l2_norm = NA_real_,
    nb_donor_effect_l2_norm = NA_real_
  )
}

empty_vi_diagnostics <- function() {
  list(
    vi_best_elbo = NA_real_,
    vi_interval_coverage = NA_real_,
    vi_mean_interval_width = NA_real_,
    vi_median_interval_width = NA_real_,
    vi_posterior_mae = NA_real_,
    vi_n_interval_rows = NA_integer_,
    vi_calibrated_interval_scale = NA_real_,
    vi_calibrated_interval_coverage = NA_real_,
    vi_calibrated_mean_interval_width = NA_real_
  )
}

empty_downstream_diagnostics <- function(status = "not_run",
                                         label_col = NA_character_,
                                         group_col = NA_character_) {
  list(
    downstream_status = status,
    downstream_label_col = label_col,
    downstream_group_col = group_col,
    downstream_n_observations = NA_integer_,
    downstream_n_evaluated = NA_integer_,
    downstream_n_classes = NA_integer_,
    downstream_n_folds = NA_integer_,
    downstream_accuracy = NA_real_,
    downstream_balanced_accuracy = NA_real_,
    downstream_macro_f1 = NA_real_,
    downstream_macro_auroc = NA_real_
  )
}

benchmark_downstream_diagnostics <- function(sim,
                                             z_pred,
                                             enabled,
                                             label_col,
                                             group_col,
                                             n_folds,
                                             seed = NULL) {
  if (!enabled) {
    return(empty_downstream_diagnostics(status = "disabled", label_col = label_col, group_col = group_col))
  }
  metadata <- sim$cell_metadata
  if (!label_col %in% colnames(metadata)) {
    return(empty_downstream_diagnostics(status = "missing_label", label_col = label_col, group_col = group_col))
  }
  groups <- NULL
  group_label <- NA_character_
  if (!is.null(group_col) && !is.na(group_col) && nzchar(group_col) && group_col %in% colnames(metadata)) {
    groups <- metadata[[group_col]]
    group_label <- group_col
  }
  metrics <- tryCatch(
    evaluate_downstream_classification(
      features = z_pred,
      labels = metadata[[label_col]],
      groups = groups,
      n_folds = n_folds,
      seed = seed
    ),
    error = function(e) {
      structure(
        empty_classification_metrics(
          status = paste0("error:", conditionMessage(e)),
          n_observations = nrow(z_pred),
          n_classes = length(unique(metadata[[label_col]]))
        ),
        class = "downstream_error"
      )
    }
  )
  list(
    downstream_status = metrics$status,
    downstream_label_col = label_col,
    downstream_group_col = group_label,
    downstream_n_observations = metrics$n_observations,
    downstream_n_evaluated = metrics$n_evaluated,
    downstream_n_classes = metrics$n_classes,
    downstream_n_folds = metrics$n_folds,
    downstream_accuracy = metrics$accuracy,
    downstream_balanced_accuracy = metrics$balanced_accuracy,
    downstream_macro_f1 = metrics$macro_f1,
    downstream_macro_auroc = metrics$macro_auroc
  )
}

vi_posterior_diagnostics <- function(vi_fit, posterior_metrics, posterior_calibration = NULL) {
  if (is.null(posterior_calibration)) {
    posterior_calibration <- list(
      selected_scale = NA_real_,
      interval_coverage = NA_real_,
      mean_interval_width = NA_real_
    )
  }
  list(
    vi_best_elbo = vi_fit$best_elbo,
    vi_interval_coverage = posterior_metrics$interval_coverage,
    vi_mean_interval_width = posterior_metrics$mean_interval_width,
    vi_median_interval_width = posterior_metrics$median_interval_width,
    vi_posterior_mae = posterior_metrics$posterior_mean_absolute_error,
    vi_n_interval_rows = posterior_metrics$n_interval_rows,
    vi_calibrated_interval_scale = posterior_calibration$selected_scale,
    vi_calibrated_interval_coverage = posterior_calibration$interval_coverage,
    vi_calibrated_mean_interval_width = posterior_calibration$mean_interval_width
  )
}

nb_optimizer_diagnostics <- function(fit) {
  trace <- fit$nb_objective_trace
  rolled_back <- fit$convergence$rolled_back
  list(
    nb_initial_objective = trace[1L],
    nb_final_objective = trace[length(trace)],
    nb_best_objective = fit$best_objective,
    nb_best_iteration = fit$best_iteration,
    nb_executed_iterations = fit$executed_iterations,
    nb_stop_reason = fit$stop_reason,
    nb_any_rollback = any(rolled_back),
    nb_rollback_count = sum(rolled_back),
    nb_study_effect_l2_norm = if (is.null(fit$study_effect)) NA_real_ else sqrt(sum(fit$study_effect^2)),
    nb_donor_effect_l2_norm = if (is.null(fit$donor_effect)) NA_real_ else sqrt(sum(fit$donor_effect^2))
  )
}

fit_seurat_cluster_baseline <- function(counts,
                                        marker_index,
                                        library_size = NULL,
                                        n_pcs = 20,
                                        resolution = 0.4,
                                        seed = NULL,
                                        assay = "RNA") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("`method = 'seurat_cluster'` requires the Seurat package.", call. = FALSE)
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("`method = 'seurat_cluster'` requires the SeuratObject package.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is_matrix_like(counts) || !matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must be a non-negative integer-like matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  object <- Seurat::CreateSeuratObject(counts = counts, assay = assay)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = min(2000L, nrow(counts)), verbose = FALSE)
  object <- Seurat::ScaleData(object, verbose = FALSE)
  n_pcs <- min(as.integer(n_pcs), ncol(object) - 1L, nrow(object) - 1L)
  if (n_pcs < 2L) {
    stop("Seurat cluster baseline requires at least two PCs.", call. = FALSE)
  }
  object <- Seurat::RunPCA(object, npcs = n_pcs, verbose = FALSE)
  object <- Seurat::FindNeighbors(object, dims = seq_len(n_pcs), verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = FALSE)

  marker_scores <- score_marker_baseline(
    counts = counts,
    marker_index = marker_index,
    library_size = library_size
  )$scores
  clusters <- as.character(object$seurat_clusters)
  state_names <- colnames(marker_scores)
  cluster_levels <- sort(unique(clusters))
  cluster_state_score <- matrix(
    NA_real_,
    nrow = length(cluster_levels),
    ncol = length(state_names),
    dimnames = list(cluster_levels, state_names)
  )
  for (cluster in cluster_levels) {
    keep <- clusters == cluster
    cluster_state_score[cluster, ] <- colMeans(marker_scores[keep, , drop = FALSE])
  }
  cluster_state <- state_names[max.col(cluster_state_score, ties.method = "first")]
  names(cluster_state) <- rownames(cluster_state_score)
  z_pred <- matrix(
    0,
    nrow = ncol(counts),
    ncol = length(state_names),
    dimnames = list(colnames(counts), state_names)
  )
  assigned_state <- cluster_state[clusters]
  z_pred[cbind(seq_len(nrow(z_pred)), match(assigned_state, state_names))] <- 1

  list(
    z_pred = z_pred,
    clusters = clusters,
    cluster_state = cluster_state,
    cluster_state_score = cluster_state_score,
    object = object
  )
}

scenario_seed_offset <- function(scenario) {
  switch(
    scenario,
    continuous = 0L,
    discrete = 10000L,
    batch_confounding = 20000L,
    rare_transition = 30000L
  )
}
