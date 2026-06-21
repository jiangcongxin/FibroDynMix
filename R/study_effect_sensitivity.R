#' Run study-effect penalty sensitivity benchmarks
#'
#' Evaluates how ridge penalties on study effects and marker orientation affect
#' FibroDynMix NB performance in a batch-confounded simulation.
#'
#' @param study_l2_grid Numeric vector of study-effect ridge penalties.
#' @param marker_l2_grid Numeric vector of marker-orientation penalties.
#' @param n_replicates Number of simulation replicates per parameter setting.
#' @param seed Optional base seed.
#' @param simulation_args Named list of arguments passed to
#'   `simulate_fibrodynmix()`.
#' @param nb_args Named list of additional arguments passed to
#'   `fit_fibrodynmix_nb()`.
#' @param include_no_study Whether to include a non-study-effect NB baseline for
#'   each replicate.
#'
#' @return A data frame with one row per replicate and parameter setting.
#' @export
run_study_effect_sensitivity <- function(study_l2_grid = c(0.01, 0.05, 0.1, 0.5, 1, 5),
                                         marker_l2_grid = c(0.01, 0.05, 0.1),
                                         n_replicates = 2,
                                         seed = 1,
                                         simulation_args = list(),
                                         nb_args = list(),
                                         include_no_study = TRUE) {
  validate_nonnegative_grid(study_l2_grid, "study_l2_grid")
  validate_nonnegative_grid(marker_l2_grid, "marker_l2_grid")
  assert_positive_integer(n_replicates, "n_replicates")

  rows <- list()
  row_index <- 1L

  for (replicate_id in seq_len(n_replicates)) {
    replicate_seed <- if (is.null(seed)) NULL else as.integer(seed + replicate_id - 1L)
    sim_args <- utils::modifyList(
      list(
        scenario = "batch_confounding",
        seed = replicate_seed
      ),
      simulation_args
    )
    sim_args$scenario <- "batch_confounding"
    if (!is.null(seed)) {
      sim_args$seed <- replicate_seed
    }
    sim <- do.call(simulate_fibrodynmix, sim_args)

    if (isTRUE(include_no_study)) {
      fit_args <- utils::modifyList(
        list(
          counts = sim$counts,
          marker_index = sim$parameters$marker_index,
          library_size = sim$cell_metadata$library_size,
          fit_study_effect = FALSE
        ),
        nb_args
      )
      fit <- do.call(fit_fibrodynmix_nb, fit_args)
      metrics <- evaluate_state_weights(sim$z, fit$z_hat)
      rows[[row_index]] <- sensitivity_row(
        replicate_id = replicate_id,
        seed = replicate_seed,
        method = "fibrodynmix_nb",
        study_l2 = NA_real_,
        marker_l2 = fit_args$marker_l2 %||% 0.05,
        metrics = metrics,
        fit = fit
      )
      row_index <- row_index + 1L
    }

    for (study_l2 in study_l2_grid) {
      for (marker_l2 in marker_l2_grid) {
        fit_args <- utils::modifyList(
          list(
            counts = sim$counts,
            marker_index = sim$parameters$marker_index,
            library_size = sim$cell_metadata$library_size,
            study_id = sim$cell_metadata$study_id,
            fit_study_effect = TRUE,
            study_l2 = study_l2,
            marker_l2 = marker_l2
          ),
          nb_args
        )
        fit_args$fit_study_effect <- TRUE
        fit_args$study_id <- sim$cell_metadata$study_id
        fit_args$study_l2 <- study_l2
        fit_args$marker_l2 <- marker_l2

        fit <- do.call(fit_fibrodynmix_nb, fit_args)
        metrics <- evaluate_state_weights(sim$z, fit$z_hat)
        rows[[row_index]] <- sensitivity_row(
          replicate_id = replicate_id,
          seed = replicate_seed,
          method = "fibrodynmix_nb_study",
          study_l2 = study_l2,
          marker_l2 = marker_l2,
          metrics = metrics,
          fit = fit
        )
        row_index <- row_index + 1L
      }
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Select study-effect and marker penalties from a sensitivity benchmark
#'
#' Chooses a simulation-calibrated default penalty pair from
#' `run_study_effect_sensitivity()` output. The selection rule prioritizes low
#' RMSE while penalizing objective degradation, large study effects, and
#' optimizer rollbacks.
#'
#' @param sensitivity_results Data frame returned by
#'   `run_study_effect_sensitivity()`.
#' @param rmse_tolerance Relative tolerance around the best RMSE. Candidates
#'   within `best_rmse * (1 + rmse_tolerance)` are preferred.
#' @param objective_tolerance Relative tolerance around the best NB objective.
#' @param weights Named numeric vector with weights for `rmse`, `objective`,
#'   `study_effect`, and `rollback`.
#'
#' @return A list with `recommended_study_l2`, `recommended_marker_l2`,
#'   `selection_reason`, and `tradeoff_table`.
#' @export
select_study_effect_penalty <- function(sensitivity_results,
                                        rmse_tolerance = 0.05,
                                        objective_tolerance = 0.03,
                                        weights = c(
                                          rmse = 0.45,
                                          objective = 0.25,
                                          study_effect = 0.2,
                                          rollback = 0.1
                                        )) {
  required <- c(
    "method",
    "study_l2",
    "marker_l2",
    "rmse",
    "nb_best_objective",
    "study_effect_l2_norm",
    "nb_any_rollback"
  )
  missing <- setdiff(required, colnames(sensitivity_results))
  if (length(missing) > 0L) {
    stop(
      sprintf("`sensitivity_results` is missing columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (length(rmse_tolerance) != 1L || is.na(rmse_tolerance) || rmse_tolerance < 0) {
    stop("`rmse_tolerance` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(objective_tolerance) != 1L || is.na(objective_tolerance) || objective_tolerance < 0) {
    stop("`objective_tolerance` must be a non-negative numeric scalar.", call. = FALSE)
  }
  weights <- validate_selection_weights(weights)

  study_rows <- sensitivity_results[sensitivity_results$method == "fibrodynmix_nb_study", , drop = FALSE]
  if (nrow(study_rows) == 0L) {
    stop("`sensitivity_results` must contain `method == 'fibrodynmix_nb_study'` rows.", call. = FALSE)
  }

  groups <- unique(study_rows[, c("study_l2", "marker_l2"), drop = FALSE])
  tradeoff_rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    keep <- study_rows$study_l2 == groups$study_l2[i] &
      study_rows$marker_l2 == groups$marker_l2[i]
    subset <- study_rows[keep, , drop = FALSE]
    tradeoff_rows[[i]] <- data.frame(
      study_l2 = groups$study_l2[i],
      marker_l2 = groups$marker_l2[i],
      n_replicates = nrow(subset),
      rmse_mean = mean(subset$rmse, na.rm = TRUE),
      rmse_sd = safe_sd(subset$rmse),
      dominant_accuracy_mean = mean(subset$dominant_accuracy, na.rm = TRUE),
      nb_best_objective_mean = mean(subset$nb_best_objective, na.rm = TRUE),
      study_effect_l2_norm_mean = mean(subset$study_effect_l2_norm, na.rm = TRUE),
      study_effect_mean_abs_mean = mean(subset$study_effect_mean_abs, na.rm = TRUE),
      rollback_rate = mean(subset$nb_any_rollback, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  tradeoff_table <- do.call(rbind, tradeoff_rows)
  rownames(tradeoff_table) <- NULL

  best_rmse <- min(tradeoff_table$rmse_mean, na.rm = TRUE)
  best_objective <- min(tradeoff_table$nb_best_objective_mean, na.rm = TRUE)
  tradeoff_table$within_rmse_tolerance <- tradeoff_table$rmse_mean <= best_rmse * (1 + rmse_tolerance)
  tradeoff_table$within_objective_tolerance <- tradeoff_table$nb_best_objective_mean <= best_objective * (1 + objective_tolerance)

  tradeoff_table$rmse_scaled <- minmax_scale(tradeoff_table$rmse_mean)
  tradeoff_table$objective_scaled <- minmax_scale(tradeoff_table$nb_best_objective_mean)
  tradeoff_table$study_effect_scaled <- minmax_scale(tradeoff_table$study_effect_l2_norm_mean)
  tradeoff_table$rollback_scaled <- minmax_scale(tradeoff_table$rollback_rate)
  tradeoff_table$selection_score <-
    weights["rmse"] * tradeoff_table$rmse_scaled +
    weights["objective"] * tradeoff_table$objective_scaled +
    weights["study_effect"] * tradeoff_table$study_effect_scaled +
    weights["rollback"] * tradeoff_table$rollback_scaled

  preferred <- tradeoff_table[
    tradeoff_table$within_rmse_tolerance & tradeoff_table$within_objective_tolerance,
    ,
    drop = FALSE
  ]
  if (nrow(preferred) == 0L) {
    preferred <- tradeoff_table
  }
  preferred <- preferred[order(preferred$selection_score, preferred$rmse_mean, preferred$nb_best_objective_mean), , drop = FALSE]
  selected <- preferred[1L, , drop = FALSE]

  selection_reason <- sprintf(
    "Selected study_l2=%s and marker_l2=%s using weighted tradeoff score %.3f; RMSE %.4f is %s the %.1f%% tolerance around best RMSE %.4f, and NB objective %.4f is %s the %.1f%% tolerance around best objective %.4f.",
    format(selected$study_l2),
    format(selected$marker_l2),
    selected$selection_score,
    selected$rmse_mean,
    if (selected$within_rmse_tolerance) "within" else "outside",
    100 * rmse_tolerance,
    best_rmse,
    selected$nb_best_objective_mean,
    if (selected$within_objective_tolerance) "within" else "outside",
    100 * objective_tolerance,
    best_objective
  )

  list(
    recommended_study_l2 = selected$study_l2,
    recommended_marker_l2 = selected$marker_l2,
    selection_reason = selection_reason,
    tradeoff_table = tradeoff_table[order(tradeoff_table$selection_score), , drop = FALSE]
  )
}

sensitivity_row <- function(replicate_id, seed, method, study_l2, marker_l2, metrics, fit) {
  diagnostics <- nb_optimizer_diagnostics(fit)
  data.frame(
    scenario = "batch_confounding",
    replicate = replicate_id,
    seed = if (is.null(seed)) NA_integer_ else seed,
    method = method,
    study_l2 = study_l2,
    marker_l2 = marker_l2,
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
    study_effect_l2_norm = study_effect_l2_norm(fit$study_effect),
    study_effect_mean_abs = study_effect_mean_abs(fit$study_effect),
    stringsAsFactors = FALSE
  )
}

study_effect_l2_norm <- function(study_effect) {
  if (is.null(study_effect)) {
    return(NA_real_)
  }
  sqrt(sum(study_effect^2))
}

study_effect_mean_abs <- function(study_effect) {
  if (is.null(study_effect)) {
    return(NA_real_)
  }
  mean(abs(study_effect))
}

validate_nonnegative_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(x < 0)) {
    stop(sprintf("`%s` must be a non-empty non-negative numeric vector.", name), call. = FALSE)
  }
}

validate_selection_weights <- function(weights) {
  required <- c("rmse", "objective", "study_effect", "rollback")
  if (!is.numeric(weights) || anyNA(weights) || any(weights < 0) || !all(required %in% names(weights))) {
    stop(
      "`weights` must be a named non-negative numeric vector containing rmse, objective, study_effect, and rollback.",
      call. = FALSE
    )
  }
  weights <- weights[required]
  if (sum(weights) == 0) {
    stop("At least one value in `weights` must be positive.", call. = FALSE)
  }
  weights / sum(weights)
}

minmax_scale <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  range_x <- range(x, na.rm = TRUE)
  if (!is.finite(range_x[1]) || !is.finite(range_x[2]) || abs(diff(range_x)) < .Machine$double.eps) {
    return(rep(0, length(x)))
  }
  (x - range_x[1]) / diff(range_x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
