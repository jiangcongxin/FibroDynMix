#' Select a FibroDynMix NB outer-iteration budget with validation-aware criteria
#'
#' Fits candidate NB models on a training cell split, transfers each fitted
#' program to held-out cells, and ranks candidates by validation metrics rather
#' than by the training objective alone.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of marker genes for each state.
#' @param library_size Optional vector of per-cell library sizes.
#' @param candidate_n_outer Positive integer candidate outer-iteration budgets.
#' @param labels Optional cell labels used for downstream validation.
#' @param groups Optional grouping variable for the holdout split and grouped
#'   downstream cross-validation.
#' @param holdout_fraction Fraction of cells or groups held out.
#' @param seed Optional split seed.
#' @param nb_args Named list passed to `fit_fibrodynmix_nb()`.
#' @param transfer_args Named list passed to `fit_fibrodynmix_transfer()`.
#' @param selection_weights Named numeric vector for validation metrics.
#' @param truth_z Optional simulated truth used only for audit metrics, not
#'   selection.
#'
#' @return A list with `selected_n_outer`, `candidate_scores`,
#'   `selection_components`, and split metadata.
#' @export
select_fibrodynmix_nb_model <- function(counts,
                                        marker_index,
                                        library_size = NULL,
                                        candidate_n_outer = c(2L, 5L, 10L, 20L),
                                        labels = NULL,
                                        groups = NULL,
                                        holdout_fraction = 0.25,
                                        seed = NULL,
                                        nb_args = list(),
                                        transfer_args = list(),
                                        selection_weights = c(
                                          heldout_nb_objective = 0.35,
                                          z_stability = 0.25,
                                          marker_gradient = 0.20,
                                          downstream_balanced_accuracy = 0.20
                                        ),
                                        truth_z = NULL) {
  if (!is_matrix_like(counts) || !matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must be a non-negative integer-like numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (!is.list(marker_index) || is.null(names(marker_index))) {
    stop("`marker_index` must be a named list.", call. = FALSE)
  }
  candidate_n_outer <- sort(unique(as.integer(candidate_n_outer)))
  if (length(candidate_n_outer) == 0L || anyNA(candidate_n_outer) || any(candidate_n_outer < 1L)) {
    stop("`candidate_n_outer` must contain positive integers.", call. = FALSE)
  }
  if (length(holdout_fraction) != 1L || is.na(holdout_fraction) || holdout_fraction <= 0 || holdout_fraction >= 1) {
    stop("`holdout_fraction` must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  if (!is.null(labels) && length(labels) != ncol(counts)) {
    stop("`labels` must contain one value per cell.", call. = FALSE)
  }
  if (!is.null(groups) && length(groups) != ncol(counts)) {
    stop("`groups` must contain one value per cell.", call. = FALSE)
  }
  if (!is.null(truth_z)) {
    truth_z <- as.matrix(truth_z)
    if (nrow(truth_z) != ncol(counts) || ncol(truth_z) != length(marker_index)) {
      stop("`truth_z` must have one row per cell and one column per state.", call. = FALSE)
    }
    if (is.null(rownames(truth_z))) {
      rownames(truth_z) <- colnames(counts)
    }
    truth_z <- truth_z[colnames(counts), , drop = FALSE]
  }

  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }
  names(library_size) <- colnames(counts)

  split <- make_validation_split(
    cell_names = colnames(counts),
    groups = groups,
    holdout_fraction = holdout_fraction,
    seed = seed
  )
  marker_scores <- score_marker_baseline(
    counts = counts,
    marker_index = marker_index,
    library_size = library_size
  )$scores

  candidate_rows <- list()
  component_rows <- list()
  z_by_candidate <- list()

  for (idx in seq_along(candidate_n_outer)) {
    n_outer <- candidate_n_outer[[idx]]
    fit_nb_args <- subset_nb_selection_args(nb_args, split$train, n_cells = ncol(counts))
    fit_args <- modifyList(
      list(
        counts = counts[, split$train, drop = FALSE],
        marker_index = marker_index,
        library_size = library_size[split$train],
        n_outer = n_outer
      ),
      fit_nb_args
    )
    fit <- do.call(fit_fibrodynmix_nb, fit_args)

    transfer_call <- modifyList(
      list(
        counts = counts[, split$holdout, drop = FALSE],
        fit = fit,
        library_size = library_size[split$holdout],
        z_init = marker_scores_to_simplex(
          marker_scores[split$holdout, , drop = FALSE],
          state_names = colnames(fit$z_hat)
        )
      ),
      transfer_args
    )
    transfer <- do.call(fit_fibrodynmix_transfer, transfer_call)
    z_all <- matrix(
      NA_real_,
      nrow = ncol(counts),
      ncol = ncol(fit$z_hat),
      dimnames = list(colnames(counts), colnames(fit$z_hat))
    )
    z_all[rownames(fit$z_hat), ] <- fit$z_hat
    z_all[rownames(transfer$z_hat), ] <- transfer$z_hat
    z_by_candidate[[as.character(n_outer)]] <- z_all

    stability_delta <- NA_real_
    if (idx > 1L) {
      previous <- z_by_candidate[[as.character(candidate_n_outer[[idx - 1L]])]]
      stability_delta <- mean(abs(z_all - previous), na.rm = TRUE)
    }
    gradient <- marker_gradient_preservation(marker_scores, z_all)
    downstream <- if (!is.null(labels)) {
      evaluate_downstream_classification(
        features = z_all,
        labels = labels,
        groups = groups,
        n_folds = min(3L, length(unique(labels))),
        seed = if (is.null(seed)) NULL else seed + n_outer
      )
    } else {
      empty_classification_metrics("not_requested", nrow(z_all), NA_integer_)
    }
    truth_metrics <- if (!is.null(truth_z)) {
      evaluate_state_weights(truth_z, z_all)
    } else {
      NULL
    }

    candidate_rows[[idx]] <- data.frame(
      n_outer = n_outer,
      heldout_nb_objective = transfer$nb_objective,
      heldout_nb_loglik = transfer$nb_loglik,
      training_best_objective = fit$best_objective,
      training_initial_objective = fit$nb_objective_trace[1L],
      training_final_objective = fit$nb_objective_trace[length(fit$nb_objective_trace)],
      objective_improvement = fit$nb_objective_trace[1L] - fit$best_objective,
      z_stability_delta_vs_previous = stability_delta,
      marker_gradient_mean_spearman = gradient$mean_spearman,
      downstream_status = downstream$status,
      downstream_balanced_accuracy = downstream$balanced_accuracy,
      downstream_macro_f1 = downstream$macro_f1,
      downstream_macro_auroc = downstream$macro_auroc,
      truth_rmse = if (is.null(truth_metrics)) NA_real_ else truth_metrics$rmse,
      truth_dominant_accuracy = if (is.null(truth_metrics)) NA_real_ else truth_metrics$dominant_accuracy,
      transfer_z_convergence_rate = transfer$z_convergence_rate,
      train_n_cells = sum(split$train),
      holdout_n_cells = sum(split$holdout),
      stop_reason = fit$stop_reason,
      stringsAsFactors = FALSE
    )
    component_rows[[idx]] <- data.frame(
      n_outer = n_outer,
      state = names(gradient$per_state_spearman),
      marker_gradient_spearman = as.numeric(gradient$per_state_spearman),
      stringsAsFactors = FALSE
    )
  }

  candidates <- do.call(rbind, candidate_rows)
  components <- do.call(rbind, component_rows)
  candidates <- add_validation_selection_score(candidates, selection_weights)
  candidates <- candidates[order(-candidates$selection_score, candidates$n_outer), , drop = FALSE]
  rownames(candidates) <- NULL

  list(
    selected_n_outer = candidates$n_outer[[1L]],
    candidate_scores = candidates,
    selection_components = components,
    split = data.frame(
      cell_id = colnames(counts),
      split = ifelse(split$holdout, "holdout", "train"),
      group = if (is.null(groups)) NA_character_ else as.character(groups),
      stringsAsFactors = FALSE
    ),
    selection_weights = selection_weights,
    selection_rule = "max weighted normalized held-out NB objective, z stability, marker-gradient preservation, and downstream validation; training objective is reported but not used"
  )
}

subset_nb_selection_args <- function(nb_args, train, n_cells) {
  out <- nb_args
  for (field in c("study_id", "donor_id")) {
    if (!is.null(out[[field]]) && length(out[[field]]) == n_cells) {
      out[[field]] <- out[[field]][train]
    }
  }
  out
}

make_validation_split <- function(cell_names, groups = NULL, holdout_fraction = 0.25, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- preserve_seed()
    on.exit(restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  n_cells <- length(cell_names)
  if (n_cells < 4L) {
    stop("Validation-aware selection requires at least four cells.", call. = FALSE)
  }
  if (!is.null(groups) && length(unique(groups)) >= 2L) {
    group_levels <- unique(as.character(groups))
    n_holdout_groups <- max(1L, round(length(group_levels) * holdout_fraction))
    n_holdout_groups <- min(n_holdout_groups, length(group_levels) - 1L)
    holdout_groups <- sample(group_levels, n_holdout_groups)
    holdout <- as.character(groups) %in% holdout_groups
  } else {
    n_holdout <- max(1L, round(n_cells * holdout_fraction))
    n_holdout <- min(n_holdout, n_cells - 1L)
    holdout_idx <- sample(seq_len(n_cells), n_holdout)
    holdout <- seq_len(n_cells) %in% holdout_idx
  }
  names(holdout) <- cell_names
  list(train = !holdout, holdout = holdout)
}

marker_scores_to_simplex <- function(marker_scores, state_names) {
  marker_scores <- as.matrix(marker_scores)
  marker_scores <- marker_scores[, state_names, drop = FALSE]
  centered <- sweep(marker_scores, 1L, rowMeans(marker_scores), "-")
  z <- softmax_rows(centered)
  colnames(z) <- state_names
  rownames(z) <- rownames(marker_scores)
  z
}

marker_gradient_preservation <- function(marker_scores, z) {
  marker_scores <- as.matrix(marker_scores)
  z <- as.matrix(z)
  shared_states <- intersect(colnames(marker_scores), colnames(z))
  if (length(shared_states) == 0L) {
    return(list(mean_spearman = NA_real_, per_state_spearman = numeric(0)))
  }
  per_state <- vapply(shared_states, function(state) {
    stats::cor(marker_scores[, state], z[, state], method = "spearman", use = "pairwise.complete.obs")
  }, numeric(1))
  list(
    mean_spearman = mean(per_state, na.rm = TRUE),
    per_state_spearman = per_state
  )
}

add_validation_selection_score <- function(candidates, selection_weights) {
  metric_map <- list(
    heldout_nb_objective = list(column = "heldout_nb_objective", direction = "lower"),
    z_stability = list(column = "z_stability_delta_vs_previous", direction = "lower"),
    marker_gradient = list(column = "marker_gradient_mean_spearman", direction = "higher"),
    downstream_balanced_accuracy = list(column = "downstream_balanced_accuracy", direction = "higher")
  )
  selection_weights <- selection_weights[names(selection_weights) %in% names(metric_map)]
  if (length(selection_weights) == 0L || sum(selection_weights) <= 0) {
    stop("`selection_weights` must contain at least one positive supported metric weight.", call. = FALSE)
  }
  selection_weights <- selection_weights / sum(selection_weights)
  candidates$selection_score <- 0
  candidates$selection_score_components <- ""

  component_strings <- vector("list", nrow(candidates))
  for (metric in names(selection_weights)) {
    spec <- metric_map[[metric]]
    values <- candidates[[spec$column]]
    finite_values <- values[is.finite(values)]
    if (length(finite_values) == 0L) {
      scaled <- rep(0.5, length(values))
    } else {
      values[!is.finite(values)] <- stats::median(finite_values)
      scaled <- minmax_metric(values, direction = spec$direction)
    }
    weighted <- selection_weights[[metric]] * scaled
    candidates$selection_score <- candidates$selection_score + weighted
    for (i in seq_along(weighted)) {
      component_strings[[i]] <- c(component_strings[[i]], sprintf("%s=%.4f", metric, weighted[[i]]))
    }
  }
  candidates$selection_score_components <- vapply(component_strings, paste, character(1), collapse = ";")
  candidates
}

minmax_metric <- function(values, direction = c("lower", "higher")) {
  direction <- match.arg(direction)
  values <- as.numeric(values)
  rng <- range(values, na.rm = TRUE)
  if (!all(is.finite(rng)) || abs(diff(rng)) <= .Machine$double.eps) {
    return(rep(0.5, length(values)))
  }
  scaled <- (values - rng[1L]) / diff(rng)
  if (direction == "lower") {
    scaled <- 1 - scaled
  }
  scaled
}
