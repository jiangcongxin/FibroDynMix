#' Evaluate inferred fibroblast state weights
#'
#' Computes state-weight recovery metrics against the true latent simplex.
#'
#' @param z_true Numeric matrix of true cell-by-state weights.
#' @param z_pred Numeric matrix of predicted cell-by-state weights.
#'
#' @return A list with scalar metrics, per-state RMSE, and a confusion matrix
#'   for dominant states.
#' @export
evaluate_state_weights <- function(z_true, z_pred) {
  z_true <- as.matrix(z_true)
  z_pred <- as.matrix(z_pred)
  validate_matching_state_matrices(z_true, z_pred)

  state_names <- colnames(z_true)
  error <- z_pred - z_true
  per_state_rmse <- sqrt(colMeans(error^2))
  names(per_state_rmse) <- state_names

  true_state <- state_names[max.col(z_true, ties.method = "first")]
  pred_state <- state_names[max.col(z_pred, ties.method = "first")]
  confusion <- table(
    true = factor(true_state, levels = state_names),
    predicted = factor(pred_state, levels = state_names)
  )

  state_tp <- as.numeric(diag(confusion))
  state_fp <- as.numeric(colSums(confusion)) - state_tp
  state_fn <- as.numeric(rowSums(confusion)) - state_tp
  per_state_precision <- safe_divide_vector(state_tp, state_tp + state_fp)
  per_state_recall <- safe_divide_vector(state_tp, state_tp + state_fn)
  per_state_f1 <- safe_divide_vector(2 * per_state_precision * per_state_recall, per_state_precision + per_state_recall)
  names(state_tp) <- state_names
  names(state_fp) <- state_names
  names(state_fn) <- state_names
  names(per_state_precision) <- state_names
  names(per_state_recall) <- state_names
  names(per_state_f1) <- state_names

  list(
    rmse = sqrt(mean(error^2)),
    mean_absolute_error = mean(abs(error)),
    dominant_accuracy = mean(true_state == pred_state),
    dominant_precision = mean(true_state == pred_state),
    dominant_recall = mean(true_state == pred_state),
    dominant_f1 = mean(true_state == pred_state),
    macro_dominant_precision = mean(per_state_precision, na.rm = TRUE),
    macro_dominant_recall = mean(per_state_recall, na.rm = TRUE),
    macro_dominant_f1 = mean(per_state_f1, na.rm = TRUE),
    mean_entropy_true = mean(simplex_entropy(z_true)),
    mean_entropy_pred = mean(simplex_entropy(z_pred)),
    per_state_rmse = per_state_rmse,
    per_state_precision = per_state_precision,
    per_state_recall = per_state_recall,
    per_state_f1 = per_state_f1,
    dominant_confusion = confusion,
    per_state_true_positives = state_tp,
    per_state_false_positives = state_fp,
    per_state_false_negatives = state_fn
  )
}

#' Evaluate dominant-state classification among confident cells
#'
#' Computes coverage and dominant accuracy only on cells whose maximum predicted
#' state probability passes a confidence threshold.
#'
#' @param z_true Numeric matrix of true cell-by-state weights.
#' @param z_pred Numeric matrix of predicted cell-by-state weights.
#' @param confidence_threshold Confidence cutoff on predicted dominant state weight.
#'
#' @return A list with coverage and confident-cell dominant accuracy.
#' @export
evaluate_dominant_state_coverage <- function(z_true, z_pred, confidence_threshold = 0.8) {
  if (length(confidence_threshold) != 1L || is.na(confidence_threshold) || confidence_threshold < 0 || confidence_threshold > 1) {
    stop("`confidence_threshold` must be in [0, 1].", call. = FALSE)
  }
  z_true <- as.matrix(z_true)
  z_pred <- as.matrix(z_pred)
  validate_matching_state_matrices(z_true, z_pred)

  state_names <- colnames(z_true)
  true_state <- state_names[max.col(z_true, ties.method = "first")]
  pred_state <- state_names[max.col(z_pred, ties.method = "first")]
  conf <- apply(z_pred, 1L, max)
  kept <- conf >= confidence_threshold
  n_kept <- sum(kept)
  list(
    confidence_threshold = confidence_threshold,
    coverage = n_kept / nrow(z_true),
    confident_cell_count = n_kept,
    confident_dominant_accuracy = if (n_kept == 0L) NA_real_ else mean(true_state[kept] == pred_state[kept]),
    confident_dominant_mistake_rate = if (n_kept == 0L) NA_real_ else 1 - mean(true_state[kept] == pred_state[kept]),
    confident_confusion = if (n_kept == 0L) NULL else table(
      true = factor(true_state[kept], levels = state_names),
      predicted = factor(pred_state[kept], levels = state_names)
    )
  )
}

#' Evaluate recovery of state-specific marker programs
#'
#' Computes state-wise and macro-average AUPRC for predicted gene-state marker
#' scores against true marker assignments.
#'
#' @param marker_truth Logical or 0/1 matrix with states in rows and genes in
#'   columns. Non-zero entries indicate true state markers.
#' @param marker_scores Numeric matrix with the same dimensions as
#'   `marker_truth`; larger values indicate stronger predicted marker support.
#'
#' @return A list with macro AUPRC and per-state AUPRC.
#' @export
evaluate_marker_recovery <- function(marker_truth, marker_scores) {
  marker_truth <- as.matrix(marker_truth)
  marker_scores <- as.matrix(marker_scores)
  if (!identical(dim(marker_truth), dim(marker_scores))) {
    stop("`marker_truth` and `marker_scores` must have identical dimensions.", call. = FALSE)
  }
  if (is.null(rownames(marker_truth))) {
    rownames(marker_truth) <- sprintf("state_%d", seq_len(nrow(marker_truth)))
  }
  if (is.null(rownames(marker_scores))) {
    rownames(marker_scores) <- rownames(marker_truth)
  }
  if (!identical(rownames(marker_truth), rownames(marker_scores))) {
    stop("`marker_truth` and `marker_scores` must have identical row names.", call. = FALSE)
  }

  truth <- marker_truth != 0
  per_state_auprc <- vapply(
    seq_len(nrow(truth)),
    function(k) binary_auprc(truth[k, ], marker_scores[k, ]),
    numeric(1)
  )
  names(per_state_auprc) <- rownames(truth)

  list(
    macro_auprc = mean(per_state_auprc, na.rm = TRUE),
    per_state_auprc = per_state_auprc
  )
}

#' Evaluate rare transition detection
#'
#' Computes threshold-free and thresholded detection metrics for rare
#' transition-like cells.
#'
#' @param transition_truth Logical vector indicating true transition cells.
#' @param transition_score Numeric vector where larger values indicate higher
#'   transition evidence.
#' @param threshold Optional score threshold. If omitted, the top `sum(truth)`
#'   cells are called positive.
#'
#' @return A list with AUROC, AUPRC, threshold, precision, recall, and F1.
#' @export
evaluate_transition_detection <- function(transition_truth,
                                          transition_score,
                                          threshold = NULL) {
  truth <- as.logical(transition_truth)
  score <- as.numeric(transition_score)
  if (length(truth) != length(score)) {
    stop("`transition_truth` and `transition_score` must have the same length.", call. = FALSE)
  }
  if (anyNA(truth) || anyNA(score)) {
    stop("`transition_truth` and `transition_score` cannot contain NA values.", call. = FALSE)
  }
  if (!any(truth) || all(truth)) {
    stop("`transition_truth` must contain both positive and negative cells.", call. = FALSE)
  }

  if (is.null(threshold)) {
    positive_count <- sum(truth)
    threshold <- sort(score, decreasing = TRUE)[positive_count]
  }

  predicted <- score >= threshold
  tp <- sum(predicted & truth)
  fp <- sum(predicted & !truth)
  fn <- sum(!predicted & truth)
  precision <- safe_divide(tp, tp + fp)
  recall <- safe_divide(tp, tp + fn)
  f1 <- safe_divide(2 * precision * recall, precision + recall)

  list(
    auroc = binary_auroc(truth, score),
    auprc = binary_auprc(truth, score),
    threshold = threshold,
    precision = precision,
    recall = recall,
    f1 = f1
  )
}

#' Evaluate downstream classification from inferred state features
#'
#' Runs a lightweight nearest-centroid classifier with stratified
#' cross-validation. The classifier is intentionally simple so benchmark
#' differences primarily reflect the supplied feature representation rather than
#' a heavily tuned downstream model.
#'
#' @param features Numeric observation-by-feature matrix.
#' @param labels Class labels, one per observation.
#' @param groups Optional grouping variable. When supplied, cross-validation
#'   folds are assigned at the group level to avoid train/test leakage.
#' @param n_folds Requested number of cross-validation folds.
#' @param seed Optional seed for fold assignment.
#'
#' @return A list with cross-validated accuracy, balanced accuracy, macro-F1,
#'   macro one-vs-rest AUROC, and fold metadata.
#' @export
evaluate_downstream_classification <- function(features,
                                               labels,
                                               groups = NULL,
                                               n_folds = 5,
                                               seed = NULL) {
  features <- as.matrix(features)
  if (!is.numeric(features) || length(dim(features)) != 2L) {
    stop("`features` must be a numeric matrix.", call. = FALSE)
  }
  labels <- as.character(labels)
  if (length(labels) != nrow(features)) {
    stop("`labels` must have one value per feature row.", call. = FALSE)
  }
  if (anyNA(features) || any(!is.finite(features))) {
    stop("`features` must contain only finite non-missing values.", call. = FALSE)
  }
  if (anyNA(labels) || any(labels == "")) {
    stop("`labels` cannot contain missing or empty values.", call. = FALSE)
  }
  assert_positive_integer(n_folds, "n_folds")

  class_levels <- sort(unique(labels))
  if (length(class_levels) < 2L) {
    return(empty_classification_metrics(
      status = "single_class",
      n_observations = nrow(features),
      n_classes = length(class_levels)
    ))
  }

  if (!is.null(groups)) {
    groups <- as.character(groups)
    if (length(groups) != nrow(features)) {
      stop("`groups` must have one value per feature row.", call. = FALSE)
    }
    if (anyNA(groups) || any(groups == "")) {
      stop("`groups` cannot contain missing or empty values.", call. = FALSE)
    }
    fold_id <- make_group_stratified_folds(labels, groups, n_folds, seed)
  } else {
    fold_id <- make_stratified_folds(labels, n_folds, seed)
  }

  if (length(unique(fold_id[!is.na(fold_id)])) < 2L) {
    return(empty_classification_metrics(
      status = "insufficient_class_replicates",
      n_observations = nrow(features),
      n_classes = length(class_levels)
    ))
  }

  predicted <- rep(NA_character_, length(labels))
  class_scores <- matrix(
    NA_real_,
    nrow = nrow(features),
    ncol = length(class_levels),
    dimnames = list(rownames(features), class_levels)
  )
  for (fold in sort(unique(fold_id))) {
    test <- fold_id == fold
    train <- !test
    train_classes <- sort(unique(labels[train]))
    if (length(train_classes) < 2L || any(!class_levels %in% train_classes)) {
      next
    }
    fit <- fit_nearest_centroid(features[train, , drop = FALSE], labels[train], class_levels)
    pred <- predict_nearest_centroid(fit, features[test, , drop = FALSE])
    predicted[test] <- pred$predicted
    class_scores[test, ] <- pred$scores
  }

  keep <- !is.na(predicted)
  if (!any(keep)) {
    return(empty_classification_metrics(
      status = "no_valid_folds",
      n_observations = nrow(features),
      n_classes = length(class_levels)
    ))
  }
  metrics <- classification_metrics(labels[keep], predicted[keep], class_levels)
  metrics$status <- if (all(keep)) "ok" else "partial_folds"
  metrics$n_observations <- nrow(features)
  metrics$n_evaluated <- sum(keep)
  metrics$n_classes <- length(class_levels)
  metrics$n_folds <- length(unique(fold_id))
  metrics$class_levels <- paste(class_levels, collapse = ";")
  metrics$macro_auroc <- macro_one_vs_rest_auroc(labels[keep], class_scores[keep, , drop = FALSE], class_levels)
  metrics
}

validate_matching_state_matrices <- function(z_true, z_pred) {
  if (!identical(dim(z_true), dim(z_pred))) {
    stop("`z_true` and `z_pred` must have identical dimensions.", call. = FALSE)
  }
  if (is.null(colnames(z_true))) {
    colnames(z_true) <- sprintf("state_%d", seq_len(ncol(z_true)))
  }
  if (is.null(colnames(z_pred))) {
    colnames(z_pred) <- colnames(z_true)
  }
  if (!identical(colnames(z_true), colnames(z_pred))) {
    stop("`z_true` and `z_pred` must have identical column names.", call. = FALSE)
  }
  invisible(TRUE)
}

simplex_entropy <- function(z) {
  z_safe <- pmax(z, .Machine$double.eps)
  -rowSums(z_safe * log(z_safe))
}

binary_auprc <- function(truth, score) {
  truth <- as.logical(truth)
  score <- as.numeric(score)
  if (!any(truth)) {
    return(NA_real_)
  }
  ord <- order(score, decreasing = TRUE)
  truth <- truth[ord]
  tp <- cumsum(truth)
  fp <- cumsum(!truth)
  precision <- tp / (tp + fp)
  recall <- tp / sum(truth)
  recall_prev <- c(0, recall[-length(recall)])
  sum((recall - recall_prev) * precision)
}

binary_auroc <- function(truth, score) {
  truth <- as.logical(truth)
  score <- as.numeric(score)
  positives <- sum(truth)
  negatives <- sum(!truth)
  if (positives == 0 || negatives == 0) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - positives * (positives + 1) / 2) / (positives * negatives)
}

safe_divide <- function(numerator, denominator) {
  if (denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

safe_divide_vector <- function(numerator, denominator) {
  out <- rep(NA_real_, length(numerator))
  keep <- !is.na(denominator) & denominator != 0
  out[keep] <- numerator[keep] / denominator[keep]
  out
}

make_stratified_folds <- function(labels, n_folds, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- preserve_seed()
    on.exit(restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  labels <- as.character(labels)
  class_counts <- table(labels)
  effective_folds <- min(as.integer(n_folds), min(class_counts))
  if (effective_folds < 2L) {
    return(rep(NA_integer_, length(labels)))
  }
  fold_id <- rep(NA_integer_, length(labels))
  for (class in names(class_counts)) {
    idx <- which(labels == class)
    idx <- sample(idx, length(idx))
    fold_id[idx] <- rep(seq_len(effective_folds), length.out = length(idx))
  }
  fold_id
}

make_group_stratified_folds <- function(labels, groups, n_folds, seed = NULL) {
  group_labels <- tapply(labels, groups, function(x) {
    tab <- sort(table(x), decreasing = TRUE)
    names(tab)[1L]
  })
  group_fold <- make_stratified_folds(as.character(group_labels), n_folds, seed)
  names(group_fold) <- names(group_labels)
  unname(group_fold[groups])
}

fit_nearest_centroid <- function(features, labels, class_levels) {
  centers <- matrix(
    NA_real_,
    nrow = length(class_levels),
    ncol = ncol(features),
    dimnames = list(class_levels, colnames(features))
  )
  for (class in class_levels) {
    centers[class, ] <- colMeans(features[labels == class, , drop = FALSE])
  }
  list(centers = centers, class_levels = class_levels)
}

predict_nearest_centroid <- function(fit, features) {
  distances <- matrix(
    NA_real_,
    nrow = nrow(features),
    ncol = length(fit$class_levels),
    dimnames = list(rownames(features), fit$class_levels)
  )
  for (class in fit$class_levels) {
    centered <- sweep(features, 2, fit$centers[class, ], "-")
    distances[, class] <- sqrt(rowSums(centered^2))
  }
  score_scale <- stats::median(distances[is.finite(distances)])
  if (is.na(score_scale) || score_scale <= 0) {
    score_scale <- 1
  }
  scores <- -distances / score_scale
  predicted <- fit$class_levels[max.col(scores, ties.method = "first")]
  list(predicted = predicted, scores = scores)
}

classification_metrics <- function(labels, predicted, class_levels) {
  confusion <- table(
    truth = factor(labels, levels = class_levels),
    predicted = factor(predicted, levels = class_levels)
  )
  tp <- as.numeric(diag(confusion))
  fp <- as.numeric(colSums(confusion)) - tp
  fn <- as.numeric(rowSums(confusion)) - tp
  precision <- safe_divide_vector(tp, tp + fp)
  recall <- safe_divide_vector(tp, tp + fn)
  f1 <- safe_divide_vector(2 * precision * recall, precision + recall)
  list(
    accuracy = mean(labels == predicted),
    balanced_accuracy = mean(recall, na.rm = TRUE),
    macro_precision = mean(precision, na.rm = TRUE),
    macro_recall = mean(recall, na.rm = TRUE),
    macro_f1 = mean(f1, na.rm = TRUE),
    confusion = confusion
  )
}

macro_one_vs_rest_auroc <- function(labels, scores, class_levels) {
  aurocs <- vapply(class_levels, function(class) {
    truth <- labels == class
    if (!any(truth) || all(truth)) {
      return(NA_real_)
    }
    binary_auroc(truth, scores[, class])
  }, numeric(1))
  mean(aurocs, na.rm = TRUE)
}

empty_classification_metrics <- function(status, n_observations, n_classes) {
  list(
    status = status,
    accuracy = NA_real_,
    balanced_accuracy = NA_real_,
    macro_precision = NA_real_,
    macro_recall = NA_real_,
    macro_f1 = NA_real_,
    macro_auroc = NA_real_,
    n_observations = n_observations,
    n_evaluated = 0L,
    n_classes = n_classes,
    n_folds = 0L,
    class_levels = NA_character_,
    confusion = NULL
  )
}

preserve_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

restore_seed <- function(seed) {
  if (is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", seed, envir = .GlobalEnv)
  }
}
