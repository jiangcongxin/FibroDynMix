#' Score cells with a marker-based baseline
#'
#' Computes a simple gene-set scoring baseline from raw counts. This function is
#' intended as a comparator for simulation benchmarks, not as the FibroDynMix
#' generative model.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of marker genes for each state. Entries can be
#'   gene names or integer row indices.
#' @param marker_weight Optional per-state marker-gene reliability weights. Supported
#'   formats are a state-by-gene matrix, a named list keyed by state, or a
#'   numeric vector.
#' @param library_size Optional vector of per-cell library sizes. If omitted,
#'   column sums of `counts` are used.
#' @param scale_factor Library-size normalization scale factor.
#' @param temperature Softmax temperature applied to centered state scores.
#'
#' @return A list with `z_pred`, `scores`, and `normalized_expression`.
#' @export
score_marker_baseline <- function(counts,
                                  marker_index,
                                  marker_weight = NULL,
                                  library_size = NULL,
                                  scale_factor = 10000,
                                  temperature = 1) {
  if (!is_matrix_like(counts)) {
    stop("`counts` must be numeric.", call. = FALSE)
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
  if (temperature <= 0 || length(temperature) != 1L || is.na(temperature)) {
    stop("`temperature` must be a positive numeric scalar.", call. = FALSE)
  }

  if (is.null(library_size) &&
      is.numeric(marker_weight) &&
      is.vector(marker_weight) &&
      length(marker_weight) == ncol(counts) &&
      length(marker_weight) != nrow(counts) &&
      length(marker_weight) != length(marker_index)) {
    library_size <- marker_weight
    marker_weight <- NULL
  }

  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }

  normalized <- log_normalize_counts(counts, library_size, scale_factor)
  state_names <- names(marker_index)
  scores <- matrix(
    NA_real_,
    nrow = ncol(counts),
    ncol = length(marker_index),
    dimnames = list(colnames(counts), state_names)
  )

  for (k in seq_along(marker_index)) {
    state_name <- state_names[k]
    marker_rows <- resolve_marker_rows(marker_index[[k]], rownames(counts))
    if (length(marker_rows) == 0L) {
      stop(sprintf("State `%s` has no markers present in `counts`.", state_names[k]), call. = FALSE)
    }
    state_weights <- resolve_marker_state_weight(
      marker_weight = marker_weight,
      state_name = state_name,
      marker_rows = marker_rows,
      gene_names = rownames(counts),
      state_names = state_names
    )
    scores[, k] <- weighted_matrix_col_means(
      x = normalized[marker_rows, , drop = FALSE],
      weights = state_weights
    )
  }

  centered_scores <- sweep(scores, 1, rowMeans(scores), "-") / temperature
  z_pred <- softmax_rows(centered_scores)
  rownames(z_pred) <- rownames(scores)
  colnames(z_pred) <- colnames(scores)

  list(
    z_pred = z_pred,
    scores = scores,
    normalized_expression = normalized
  )
}

#' Run a marker-scoring simulation benchmark
#'
#' Simulates one FibroDynMix scenario, applies the marker-based baseline, and
#' evaluates state-weight recovery.
#'
#' @param ... Arguments passed to `simulate_fibrodynmix()`.
#' @param temperature Softmax temperature for `score_marker_baseline()`.
#' @param marker_weight Optional `score_marker_baseline()` marker weights.
#'
#' @return A list with `simulation`, `baseline`, and `metrics`.
#' @export
run_marker_scoring_benchmark <- function(..., temperature = 1, marker_weight = NULL) {
  sim <- simulate_fibrodynmix(...)
  baseline <- score_marker_baseline(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    marker_weight = marker_weight,
    library_size = sim$cell_metadata$library_size,
    temperature = temperature
  )
  metrics <- evaluate_state_weights(sim$z, baseline$z_pred)

  list(
    simulation = sim,
    baseline = baseline,
    metrics = metrics
  )
}

weighted_matrix_col_means <- function(x, weights) {
  x <- as.matrix(x)
  weights <- as.numeric(weights)
  if (length(weights) != nrow(x)) {
    stop("`weights` length must match the number of marker rows in `x`.", call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("`marker_weight` must contain only finite numeric values.", call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("`marker_weight` values must be non-negative.", call. = FALSE)
  }
  if (all(weights == 0)) {
    return(rep(NA_real_, ncol(x)))
  }
  denom <- sum(weights)
  if (denom <= 0 || !is.finite(denom)) {
    stop("`weights` must contain at least one positive finite value.", call. = FALSE)
  }
  as.numeric((weights %*% x) / denom)
}

resolve_marker_state_weight <- function(marker_weight,
                                       state_name,
                                       marker_rows,
                                       gene_names,
                                       state_names) {
  n_markers <- length(marker_rows)
  if (is.null(marker_weight)) {
    return(rep(1, n_markers))
  }

  if (is.list(marker_weight)) {
    if (!state_name %in% names(marker_weight)) {
      return(rep(1, n_markers))
    }
    weights <- as.numeric(marker_weight[[state_name]])
    if (length(weights) == 1L) {
      return(rep(weights, n_markers))
    }
    if (length(weights) != n_markers) {
      stop("`marker_weight` list entries must be scalar or match per-state marker count.", call. = FALSE)
    }
    return(weights)
  }

  if (is.matrix(marker_weight) || is.data.frame(marker_weight)) {
    weights <- as.matrix(marker_weight)
    if (nrow(weights) != length(state_names) || ncol(weights) != length(gene_names)) {
      stop("`marker_weight` matrix must be state-by-gene and match marker names.", call. = FALSE)
    }
    if (!is.null(rownames(weights))) {
      state_match <- match(state_name, rownames(weights))
      if (is.na(state_match)) {
        stop("`marker_weight` matrix rows must include all states.", call. = FALSE)
      }
      weights <- weights[state_match, , drop = FALSE]
    }
    if (!is.null(colnames(weights))) {
      gene_match <- match(gene_names, colnames(weights))
      if (anyNA(gene_match)) {
        stop("`marker_weight` matrix columns must include all genes.", call. = FALSE)
      }
      weights <- weights[, gene_match, drop = FALSE]
    }
    return(as.numeric(weights[1L, marker_rows]))
  }

  if (is.numeric(marker_weight) && is.vector(marker_weight)) {
    weights <- as.numeric(marker_weight)
    if (length(weights) == 1L) {
      return(rep(weights, n_markers))
    }
    if (length(weights) == n_markers) {
      return(weights)
    }
    if (length(weights) == length(gene_names)) {
      if (!is.null(names(weights))) {
        gene_match <- match(gene_names, names(weights))
        if (anyNA(gene_match)) {
          stop("`marker_weight` names must cover all genes.", call. = FALSE)
        }
        weights <- weights[gene_match]
      }
      return(weights[marker_rows])
    }
    if (!is.null(names(weights)) && length(weights) == length(state_names)) {
      state_match <- match(state_name, names(weights))
      if (is.na(state_match)) {
        stop("`marker_weight` names must cover all states.", call. = FALSE)
      }
      return(rep(weights[state_match], n_markers))
    }
  }

  stop("Unsupported `marker_weight` format. Use NULL, numeric vector, matrix, or state-keyed list.", call. = FALSE)
}

resolve_marker_rows <- function(markers, gene_names) {
  if (is.numeric(markers)) {
    markers <- markers[markers >= 1 & markers <= length(gene_names)]
    return(as.integer(markers))
  }

  markers <- as.character(markers)
  which(gene_names %in% markers)
}
