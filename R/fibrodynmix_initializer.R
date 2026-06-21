#' Fit a count-based FibroDynMix initializer
#'
#' Learns an interpretable cell-state simplex and state-gene program from raw
#' counts using weak marker priors and alternating least-squares updates on
#' log-normalized expression. This is an initializer for later probabilistic
#' FibroDynMix inference, not the final Bayesian model.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of weak prior marker genes for each state.
#'   Entries can be gene names or integer row indices.
#' @param library_size Optional vector of per-cell library sizes. If omitted,
#'   column sums of `counts` are used.
#' @param n_iter Number of alternating update iterations.
#' @param scale_factor Library-size normalization scale factor.
#' @param marker_boost Additive initialization boost for marker genes in their
#'   corresponding state program.
#' @param ridge Ridge penalty used in least-squares program updates.
#' @param temperature Softmax temperature for cell-state updates.
#' @param prior_weight Weight of the weak marker-prior logits retained during
#'   cell-state updates.
#' @param program_weight Weight of learned state-program reconstruction logits
#'   during cell-state updates.
#' @param return_normalized_expression Whether to retain the normalized
#'   gene-by-cell expression matrix in the returned object.
#' @param verbose Whether to print objective values.
#'
#' @return A list with `z_hat`, `beta_hat`, `alpha_hat`, `objective`,
#'   `marker_index`, and `normalized_expression`.
#' @export
fit_fibrodynmix_initializer <- function(counts,
                                        marker_index,
                                        library_size = NULL,
                                        n_iter = 25,
                                        scale_factor = 10000,
                                        marker_boost = 1,
                                        ridge = 0.1,
                                        temperature = 1,
                                        prior_weight = 1,
                                        program_weight = 0.1,
                                        return_normalized_expression = FALSE,
                                        verbose = FALSE) {
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
  assert_positive_integer(n_iter, "n_iter")
  if (temperature <= 0 || length(temperature) != 1L || is.na(temperature)) {
    stop("`temperature` must be a positive numeric scalar.", call. = FALSE)
  }
  if (ridge < 0 || length(ridge) != 1L || is.na(ridge)) {
    stop("`ridge` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (prior_weight < 0 || length(prior_weight) != 1L || is.na(prior_weight)) {
    stop("`prior_weight` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (program_weight < 0 || length(program_weight) != 1L || is.na(program_weight)) {
    stop("`program_weight` must be a non-negative numeric scalar.", call. = FALSE)
  }

  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }

  normalized <- log_normalize_counts(counts, library_size, scale_factor)
  alpha_hat <- matrix_row_means_core(normalized)

  state_names <- names(marker_index)
  marker_rows <- lapply(marker_index, resolve_marker_rows, gene_names = rownames(counts))
  empty_states <- state_names[vapply(marker_rows, length, integer(1)) == 0L]
  if (length(empty_states) > 0L) {
    stop(
      sprintf("States have no markers present in `counts`: %s", paste(empty_states, collapse = ", ")),
      call. = FALSE
    )
  }

  beta_hat <- initialize_beta_from_markers(normalized, alpha_hat, marker_rows, state_names, marker_boost)
  prior_scores <- marker_prior_scores(normalized, marker_rows, state_names)
  z_hat <- update_z_from_beta(
    normalized,
    alpha_hat,
    beta_hat,
    temperature,
    prior_scores,
    prior_weight,
    program_weight
  )
  objective <- numeric(n_iter + 1L)
  objective[1L] <- reconstruction_mse(normalized, alpha_hat, z_hat, beta_hat)

  for (iter in seq_len(n_iter)) {
    beta_hat <- update_beta_from_z(normalized, alpha_hat, z_hat, ridge)
    beta_hat <- orient_beta_with_markers(beta_hat, marker_rows, marker_boost = marker_boost / 10)
    z_hat <- update_z_from_beta(
      normalized,
      alpha_hat,
      beta_hat,
      temperature,
      prior_scores,
      prior_weight,
      program_weight
    )
    objective[iter + 1L] <- reconstruction_mse(normalized, alpha_hat, z_hat, beta_hat)
    if (isTRUE(verbose)) {
      message(sprintf("iter %d objective %.6f", iter, objective[iter + 1L]))
    }
  }

  rownames(z_hat) <- colnames(counts)
  colnames(z_hat) <- state_names
  rownames(beta_hat) <- state_names
  colnames(beta_hat) <- rownames(counts)

  out <- list(
    z_hat = z_hat,
    beta_hat = beta_hat,
    alpha_hat = alpha_hat,
    objective = objective,
    marker_index = marker_index,
    prior_scores = prior_scores
  )
  if (isTRUE(return_normalized_expression)) {
    out$normalized_expression <- normalized
  }
  out
}

initialize_beta_from_markers <- function(normalized, alpha_hat, marker_rows, state_names, marker_boost) {
  n_states <- length(marker_rows)
  n_genes <- nrow(normalized)
  beta_hat <- matrix(0, nrow = n_states, ncol = n_genes)

  gene_sd <- centered_gene_sd(normalized, alpha_hat)
  gene_sd[is.na(gene_sd)] <- 0
  for (k in seq_along(marker_rows)) {
    beta_hat[k, ] <- rnorm(n_genes, 0, 0.01)
    beta_hat[k, marker_rows[[k]]] <- gene_sd[marker_rows[[k]]] + marker_boost
  }

  rownames(beta_hat) <- state_names
  colnames(beta_hat) <- rownames(normalized)
  beta_hat
}

marker_prior_scores <- function(normalized, marker_rows, state_names) {
  scores <- matrix(
    NA_real_,
    nrow = ncol(normalized),
    ncol = length(marker_rows),
    dimnames = list(colnames(normalized), state_names)
  )
  for (k in seq_along(marker_rows)) {
    scores[, k] <- matrix_col_means(normalized[marker_rows[[k]], , drop = FALSE])
  }
  sweep(scores, 1, rowMeans(scores), "-")
}

update_z_from_beta <- function(normalized,
                               alpha_hat,
                               beta_hat,
                               temperature,
                               prior_scores = NULL,
                               prior_weight = 0,
                               program_weight = 1) {
  scores <- as.matrix(Matrix::t(normalized) %*% t(beta_hat))
  center_offset <- as.numeric(beta_hat %*% alpha_hat)
  scores <- sweep(scores, 2L, center_offset, "-")
  beta_norm <- sqrt(rowSums(beta_hat^2))
  beta_norm[beta_norm == 0] <- 1
  scores <- sweep(scores, 2, beta_norm, "/")
  scores <- program_weight * scores
  if (!is.null(prior_scores) && prior_weight > 0) {
    scores <- scores + prior_weight * prior_scores
  }
  scores <- sweep(scores, 1, rowMeans(scores), "-") / temperature
  softmax_rows(scores)
}

update_beta_from_z <- function(normalized, alpha_hat, z_hat, ridge) {
  ztz <- crossprod(z_hat)
  penalty <- diag(ridge, nrow = ncol(z_hat), ncol = ncol(z_hat))
  rhs <- t(z_hat) %*% as.matrix(Matrix::t(normalized))
  rhs <- rhs - outer(colSums(z_hat), alpha_hat)
  solve(ztz + penalty, rhs)
}

orient_beta_with_markers <- function(beta_hat, marker_rows, marker_boost) {
  for (k in seq_along(marker_rows)) {
    beta_hat[k, marker_rows[[k]]] <- beta_hat[k, marker_rows[[k]]] + marker_boost
  }
  beta_hat
}

reconstruction_mse <- function(normalized, alpha_hat, z_hat, beta_hat) {
  centered <- as.matrix(normalized)
  centered <- sweep(centered, 1L, alpha_hat, "-")
  fitted <- z_hat %*% beta_hat
  mean((t(centered) - fitted)^2)
}

centered_gene_sd <- function(normalized, alpha_hat) {
  n_cells <- ncol(normalized)
  if (n_cells <= 1L) {
    return(rep(0, nrow(normalized)))
  }
  row_sum_sq <- matrix_row_sums(normalized * normalized)
  variance <- (row_sum_sq - n_cells * alpha_hat^2) / (n_cells - 1L)
  sqrt(pmax(variance, 0))
}
