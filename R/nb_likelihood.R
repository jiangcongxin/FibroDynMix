#' FibroDynMix negative-binomial log-likelihood
#'
#' Evaluates the raw-count negative-binomial likelihood under the FibroDynMix
#' mean model:
#'
#' `log(mu_ig) = log(library_i) + alpha_g + sum_k z_ik beta_kg + effects_ig`.
#'
#' Optional study and donor effects are added on the linear predictor scale when
#' both effect matrices and matching cell-level identifiers are supplied.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param z Cell-by-state simplex matrix.
#' @param beta State-by-gene expression program matrix.
#' @param alpha Gene-level baseline vector.
#' @param phi Gene-specific negative-binomial size/overdispersion vector.
#' @param library_size Positive vector with one library size per cell.
#' @param study_effect Optional study-by-gene effect matrix.
#' @param donor_effect Optional donor-by-gene effect matrix.
#' @param study_id Optional cell-level study identifiers matching rows of
#'   `study_effect`.
#' @param donor_id Optional cell-level donor identifiers matching rows of
#'   `donor_effect`.
#' @param return_matrix If `TRUE`, return a gene-by-cell log-likelihood matrix;
#'   otherwise return the scalar sum.
#'
#' @return Numeric scalar or gene-by-cell matrix.
#' @export
fibrodynmix_nb_loglik <- function(counts,
                                  z,
                                  beta,
                                  alpha,
                                  phi,
                                  library_size,
                                  study_effect = NULL,
                                  donor_effect = NULL,
                                  study_id = NULL,
                                  donor_id = NULL,
                                  return_matrix = FALSE) {
  parts <- prepare_nb_inputs(
    counts = counts,
    z = z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id
  )

  if (isTRUE(return_matrix)) {
    mu <- fibrodynmix_nb_mu_from_prepared(parts)
    loglik <- stats::dnbinom(
      x = as.vector(as.matrix(parts$counts)),
      size = rep(parts$phi, times = ncol(parts$counts)),
      mu = as.vector(mu),
      log = TRUE
    )
    loglik <- matrix(
      loglik,
      nrow = nrow(parts$counts),
      ncol = ncol(parts$counts),
      dimnames = dimnames(parts$counts)
    )
    return(loglik)
  }
  fibrodynmix_nb_loglik_scalar(parts)
}

#' FibroDynMix negative-binomial deviance
#'
#' Computes twice the difference between the saturated negative-binomial
#' log-likelihood and the fitted FibroDynMix negative-binomial log-likelihood,
#' using the supplied gene-specific `phi`.
#'
#' @inheritParams fibrodynmix_nb_loglik
#'
#' @return Numeric scalar deviance.
#' @export
fibrodynmix_nb_deviance <- function(counts,
                                    z,
                                    beta,
                                    alpha,
                                    phi,
                                    library_size,
                                    study_effect = NULL,
                                    donor_effect = NULL,
                                    study_id = NULL,
                                    donor_id = NULL) {
  fitted_loglik <- fibrodynmix_nb_loglik(
    counts = counts,
    z = z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id,
    return_matrix = FALSE
  )
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  phi <- validate_gene_vector(phi, nrow(counts), rownames(counts), "phi")
  saturated_loglik <- sum(vapply(seq_len(nrow(counts)), function(g) {
    y <- matrix_gene_vector(counts, g)
    mu <- pmax(y, .Machine$double.eps)
    sum(stats::dnbinom(y, size = phi[g], mu = mu, log = TRUE))
  }, numeric(1)))

  2 * (saturated_loglik - fitted_loglik)
}

#' FibroDynMix negative-binomial optimization objective
#'
#' Returns a minimization objective based on the negative raw-count
#' negative-binomial log-likelihood. Optional ridge penalties can be applied to
#' `beta`, `study_effect`, and `donor_effect`.
#'
#' @inheritParams fibrodynmix_nb_loglik
#' @param beta_l2 Non-negative L2 penalty for `beta`.
#' @param marker_target Optional state-by-gene matrix used as a weak marker
#'   orientation target for `beta`.
#' @param marker_l2 Non-negative L2 penalty around `marker_target`.
#' @param effect_l2 Non-negative L2 penalty for study effects and, by default,
#'   donor effects.
#' @param donor_effect_l2 Optional non-negative L2 penalty for donor effects. If
#'   `NULL`, `effect_l2` is used.
#' @param average If `TRUE`, divide the objective by the number of observed
#'   count entries.
#'
#' @return Numeric scalar objective to minimize.
#' @export
fibrodynmix_nb_objective <- function(counts,
                                     z,
                                     beta,
                                     alpha,
                                     phi,
                                     library_size,
                                     study_effect = NULL,
                                     donor_effect = NULL,
                                     study_id = NULL,
                                     donor_id = NULL,
                                     beta_l2 = 0,
                                     marker_target = NULL,
                                     marker_l2 = 0,
                                     effect_l2 = 0,
                                     donor_effect_l2 = NULL,
                                     average = FALSE) {
  if (length(beta_l2) != 1L || is.na(beta_l2) || beta_l2 < 0) {
    stop("`beta_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(effect_l2) != 1L || is.na(effect_l2) || effect_l2 < 0) {
    stop("`effect_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (is.null(donor_effect_l2)) {
    donor_effect_l2 <- effect_l2
  }
  if (length(donor_effect_l2) != 1L || is.na(donor_effect_l2) || donor_effect_l2 < 0) {
    stop("`donor_effect_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(marker_l2) != 1L || is.na(marker_l2) || marker_l2 < 0) {
    stop("`marker_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }

  loglik <- fibrodynmix_nb_loglik(
    counts = counts,
    z = z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id,
    return_matrix = FALSE
  )

  penalty <- beta_l2 * sum(as.matrix(beta)^2)
  if (!is.null(marker_target)) {
    marker_target <- as.matrix(marker_target)
    if (!identical(dim(marker_target), dim(as.matrix(beta)))) {
      stop("`marker_target` must have the same dimensions as `beta`.", call. = FALSE)
    }
    penalty <- penalty + marker_l2 * sum((as.matrix(beta) - marker_target)^2)
  }
  if (!is.null(study_effect)) {
    penalty <- penalty + effect_l2 * sum(as.matrix(study_effect)^2)
  }
  if (!is.null(donor_effect)) {
    penalty <- penalty + donor_effect_l2 * sum(as.matrix(donor_effect)^2)
  }

  objective <- -loglik + penalty
  if (isTRUE(average)) {
    objective <- objective / matrix_n_entries(counts)
  }
  objective
}

prepare_nb_inputs <- function(counts,
                              z,
                              beta,
                              alpha,
                              phi,
                              library_size,
                              study_effect,
                              donor_effect,
                              study_id,
                              donor_id) {
  z <- as.matrix(z)
  beta <- as.matrix(beta)

  if (!is_matrix_like(counts) || !matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must be a non-negative integer-like numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(z) || anyNA(z) || any(z < 0)) {
    stop("`z` must be a non-negative numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be a numeric matrix without NA values.", call. = FALSE)
  }
  if (nrow(z) != ncol(counts)) {
    stop("`z` must have one row per cell/column in `counts`.", call. = FALSE)
  }
  if (ncol(z) != nrow(beta)) {
    stop("`ncol(z)` must equal `nrow(beta)`.", call. = FALSE)
  }
  if (ncol(beta) != nrow(counts)) {
    stop("`beta` must have one column per gene/row in `counts`.", call. = FALSE)
  }

  row_sums <- rowSums(z)
  if (any(abs(row_sums - 1) > 1e-4)) {
    stop("Rows of `z` must sum to one.", call. = FALSE)
  }

  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (is.null(rownames(z))) {
    rownames(z) <- colnames(counts)
  }
  if (is.null(colnames(beta))) {
    colnames(beta) <- rownames(counts)
  }

  if (!identical(rownames(z), colnames(counts))) {
    stop("`rownames(z)` must match `colnames(counts)`.", call. = FALSE)
  }
  if (!identical(colnames(beta), rownames(counts))) {
    stop("`colnames(beta)` must match `rownames(counts)`.", call. = FALSE)
  }

  alpha <- validate_gene_vector(alpha, nrow(counts), rownames(counts), "alpha")
  phi <- validate_gene_vector(phi, nrow(counts), rownames(counts), "phi")
  if (any(phi <= 0)) {
    stop("`phi` must contain positive values.", call. = FALSE)
  }

  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }
  names(library_size) <- colnames(counts)

  study_effect <- validate_effect_matrix(study_effect, rownames(counts), "study_effect")
  donor_effect <- validate_effect_matrix(donor_effect, rownames(counts), "donor_effect")

  if (!is.null(study_effect)) {
    if (is.null(study_id) || length(study_id) != ncol(counts)) {
      stop("`study_id` must be supplied with one value per cell when `study_effect` is used.", call. = FALSE)
    }
    if (any(!study_id %in% rownames(study_effect))) {
      stop("All `study_id` values must appear in `rownames(study_effect)`.", call. = FALSE)
    }
  }
  if (!is.null(donor_effect)) {
    if (is.null(donor_id) || length(donor_id) != ncol(counts)) {
      stop("`donor_id` must be supplied with one value per cell when `donor_effect` is used.", call. = FALSE)
    }
    if (any(!donor_id %in% rownames(donor_effect))) {
      stop("All `donor_id` values must appear in `rownames(donor_effect)`.", call. = FALSE)
    }
  }

  list(
    counts = counts,
    z = z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id
  )
}

fibrodynmix_nb_mu_from_prepared <- function(parts) {
  eta <- parts$z %*% parts$beta
  eta <- sweep(eta, 2, parts$alpha, "+")

  if (!is.null(parts$study_effect)) {
    eta <- eta + parts$study_effect[parts$study_id, , drop = FALSE]
  }
  if (!is.null(parts$donor_effect)) {
    eta <- eta + parts$donor_effect[parts$donor_id, , drop = FALSE]
  }

  log_mu <- sweep(eta, 1, log(parts$library_size), "+")
  log_mu <- pmin(pmax(log_mu, -745), 700)
  mu_cell_gene <- exp(log_mu)
  mu_gene_cell <- t(mu_cell_gene)
  dimnames(mu_gene_cell) <- dimnames(parts$counts)
  pmax(mu_gene_cell, .Machine$double.eps)
}

fibrodynmix_nb_loglik_scalar <- function(parts) {
  eta <- parts$z %*% parts$beta
  eta <- sweep(eta, 2, parts$alpha, "+")
  if (!is.null(parts$study_effect)) {
    eta <- eta + parts$study_effect[parts$study_id, , drop = FALSE]
  }
  if (!is.null(parts$donor_effect)) {
    eta <- eta + parts$donor_effect[parts$donor_id, , drop = FALSE]
  }
  log_mu_cell_gene <- sweep(eta, 1, log(parts$library_size), "+")
  log_mu_cell_gene <- pmin(pmax(log_mu_cell_gene, -745), 700)

  sum(vapply(seq_len(nrow(parts$counts)), function(g) {
    y <- matrix_gene_vector(parts$counts, g)
    mu <- pmax(exp(log_mu_cell_gene[, g]), .Machine$double.eps)
    sum(stats::dnbinom(y, size = parts$phi[g], mu = mu, log = TRUE))
  }, numeric(1)))
}

validate_gene_vector <- function(x, n_genes, gene_names, name) {
  x <- as.numeric(x)
  if (length(x) != n_genes || anyNA(x)) {
    stop(sprintf("`%s` must contain one numeric value per gene.", name), call. = FALSE)
  }
  names(x) <- gene_names
  x
}

validate_effect_matrix <- function(effect, gene_names, name) {
  if (is.null(effect)) {
    return(NULL)
  }
  effect <- as.matrix(effect)
  if (!is.numeric(effect) || anyNA(effect)) {
    stop(sprintf("`%s` must be a numeric matrix without NA values.", name), call. = FALSE)
  }
  if (ncol(effect) != length(gene_names)) {
    stop(sprintf("`%s` must have one column per gene.", name), call. = FALSE)
  }
  if (is.null(rownames(effect))) {
    stop(sprintf("`%s` must have row names.", name), call. = FALSE)
  }
  if (is.null(colnames(effect))) {
    colnames(effect) <- gene_names
  }
  if (!identical(colnames(effect), gene_names)) {
    stop(sprintf("`colnames(%s)` must match `rownames(counts)`.", name), call. = FALSE)
  }
  effect
}
