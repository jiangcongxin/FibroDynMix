#' Prepare a fixed logistic-normal prior for FibroDynMix state weights
#'
#' A K-state simplex is represented by K - 1 reference logits. This function
#' validates a Gaussian prior on those logits and aligns its cells and state
#' coordinates by name. The resulting object can be supplied as `z_prior` to
#' [fit_fibrodynmix_nb()] or [fibrodynmix_nb_objective()].
#'
#' @param eta_mean A length-`n_logits` mean vector or an `n_cells` by
#'   `n_logits` matrix of per-cell means.
#' @param n_cells Number of cells represented by the prior.
#' @param n_logits Number of reference logits, equal to the number of states
#'   minus one.
#' @param covariance Positive-definite covariance matrix for reference logits.
#'   Supply exactly one of `covariance` or `precision`.
#' @param precision Positive-definite precision matrix for reference logits.
#' @param cell_names Optional unique cell names, in the order used by the fit.
#' @param logit_names Optional names of non-reference states, in the order used
#'   by the fit.
#' @param covariance_ridge Non-negative diagonal ridge added to a supplied
#'   covariance before inversion.
#'
#' @return An object of class `fibrodynmix_logistic_normal_prior`.
#' @export
prepare_logistic_normal_prior <- function(eta_mean,
                                           n_cells,
                                           n_logits,
                                           covariance = NULL,
                                           precision = NULL,
                                           cell_names = NULL,
                                           logit_names = NULL,
                                           covariance_ridge = 0) {
  lnp_assert_positive_integer(n_cells, "n_cells")
  lnp_assert_positive_integer(n_logits, "n_logits")
  cell_names <- lnp_normalize_names(cell_names, n_cells, "cell", "cell_names")
  logit_names <- lnp_normalize_names(logit_names, n_logits, "eta", "logit_names")
  covariance_ridge <- lnp_assert_nonnegative_scalar(
    covariance_ridge,
    "covariance_ridge"
  )

  normalized_mean <- normalize_logistic_normal_prior_mean(
    eta_mean = eta_mean,
    n_cells = n_cells,
    n_logits = n_logits,
    cell_names = cell_names,
    logit_names = logit_names
  )
  normalized_scale <- normalize_logistic_normal_prior_scale(
    covariance = covariance,
    precision = precision,
    n_logits = n_logits,
    logit_names = logit_names,
    covariance_ridge = covariance_ridge
  )

  structure(
    c(
      list(
        mean = normalized_mean,
        n_cells = n_cells,
        n_logits = n_logits,
        cell_names = cell_names,
        logit_names = logit_names
      ),
      normalized_scale
    ),
    class = c("fibrodynmix_logistic_normal_prior", "list")
  )
}

normalize_logistic_normal_prior_mean <- function(eta_mean,
                                                  n_cells,
                                                  n_logits,
                                                  cell_names = NULL,
                                                  logit_names = NULL) {
  lnp_assert_positive_integer(n_cells, "n_cells")
  lnp_assert_positive_integer(n_logits, "n_logits")
  cell_names <- lnp_normalize_names(cell_names, n_cells, "cell", "cell_names")
  logit_names <- lnp_normalize_names(logit_names, n_logits, "eta", "logit_names")

  if (is.data.frame(eta_mean)) {
    eta_mean <- as.matrix(eta_mean)
  }

  if (is.matrix(eta_mean)) {
    if (!is.numeric(eta_mean) ||
      nrow(eta_mean) != n_cells ||
      ncol(eta_mean) != n_logits) {
      stop(
        "`eta_mean` must be a numeric matrix with `n_cells` rows and `n_logits` columns.",
        call. = FALSE
      )
    }
    eta_mean <- lnp_reorder_rectangular_matrix(
      x = eta_mean,
      row_names = cell_names,
      col_names = logit_names,
      argument = "eta_mean"
    )
  } else {
    if (!is.numeric(eta_mean) || length(eta_mean) != n_logits) {
      stop(
        "`eta_mean` must be a length-`n_logits` numeric vector or an `n_cells` by `n_logits` matrix.",
        call. = FALSE
      )
    }
    if (!is.null(names(eta_mean))) {
      eta_mean <- lnp_reorder_named_vector(
        x = eta_mean,
        target_names = logit_names,
        argument = "eta_mean"
      )
    }
    eta_mean <- matrix(
      rep(as.numeric(eta_mean), each = n_cells),
      nrow = n_cells,
      ncol = n_logits,
      dimnames = list(cell_names, logit_names)
    )
  }

  if (anyNA(eta_mean) || any(!is.finite(eta_mean))) {
    stop("`eta_mean` must contain only finite values.", call. = FALSE)
  }
  storage.mode(eta_mean) <- "double"
  eta_mean
}

normalize_logistic_normal_prior_scale <- function(covariance = NULL,
                                                   precision = NULL,
                                                   n_logits,
                                                   logit_names = NULL,
                                                   covariance_ridge = 0) {
  lnp_assert_positive_integer(n_logits, "n_logits")
  logit_names <- lnp_normalize_names(logit_names, n_logits, "eta", "logit_names")
  covariance_ridge <- lnp_assert_nonnegative_scalar(
    covariance_ridge,
    "covariance_ridge"
  )

  has_covariance <- !is.null(covariance)
  has_precision <- !is.null(precision)
  if (identical(has_covariance, has_precision)) {
    stop("Supply exactly one of `covariance` or `precision`.", call. = FALSE)
  }
  if (has_precision && covariance_ridge > 0) {
    stop(
      "`covariance_ridge` can only be used when `covariance` is supplied.",
      call. = FALSE
    )
  }

  if (has_covariance) {
    covariance <- lnp_validate_pd_matrix(
      x = covariance,
      n_logits = n_logits,
      logit_names = logit_names,
      argument = "covariance",
      check_pd = FALSE
    )
    if (covariance_ridge > 0) {
      covariance <- covariance + diag(covariance_ridge, n_logits)
      dimnames(covariance) <- list(logit_names, logit_names)
    }
    covariance <- lnp_validate_pd_matrix(
      x = covariance,
      n_logits = n_logits,
      logit_names = logit_names,
      argument = "covariance",
      check_pd = TRUE
    )
    covariance_chol <- chol(covariance)
    precision <- chol2inv(covariance_chol)
    dimnames(precision) <- list(logit_names, logit_names)
    matrix_source <- "covariance"
  } else {
    precision <- lnp_validate_pd_matrix(
      x = precision,
      n_logits = n_logits,
      logit_names = logit_names,
      argument = "precision",
      check_pd = TRUE
    )
    precision_chol <- chol(precision)
    covariance <- chol2inv(precision_chol)
    dimnames(covariance) <- list(logit_names, logit_names)
    matrix_source <- "precision"
  }

  covariance_chol <- chol(covariance)
  logdet_covariance <- 2 * sum(log(diag(covariance_chol)))
  list(
    covariance = covariance,
    precision = precision,
    logdet_covariance = logdet_covariance,
    logdet_precision = -logdet_covariance,
    covariance_ridge = covariance_ridge,
    matrix_source = matrix_source
  )
}

validate_logistic_normal_prior <- function(prior,
                                            n_cells = NULL,
                                            n_logits = NULL,
                                            cell_names = NULL,
                                            logit_names = NULL) {
  if (!inherits(prior, "fibrodynmix_logistic_normal_prior")) {
    stop(
      "`prior` must be created by `prepare_logistic_normal_prior()` or `estimate_group_logistic_normal_prior()`.",
      call. = FALSE
    )
  }
  required <- c("mean", "covariance", "n_cells", "n_logits", "cell_names", "logit_names")
  missing <- setdiff(required, names(prior))
  if (length(missing) > 0L) {
    stop(
      sprintf("`prior` is missing required fields: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.null(n_cells) && (!is.numeric(n_cells) || length(n_cells) != 1L || n_cells != prior$n_cells)) {
    stop("`prior` has an incompatible number of cells.", call. = FALSE)
  }
  if (!is.null(n_logits) && (!is.numeric(n_logits) || length(n_logits) != 1L || n_logits != prior$n_logits)) {
    stop("`prior` has an incompatible number of reference logits.", call. = FALSE)
  }
  if (!is.null(cell_names)) {
    cell_names <- lnp_normalize_names(cell_names, prior$n_cells, "cell", "cell_names")
    if (!identical(cell_names, prior$cell_names)) {
      stop("`prior` has incompatible cell names.", call. = FALSE)
    }
  }
  if (!is.null(logit_names)) {
    logit_names <- lnp_normalize_names(logit_names, prior$n_logits, "eta", "logit_names")
    if (!identical(logit_names, prior$logit_names)) {
      stop("`prior` has incompatible reference-logit names.", call. = FALSE)
    }
  }

  normalized <- prepare_logistic_normal_prior(
    eta_mean = prior$mean,
    n_cells = prior$n_cells,
    n_logits = prior$n_logits,
    covariance = prior$covariance,
    cell_names = prior$cell_names,
    logit_names = prior$logit_names
  )
  extra_fields <- setdiff(names(prior), names(normalized))
  for (field in extra_fields) {
    normalized[[field]] <- prior[[field]]
  }
  normalized
}

# Aligns and validates a prepared logistic-normal prior against a simplex
# matrix.  The final state is the reference state, matching
# `simplex_to_logits()` and the NB optimizer.
validate_logistic_normal_prior_for_simplex <- function(prior, z) {
  if (is.null(prior)) {
    return(NULL)
  }
  z <- lnp_normalize_simplex_matrix(z)
  validate_logistic_normal_prior(
    prior = prior,
    n_cells = nrow(z),
    n_logits = ncol(z) - 1L,
    cell_names = rownames(z),
    logit_names = colnames(z)[seq_len(ncol(z) - 1L)]
  )
}

# Builds a cell-specific marker-logit anchor from an existing simplex
# coordinate.  This is a proximal regularizer, not an estimate of a fully
# hierarchical posterior: its role is to keep NB refinement in the
# user-specified marker coordinate while beta and nuisance terms are updated.
prepare_initializer_logit_anchor <- function(z,
                                             logit_sd = 0.1,
                                             reference_state = ncol(z)) {
  z <- lnp_normalize_simplex_matrix(z)
  if (length(logit_sd) != 1L || !is.numeric(logit_sd) || is.na(logit_sd) ||
    !is.finite(logit_sd) || logit_sd <= 0) {
    stop("`logit_sd` must be a positive finite numeric scalar.", call. = FALSE)
  }
  reference_index <- lnp_resolve_reference_state(reference_state, colnames(z))
  eta_mean <- simplex_matrix_to_reference_logits(
    z = z,
    reference_state = reference_index
  )
  prior <- prepare_logistic_normal_prior(
    eta_mean = eta_mean,
    n_cells = nrow(z),
    n_logits = ncol(z) - 1L,
    covariance = diag(logit_sd^2, ncol(z) - 1L),
    cell_names = rownames(z),
    logit_names = colnames(eta_mean)
  )
  prior$prior_origin <- "cell_specific_marker_initializer_logit_anchor"
  prior$reference_state <- colnames(z)[reference_index]
  prior$reference_index <- reference_index
  prior$anchor_logit_sd <- logit_sd
  prior
}

logistic_normal_prior_quadratic <- function(eta,
                                             prior,
                                             cell_index = NULL,
                                             average = FALSE) {
  lnp_assert_prepared_logistic_normal_prior(prior)
  eta_parts <- lnp_normalize_eta_for_prior(
    eta = eta,
    prior = prior,
    cell_index = cell_index
  )
  centered <- eta_parts$eta - eta_parts$mean
  quadratic <- 0.5 * sum((centered %*% prior$precision) * centered)
  if (isTRUE(average)) {
    quadratic <- quadratic / nrow(centered)
  }
  as.numeric(quadratic)
}

#' Estimate an empirical group logistic-normal prior from state weights
#'
#' Estimates group-specific reference-logit means and a pooled within-group
#' covariance from a simplex matrix. This is an empirical-Bayes anchor: it is
#' intended to be estimated from an external or marker-guided initialization,
#' then held fixed during negative-binomial optimization. It is not a full
#' posterior hierarchical model.
#'
#' @param z Cell-by-state simplex matrix, usually a marker-guided
#'   initialization.
#' @param group_id One group identifier per cell, such as donor or study.
#' @param covariance_ridge Non-negative diagonal ridge used to stabilize the
#'   pooled covariance estimate.
#' @param reference_state State name or index used as the reference logit.
#' @param epsilon Positive lower bound used before log-ratio transformation.
#'
#' @return A `fibrodynmix_logistic_normal_prior` object with empirical group
#'   means and pooled covariance metadata.
#' @export
estimate_group_logistic_normal_prior <- function(z,
                                                  group_id,
                                                  covariance_ridge = 1e-4,
                                                  reference_state = ncol(z),
                                                  epsilon = .Machine$double.eps) {
  z <- lnp_normalize_simplex_matrix(z)
  n_cells <- nrow(z)
  n_states <- ncol(z)
  cell_names <- rownames(z)
  state_names <- colnames(z)
  reference_index <- lnp_resolve_reference_state(reference_state, state_names)
  logit_names <- state_names[-reference_index]
  group_id <- lnp_normalize_group_id(group_id, cell_names)
  covariance_ridge <- lnp_assert_nonnegative_scalar(
    covariance_ridge,
    "covariance_ridge"
  )
  if (length(epsilon) != 1L || !is.numeric(epsilon) || is.na(epsilon) ||
    !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("`epsilon` must be a finite scalar in (0, 1).", call. = FALSE)
  }

  eta <- simplex_matrix_to_reference_logits(
    z = z,
    reference_state = reference_index,
    epsilon = epsilon
  )
  group_levels <- sort(unique(group_id))
  group_index <- match(group_id, group_levels)
  group_means <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = ncol(eta),
    dimnames = list(group_levels, colnames(eta))
  )
  for (g in seq_along(group_levels)) {
    group_means[g, ] <- colMeans(eta[group_index == g, , drop = FALSE])
  }
  eta_mean <- group_means[group_index, , drop = FALSE]
  rownames(eta_mean) <- cell_names
  residuals <- eta - eta_mean
  residual_df <- n_cells - length(group_levels)
  denominator <- max(residual_df, 1L)
  empirical_covariance <- crossprod(residuals) / denominator
  dimnames(empirical_covariance) <- list(logit_names, logit_names)

  prior <- prepare_logistic_normal_prior(
    eta_mean = eta_mean,
    n_cells = n_cells,
    n_logits = n_states - 1L,
    covariance = empirical_covariance,
    cell_names = cell_names,
    logit_names = logit_names,
    covariance_ridge = covariance_ridge
  )
  prior$group_id <- group_id
  prior$group_levels <- group_levels
  prior$group_means <- group_means
  prior$group_sizes <- as.integer(tabulate(group_index, nbins = length(group_levels)))
  names(prior$group_sizes) <- group_levels
  prior$empirical_covariance <- empirical_covariance
  prior$residual_df <- residual_df
  prior$reference_state <- state_names[reference_index]
  prior$reference_index <- reference_index
  prior$state_names <- state_names
  prior$epsilon <- epsilon
  prior
}

simplex_matrix_to_reference_logits <- function(z,
                                               reference_state = ncol(z),
                                               epsilon = .Machine$double.eps) {
  z <- lnp_normalize_simplex_matrix(z)
  state_names <- colnames(z)
  reference_index <- lnp_resolve_reference_state(reference_state, state_names)
  if (length(epsilon) != 1L || !is.numeric(epsilon) || is.na(epsilon) ||
    !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("`epsilon` must be a finite scalar in (0, 1).", call. = FALSE)
  }
  z <- pmax(z, epsilon)
  reference <- z[, reference_index]
  eta <- log(sweep(z[, -reference_index, drop = FALSE], 1L, reference, "/"))
  rownames(eta) <- rownames(z)
  colnames(eta) <- state_names[-reference_index]
  eta
}

lnp_normalize_eta_for_prior <- function(eta, prior, cell_index = NULL) {
  if (is.data.frame(eta)) {
    eta <- as.matrix(eta)
  }
  if (is.matrix(eta)) {
    if (!is.numeric(eta) || ncol(eta) != prior$n_logits) {
      stop("`eta` must have one column per reference logit.", call. = FALSE)
    }
    eta <- lnp_reorder_eta_columns(eta, prior$logit_names)
    if (is.null(cell_index)) {
      if (nrow(eta) != prior$n_cells) {
        stop("A matrix `eta` must have one row per prior cell when `cell_index` is omitted.", call. = FALSE)
      }
      eta <- lnp_reorder_eta_rows(eta, prior$cell_names)
      mean <- prior$mean
    } else {
      indices <- lnp_resolve_cell_index(cell_index, prior)
      if (nrow(eta) != length(indices)) {
        stop("`eta` rows must match the length of `cell_index`.", call. = FALSE)
      }
      expected_cells <- prior$cell_names[indices]
      if (!is.null(rownames(eta))) {
        eta <- lnp_reorder_eta_rows(eta, expected_cells)
      } else {
        rownames(eta) <- expected_cells
      }
      mean <- prior$mean[indices, , drop = FALSE]
    }
  } else {
    if (!is.numeric(eta) || length(eta) != prior$n_logits) {
      stop("A vector `eta` must have one value per reference logit.", call. = FALSE)
    }
    if (is.null(cell_index) && prior$n_cells != 1L) {
      stop("`cell_index` is required for a vector `eta` with a per-cell prior.", call. = FALSE)
    }
    indices <- if (is.null(cell_index)) 1L else lnp_resolve_cell_index(cell_index, prior)
    if (length(indices) != 1L) {
      stop("A vector `eta` requires a single `cell_index`.", call. = FALSE)
    }
    if (!is.null(names(eta))) {
      eta <- lnp_reorder_named_vector(eta, prior$logit_names, "eta")
    }
    eta <- matrix(
      as.numeric(eta),
      nrow = 1L,
      dimnames = list(prior$cell_names[indices], prior$logit_names)
    )
    mean <- prior$mean[indices, , drop = FALSE]
  }
  if (anyNA(eta) || any(!is.finite(eta))) {
    stop("`eta` must contain only finite values.", call. = FALSE)
  }
  list(eta = eta, mean = mean)
}

lnp_assert_prepared_logistic_normal_prior <- function(prior) {
  if (!inherits(prior, "fibrodynmix_logistic_normal_prior")) {
    stop(
      "`prior` must be created by `prepare_logistic_normal_prior()` or `estimate_group_logistic_normal_prior()`.",
      call. = FALSE
    )
  }
  required <- c("mean", "precision", "n_cells", "n_logits", "cell_names", "logit_names")
  missing <- setdiff(required, names(prior))
  if (length(missing) > 0L || !is.matrix(prior$mean) || !is.matrix(prior$precision) ||
    nrow(prior$mean) != prior$n_cells || ncol(prior$mean) != prior$n_logits ||
    nrow(prior$precision) != prior$n_logits || ncol(prior$precision) != prior$n_logits ||
    anyNA(prior$mean) || any(!is.finite(prior$mean)) ||
    anyNA(prior$precision) || any(!is.finite(prior$precision))) {
    stop("`prior` is not a valid prepared logistic-normal prior.", call. = FALSE)
  }
  invisible(prior)
}

lnp_normalize_simplex_matrix <- function(z) {
  if (is.data.frame(z)) {
    z <- as.matrix(z)
  }
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) < 1L || ncol(z) < 2L) {
    stop("`z` must be a numeric matrix with at least two state columns.", call. = FALSE)
  }
  if (anyNA(z) || any(!is.finite(z)) || any(z < 0)) {
    stop("`z` must contain finite, non-negative simplex weights.", call. = FALSE)
  }
  row_total <- rowSums(z)
  if (any(row_total <= 0)) {
    stop("Every row of `z` must have a positive sum.", call. = FALSE)
  }
  z <- sweep(z, 1L, row_total, "/")
  if (is.null(rownames(z))) {
    rownames(z) <- paste0("cell_", seq_len(nrow(z)))
  }
  if (is.null(colnames(z))) {
    colnames(z) <- paste0("state_", seq_len(ncol(z)))
  }
  if (anyNA(rownames(z)) || any(!nzchar(rownames(z))) || anyDuplicated(rownames(z))) {
    stop("`z` row names must be non-missing and unique when supplied.", call. = FALSE)
  }
  if (anyNA(colnames(z)) || any(!nzchar(colnames(z))) || anyDuplicated(colnames(z))) {
    stop("`z` column names must be non-missing and unique when supplied.", call. = FALSE)
  }
  storage.mode(z) <- "double"
  z
}

lnp_normalize_group_id <- function(group_id, cell_names) {
  if (length(group_id) != length(cell_names) || anyNA(group_id)) {
    stop("`group_id` must have one non-missing value per cell.", call. = FALSE)
  }
  if (!is.null(names(group_id))) {
    group_id <- lnp_reorder_named_vector(
      x = group_id,
      target_names = cell_names,
      argument = "group_id"
    )
  }
  group_id <- as.character(group_id)
  if (any(!nzchar(group_id))) {
    stop("`group_id` values must be non-empty.", call. = FALSE)
  }
  unname(group_id)
}

lnp_resolve_reference_state <- function(reference_state, state_names) {
  if (is.character(reference_state)) {
    if (length(reference_state) != 1L || is.na(reference_state)) {
      stop("`reference_state` must name one state or give one state index.", call. = FALSE)
    }
    reference_index <- match(reference_state, state_names)
    if (is.na(reference_index)) {
      stop("`reference_state` does not match a state name.", call. = FALSE)
    }
  } else if (is.numeric(reference_state) && length(reference_state) == 1L &&
    !is.na(reference_state) && reference_state == as.integer(reference_state)) {
    reference_index <- as.integer(reference_state)
    if (reference_index < 1L || reference_index > length(state_names)) {
      stop("`reference_state` is outside the state columns of `z`.", call. = FALSE)
    }
  } else {
    stop("`reference_state` must name one state or give one state index.", call. = FALSE)
  }
  reference_index
}

lnp_resolve_cell_index <- function(cell_index, prior) {
  if (is.character(cell_index)) {
    if (anyNA(cell_index)) {
      stop("`cell_index` cannot contain missing values.", call. = FALSE)
    }
    indices <- match(cell_index, prior$cell_names)
    if (anyNA(indices)) {
      stop("`cell_index` contains unknown cell names.", call. = FALSE)
    }
  } else if (is.numeric(cell_index) && all(!is.na(cell_index)) &&
    all(cell_index == as.integer(cell_index))) {
    indices <- as.integer(cell_index)
    if (any(indices < 1L | indices > prior$n_cells)) {
      stop("`cell_index` is outside the prior cell range.", call. = FALSE)
    }
  } else {
    stop("`cell_index` must contain cell names or positive integer indices.", call. = FALSE)
  }
  if (anyDuplicated(indices)) {
    stop("`cell_index` cannot contain duplicate cells.", call. = FALSE)
  }
  indices
}

lnp_reorder_rectangular_matrix <- function(x, row_names, col_names, argument) {
  if (!is.null(rownames(x))) {
    lnp_validate_matrix_names(rownames(x), row_names, argument, "row")
    x <- x[row_names, , drop = FALSE]
  } else {
    rownames(x) <- row_names
  }
  if (!is.null(colnames(x))) {
    lnp_validate_matrix_names(colnames(x), col_names, argument, "column")
    x <- x[, col_names, drop = FALSE]
  } else {
    colnames(x) <- col_names
  }
  x
}

lnp_reorder_eta_columns <- function(eta, logit_names) {
  if (!is.null(colnames(eta))) {
    lnp_validate_matrix_names(colnames(eta), logit_names, "eta", "column")
    eta <- eta[, logit_names, drop = FALSE]
  } else {
    colnames(eta) <- logit_names
  }
  eta
}

lnp_reorder_eta_rows <- function(eta, cell_names) {
  if (!is.null(rownames(eta))) {
    lnp_validate_matrix_names(rownames(eta), cell_names, "eta", "row")
    eta <- eta[cell_names, , drop = FALSE]
  } else {
    rownames(eta) <- cell_names
  }
  eta
}

lnp_validate_pd_matrix <- function(x,
                                   n_logits,
                                   logit_names,
                                   argument,
                                   check_pd = TRUE,
                                   tolerance = sqrt(.Machine$double.eps)) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) ||
    nrow(x) != n_logits || ncol(x) != n_logits) {
    stop(
      sprintf("`%s` must be a numeric %d by %d matrix.", argument, n_logits, n_logits),
      call. = FALSE
    )
  }
  x <- lnp_reorder_rectangular_matrix(
    x = x,
    row_names = logit_names,
    col_names = logit_names,
    argument = argument
  )
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", argument), call. = FALSE)
  }
  max_asymmetry <- max(abs(x - t(x)))
  scale <- max(1, max(abs(x)))
  if (max_asymmetry > tolerance * scale) {
    stop(sprintf("`%s` must be symmetric.", argument), call. = FALSE)
  }
  x <- 0.5 * (x + t(x))
  dimnames(x) <- list(logit_names, logit_names)
  if (isTRUE(check_pd)) {
    chol_result <- tryCatch(chol(x), error = function(e) NULL)
    if (is.null(chol_result)) {
      stop(sprintf("`%s` must be positive definite.", argument), call. = FALSE)
    }
  }
  storage.mode(x) <- "double"
  x
}

lnp_reorder_named_vector <- function(x, target_names, argument) {
  if (is.null(names(x))) {
    return(x)
  }
  lnp_validate_matrix_names(names(x), target_names, argument, "names")
  unname(x[target_names])
}

lnp_validate_matrix_names <- function(observed, expected, argument, margin) {
  if (length(observed) != length(expected) || anyNA(observed) ||
    any(!nzchar(observed)) || anyDuplicated(observed) ||
    !setequal(observed, expected)) {
    stop(
      sprintf("`%s` %s names must uniquely match the expected names.", argument, margin),
      call. = FALSE
    )
  }
}

lnp_normalize_names <- function(x, n, prefix, argument) {
  if (is.null(x)) {
    return(paste0(prefix, "_", seq_len(n)))
  }
  if (!is.character(x) || length(x) != n || anyNA(x) ||
    any(!nzchar(x)) || anyDuplicated(x)) {
    stop(
      sprintf("`%s` must be a unique, non-missing character vector of length %d.", argument, n),
      call. = FALSE
    )
  }
  unname(x)
}

lnp_assert_positive_integer <- function(x, argument) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
    x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive integer scalar.", argument), call. = FALSE)
  }
  invisible(as.integer(x))
}

lnp_assert_nonnegative_scalar <- function(x, argument) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be a finite, non-negative numeric scalar.", argument), call. = FALSE)
  }
  as.numeric(x)
}
