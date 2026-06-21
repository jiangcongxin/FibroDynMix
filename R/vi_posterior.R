#' Fit a lightweight logistic-normal variational posterior for FibroDynMix
#'
#' Fits the current negative-binomial FibroDynMix optimizer, then constructs a
#' mean-field logistic-normal variational posterior over each cell's latent
#' simplex weights. This is a posterior layer around the implemented raw-count
#' NB program: `alpha`, `beta`, `phi`, study effects, and donor effects are fixed
#' at the fitted mode while `q(z_i)` is represented in cell-logit space.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of weak prior marker genes for each state.
#' @param library_size Optional vector of per-cell library sizes. If omitted,
#'   column sums of `counts` are used.
#' @param nb_args Named list passed to `fit_fibrodynmix_nb()`.
#' @param n_draws Number of posterior draws used for summaries.
#' @param n_elbo_draws Number of Monte Carlo draws used per ELBO evaluation.
#' @param n_vi_iter Number of scale-refinement iterations for the diagonal
#'   logistic-normal posterior.
#' @param posterior_scale Initial logit posterior standard deviation.
#' @param scale_grid Multipliers evaluated around the current posterior scale.
#' @param prior_sd Logistic-normal prior standard deviation for cell logits.
#' @param probs Lower and upper probabilities for credible intervals.
#' @param seed Optional random seed.
#' @param cell_metadata Optional cell metadata for sample-level composition
#'   posterior summaries.
#' @param sample_col Optional column in `cell_metadata` used for sample-level
#'   composition summaries.
#' @param keep_draws Whether to retain full posterior cell-level draws.
#'
#' @return A list with the NB fit, variational parameters, posterior draws or
#'   summaries, credible intervals, and ELBO-like trace.
#' @export
fit_fibrodynmix_vi <- function(counts,
                               marker_index,
                               library_size = NULL,
                               nb_args = list(),
                               n_draws = 50,
                               n_elbo_draws = 10,
                               n_vi_iter = 3,
                               posterior_scale = 0.25,
                               scale_grid = c(0.5, 1, 1.5),
                               prior_sd = 2,
                               probs = c(0.025, 0.975),
                               seed = 1,
                               cell_metadata = NULL,
                               sample_col = NULL,
                               keep_draws = TRUE) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || anyNA(counts) || any(counts < 0) || any(counts != round(counts))) {
    stop("`counts` must be a non-negative integer-like numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (is.null(library_size)) {
    library_size <- colSums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }
  assert_positive_integer(n_draws, "n_draws")
  assert_positive_integer(n_elbo_draws, "n_elbo_draws")
  assert_positive_integer(n_vi_iter, "n_vi_iter")
  if (length(posterior_scale) != 1L || is.na(posterior_scale) || posterior_scale <= 0) {
    stop("`posterior_scale` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(scale_grid) || length(scale_grid) == 0L || anyNA(scale_grid) || any(scale_grid <= 0)) {
    stop("`scale_grid` must contain positive numeric multipliers.", call. = FALSE)
  }
  if (length(prior_sd) != 1L || is.na(prior_sd) || prior_sd <= 0) {
    stop("`prior_sd` must be a positive numeric scalar.", call. = FALSE)
  }
  if (length(probs) != 2L || anyNA(probs) || any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must contain two increasing probabilities in [0, 1].", call. = FALSE)
  }
  if (!is.null(cell_metadata)) {
    cell_metadata <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
    if (nrow(cell_metadata) != ncol(counts)) {
      stop("`cell_metadata` must contain one row per cell.", call. = FALSE)
    }
    if (is.null(rownames(cell_metadata)) || any(!colnames(counts) %in% rownames(cell_metadata))) {
      rownames(cell_metadata) <- colnames(counts)
    }
    if (is.null(sample_col)) {
      cell_metadata$vi_sample <- "all"
      sample_col <- "vi_sample"
    }
    if (!sample_col %in% colnames(cell_metadata)) {
      stop("`sample_col` must name a column in `cell_metadata`.", call. = FALSE)
    }
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  nb_call <- utils::modifyList(
    list(counts = counts, marker_index = marker_index, library_size = library_size),
    nb_args
  )
  nb_fit <- do.call(fit_fibrodynmix_nb, nb_call)
  eta_mean <- t(apply(nb_fit$z_hat, 1L, simplex_to_logits))
  rownames(eta_mean) <- rownames(nb_fit$z_hat)
  colnames(eta_mean) <- paste0("eta_", seq_len(ncol(eta_mean)))

  current_sd <- matrix(
    posterior_scale,
    nrow = nrow(eta_mean),
    ncol = ncol(eta_mean),
    dimnames = dimnames(eta_mean)
  )
  elbo_rows <- list()
  row_index <- 1L
  best_elbo <- -Inf
  best_sd <- current_sd

  for (iter in seq_len(n_vi_iter)) {
    candidate_scales <- unique(as.numeric(scale_grid) * current_sd[1L, 1L])
    candidate_scales <- candidate_scales[candidate_scales > 0]
    for (scale in candidate_scales) {
      eta_sd <- matrix(scale, nrow = nrow(eta_mean), ncol = ncol(eta_mean), dimnames = dimnames(eta_mean))
      elbo <- estimate_logistic_normal_elbo(
        counts = counts,
        eta_mean = eta_mean,
        eta_sd = eta_sd,
        fit = nb_fit,
        library_size = library_size,
        prior_sd = prior_sd,
        n_draws = n_elbo_draws
      )
      elbo_rows[[row_index]] <- data.frame(
        iteration = iter,
        posterior_scale = scale,
        expected_loglik = elbo$expected_loglik,
        log_prior = elbo$log_prior,
        entropy = elbo$entropy,
        elbo = elbo$elbo,
        stringsAsFactors = FALSE
      )
      if (elbo$elbo > best_elbo) {
        best_elbo <- elbo$elbo
        best_sd <- eta_sd
      }
      row_index <- row_index + 1L
    }
    current_sd <- best_sd
  }

  elbo_trace <- do.call(rbind, elbo_rows)
  rownames(elbo_trace) <- NULL
  posterior_draws <- draw_logistic_normal_posterior(
    eta_mean = eta_mean,
    eta_sd = best_sd,
    n_draws = n_draws,
    state_names = colnames(nb_fit$z_hat)
  )
  cell_summary <- summarize_vi_cell_draws(posterior_draws, probs = probs)
  sample_summary <- NULL
  if (!is.null(cell_metadata)) {
    sample_draws <- summarize_vi_sample_draws(
      posterior_draws = posterior_draws,
      cell_metadata = cell_metadata,
      sample_col = sample_col
    )
    sample_summary <- summarize_draws(sample_draws, group_cols = c("sample_id", "state"), value_col = "composition", probs = probs)
  } else {
    sample_draws <- NULL
  }

  result <- list(
    method = "logistic_normal_vi",
    nb_fit = nb_fit,
    eta_mean = eta_mean,
    eta_sd = best_sd,
    z_mean = vi_z_mean_from_draws(posterior_draws, state_names = colnames(nb_fit$z_hat), cell_names = rownames(nb_fit$z_hat)),
    cell_summary = cell_summary,
    sample_summary = sample_summary,
    sample_draws = sample_draws,
    elbo_trace = elbo_trace,
    best_elbo = best_elbo,
    prior_sd = prior_sd,
    probs = probs
  )
  if (isTRUE(keep_draws)) {
    result$cell_draws <- posterior_draws
  }
  class(result) <- c("FibroDynMixVI", class(result))
  result
}

#' Evaluate posterior interval calibration for state weights
#'
#' Compares true simulated latent state weights with posterior credible
#' intervals returned by `fit_fibrodynmix_vi()`.
#'
#' @param z_true Cell-by-state matrix of true simplex weights.
#' @param cell_summary_z Data frame containing `cell_id`, `state`, `lower`,
#'   `upper`, and optionally `mean` columns.
#'
#' @return A list with scalar interval coverage and width diagnostics.
#' @export
evaluate_posterior_intervals <- function(z_true, cell_summary_z) {
  z_true <- as.matrix(z_true)
  required <- c("cell_id", "state", "lower", "upper")
  missing <- setdiff(required, colnames(cell_summary_z))
  if (length(missing) > 0L) {
    stop(sprintf("`cell_summary_z` is missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (is.null(rownames(z_true))) {
    rownames(z_true) <- sprintf("cell_%d", seq_len(nrow(z_true)))
  }
  if (is.null(colnames(z_true))) {
    colnames(z_true) <- sprintf("state_%d", seq_len(ncol(z_true)))
  }

  keep <- cell_summary_z$cell_id %in% rownames(z_true) & cell_summary_z$state %in% colnames(z_true)
  if (!any(keep)) {
    stop("No posterior interval rows match `z_true` cell and state names.", call. = FALSE)
  }
  intervals <- cell_summary_z[keep, , drop = FALSE]
  truth <- z_true[cbind(intervals$cell_id, intervals$state)]
  covered <- truth >= intervals$lower & truth <= intervals$upper
  width <- intervals$upper - intervals$lower
  mean_error <- if ("mean" %in% colnames(intervals)) {
    mean(abs(intervals$mean - truth), na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    interval_coverage = mean(covered, na.rm = TRUE),
    mean_interval_width = mean(width, na.rm = TRUE),
    median_interval_width = stats::median(width, na.rm = TRUE),
    posterior_mean_absolute_error = mean_error,
    n_interval_rows = nrow(intervals)
  )
}

#' Calibrate posterior interval scale in simulation
#'
#' Uses known simulated latent state weights to choose a multiplicative expansion
#' of posterior intervals. This is intended for simulation calibration and
#' method diagnostics, not for real-data inference where latent truth is unknown.
#'
#' @inheritParams evaluate_posterior_intervals
#' @param scale_grid Positive interval-width multipliers to evaluate.
#' @param target_coverage Target empirical coverage.
#'
#' @return A list with the selected scale and calibrated interval diagnostics.
#' @export
calibrate_posterior_interval_scale <- function(z_true,
                                               cell_summary_z,
                                               scale_grid = c(1, 1.5, 2, 3, 4, 6, 8, 12, 24, 100),
                                               target_coverage = 0.9) {
  if (!is.numeric(scale_grid) || length(scale_grid) == 0L || anyNA(scale_grid) || any(scale_grid <= 0)) {
    stop("`scale_grid` must contain positive numeric values.", call. = FALSE)
  }
  if (length(target_coverage) != 1L || is.na(target_coverage) || target_coverage <= 0 || target_coverage > 1) {
    stop("`target_coverage` must be in (0, 1].", call. = FALSE)
  }
  if (!"mean" %in% colnames(cell_summary_z)) {
    stop("`cell_summary_z` must contain a `mean` column for interval scaling.", call. = FALSE)
  }

  rows <- vector("list", length(scale_grid))
  for (i in seq_along(scale_grid)) {
    scaled <- scale_posterior_intervals(cell_summary_z, scale_grid[i])
    metrics <- evaluate_posterior_intervals(z_true, scaled)
    rows[[i]] <- data.frame(
      interval_scale = scale_grid[i],
      interval_coverage = metrics$interval_coverage,
      mean_interval_width = metrics$mean_interval_width,
      median_interval_width = metrics$median_interval_width,
      posterior_mean_absolute_error = metrics$posterior_mean_absolute_error,
      stringsAsFactors = FALSE
    )
  }
  calibration <- do.call(rbind, rows)
  candidates <- calibration[calibration$interval_coverage >= target_coverage, , drop = FALSE]
  if (nrow(candidates) > 0L) {
    selected <- candidates[which.min(candidates$mean_interval_width), , drop = FALSE]
  } else {
    selected <- calibration[which.max(calibration$interval_coverage), , drop = FALSE]
  }
  list(
    selected_scale = selected$interval_scale[1],
    target_coverage = target_coverage,
    interval_coverage = selected$interval_coverage[1],
    mean_interval_width = selected$mean_interval_width[1],
    median_interval_width = selected$median_interval_width[1],
    posterior_mean_absolute_error = selected$posterior_mean_absolute_error[1],
    calibration_table = calibration
  )
}

scale_posterior_intervals <- function(cell_summary_z, interval_scale) {
  scaled <- cell_summary_z
  lower_width <- pmax(scaled$mean - scaled$lower, 0)
  upper_width <- pmax(scaled$upper - scaled$mean, 0)
  half_width <- pmax(lower_width, upper_width, 1e-6)
  scaled$lower <- pmax(0, scaled$mean - interval_scale * half_width)
  scaled$upper <- pmin(1, scaled$mean + interval_scale * half_width)
  scaled
}

estimate_logistic_normal_elbo <- function(counts,
                                          eta_mean,
                                          eta_sd,
                                          fit,
                                          library_size,
                                          prior_sd,
                                          n_draws) {
  loglik <- numeric(n_draws)
  for (draw in seq_len(n_draws)) {
    z_draw <- eta_draw_to_z(draw_eta_matrix(eta_mean, eta_sd), state_names = rownames(fit$beta_hat))
    loglik[draw] <- fibrodynmix_nb_loglik(
      counts = counts,
      z = z_draw,
      beta = fit$beta_hat,
      alpha = fit$alpha_hat,
      phi = fit$phi_hat,
      library_size = library_size,
      study_effect = fit$study_effect,
      donor_effect = fit$donor_effect,
      study_id = fit$study_id,
      donor_id = fit$donor_id
    )
  }
  log_prior <- sum(stats::dnorm(as.vector(eta_mean), mean = 0, sd = prior_sd, log = TRUE))
  entropy <- sum(0.5 * log(2 * pi * exp(1) * as.vector(eta_sd)^2))
  expected_loglik <- mean(loglik)
  list(
    expected_loglik = expected_loglik,
    log_prior = log_prior,
    entropy = entropy,
    elbo = expected_loglik + log_prior + entropy
  )
}

draw_logistic_normal_posterior <- function(eta_mean, eta_sd, n_draws, state_names) {
  draws <- vector("list", n_draws)
  for (draw in seq_len(n_draws)) {
    z <- eta_draw_to_z(draw_eta_matrix(eta_mean, eta_sd), state_names = state_names)
    draws[[draw]] <- vi_cell_draw_rows(z = z, draw = draw)
  }
  result <- do.call(rbind, draws)
  rownames(result) <- NULL
  result
}

draw_eta_matrix <- function(eta_mean, eta_sd) {
  eta <- stats::rnorm(length(eta_mean), mean = as.vector(eta_mean), sd = as.vector(eta_sd))
  matrix(eta, nrow = nrow(eta_mean), ncol = ncol(eta_mean), dimnames = dimnames(eta_mean))
}

eta_draw_to_z <- function(eta, state_names) {
  z <- t(apply(eta, 1L, logits_to_simplex))
  rownames(z) <- rownames(eta)
  colnames(z) <- state_names
  z
}

vi_cell_draw_rows <- function(z, draw) {
  entropy <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))
  rows <- vector("list", ncol(z))
  for (k in seq_len(ncol(z))) {
    rows[[k]] <- data.frame(
      draw = draw,
      cell_id = rownames(z),
      state = colnames(z)[k],
      z = z[, k],
      entropy = entropy,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_vi_cell_draws <- function(posterior_draws, probs) {
  z_summary <- summarize_draws(posterior_draws, group_cols = c("cell_id", "state"), value_col = "z", probs = probs)
  entropy_draws <- unique(posterior_draws[, c("draw", "cell_id", "entropy"), drop = FALSE])
  entropy_summary <- summarize_draws(
    transform(entropy_draws, state = "entropy"),
    group_cols = c("cell_id", "state"),
    value_col = "entropy",
    probs = probs
  )
  colnames(entropy_summary)[colnames(entropy_summary) == "state"] <- "quantity"
  list(z = z_summary, entropy = entropy_summary)
}

summarize_vi_sample_draws <- function(posterior_draws, cell_metadata, sample_col) {
  cell_order <- unique(posterior_draws$cell_id)
  sample_id <- as.character(cell_metadata[cell_order, sample_col])
  if (anyNA(sample_id)) {
    stop("Could not align `cell_metadata` to posterior draw cell IDs.", call. = FALSE)
  }
  names(sample_id) <- cell_order

  merged <- posterior_draws
  merged$sample_id <- sample_id[merged$cell_id]
  result <- stats::aggregate(z ~ draw + sample_id + state, data = merged, FUN = mean)
  colnames(result)[colnames(result) == "z"] <- "composition"
  result
}

vi_z_mean_from_draws <- function(posterior_draws, state_names, cell_names) {
  z_summary <- stats::aggregate(z ~ cell_id + state, data = posterior_draws, FUN = mean)
  z <- matrix(
    NA_real_,
    nrow = length(cell_names),
    ncol = length(state_names),
    dimnames = list(cell_names, state_names)
  )
  for (i in seq_len(nrow(z_summary))) {
    z[z_summary$cell_id[i], z_summary$state[i]] <- z_summary$z[i]
  }
  z
}
