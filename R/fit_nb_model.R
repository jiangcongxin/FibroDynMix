#' Fit the first FibroDynMix negative-binomial model
#'
#' Fits a lightweight raw-count negative-binomial FibroDynMix model by
#' alternating optimization. By default this optimizer estimates the
#' library-offset count model:
#'
#' `log(mu_ig) = log(library_i) + alpha_g + sum_k z_ik beta_kg`.
#'
#' With `fit_study_effect = TRUE` or `fit_donor_effect = TRUE`, it additionally
#' fits ridge-penalized study-by-gene and donor-by-gene effects.
#'
#' It is intended as the bridge from the log-normalized initializer toward the
#' full hierarchical Bayesian/VI implementation.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of weak prior marker genes for each state.
#' @param library_size Optional vector of per-cell library sizes. If omitted,
#'   column sums of `counts` are used.
#' @param study_id Optional cell-level study identifiers. When supplied with
#'   `fit_study_effect = TRUE`, a ridge-penalized study-by-gene effect is fit.
#' @param donor_id Optional cell-level donor identifiers. When supplied with
#'   `fit_donor_effect = TRUE`, a ridge-penalized donor-by-gene effect is fit.
#' @param n_outer Number of alternating NB optimization iterations.
#' @param initializer_args Named list of arguments passed to
#'   `fit_fibrodynmix_initializer()`.
#' @param estimate_phi Whether to update gene-specific `phi` by
#'   method-of-moments after each outer iteration.
#' @param phi_init Optional initial gene-specific NB size vector.
#' @param beta_l2 Non-negative L2 penalty for state-gene program coefficients.
#' @param marker_l2 Non-negative L2 penalty that keeps prior marker coefficients
#'   oriented toward positive state-specific programs.
#' @param marker_weight Optional per-state, per-gene marker reliability weights.
#' @param fit_study_effect Whether to fit ridge-penalized study-by-gene effects.
#' @param study_l2 Non-negative L2 penalty for study effects.
#' @param fit_donor_effect Whether to fit ridge-penalized donor-by-gene effects.
#' @param donor_l2 Non-negative L2 penalty for donor effects.
#' @param z_l2 Non-negative L2 penalty for cell logits in simplex updates.
#' @param optimizer Optimizer for NB subproblems. One of `"BFGS"` or `"L-BFGS-B"`.
#' @param optimizer_control Control list for `stats::optim()`, merged with
#'   `maxit` defaults.
#' @param maxit_beta Maximum `optim()` iterations for each gene-level
#'   `alpha/beta` update.
#' @param maxit_z Maximum `optim()` iterations for each cell-level `z` update.
#' @param early_stopping Whether to stop when the best objective has not
#'   improved for `patience` iterations.
#' @param patience Number of non-improving iterations tolerated when
#'   `early_stopping = TRUE`.
#' @param min_delta Minimum objective decrease required to count as an
#'   improvement.
#' @param objective_rel_tol Relative objective-tolerance threshold for early
#'   stopping when improvements are flat across `stagnation_window`.
#' @param stagnation_window Number of recent iterations checked for objective
#'   stagnation.
#' @param objective_abs_tol Absolute objective-tolerance threshold for flat
#'   objective checks.
#' @param rollback_to_best Whether to return the best-seen parameters instead
#'   of the last parameters.
#' @param verbose Whether to print objective values.
#'
#' @return A list with fitted `z_hat`, `beta_hat`, `alpha_hat`, `phi_hat`,
#'   `nb_objective_trace`, `converged`, and initializer/optimizer metadata.
#' @export
fit_fibrodynmix_nb <- function(counts,
                               marker_index,
                               library_size = NULL,
                               study_id = NULL,
                               n_outer = 5,
                               initializer_args = list(),
                               estimate_phi = TRUE,
                               phi_init = NULL,
                               beta_l2 = 0.01,
                               marker_l2 = 0.05,
                               marker_weight = NULL,
                               fit_study_effect = FALSE,
                               study_l2 = 0.1,
                               donor_id = NULL,
                               fit_donor_effect = FALSE,
                               donor_l2 = 0.1,
                               z_l2 = 0.001,
                               optimizer = c("BFGS", "L-BFGS-B"),
                               optimizer_control = list(),
                               maxit_beta = 50,
                               maxit_z = 35,
                               early_stopping = TRUE,
                               patience = 2,
                               min_delta = 1e-5,
                               objective_rel_tol = 1e-6,
                               stagnation_window = 5L,
                               objective_abs_tol = 1e-8,
                               rollback_to_best = TRUE,
                               verbose = FALSE) {
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
  assert_positive_integer(n_outer, "n_outer")
  if (length(beta_l2) != 1L || is.na(beta_l2) || beta_l2 < 0) {
    stop("`beta_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(marker_l2) != 1L || is.na(marker_l2) || marker_l2 < 0) {
    stop("`marker_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(study_l2) != 1L || is.na(study_l2) || study_l2 < 0) {
    stop("`study_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(donor_l2) != 1L || is.na(donor_l2) || donor_l2 < 0) {
    stop("`donor_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(z_l2) != 1L || is.na(z_l2) || z_l2 < 0) {
    stop("`z_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  assert_positive_integer(maxit_beta, "maxit_beta")
  assert_positive_integer(maxit_z, "maxit_z")
  assert_positive_integer(patience, "patience")
  optimizer <- match.arg(optimizer)
  if (!is.list(optimizer_control)) {
    stop("`optimizer_control` must be a list.", call. = FALSE)
  }
  if (length(min_delta) != 1L || is.na(min_delta) || min_delta < 0) {
    stop("`min_delta` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (length(objective_rel_tol) != 1L || is.na(objective_rel_tol) || objective_rel_tol < 0) {
    stop("`objective_rel_tol` must be a non-negative numeric scalar.", call. = FALSE)
  }
  assert_positive_integer(stagnation_window, "stagnation_window")
  if (length(objective_abs_tol) != 1L || is.na(objective_abs_tol) || objective_abs_tol < 0) {
    stop("`objective_abs_tol` must be a non-negative numeric scalar.", call. = FALSE)
  }

  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }
  if (isTRUE(fit_study_effect)) {
    if (is.null(study_id) || length(study_id) != ncol(counts) || anyNA(study_id)) {
      stop("`study_id` must contain one non-NA value per cell when `fit_study_effect = TRUE`.", call. = FALSE)
    }
    study_id <- as.character(study_id)
  }
  if (isTRUE(fit_donor_effect)) {
    if (is.null(donor_id) || length(donor_id) != ncol(counts) || anyNA(donor_id)) {
      stop("`donor_id` must contain one non-NA value per cell when `fit_donor_effect = TRUE`.", call. = FALSE)
    }
    donor_id <- as.character(donor_id)
  }

  initializer_call <- utils::modifyList(
    list(
      counts = counts,
      marker_index = marker_index,
      library_size = library_size
    ),
    initializer_args
  )
  init <- do.call(fit_fibrodynmix_initializer, initializer_call)
  z_hat <- init$z_hat
  beta_hat <- init$beta_hat
  marker_target <- build_marker_beta_target(
    marker_index = marker_index,
    gene_names = rownames(counts),
    state_names = colnames(z_hat),
    marker_weight = marker_weight
  )

  alpha_hat <- initialize_nb_alpha(counts, library_size)
  if (is.null(phi_init)) {
    phi_hat <- estimate_phi_moments(counts, library_size)
  } else {
    phi_hat <- validate_gene_vector(phi_init, nrow(counts), rownames(counts), "phi_init")
    if (any(phi_hat <= 0)) {
      stop("`phi_init` must contain positive values.", call. = FALSE)
    }
  }
  study_effect <- if (isTRUE(fit_study_effect)) {
    initialize_group_effect(study_id, rownames(counts))
  } else {
    NULL
  }
  donor_effect <- if (isTRUE(fit_donor_effect)) {
    initialize_group_effect(donor_id, rownames(counts))
  } else {
    NULL
  }

  objective_trace <- numeric(n_outer + 1L)
  objective_trace[1L] <- fibrodynmix_nb_objective(
    counts = counts,
    z = z_hat,
    beta = beta_hat,
    alpha = alpha_hat,
    phi = phi_hat,
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id,
    beta_l2 = beta_l2,
    marker_target = marker_target,
    marker_l2 = marker_l2,
    effect_l2 = study_l2,
    donor_effect_l2 = donor_l2,
    average = TRUE
  )

  converged_beta <- logical(0)
  converged_z <- logical(0)
  converged_study <- logical(0)
  converged_donor <- logical(0)
  rollback_flags <- logical(0)
  best_objective <- objective_trace[1L]
  best_state <- list(
    z = z_hat,
    beta = beta_hat,
    alpha = alpha_hat,
    phi = phi_hat,
    study_effect = study_effect,
    donor_effect = donor_effect,
    iteration = 0L
  )
  stale_iterations <- 0L
  executed_iterations <- 0L
  stop_reason <- "max_iterations"

  for (iter in seq_len(n_outer)) {
    previous_state <- list(z = z_hat, beta = beta_hat, alpha = alpha_hat, phi = phi_hat, study_effect = study_effect, donor_effect = donor_effect)
    previous_objective <- objective_trace[iter]

    beta_update <- update_alpha_beta_nb(
      counts = counts,
      z = z_hat,
      beta = beta_hat,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size,
      study_effect = study_effect,
      study_id = study_id,
      donor_effect = donor_effect,
      donor_id = donor_id,
      beta_l2 = beta_l2,
      marker_target = marker_target,
      marker_l2 = marker_l2,
      optimizer = optimizer,
      optimizer_control = optimizer_control,
      maxit = maxit_beta
    )
    alpha_hat <- beta_update$alpha
    beta_hat <- beta_update$beta
    converged_beta[iter] <- beta_update$converged

    if (isTRUE(fit_study_effect)) {
      study_update <- update_group_effect_nb(
        counts = counts,
        z = z_hat,
        beta = beta_hat,
        alpha = alpha_hat,
        phi = phi_hat,
        library_size = library_size,
        effect_id = study_id,
        effect = study_effect,
        effect_l2 = study_l2,
        other_effect = donor_effect,
        other_id = donor_id,
        optimizer = optimizer,
        optimizer_control = optimizer_control,
        maxit = maxit_beta
      )
      study_effect <- study_update$effect
      converged_study[iter] <- study_update$converged
    }

    if (isTRUE(fit_donor_effect)) {
      donor_update <- update_group_effect_nb(
        counts = counts,
        z = z_hat,
        beta = beta_hat,
        alpha = alpha_hat,
        phi = phi_hat,
        library_size = library_size,
        effect_id = donor_id,
        effect = donor_effect,
        effect_l2 = donor_l2,
        other_effect = study_effect,
        other_id = study_id,
        optimizer = optimizer,
        optimizer_control = optimizer_control,
        maxit = maxit_beta
      )
      donor_effect <- donor_update$effect
      converged_donor[iter] <- donor_update$converged
    }

    z_update <- update_z_nb(
      counts = counts,
      z = z_hat,
      beta = beta_hat,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size,
      study_effect = study_effect,
      study_id = study_id,
      donor_effect = donor_effect,
      donor_id = donor_id,
      z_l2 = z_l2,
      optimizer = optimizer,
      optimizer_control = optimizer_control,
      maxit = maxit_z
    )
    z_hat <- z_update$z
    converged_z[iter] <- z_update$converged

    if (isTRUE(estimate_phi)) {
      mu <- fibrodynmix_nb_mu_public(counts, z_hat, beta_hat, alpha_hat, library_size, study_effect, study_id, donor_effect, donor_id)
      phi_hat <- estimate_phi_moments(counts, library_size, mu = mu)
    }

    objective_trace[iter + 1L] <- fibrodynmix_nb_objective(
      counts = counts,
      z = z_hat,
      beta = beta_hat,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size,
      study_effect = study_effect,
      donor_effect = donor_effect,
      study_id = study_id,
      donor_id = donor_id,
      beta_l2 = beta_l2,
      marker_target = marker_target,
      marker_l2 = marker_l2,
      effect_l2 = study_l2,
      donor_effect_l2 = donor_l2,
      average = TRUE
    )

    executed_iterations <- iter
    rollback_flags[iter] <- FALSE
    if (!is.finite(objective_trace[iter + 1L])) {
      z_hat <- previous_state$z
      beta_hat <- previous_state$beta
      alpha_hat <- previous_state$alpha
      phi_hat <- previous_state$phi
      study_effect <- previous_state$study_effect
      donor_effect <- previous_state$donor_effect
      objective_trace[iter + 1L] <- previous_objective
      rollback_flags[iter] <- TRUE
      stop_reason <- "non_finite_objective"
      break
    }

    if (objective_trace[iter + 1L] > previous_objective + min_delta) {
      z_hat <- previous_state$z
      beta_hat <- previous_state$beta
      alpha_hat <- previous_state$alpha
      phi_hat <- previous_state$phi
      study_effect <- previous_state$study_effect
      donor_effect <- previous_state$donor_effect
      objective_trace[iter + 1L] <- previous_objective
      rollback_flags[iter] <- TRUE
    }

    if (objective_trace[iter + 1L] < best_objective - min_delta) {
      best_objective <- objective_trace[iter + 1L]
      best_state <- list(
        z = z_hat,
        beta = beta_hat,
        alpha = alpha_hat,
        phi = phi_hat,
        study_effect = study_effect,
        donor_effect = donor_effect,
        iteration = iter
      )
      stale_iterations <- 0L
    } else {
      stale_iterations <- stale_iterations + 1L
    }

    if (isTRUE(verbose)) {
      message(sprintf(
        "outer %d NB objective %.6f%s",
        iter,
        objective_trace[iter + 1L],
        if (rollback_flags[iter]) " (rolled back)" else ""
      ))
    }

    if (isTRUE(early_stopping) && iter >= stagnation_window) {
      window_start <- iter + 1L - stagnation_window
      window <- objective_trace[window_start:(iter + 1L)]
      window_gain <- window[1L] - window[length(window)]
      window_rel <- window_gain / max(abs(window[1L]), .Machine$double.eps)
      if (window_gain <= objective_abs_tol || window_rel <= objective_rel_tol) {
        stale_iterations <- stale_iterations + 1L
      } else {
        stale_iterations <- 0L
      }
    }
    if (isTRUE(early_stopping) && stale_iterations >= patience) {
      stop_reason <- "early_stopping"
      break
    }
  }

  objective_trace <- objective_trace[seq_len(executed_iterations + 1L)]
  converged_beta <- converged_beta[seq_len(executed_iterations)]
  converged_z <- converged_z[seq_len(executed_iterations)]
  converged_study <- converged_study[seq_len(executed_iterations)]
  converged_donor <- converged_donor[seq_len(executed_iterations)]
  rollback_flags <- rollback_flags[seq_len(executed_iterations)]

  if (isTRUE(rollback_to_best)) {
    z_hat <- best_state$z
    beta_hat <- best_state$beta
    alpha_hat <- best_state$alpha
    phi_hat <- best_state$phi
    study_effect <- best_state$study_effect
    donor_effect <- best_state$donor_effect
  }

  list(
    z_hat = z_hat,
    beta_hat = beta_hat,
    alpha_hat = alpha_hat,
    phi_hat = phi_hat,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = if (isTRUE(fit_study_effect)) study_id else NULL,
    donor_id = if (isTRUE(fit_donor_effect)) donor_id else NULL,
    nb_objective_trace = objective_trace,
    best_objective = best_objective,
    best_iteration = best_state$iteration,
    executed_iterations = executed_iterations,
    stop_reason = stop_reason,
    converged = length(converged_beta) > 0L && all(converged_beta) && all(converged_z) &&
      (!isTRUE(fit_study_effect) || all(converged_study)) &&
      (!isTRUE(fit_donor_effect) || all(converged_donor)),
    convergence = list(beta = converged_beta, z = converged_z, study = converged_study, donor = converged_donor, rolled_back = rollback_flags),
    marker_target = marker_target,
    initializer = init,
    call = match.call()
  )
}

update_alpha_beta_nb <- function(counts,
                                 z,
                                 beta,
                                 alpha,
                                 phi,
                                 library_size,
                                 study_effect,
                                 study_id,
                                 donor_effect,
                                 donor_id,
                                 beta_l2,
                                 marker_target,
                                 marker_l2,
                                 optimizer,
                                 optimizer_control,
                                 maxit) {
  n_genes <- nrow(counts)
  n_states <- ncol(z)
  new_beta <- beta
  new_alpha <- alpha
  converged <- logical(n_genes)

  for (g in seq_len(n_genes)) {
    y <- matrix_gene_vector(counts, g)
    phi_g <- phi[g]
    start <- c(alpha[g], beta[, g])
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        eta <- par[1L] + as.vector(z %*% par[-1L])
        if (!is.null(study_effect)) {
          eta <- eta + study_effect[study_id, g]
        }
        if (!is.null(donor_effect)) {
          eta <- eta + donor_effect[donor_id, g]
        }
        log_mu <- log(library_size) + eta
        log_mu <- pmin(pmax(log_mu, -745), 700)
        nll <- -sum(stats::dnbinom(y, size = phi_g, mu = exp(log_mu), log = TRUE))
        marker_penalty <- marker_l2 * sum((par[-1L] - marker_target[, g])^2)
        nll + beta_l2 * sum(par[-1L]^2) + marker_penalty
      },
      method = optimizer,
      control = modifyList(
        list(maxit = maxit),
        optimizer_control
      )
    )
    new_alpha[g] <- fit$par[1L]
    new_beta[, g] <- fit$par[-1L]
    converged[g] <- fit$convergence == 0
  }

  colnames(new_beta) <- rownames(counts)
  rownames(new_beta) <- colnames(z)
  names(new_alpha) <- rownames(counts)
  list(alpha = new_alpha, beta = new_beta, converged = mean(converged) >= 0.95)
}

update_z_nb <- function(counts,
                        z,
                        beta,
                        alpha,
                        phi,
                        library_size,
                        study_effect,
                        study_id,
                        donor_effect,
                        donor_id,
                        z_l2,
                        optimizer,
                        optimizer_control,
                        maxit,
                        chunk_size = NULL,
                        parallel = FALSE,
                        n_workers = 2) {
  n_cells <- ncol(counts)
  n_states <- nrow(beta)
  new_z <- z
  if (is.null(chunk_size)) {
    chunk_size <- n_cells
  }
  assert_positive_integer(chunk_size, "chunk_size")
  assert_positive_integer(n_workers, "n_workers")

  chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))
  optimizer_results <- vector("list", n_cells)
  use_parallel <- isTRUE(parallel) && n_workers > 1L && n_cells > 1L

  for (chunk in chunks) {
    run_one <- function(i) {
      optimize_one_z_nb(
        i = i,
        counts = counts,
        z_start = z,
        beta = beta,
        alpha = alpha,
        phi = phi,
        library_size = library_size,
        study_effect = study_effect,
        study_id = study_id,
        donor_effect = donor_effect,
        donor_id = donor_id,
        z_l2 = z_l2,
        optimizer = optimizer,
        optimizer_control = optimizer_control,
        maxit = maxit
      )
    }
    chunk_results <- if (use_parallel) {
      parallel::mclapply(chunk, run_one, mc.cores = n_workers)
    } else {
      lapply(chunk, run_one)
    }
    optimizer_results[chunk] <- chunk_results
  }

  converged <- vapply(optimizer_results, function(x) x$converged, logical(1))
  for (i in seq_len(n_cells)) {
    new_z[i, ] <- optimizer_results[[i]]$z
  }
  rownames(new_z) <- colnames(counts)
  colnames(new_z) <- rownames(beta)
  convergence_rate <- mean(converged)
  cell_diagnostics <- data.frame(
    cell_id = colnames(counts),
    converged = converged,
    convergence_code = vapply(optimizer_results, function(x) x$convergence_code, integer(1)),
    objective = vapply(optimizer_results, function(x) x$objective, numeric(1)),
    function_evaluations = vapply(optimizer_results, function(x) x$function_evaluations, integer(1)),
    gradient_evaluations = vapply(optimizer_results, function(x) x$gradient_evaluations, integer(1)),
    max_state_weight = apply(new_z, 1L, max),
    entropy = -rowSums(pmax(new_z, .Machine$double.eps) * log(pmax(new_z, .Machine$double.eps))),
    stringsAsFactors = FALSE
  )
  list(
    z = new_z,
    converged = convergence_rate >= 0.95,
    convergence_rate = convergence_rate,
    n_nonconverged = sum(!converged),
    cell_converged = converged,
    cell_diagnostics = cell_diagnostics,
    chunk_size = chunk_size,
    n_chunks = length(chunks)
  )
}

optimize_one_z_nb <- function(i,
                              counts,
                              z_start,
                              beta,
                              alpha,
                              phi,
                              library_size,
                              study_effect,
                              study_id,
                              donor_effect,
                              donor_id,
                              z_l2,
                              optimizer,
                              optimizer_control,
                              maxit) {
  y <- matrix_cell_vector(counts, i)
  start <- simplex_to_logits(z_start[i, ])
  fit <- stats::optim(
    par = start,
    fn = function(logits) {
      z_i <- logits_to_simplex(logits)
      eta <- as.vector(z_i %*% beta) + alpha
      if (!is.null(study_effect)) {
        eta <- eta + study_effect[study_id[i], ]
      }
      if (!is.null(donor_effect)) {
        eta <- eta + donor_effect[donor_id[i], ]
      }
      log_mu <- log(library_size[i]) + eta
      log_mu <- pmin(pmax(log_mu, -745), 700)
      nll <- -sum(stats::dnbinom(y, size = phi, mu = exp(log_mu), log = TRUE))
      nll + z_l2 * sum(logits^2)
    },
    method = optimizer,
    control = modifyList(
      list(maxit = maxit),
      optimizer_control
    )
  )
  list(
    z = logits_to_simplex(fit$par),
    converged = fit$convergence == 0,
    convergence_code = as.integer(fit$convergence),
    objective = as.numeric(fit$value),
    function_evaluations = as.integer(fit$counts[["function"]]),
    gradient_evaluations = as.integer(fit$counts[["gradient"]])
  )
}

update_group_effect_nb <- function(counts,
                                   z,
                                   beta,
                                   alpha,
                                   phi,
                                   library_size,
                                   effect_id,
                                   effect,
                                   effect_l2,
                                   other_effect = NULL,
                                   other_id = NULL,
                                   optimizer,
                                   optimizer_control,
                                   maxit) {
  effect_levels <- rownames(effect)
  n_genes <- nrow(counts)
  new_effect <- effect
  converged <- logical(n_genes)
  base_eta <- z %*% beta
  base_eta <- sweep(base_eta, 2, alpha, "+")
  if (!is.null(other_effect)) {
    base_eta <- base_eta + other_effect[other_id, , drop = FALSE]
  }

  for (g in seq_len(n_genes)) {
    y <- matrix_gene_vector(counts, g)
    phi_g <- phi[g]
    start <- effect[, g]
    fit <- stats::optim(
      par = start,
      fn = function(effect_g) {
        names(effect_g) <- effect_levels
        eta <- base_eta[, g] + effect_g[effect_id]
        log_mu <- log(library_size) + eta
        log_mu <- pmin(pmax(log_mu, -745), 700)
        nll <- -sum(stats::dnbinom(y, size = phi_g, mu = exp(log_mu), log = TRUE))
        nll + effect_l2 * sum(effect_g^2)
      },
      method = optimizer,
      control = modifyList(
        list(maxit = maxit),
        optimizer_control
      )
    )
    new_effect[, g] <- fit$par
    converged[g] <- fit$convergence == 0
  }

  rownames(new_effect) <- effect_levels
  colnames(new_effect) <- rownames(counts)
  list(effect = new_effect, converged = mean(converged) >= 0.95)
}

initialize_nb_alpha <- function(counts, library_size) {
  alpha <- log((matrix_row_means_core(counts) + 0.1) / mean(library_size))
  names(alpha) <- rownames(counts)
  alpha
}

estimate_phi_moments <- function(counts, library_size, mu = NULL, min_phi = 0.1, max_phi = 1000) {
  if (is.null(mu)) {
    gene_rate <- (matrix_row_sums(counts) + 0.1) / sum(library_size)
    mu <- outer(gene_rate, library_size)
    dimnames(mu) <- dimnames(counts)
  }
  mean_mu <- rowMeans(mu)
  residual_var <- vapply(seq_len(nrow(counts)), function(g) {
    mean((matrix_gene_vector(counts, g) - mu[g, ])^2)
  }, numeric(1))
  denom <- pmax(residual_var - mean_mu, .Machine$double.eps)
  phi <- mean_mu^2 / denom
  phi[!is.finite(phi)] <- max_phi
  phi <- pmin(pmax(phi, min_phi), max_phi)
  names(phi) <- rownames(counts)
  phi
}

build_marker_beta_target <- function(marker_index,
                                     gene_names,
                                     state_names,
                                     marker_value = 1,
                                     background_value = 0,
                                     marker_weight = NULL) {
  marker_target <- matrix(
    background_value,
    nrow = length(state_names),
    ncol = length(gene_names),
    dimnames = list(state_names, gene_names)
  )
  for (state in state_names) {
    if (!state %in% names(marker_index)) {
      next
    }
    marker_rows <- resolve_marker_rows(marker_index[[state]], gene_names)
    if (length(marker_rows) > 0L) {
      marker_target[state, marker_rows] <- marker_value
    }
  }

  if (!is.null(marker_weight)) {
    marker_weight <- expand_marker_weight_matrix(
      marker_weight = marker_weight,
      state_names = state_names,
      gene_names = gene_names,
      marker_target = marker_target
    )
    marker_target <- marker_target * marker_weight
  }

  marker_target
}

expand_marker_weight_matrix <- function(marker_weight,
                                        state_names,
                                        gene_names,
                                        marker_target) {
  n_states <- length(state_names)
  n_genes <- length(gene_names)

  if (is.matrix(marker_weight) || is.data.frame(marker_weight)) {
    marker_weight <- as.matrix(marker_weight)
    if (nrow(marker_weight) != n_states || ncol(marker_weight) != n_genes) {
      stop("`marker_weight` matrix must have one row per state and one column per gene.", call. = FALSE)
    }
    if (!is.null(rownames(marker_weight))) {
      row_match <- match(state_names, rownames(marker_weight))
      if (anyNA(row_match)) {
        stop("`marker_weight` row names must include all state names.", call. = FALSE)
      }
      marker_weight <- marker_weight[row_match, , drop = FALSE]
    }
    if (!is.null(colnames(marker_weight))) {
      col_match <- match(gene_names, colnames(marker_weight))
      if (anyNA(col_match)) {
        stop("`marker_weight` column names must include all gene names.", call. = FALSE)
      }
      marker_weight <- marker_weight[, col_match, drop = FALSE]
    }
  } else if (is.numeric(marker_weight) && is.vector(marker_weight)) {
    weights <- as.numeric(marker_weight)
    if (length(weights) == 1L) {
      marker_weight <- matrix(weights, nrow = n_states, ncol = n_genes)
    } else if (length(weights) == n_states) {
      if (!is.null(names(weights))) {
        row_match <- match(state_names, names(weights))
        if (anyNA(row_match)) {
          stop("`marker_weight` names must cover all state names when length equals `length(state_names)`.", call. = FALSE)
        }
        weights <- weights[row_match]
      }
      marker_weight <- matrix(rep(weights, times = n_genes), nrow = n_states, ncol = n_genes)
    } else if (length(weights) == n_genes) {
      if (!is.null(names(weights))) {
        col_match <- match(gene_names, names(weights))
        if (anyNA(col_match)) {
          stop("`marker_weight` names must cover all gene names when length equals `length(gene_names)`.", call. = FALSE)
        }
        weights <- weights[col_match]
      }
      marker_weight <- matrix(rep(weights, times = n_states), nrow = n_states, ncol = n_genes, byrow = TRUE)
    } else if (length(weights) == n_states * n_genes) {
      marker_weight <- matrix(weights, nrow = n_states, ncol = n_genes, byrow = TRUE)
    } else {
      stop("`marker_weight` must be a matrix or a numeric vector of length 1, `length(state_names)`, `length(gene_names)`, or `length(state_names) * length(gene_names)`.", call. = FALSE)
    }
  } else {
    stop("`marker_weight` must be a numeric matrix/data frame or numeric vector.", call. = FALSE)
  }

  if (anyNA(marker_weight) || any(!is.finite(marker_weight))) {
    stop("`marker_weight` must contain only finite values.", call. = FALSE)
  }
  if (any(marker_weight < 0)) {
    stop("`marker_weight` values must be non-negative.", call. = FALSE)
  }

  dimnames(marker_weight) <- dimnames(marker_target)
  marker_weight
}

simplex_to_logits <- function(z_i) {
  z_i <- pmax(as.numeric(z_i), .Machine$double.eps)
  log(z_i[-length(z_i)] / z_i[length(z_i)])
}

logits_to_simplex <- function(logits) {
  logits <- c(as.numeric(logits), 0)
  logits <- logits - max(logits)
  exp_logits <- exp(logits)
  exp_logits / sum(exp_logits)
}

initialize_group_effect <- function(group_id, gene_names) {
  group_levels <- sort(unique(as.character(group_id)))
  matrix(
    0,
    nrow = length(group_levels),
    ncol = length(gene_names),
    dimnames = list(group_levels, gene_names)
  )
}

fibrodynmix_nb_mu_public <- function(counts,
                                     z,
                                     beta,
                                     alpha,
                                     library_size,
                                     study_effect = NULL,
                                     study_id = NULL,
                                     donor_effect = NULL,
                                     donor_id = NULL) {
  parts <- prepare_nb_inputs(
    counts = counts,
    z = z,
    beta = beta,
    alpha = alpha,
    phi = rep(1, nrow(counts)),
    library_size = library_size,
    study_effect = study_effect,
    donor_effect = donor_effect,
    study_id = study_id,
    donor_id = donor_id
  )
  fibrodynmix_nb_mu_from_prepared(parts)
}
