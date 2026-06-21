#' Transfer a fitted FibroDynMix program to held-out cells
#'
#' Estimates held-out cell state mixtures under a fixed fitted FibroDynMix
#' program. The transfer step keeps `alpha`, `beta`, and `phi` fixed and
#' optimizes only held-out cell simplex weights under the raw-count
#' negative-binomial likelihood.
#'
#' @param counts Gene-by-cell raw count matrix for held-out cells.
#' @param fit A fitted object returned by `fit_fibrodynmix_nb()`.
#' @param library_size Optional held-out library sizes. If omitted, column sums
#'   of `counts` are used.
#' @param z_l2 Non-negative L2 penalty for held-out cell logits.
#' @param maxit_z Maximum `optim()` iterations per held-out cell.
#' @param optimizer Optimizer used by held-out `z` updates.
#' @param optimizer_control Control list passed to `stats::optim()`.
#' @param z_init Optional warm-start state-weight matrix for held-out cells.
#'   Rows must match cells and columns must match fitted states.
#' @param chunk_size Optional number of cells to optimize per chunk. Chunking
#'   bounds memory use and makes long transfers easier to monitor.
#' @param parallel Whether to optimize cells within each chunk in parallel.
#' @param n_workers Number of parallel workers when `parallel = TRUE`.
#' @param return_cell_diagnostics Whether to include a per-cell diagnostics
#'   table with convergence code, iteration counts, and objective values.
#'
#' @return A list with `z_hat`, `nb_objective`, `nb_loglik`, and transfer
#'   diagnostics, including held-out cell simplex convergence rate.
#' @export
fit_fibrodynmix_transfer <- function(counts,
                                     fit,
                                     library_size = NULL,
                                     z_l2 = 0.001,
                                     maxit_z = 35,
                                     optimizer = c("BFGS", "L-BFGS-B"),
                                     optimizer_control = list(),
                                     z_init = NULL,
                                     chunk_size = NULL,
                                     parallel = FALSE,
                                     n_workers = 2,
                                     return_cell_diagnostics = TRUE) {
  if (!is.list(fit) || is.null(fit$beta_hat) || is.null(fit$alpha_hat) || is.null(fit$phi_hat)) {
    stop("`fit` must be returned by `fit_fibrodynmix_nb()`.", call. = FALSE)
  }
  if (!is_matrix_like(counts) || !matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must be a non-negative integer-like numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    stop("`counts` must have gene rownames for transfer.", call. = FALSE)
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  common_genes <- intersect(colnames(fit$beta_hat), rownames(counts))
  if (length(common_genes) < 2L) {
    stop("Transfer requires at least two genes shared with the fitted program.", call. = FALSE)
  }
  counts <- counts[common_genes, , drop = FALSE]
  beta <- fit$beta_hat[, common_genes, drop = FALSE]
  alpha <- fit$alpha_hat[common_genes]
  phi <- fit$phi_hat[common_genes]

  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per held-out cell.", call. = FALSE)
  }
  if (length(z_l2) != 1L || is.na(z_l2) || z_l2 < 0) {
    stop("`z_l2` must be a non-negative numeric scalar.", call. = FALSE)
  }
  assert_positive_integer(maxit_z, "maxit_z")
  if (!is.null(chunk_size)) {
    assert_positive_integer(chunk_size, "chunk_size")
  }
  optimizer <- match.arg(optimizer)
  if (!is.list(optimizer_control)) {
    stop("`optimizer_control` must be a list.", call. = FALSE)
  }
  if (length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  assert_positive_integer(n_workers, "n_workers")

  n_cells <- ncol(counts)
  n_states <- nrow(beta)
  start_z <- resolve_transfer_warm_start(z_init, cell_names = colnames(counts), state_names = rownames(beta))
  transfer <- update_z_nb(
    counts = counts,
    z = start_z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    study_effect = NULL,
    study_id = NULL,
    donor_effect = NULL,
    donor_id = NULL,
    z_l2 = z_l2,
    optimizer = optimizer,
    optimizer_control = optimizer_control,
    maxit = maxit_z,
    chunk_size = chunk_size,
    parallel = parallel,
    n_workers = n_workers
  )

  objective <- fibrodynmix_nb_objective(
    counts = counts,
    z = transfer$z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size,
    average = TRUE
  )
  loglik <- fibrodynmix_nb_loglik(
    counts = counts,
    z = transfer$z,
    beta = beta,
    alpha = alpha,
    phi = phi,
    library_size = library_size
  )

  list(
    z_hat = transfer$z,
    nb_objective = objective,
    nb_loglik = loglik,
    n_cells = n_cells,
    n_genes = nrow(counts),
    shared_genes = common_genes,
    converged = transfer$converged,
    z_convergence_rate = transfer$convergence_rate,
    n_nonconverged_cells = transfer$n_nonconverged,
    cell_converged = transfer$cell_converged,
    cell_diagnostics = if (isTRUE(return_cell_diagnostics)) transfer$cell_diagnostics else NULL,
    chunk_size = transfer$chunk_size,
    n_chunks = transfer$n_chunks,
    parallel = isTRUE(parallel),
    n_workers = if (isTRUE(parallel)) n_workers else 1L
  )
}

resolve_transfer_warm_start <- function(z_init, cell_names, state_names) {
  n_cells <- length(cell_names)
  n_states <- length(state_names)
  if (is.null(z_init)) {
    return(matrix(
      1 / n_states,
      nrow = n_cells,
      ncol = n_states,
      dimnames = list(cell_names, state_names)
    ))
  }
  z_init <- as.matrix(z_init)
  storage.mode(z_init) <- "double"
  if (nrow(z_init) != n_cells || ncol(z_init) != n_states) {
    stop("`z_init` must have one row per held-out cell and one column per fitted state.", call. = FALSE)
  }
  if (!is.null(rownames(z_init))) {
    if (!all(cell_names %in% rownames(z_init))) {
      stop("`z_init` row names must contain all held-out cell names.", call. = FALSE)
    }
    z_init <- z_init[cell_names, , drop = FALSE]
  }
  if (!is.null(colnames(z_init))) {
    if (!all(state_names %in% colnames(z_init))) {
      stop("`z_init` column names must contain all fitted state names.", call. = FALSE)
    }
    z_init <- z_init[, state_names, drop = FALSE]
  }
  if (anyNA(z_init) || any(!is.finite(z_init)) || any(z_init < 0) || any(rowSums(z_init) <= 0)) {
    stop("`z_init` must contain non-negative finite state weights with positive row sums.", call. = FALSE)
  }
  z_init <- z_init / rowSums(z_init)
  dimnames(z_init) <- list(cell_names, state_names)
  z_init
}

#' Run leave-study-out FibroDynMix transfer benchmark
#'
#' Simulates a multi-study dataset, fits FibroDynMix on all but one study, and
#' transfers the fitted program to the held-out study by optimizing held-out
#' cell state mixtures under a fixed raw-count NB program.
#'
#' @param n_replicates Number of simulation replicates.
#' @param seed Optional base seed.
#' @param simulation_args Named list passed to `simulate_fibrodynmix()`.
#' @param train_fit_args Named list passed to `fit_fibrodynmix_nb()` for the
#'   training studies.
#' @param transfer_args Named list passed to `fit_fibrodynmix_transfer()`.
#' @param holdout_studies Optional study identifiers to hold out. If `NULL`,
#'   every simulated study is held out once.
#' @param keep_results Whether to retain full simulation, train fit, and
#'   transfer fit objects as an attribute.
#'
#' @return A data frame with transfer metrics per replicate and held-out study.
#' @export
run_cross_cohort_transfer_benchmark <- function(n_replicates = 2,
                                                seed = 1,
                                                simulation_args = list(),
                                                train_fit_args = list(),
                                                transfer_args = list(),
                                                holdout_studies = NULL,
                                                keep_results = FALSE) {
  assert_positive_integer(n_replicates, "n_replicates")
  rows <- list()
  full_results <- list()
  row_index <- 1L

  for (replicate_id in seq_len(n_replicates)) {
    replicate_seed <- if (is.null(seed)) NULL else as.integer(seed + replicate_id - 1L)
    sim_args <- utils::modifyList(
      list(
        n_studies = 3,
        donors_per_study = 2,
        cells_per_donor = 8,
        n_genes = 120,
        marker_genes_per_state = 5,
        scenario = "batch_confounding",
        seed = replicate_seed
      ),
      simulation_args
    )
    if (is.null(seed)) {
      sim_args$seed <- NULL
    }
    sim <- do.call(simulate_fibrodynmix, sim_args)
    studies <- sort(unique(sim$cell_metadata$study_id))
    holdouts <- if (is.null(holdout_studies)) studies else intersect(holdout_studies, studies)
    if (length(holdouts) == 0L) {
      stop("No requested `holdout_studies` are present in the simulation.", call. = FALSE)
    }

    for (holdout in holdouts) {
      train_cells <- sim$cell_metadata$study_id != holdout
      test_cells <- sim$cell_metadata$study_id == holdout
      train_counts <- sim$counts[, train_cells, drop = FALSE]
      test_counts <- sim$counts[, test_cells, drop = FALSE]
      train_study_id <- sim$cell_metadata$study_id[train_cells]
      train_donor_id <- sim$cell_metadata$donor_id[train_cells]

      train_args <- utils::modifyList(
        list(
          counts = train_counts,
          marker_index = sim$parameters$marker_index,
          library_size = sim$cell_metadata$library_size[train_cells],
          study_id = train_study_id,
          donor_id = train_donor_id,
          fit_study_effect = length(unique(train_study_id)) > 1L,
          fit_donor_effect = length(unique(train_donor_id)) > length(unique(train_study_id)),
          n_outer = 2,
          initializer_args = list(n_iter = 2),
          maxit_beta = 10,
          maxit_z = 10
        ),
        train_fit_args
      )
      train_fit <- do.call(fit_fibrodynmix_nb, train_args)

      heldout_args <- utils::modifyList(
        list(
          counts = test_counts,
          fit = train_fit,
          library_size = sim$cell_metadata$library_size[test_cells],
          maxit_z = 10
        ),
        transfer_args
      )
      transfer_fit <- do.call(fit_fibrodynmix_transfer, heldout_args)
      metrics <- evaluate_state_weights(sim$z[test_cells, , drop = FALSE], transfer_fit$z_hat)

      rows[[row_index]] <- data.frame(
        replicate = replicate_id,
        seed = if (is.null(replicate_seed)) NA_integer_ else replicate_seed,
        holdout_study = holdout,
        n_train_cells = sum(train_cells),
        n_test_cells = sum(test_cells),
        n_shared_genes = length(transfer_fit$shared_genes),
        transfer_rmse = metrics$rmse,
        transfer_mean_absolute_error = metrics$mean_absolute_error,
        transfer_dominant_accuracy = metrics$dominant_accuracy,
        transfer_mean_entropy_pred = metrics$mean_entropy_pred,
        transfer_nb_objective = transfer_fit$nb_objective,
        transfer_nb_loglik = transfer_fit$nb_loglik,
        transfer_converged = transfer_fit$converged,
        transfer_z_convergence_rate = transfer_fit$z_convergence_rate,
        transfer_n_nonconverged_cells = transfer_fit$n_nonconverged_cells,
        train_best_objective = train_fit$best_objective,
        train_stop_reason = train_fit$stop_reason,
        stringsAsFactors = FALSE
      )
      if (keep_results) {
        full_results[[row_index]] <- list(
          simulation = sim,
          train_fit = train_fit,
          transfer_fit = transfer_fit,
          metrics = metrics
        )
      }
      row_index <- row_index + 1L
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (keep_results) {
    attr(result, "results") <- full_results
  }
  result
}
