#' Bootstrap FibroDynMix state uncertainty
#'
#' Runs a lightweight cell bootstrap to estimate uncertainty in cell-state
#' mixtures, entropy, sample-level composition, and marker-program stability.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of weak prior marker genes for each state.
#' @param library_size Optional vector of per-cell library sizes. If omitted,
#'   column sums of `counts` are used.
#' @param cell_metadata Optional cell metadata. If supplied, `sample_col` can be
#'   used for sample-level composition intervals.
#' @param sample_col Column in `cell_metadata` used to aggregate sample-level
#'   composition. If `NULL`, all cells are treated as one sample.
#' @param method Fitting method: `initializer`, `nb`, or `nb_study`.
#' @param n_boot Number of bootstrap replicates.
#' @param seed Optional random seed.
#' @param fit_args Additional arguments passed to the selected fit function.
#' @param keep_fits Whether to retain full fit objects for each replicate.
#'
#' @return A list containing bootstrap draws and uncertainty summaries.
#' @export
bootstrap_fibrodynmix <- function(counts,
                                  marker_index,
                                  library_size = NULL,
                                  cell_metadata = NULL,
                                  sample_col = NULL,
  method = c("initializer", "nb", "nb_study"),
                                  n_boot = 20,
                                  seed = 1,
                                  fit_args = list(),
  keep_fits = FALSE) {
  method <- match.arg(method)
  if (!is_matrix_like(counts) || !matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must be a non-negative integer-like numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (is.null(library_size)) {
    library_size <- matrix_col_sums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size <= 0)) {
    stop("`library_size` must be a positive vector with one value per cell.", call. = FALSE)
  }
  assert_positive_integer(n_boot, "n_boot")

  if (is.null(cell_metadata)) {
    cell_metadata <- data.frame(cell_id = colnames(counts), bootstrap_sample = "all", stringsAsFactors = FALSE)
    sample_col <- "bootstrap_sample"
  } else {
    cell_metadata <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
    if (nrow(cell_metadata) != ncol(counts)) {
      stop("`cell_metadata` must contain one row per cell.", call. = FALSE)
    }
    if (is.null(sample_col)) {
      cell_metadata$bootstrap_sample <- "all"
      sample_col <- "bootstrap_sample"
    }
    if (!sample_col %in% colnames(cell_metadata)) {
      stop("`sample_col` must be a column in `cell_metadata`.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  base_fit <- fit_bootstrap_method(
    method = method,
    counts = counts,
    marker_index = marker_index,
    library_size = library_size,
    cell_metadata = cell_metadata,
    fit_args = fit_args
  )
  state_names <- colnames(base_fit$z_hat)
  base_cell_summary <- cell_uncertainty_rows(
    z = base_fit$z_hat,
    replicate = 0L,
    cell_id = rownames(base_fit$z_hat)
  )
  base_sample_summary <- sample_composition_rows(
    z = base_fit$z_hat,
    cell_metadata = cell_metadata,
    sample_col = sample_col,
    replicate = 0L
  )
  marker_rows <- marker_stability_rows(base_fit$beta_hat, replicate = 0L)

  cell_rows <- list(base_cell_summary)
  sample_rows <- list(base_sample_summary)
  fit_objects <- if (keep_fits) list(base_fit) else NULL

  for (b in seq_len(n_boot)) {
    boot_idx <- sample(seq_len(ncol(counts)), size = ncol(counts), replace = TRUE)
    boot_counts <- counts[, boot_idx, drop = FALSE]
    boot_library <- library_size[boot_idx]
    boot_metadata <- cell_metadata[boot_idx, , drop = FALSE]
    boot_cell_names <- paste0(colnames(counts)[boot_idx], "__boot", b, "_", seq_along(boot_idx))
    colnames(boot_counts) <- boot_cell_names

    fit <- fit_bootstrap_method(
      method = method,
      counts = boot_counts,
      marker_index = marker_index,
      library_size = boot_library,
      cell_metadata = boot_metadata,
      fit_args = fit_args
    )

    cell_rows[[b + 1L]] <- cell_uncertainty_rows(
      z = fit$z_hat,
      replicate = b,
      cell_id = colnames(counts)[boot_idx]
    )
    sample_rows[[b + 1L]] <- sample_composition_rows(
      z = fit$z_hat,
      cell_metadata = boot_metadata,
      sample_col = sample_col,
      replicate = b
    )
    marker_rows <- rbind(marker_rows, marker_stability_rows(fit$beta_hat, replicate = b))
    if (keep_fits) {
      fit_objects[[b + 1L]] <- fit
    }
  }

  cell_draws <- do.call(rbind, cell_rows)
  sample_draws <- do.call(rbind, sample_rows)
  rownames(cell_draws) <- NULL
  rownames(sample_draws) <- NULL
  rownames(marker_rows) <- NULL

  summaries <- summarize_bootstrap_uncertainty(
    cell_draws = cell_draws,
    sample_draws = sample_draws,
    marker_draws = marker_rows
  )

  result <- list(
    method = method,
    n_boot = n_boot,
    state_names = state_names,
    base_fit = base_fit,
    cell_draws = cell_draws,
    sample_draws = sample_draws,
    marker_draws = marker_rows,
    cell_summary = summaries$cell_summary,
    sample_summary = summaries$sample_summary,
    marker_summary = summaries$marker_summary
  )
  if (keep_fits) {
    result$fits <- fit_objects
  }
  result
}

#' Summarize FibroDynMix bootstrap uncertainty draws
#'
#' @param cell_draws Cell-level bootstrap state/entropy draws.
#' @param sample_draws Sample-level bootstrap composition draws.
#' @param marker_draws State-gene program bootstrap draws.
#' @param probs Lower and upper quantiles used for intervals.
#'
#' @return A list with cell, sample, and marker summary tables.
#' @export
summarize_bootstrap_uncertainty <- function(cell_draws,
                                            sample_draws,
                                            marker_draws,
                                            probs = c(0.025, 0.975)) {
  if (length(probs) != 2L || anyNA(probs) || any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must contain two increasing probabilities in [0, 1].", call. = FALSE)
  }
  list(
    cell_summary = summarize_draws(cell_draws, group_cols = c("cell_id", "state"), value_col = "z", probs = probs),
    sample_summary = summarize_draws(sample_draws, group_cols = c("sample_id", "state"), value_col = "composition", probs = probs),
    marker_summary = summarize_draws(marker_draws, group_cols = c("state", "gene"), value_col = "beta_abs", probs = probs)
  )
}

fit_bootstrap_method <- function(method, counts, marker_index, library_size, cell_metadata, fit_args) {
  if (method == "initializer") {
    args <- utils::modifyList(
      list(counts = counts, marker_index = marker_index, library_size = library_size),
      fit_args
    )
    return(do.call(fit_fibrodynmix_initializer, args))
  }
  if (method == "nb") {
    args <- utils::modifyList(
      list(counts = counts, marker_index = marker_index, library_size = library_size),
      fit_args
    )
    return(do.call(fit_fibrodynmix_nb, args))
  }
  if (method == "nb_study") {
    if (!"study_id" %in% colnames(cell_metadata)) {
      stop("`cell_metadata` must contain `study_id` when `method = 'nb_study'`.", call. = FALSE)
    }
    args <- utils::modifyList(
      list(
        counts = counts,
        marker_index = marker_index,
        library_size = library_size,
        study_id = cell_metadata$study_id,
        fit_study_effect = TRUE
      ),
      fit_args
    )
    args$study_id <- cell_metadata$study_id
    args$fit_study_effect <- TRUE
    return(do.call(fit_fibrodynmix_nb, args))
  }
  stop("Unsupported bootstrap method.", call. = FALSE)
}

cell_uncertainty_rows <- function(z, replicate, cell_id) {
  entropy <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))
  rows <- vector("list", ncol(z))
  for (k in seq_len(ncol(z))) {
    rows[[k]] <- data.frame(
      replicate = replicate,
      cell_id = cell_id,
      state = colnames(z)[k],
      z = z[, k],
      entropy = entropy,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

sample_composition_rows <- function(z, cell_metadata, sample_col, replicate) {
  samples <- as.character(cell_metadata[[sample_col]])
  rows <- list()
  row_index <- 1L
  for (sample_id in sort(unique(samples))) {
    keep <- samples == sample_id
    composition <- colMeans(z[keep, , drop = FALSE])
    for (state in names(composition)) {
      rows[[row_index]] <- data.frame(
        replicate = replicate,
        sample_id = sample_id,
        state = state,
        composition = composition[state],
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  do.call(rbind, rows)
}

marker_stability_rows <- function(beta_hat, replicate) {
  rows <- vector("list", nrow(beta_hat))
  for (k in seq_len(nrow(beta_hat))) {
    rows[[k]] <- data.frame(
      replicate = replicate,
      state = rownames(beta_hat)[k],
      gene = colnames(beta_hat),
      beta = beta_hat[k, ],
      beta_abs = abs(beta_hat[k, ]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_draws <- function(draws, group_cols, value_col, probs) {
  missing <- setdiff(c(group_cols, value_col), colnames(draws))
  if (length(missing) > 0L) {
    stop(sprintf("Draw table is missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  groups <- unique(draws[, group_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(draws))
    for (col in group_cols) {
      keep <- keep & draws[[col]] == groups[[col]][i]
    }
    values <- draws[[value_col]][keep]
    row <- groups[i, , drop = FALSE]
    row$mean <- mean(values, na.rm = TRUE)
    row$sd <- stats::sd(values, na.rm = TRUE)
    row$lower <- as.numeric(stats::quantile(values, probs[1], na.rm = TRUE, names = FALSE))
    row$upper <- as.numeric(stats::quantile(values, probs[2], na.rm = TRUE, names = FALSE))
    rows[[i]] <- row
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
