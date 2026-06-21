#' Prepare real single-cell count data for FibroDynMix
#'
#' Validates and aligns a raw gene-by-cell count matrix, cell metadata, and
#' weak marker priors before fitting FibroDynMix. This is the standard real-data
#' entry point for scripts and manuscript analyses.
#'
#' @param counts Gene-by-cell raw UMI count matrix.
#' @param cell_metadata Data frame with one row per cell.
#' @param marker_index Named list of weak prior marker genes for each state.
#'   Entries can be gene names or integer row indices.
#' @param cell_id_col Optional metadata column containing cell identifiers. If
#'   supplied, metadata is aligned to `colnames(counts)` by this column.
#' @param study_col Optional metadata column containing study/batch identifiers.
#' @param donor_col Optional metadata column containing donor/sample identifiers.
#' @param library_size_col Optional metadata column containing precomputed
#'   library sizes. If omitted, column sums of `counts` are used.
#' @param min_cells_per_gene Keep genes detected in at least this many cells.
#' @param min_counts_per_gene Keep genes with at least this total raw count.
#' @param drop_zero_library_cells Whether to remove cells with zero library size.
#' @param require_all_states Whether every state must retain at least one marker
#'   after filtering.
#'
#' @return A `FibroDynMixData` list with filtered `counts`, aligned
#'   `cell_metadata`, `marker_index`, `library_size`, optional `study_id` and
#'   `donor_id`, and audit summaries.
#' @export
prepare_fibrodynmix_data <- function(counts,
                                     cell_metadata,
                                     marker_index,
                                     cell_id_col = NULL,
                                     study_col = NULL,
                                     donor_col = NULL,
                                     library_size_col = NULL,
                                     min_cells_per_gene = 1,
                                     min_counts_per_gene = 1,
                                     drop_zero_library_cells = TRUE,
                                     require_all_states = TRUE) {
  validate_raw_count_matrix(counts)
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (anyDuplicated(rownames(counts))) {
    stop("`counts` rownames must be unique gene identifiers.", call. = FALSE)
  }
  if (anyDuplicated(colnames(counts))) {
    stop("`counts` colnames must be unique cell identifiers.", call. = FALSE)
  }
  if (!is.data.frame(cell_metadata)) {
    stop("`cell_metadata` must be a data frame.", call. = FALSE)
  }
  if (!is.list(marker_index) || is.null(names(marker_index)) || any(names(marker_index) == "")) {
    stop("`marker_index` must be a named list with non-empty state names.", call. = FALSE)
  }
  validate_nonnegative_integer_scalar(min_cells_per_gene, "min_cells_per_gene")
  validate_nonnegative_numeric_scalar(min_counts_per_gene, "min_counts_per_gene")

  cell_metadata <- align_cell_metadata(cell_metadata, colnames(counts), cell_id_col)

  if (is.null(library_size_col)) {
    library_size <- matrix_col_sums(counts)
  } else {
    validate_metadata_column(cell_metadata, library_size_col, "library_size_col")
    library_size <- as.numeric(cell_metadata[[library_size_col]])
  }
  if (length(library_size) != ncol(counts) || anyNA(library_size) || any(library_size < 0)) {
    stop("Library sizes must be non-negative and contain one value per cell.", call. = FALSE)
  }

  keep_cells <- rep(TRUE, ncol(counts))
  if (isTRUE(drop_zero_library_cells)) {
    keep_cells <- library_size > 0
  } else if (any(library_size <= 0)) {
    stop("Zero-library cells are present; use `drop_zero_library_cells = TRUE` or filter them first.", call. = FALSE)
  }
  if (!any(keep_cells)) {
    stop("No cells remain after library-size filtering.", call. = FALSE)
  }

  counts <- counts[, keep_cells, drop = FALSE]
  cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]
  library_size <- library_size[keep_cells]

  keep_genes <- matrix_row_sums(counts > 0) >= min_cells_per_gene & matrix_row_sums(counts) >= min_counts_per_gene
  if (!any(keep_genes)) {
    stop("No genes remain after gene filtering.", call. = FALSE)
  }
  counts <- counts[keep_genes, , drop = FALSE]

  marker_rows <- lapply(marker_index, resolve_marker_rows, gene_names = rownames(counts))
  marker_index_filtered <- lapply(marker_rows, function(idx) rownames(counts)[idx])
  names(marker_index_filtered) <- names(marker_index)
  marker_summary <- data.frame(
    state = names(marker_index),
    supplied_markers = vapply(marker_index, length, integer(1)),
    retained_markers = vapply(marker_index_filtered, length, integer(1)),
    stringsAsFactors = FALSE
  )
  empty_states <- marker_summary$state[marker_summary$retained_markers == 0L]
  if (isTRUE(require_all_states) && length(empty_states) > 0L) {
    stop(
      sprintf("States have no retained markers after filtering: %s", paste(empty_states, collapse = ", ")),
      call. = FALSE
    )
  }

  study_id <- NULL
  if (!is.null(study_col)) {
    validate_metadata_column(cell_metadata, study_col, "study_col")
    study_id <- as.character(cell_metadata[[study_col]])
    if (anyNA(study_id) || any(study_id == "")) {
      stop("`study_col` must not contain missing or empty values.", call. = FALSE)
    }
  }

  donor_id <- NULL
  if (!is.null(donor_col)) {
    validate_metadata_column(cell_metadata, donor_col, "donor_col")
    donor_id <- as.character(cell_metadata[[donor_col]])
    if (anyNA(donor_id) || any(donor_id == "")) {
      stop("`donor_col` must not contain missing or empty values.", call. = FALSE)
    }
  }

  rownames(cell_metadata) <- colnames(counts)
  filter_summary <- list(
    input_genes = length(keep_genes),
    retained_genes = nrow(counts),
    input_cells = length(keep_cells),
    retained_cells = ncol(counts),
    dropped_zero_library_cells = sum(!keep_cells),
    min_cells_per_gene = min_cells_per_gene,
    min_counts_per_gene = min_counts_per_gene
  )

  out <- list(
    counts = counts,
    cell_metadata = cell_metadata,
    marker_index = marker_index_filtered,
    library_size = as.numeric(library_size),
    study_id = study_id,
    donor_id = donor_id,
    marker_summary = marker_summary,
    filter_summary = filter_summary
  )
  class(out) <- c("FibroDynMixData", class(out))
  out
}

#' Fit FibroDynMix from a prepared real-data object
#'
#' Convenience wrapper around `fit_fibrodynmix_nb()` for objects returned by
#' `prepare_fibrodynmix_data()`.
#'
#' @param data A `FibroDynMixData` object.
#' @param fit_study_effect Whether to fit study-by-gene effects. If `NULL`, the
#'   value is inferred from whether `data$study_id` is available.
#' @param fit_donor_effect Whether to fit donor-by-gene effects. If `NULL`, the
#'   value is inferred when donor identifiers are available and represent a
#'   finer level than study identifiers.
#' @param ... Additional arguments passed to `fit_fibrodynmix_nb()`.
#'
#' @return A fitted FibroDynMix NB model.
#' @export
fit_fibrodynmix_prepared <- function(data, fit_study_effect = NULL, fit_donor_effect = NULL, ...) {
  if (!inherits(data, "FibroDynMixData")) {
    stop("`data` must be returned by `prepare_fibrodynmix_data()`.", call. = FALSE)
  }
  if (is.null(fit_study_effect)) {
    fit_study_effect <- !is.null(data$study_id)
  }
  if (is.null(fit_donor_effect)) {
    fit_donor_effect <- !is.null(data$donor_id) &&
      (is.null(data$study_id) || length(unique(data$donor_id)) > length(unique(data$study_id)))
  }
  fit_fibrodynmix_nb(
    counts = data$counts,
    marker_index = data$marker_index,
    library_size = data$library_size,
    study_id = data$study_id,
    donor_id = data$donor_id,
    fit_study_effect = fit_study_effect,
    fit_donor_effect = fit_donor_effect,
    ...
  )
}

validate_raw_count_matrix <- function(counts) {
  if (!is_matrix_like(counts)) {
    stop("`counts` must be a numeric matrix.", call. = FALSE)
  }
  if (!matrix_is_nonnegative_integerish(counts)) {
    stop("`counts` must contain non-negative integer-like raw counts.", call. = FALSE)
  }
}

align_cell_metadata <- function(cell_metadata, cell_names, cell_id_col) {
  if (is.null(cell_id_col)) {
    if (nrow(cell_metadata) != length(cell_names)) {
      stop("`cell_metadata` must have one row per cell when `cell_id_col` is not supplied.", call. = FALSE)
    }
    rownames(cell_metadata) <- cell_names
    return(cell_metadata)
  }

  validate_metadata_column(cell_metadata, cell_id_col, "cell_id_col")
  ids <- as.character(cell_metadata[[cell_id_col]])
  if (anyNA(ids) || any(ids == "")) {
    stop("`cell_id_col` must not contain missing or empty cell identifiers.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("`cell_id_col` must contain unique cell identifiers.", call. = FALSE)
  }
  if (!all(cell_names %in% ids)) {
    missing_ids <- cell_names[!cell_names %in% ids]
    stop(
      sprintf("`cell_metadata` is missing %d cells, including: %s", length(missing_ids), paste(utils::head(missing_ids, 5), collapse = ", ")),
      call. = FALSE
    )
  }
  cell_metadata[match(cell_names, ids), , drop = FALSE]
}

validate_metadata_column <- function(cell_metadata, column, argument_name) {
  if (length(column) != 1L || is.na(column) || !column %in% colnames(cell_metadata)) {
    stop(sprintf("`%s` must name a column in `cell_metadata`.", argument_name), call. = FALSE)
  }
}

validate_nonnegative_integer_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 0 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a non-negative integer scalar.", name), call. = FALSE)
  }
}

validate_nonnegative_numeric_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 0) {
    stop(sprintf("`%s` must be a non-negative numeric scalar.", name), call. = FALSE)
  }
}
