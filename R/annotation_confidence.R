#' Evaluate FibroDynMix annotation confidence
#'
#' Summarizes cell-level state-weight confidence, state support, optional
#' cluster-state agreement, optional condition-state composition, and optional
#' marker support from expression values.
#'
#' @param cell_weights Data frame or matrix with one row per cell and state
#'   weight columns.
#' @param state_cols Optional state-weight columns. If omitted, numeric columns
#'   excluding common metadata columns are used.
#' @param cell_col Column containing cell identifiers. If absent, row names or
#'   row numbers are used.
#' @param metadata Optional cell metadata data frame.
#' @param metadata_cell_col Metadata column containing cell identifiers. If
#'   `NULL`, row names are used when present.
#' @param cluster_col Optional Seurat or other cluster column used to evaluate
#'   cluster-state agreement.
#' @param condition_col Optional condition column used for state composition.
#' @param expression Optional gene-by-cell expression matrix used for marker
#'   support. Raw, log-normalized, or scaled values are accepted; values are
#'   compared within the supplied matrix.
#' @param marker_index Optional named list of marker genes per state.
#' @param high_threshold Maximum state weight required for high confidence.
#' @param moderate_threshold Maximum state weight required for moderate
#'   confidence.
#' @param min_state_cells Minimum dominant-cell count for a state to be
#'   considered supported.
#' @param min_high_confidence_fraction Minimum fraction of high-confidence cells
#'   for a state to be considered supported.
#' @param min_marker_log2_ratio Minimum marker-score log2 ratio of own-state
#'   cells versus other cells. Only used when marker support is supplied.
#'
#' @return A `FibroDynMixAnnotationEvaluation` list with `cell_summary`,
#'   `state_summary`, and optional `cluster_agreement`,
#'   `condition_composition`, and `marker_support` tables.
#' @export
evaluate_fibrodynmix_annotation <- function(cell_weights,
                                            state_cols = NULL,
                                            cell_col = "cell_id",
                                            metadata = NULL,
                                            metadata_cell_col = NULL,
                                            cluster_col = NULL,
                                            condition_col = NULL,
                                            expression = NULL,
                                            marker_index = NULL,
                                            high_threshold = 0.7,
                                            moderate_threshold = 0.4,
                                            min_state_cells = 20,
                                            min_high_confidence_fraction = 0.2,
                                            min_marker_log2_ratio = 0.25) {
  cell_weights <- as.data.frame(cell_weights, stringsAsFactors = FALSE, check.names = FALSE)
  state_cols <- infer_state_columns(cell_weights, state_cols)
  thresholds <- validate_confidence_thresholds(high_threshold, moderate_threshold)
  validate_support_thresholds(
    min_state_cells = min_state_cells,
    min_high_confidence_fraction = min_high_confidence_fraction,
    min_marker_log2_ratio = min_marker_log2_ratio
  )

  cell_id <- resolve_cell_ids(cell_weights, cell_col)
  z <- as.matrix(cell_weights[, state_cols, drop = FALSE])
  storage.mode(z) <- "double"
  validate_state_weight_matrix(z)
  rownames(z) <- cell_id

  dominant_state <- state_cols[max.col(z, ties.method = "first")]
  max_state_weight <- apply(z, 1L, max)
  entropy <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))
  normalized_entropy <- entropy / log(ncol(z))
  confidence_class <- classify_fibrodynmix_confidence(
    max_state_weight,
    high_threshold = thresholds$high,
    moderate_threshold = thresholds$moderate
  )

  cell_summary <- data.frame(
    cell_id = cell_id,
    dominant_state = dominant_state,
    max_state_weight = max_state_weight,
    entropy = entropy,
    normalized_entropy = normalized_entropy,
    confidence_class = confidence_class,
    stringsAsFactors = FALSE
  )

  metadata_aligned <- NULL
  if (!is.null(metadata)) {
    metadata_aligned <- align_evaluation_metadata(metadata, cell_id, metadata_cell_col)
    for (column in setdiff(colnames(metadata_aligned), colnames(cell_summary))) {
      cell_summary[[column]] <- metadata_aligned[[column]]
    }
  }

  state_summary <- summarize_fibrodynmix_states(cell_summary, state_cols)

  cluster_agreement <- NULL
  if (!is.null(cluster_col)) {
    check_evaluation_columns(cell_summary, cluster_col, "cell_summary")
    cluster_agreement <- summarize_cluster_state_agreement(cell_summary, cluster_col)
  }

  condition_composition <- NULL
  if (!is.null(condition_col)) {
    check_evaluation_columns(cell_summary, condition_col, "cell_summary")
    condition_composition <- summarize_condition_state_composition(cell_summary, z, condition_col)
  }

  marker_support <- NULL
  if (!is.null(expression) || !is.null(marker_index)) {
    if (is.null(expression) || is.null(marker_index)) {
      stop("`expression` and `marker_index` must be supplied together for marker support.", call. = FALSE)
    }
    marker_support <- summarize_marker_support(expression, marker_index, cell_summary)
  }
  state_summary <- add_state_support_labels(
    state_summary,
    marker_support = marker_support,
    min_state_cells = min_state_cells,
    min_high_confidence_fraction = min_high_confidence_fraction,
    min_marker_log2_ratio = min_marker_log2_ratio
  )

  out <- list(
    cell_summary = cell_summary,
    state_summary = state_summary,
    cluster_agreement = cluster_agreement,
    condition_composition = condition_composition,
    marker_support = marker_support,
    thresholds = thresholds,
    support_thresholds = list(
      min_state_cells = as.integer(min_state_cells),
      min_high_confidence_fraction = as.numeric(min_high_confidence_fraction),
      min_marker_log2_ratio = as.numeric(min_marker_log2_ratio)
    )
  )
  class(out) <- c("FibroDynMixAnnotationEvaluation", class(out))
  out
}

validate_support_thresholds <- function(min_state_cells,
                                        min_high_confidence_fraction,
                                        min_marker_log2_ratio) {
  if (length(min_state_cells) != 1L || is.na(min_state_cells) || min_state_cells < 0 || min_state_cells != as.integer(min_state_cells)) {
    stop("`min_state_cells` must be a non-negative integer.", call. = FALSE)
  }
  if (length(min_high_confidence_fraction) != 1L ||
      is.na(min_high_confidence_fraction) ||
      min_high_confidence_fraction < 0 ||
      min_high_confidence_fraction > 1) {
    stop("`min_high_confidence_fraction` must be in [0, 1].", call. = FALSE)
  }
  if (length(min_marker_log2_ratio) != 1L ||
      is.na(min_marker_log2_ratio) ||
      !is.finite(min_marker_log2_ratio)) {
    stop("`min_marker_log2_ratio` must be a finite numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

classify_fibrodynmix_confidence <- function(max_state_weight,
                                            high_threshold = 0.7,
                                            moderate_threshold = 0.4) {
  thresholds <- validate_confidence_thresholds(high_threshold, moderate_threshold)
  out <- ifelse(
    max_state_weight >= thresholds$high,
    "high",
    ifelse(max_state_weight >= thresholds$moderate, "moderate", "low")
  )
  factor(out, levels = c("high", "moderate", "low"))
}

validate_confidence_thresholds <- function(high_threshold, moderate_threshold) {
  if (length(high_threshold) != 1L || is.na(high_threshold) || high_threshold <= 0 || high_threshold > 1) {
    stop("`high_threshold` must be in (0, 1].", call. = FALSE)
  }
  if (length(moderate_threshold) != 1L || is.na(moderate_threshold) || moderate_threshold < 0 || moderate_threshold >= high_threshold) {
    stop("`moderate_threshold` must be non-negative and less than `high_threshold`.", call. = FALSE)
  }
  list(high = as.numeric(high_threshold), moderate = as.numeric(moderate_threshold))
}

resolve_cell_ids <- function(cell_weights, cell_col) {
  if (!is.null(cell_col) && cell_col %in% colnames(cell_weights)) {
    cell_id <- as.character(cell_weights[[cell_col]])
  } else if (!is.null(rownames(cell_weights)) && !identical(rownames(cell_weights), as.character(seq_len(nrow(cell_weights))))) {
    cell_id <- rownames(cell_weights)
  } else {
    cell_id <- sprintf("cell_%d", seq_len(nrow(cell_weights)))
  }
  if (anyNA(cell_id) || any(cell_id == "") || anyDuplicated(cell_id)) {
    stop("Cell identifiers must be unique and non-missing.", call. = FALSE)
  }
  cell_id
}

validate_state_weight_matrix <- function(z) {
  if (ncol(z) < 2L || nrow(z) == 0L || anyNA(z) || any(!is.finite(z)) || any(z < 0)) {
    stop("State-weight columns must be non-negative finite numeric values with at least two states.", call. = FALSE)
  }
  row_sums <- rowSums(z)
  if (any(row_sums <= 0)) {
    stop("Each cell must have positive total state weight.", call. = FALSE)
  }
  invisible(TRUE)
}

align_evaluation_metadata <- function(metadata, cell_id, metadata_cell_col) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(metadata_cell_col)) {
    check_evaluation_columns(metadata, metadata_cell_col, "metadata")
    ids <- as.character(metadata[[metadata_cell_col]])
  } else if (!is.null(rownames(metadata)) && !identical(rownames(metadata), as.character(seq_len(nrow(metadata))))) {
    ids <- rownames(metadata)
    metadata$.metadata_cell_id <- ids
  } else if ("cell_id" %in% colnames(metadata)) {
    ids <- as.character(metadata$cell_id)
  } else {
    stop("`metadata_cell_col` is required when metadata row names do not contain cell identifiers.", call. = FALSE)
  }
  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop("Metadata cell identifiers must be unique and non-missing.", call. = FALSE)
  }
  if (!all(cell_id %in% ids)) {
    stop("`metadata` is missing cells present in `cell_weights`.", call. = FALSE)
  }
  metadata[match(cell_id, ids), , drop = FALSE]
}

summarize_fibrodynmix_states <- function(cell_summary, state_cols) {
  rows <- lapply(state_cols, function(state) {
    idx <- cell_summary$dominant_state == state
    if (!any(idx)) {
      return(data.frame(
        state = state,
        n_cells = 0L,
        fraction = 0,
        median_max_weight = NA_real_,
        mean_max_weight = NA_real_,
        high_confidence_fraction = NA_real_,
        median_normalized_entropy = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      state = state,
      n_cells = sum(idx),
      fraction = mean(idx),
      median_max_weight = stats::median(cell_summary$max_state_weight[idx]),
      mean_max_weight = mean(cell_summary$max_state_weight[idx]),
      high_confidence_fraction = mean(cell_summary$confidence_class[idx] == "high"),
      median_normalized_entropy = stats::median(cell_summary$normalized_entropy[idx]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$n_cells, out$state), , drop = FALSE]
}

summarize_cluster_state_agreement <- function(cell_summary, cluster_col) {
  tab <- as.data.frame(table(
    cluster = cell_summary[[cluster_col]],
    state = cell_summary$dominant_state
  ), stringsAsFactors = FALSE)
  colnames(tab)[3L] <- "n_cells"
  totals <- stats::aggregate(n_cells ~ cluster, tab, sum)
  tab <- merge(tab, totals, by = "cluster", suffixes = c("", "_cluster_total"))
  tab$fraction_within_cluster <- ifelse(tab$n_cells_cluster_total > 0, tab$n_cells / tab$n_cells_cluster_total, NA_real_)
  tab <- tab[order(tab$cluster, -tab$fraction_within_cluster), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

summarize_condition_state_composition <- function(cell_summary, z, condition_col) {
  conditions <- unique(as.character(cell_summary[[condition_col]]))
  rows <- lapply(conditions, function(condition) {
    idx <- as.character(cell_summary[[condition_col]]) == condition
    data.frame(
      condition = condition,
      state = colnames(z),
      composition = colMeans(z[idx, , drop = FALSE]),
      dominant_fraction = as.numeric(prop.table(table(factor(cell_summary$dominant_state[idx], levels = colnames(z))))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_marker_support <- function(expression, marker_index, cell_summary) {
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("`expression` must have gene rownames and cell colnames.", call. = FALSE)
  }
  if (!is.list(marker_index) || is.null(names(marker_index))) {
    stop("`marker_index` must be a named list.", call. = FALSE)
  }
  cells <- intersect(cell_summary$cell_id, colnames(expression))
  if (length(cells) == 0L) {
    stop("`expression` column names do not overlap evaluated cells.", call. = FALSE)
  }
  cell_summary <- cell_summary[match(cells, cell_summary$cell_id), , drop = FALSE]
  expression <- expression[, cells, drop = FALSE]
  rows <- lapply(names(marker_index), function(state) {
    genes <- unique(as.character(marker_index[[state]]))
    genes <- genes[genes %in% rownames(expression)]
    if (length(genes) == 0L) {
      return(data.frame(
        state = state,
        n_marker_genes = 0L,
        retained_marker_genes = "",
        mean_score_own_state = NA_real_,
        mean_score_other_states = NA_real_,
        log2_ratio_own_vs_other = NA_real_,
        pct_own_cells_score_positive = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    score <- marker_col_means(expression, genes)
    own <- score[cell_summary$dominant_state == state]
    other <- score[cell_summary$dominant_state != state]
    data.frame(
      state = state,
      n_marker_genes = length(genes),
      retained_marker_genes = paste(genes, collapse = ","),
      mean_score_own_state = if (length(own) == 0L) NA_real_ else mean(own),
      mean_score_other_states = if (length(other) == 0L) NA_real_ else mean(other),
      log2_ratio_own_vs_other = if (length(own) == 0L || length(other) == 0L) NA_real_ else log2((mean(own) + 0.01) / (mean(other) + 0.01)),
      pct_own_cells_score_positive = if (length(own) == 0L) NA_real_ else mean(own > 0),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$log2_ratio_own_vs_other, out$state, na.last = TRUE), , drop = FALSE]
}

add_state_support_labels <- function(state_summary,
                                     marker_support = NULL,
                                     min_state_cells = 20,
                                     min_high_confidence_fraction = 0.2,
                                     min_marker_log2_ratio = 0.25) {
  out <- state_summary
  out$marker_log2_ratio_own_vs_other <- NA_real_
  out$n_marker_genes <- NA_integer_
  out$marker_support_status <- if (is.null(marker_support)) "not_evaluated" else "missing"

  if (!is.null(marker_support) && nrow(marker_support) > 0L) {
    marker_support <- as.data.frame(marker_support, stringsAsFactors = FALSE)
    idx <- match(out$state, marker_support$state)
    has_marker_row <- !is.na(idx)
    out$marker_log2_ratio_own_vs_other[has_marker_row] <- marker_support$log2_ratio_own_vs_other[idx[has_marker_row]]
    out$n_marker_genes[has_marker_row] <- marker_support$n_marker_genes[idx[has_marker_row]]
    out$marker_support_status[has_marker_row] <- ifelse(
      is.na(out$n_marker_genes[has_marker_row]) | out$n_marker_genes[has_marker_row] == 0L,
      "no_markers_retained",
      ifelse(
        !is.na(out$marker_log2_ratio_own_vs_other[has_marker_row]) &
          out$marker_log2_ratio_own_vs_other[has_marker_row] >= min_marker_log2_ratio,
        "supported",
        "weak"
      )
    )
  }

  cell_count_ok <- out$n_cells >= min_state_cells
  confidence_ok <- !is.na(out$high_confidence_fraction) &
    out$high_confidence_fraction >= min_high_confidence_fraction
  marker_ok <- out$marker_support_status %in% c("supported", "not_evaluated")
  out$formal_ready <- cell_count_ok & confidence_ok & marker_ok
  out$support_label <- ifelse(
    out$n_cells == 0L,
    "unsupported",
    ifelse(out$formal_ready, "supported", "exploratory")
  )
  out$support_reason <- vapply(seq_len(nrow(out)), function(i) {
    reasons <- character()
    if (out$n_cells[[i]] < min_state_cells) {
      reasons <- c(reasons, "low_cell_count")
    }
    if (is.na(out$high_confidence_fraction[[i]]) ||
        out$high_confidence_fraction[[i]] < min_high_confidence_fraction) {
      reasons <- c(reasons, "low_high_confidence_fraction")
    }
    if (out$marker_support_status[[i]] == "weak") {
      reasons <- c(reasons, "weak_marker_enrichment")
    } else if (out$marker_support_status[[i]] == "no_markers_retained") {
      reasons <- c(reasons, "no_marker_genes_retained")
    } else if (out$marker_support_status[[i]] == "missing") {
      reasons <- c(reasons, "marker_support_missing")
    }
    if (length(reasons) == 0L) {
      return("passes_thresholds")
    }
    paste(reasons, collapse = ";")
  }, character(1L))
  out
}

marker_col_means <- function(expression, genes) {
  subset <- expression[genes, , drop = FALSE]
  if (inherits(subset, "sparseMatrix") || inherits(subset, "Matrix")) {
    return(as.numeric(Matrix::colMeans(subset)))
  }
  colMeans(as.matrix(subset))
}

check_evaluation_columns <- function(data, columns, data_name) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0L) {
    stop(sprintf("`%s` is missing required column(s): %s.", data_name, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}
