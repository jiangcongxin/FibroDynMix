#' Prepare a Seurat object for FibroDynMix
#'
#' Extracts raw counts and cell metadata from a Seurat object and returns a
#' `FibroDynMixData` object suitable for model fitting.
#'
#' @param object A Seurat object.
#' @param marker_index Named list of weak prior marker genes for each state.
#' @param assay Assay name. If `NULL`, the Seurat default assay is used.
#' @param layer Seurat v5 assay layer to read. Defaults to `counts`.
#' @param slot Legacy Seurat assay slot fallback. If supplied, this is used when
#'   layer extraction fails.
#' @param study_col Optional metadata column containing study/batch identifiers.
#' @param donor_col Optional metadata column containing donor/sample identifiers.
#' @param library_size_col Optional metadata column containing precomputed
#'   library sizes.
#' @param ... Additional arguments passed to `prepare_fibrodynmix_data()`.
#'
#' @return A `FibroDynMixData` object.
#' @export
prepare_fibrodynmix_seurat <- function(object,
                                       marker_index,
                                       assay = NULL,
                                       layer = "counts",
                                       slot = NULL,
                                       study_col = NULL,
                                       donor_col = NULL,
                                       library_size_col = NULL,
                                       ...) {
  require_seurat_object()
  assay <- resolve_seurat_assay(object, assay)
  counts <- extract_seurat_counts(object, assay = assay, layer = layer, slot = slot)
  metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  metadata$.fibrodynmix_cell_id <- rownames(metadata)

  prepare_fibrodynmix_data(
    counts = counts,
    cell_metadata = metadata,
    marker_index = marker_index,
    cell_id_col = ".fibrodynmix_cell_id",
    study_col = study_col,
    donor_col = donor_col,
    library_size_col = library_size_col,
    ...
  )
}

#' Fit FibroDynMix inside a Seurat workflow
#'
#' Convenience wrapper that prepares a Seurat object, fits FibroDynMix, and
#' optionally attaches state weights and diagnostics back to the Seurat object.
#'
#' @param object A Seurat object.
#' @param marker_index Named list of weak prior marker genes for each state.
#' @param assay Assay name. If `NULL`, the Seurat default assay is used.
#' @param layer Seurat v5 assay layer to read. Defaults to `counts`.
#' @param slot Legacy Seurat assay slot fallback.
#' @param study_col Optional metadata column containing study/batch identifiers.
#' @param donor_col Optional metadata column containing donor/sample identifiers.
#' @param library_size_col Optional metadata column containing precomputed
#'   library sizes.
#' @param attach Whether to attach results to the returned Seurat object.
#' @param prefix Prefix for added metadata columns.
#' @param add_reduction Whether to add a `DimReduc` containing state weights.
#' @param reduction_name Name of the Seurat reduction to create.
#' @param key Key prefix for reduction dimensions.
#' @param fit_args Additional arguments passed to `fit_fibrodynmix_prepared()`.
#' @param prepare_args Additional arguments passed to
#'   `prepare_fibrodynmix_seurat()`.
#'
#' @return A `FibroDynMixSeuratFit` list containing `seurat`, `fit`, `data`, and
#'   `assay`.
#' @export
fit_fibrodynmix_seurat <- function(object,
                                   marker_index,
                                   assay = NULL,
                                   layer = "counts",
                                   slot = NULL,
                                   study_col = NULL,
                                   donor_col = NULL,
                                   library_size_col = NULL,
                                   attach = TRUE,
                                   prefix = "fibrodynmix_",
                                   add_reduction = TRUE,
                                   reduction_name = "fibrodynmix",
                                   key = "FDM_",
                                   fit_args = list(),
                                   prepare_args = list()) {
  assay <- resolve_seurat_assay(object, assay)
  prepared <- do.call(
    prepare_fibrodynmix_seurat,
    c(
      list(
        object = object,
        marker_index = marker_index,
        assay = assay,
        layer = layer,
        slot = slot,
        study_col = study_col,
        donor_col = donor_col,
        library_size_col = library_size_col
      ),
      prepare_args
    )
  )
  fit <- do.call(fit_fibrodynmix_prepared, c(list(data = prepared), fit_args))

  seurat_out <- object
  if (isTRUE(attach)) {
    seurat_out <- add_fibrodynmix_to_seurat(
      object = seurat_out,
      fit = fit,
      cells = colnames(prepared$counts),
      assay = assay,
      prefix = prefix,
      add_reduction = add_reduction,
      reduction_name = reduction_name,
      key = key
    )
  }

  out <- list(
    seurat = seurat_out,
    fit = fit,
    data = prepared,
    assay = assay,
    layer = layer,
    reduction_name = if (isTRUE(add_reduction) && isTRUE(attach)) reduction_name else NA_character_,
    metadata_prefix = prefix
  )
  class(out) <- c("FibroDynMixSeuratFit", class(out))
  out
}

#' Run an end-to-end FibroDynMix Seurat workflow
#'
#' Filters an optional target cell type, prepares raw counts, fits FibroDynMix
#' on all or sampled cells, transfers the fitted program to all selected cells,
#' attaches state weights and confidence fields to the Seurat object, and
#' returns annotation evaluation summaries.
#'
#' @param object A Seurat object.
#' @param marker_index Optional named marker list. If `NULL`, markers are loaded
#'   with `get_fibrodynmix_markers()`.
#' @param species Species passed to `get_fibrodynmix_markers()`.
#' @param marker_context Context passed to `get_fibrodynmix_markers()`.
#' @param cell_type_col Optional metadata column containing cell type labels.
#' @param target_cell_type Optional cell type label to retain.
#' @param condition_col Optional condition column.
#' @param study_col Optional study/batch column for fitting.
#' @param donor_col Optional donor/sample column.
#' @param assay Assay name. If `NULL`, the Seurat default assay is used.
#' @param layer Seurat v5 assay layer to read. Defaults to `counts`.
#' @param slot Legacy Seurat assay slot fallback.
#' @param max_fit_cells Optional maximum number of cells used for fitting.
#'   Transfer is still run on all selected cells.
#' @param sample_col Optional metadata column used for balanced fit-cell
#'   sampling. Defaults to `donor_col`, then `condition_col`, when available.
#' @param seed Optional sampling seed.
#' @param prefix Prefix for attached metadata columns.
#' @param fit_args Additional arguments passed to `fit_fibrodynmix_prepared()`.
#' @param transfer_args Additional arguments passed to
#'   `fit_fibrodynmix_transfer()`.
#' @param prepare_args Additional arguments passed to
#'   `prepare_fibrodynmix_data()`.
#'
#' @return A `FibroDynMixSeuratWorkflow` list with the annotated Seurat object,
#'   prepared data, fit, transfer result, evaluation, and selected cells.
#' @export
run_fibrodynmix_seurat_workflow <- function(object,
                                            marker_index = NULL,
                                            species = c("human", "mouse"),
                                            marker_context = c("generic", "skin", "scar", "caf"),
                                            cell_type_col = NULL,
                                            target_cell_type = NULL,
                                            condition_col = NULL,
                                            study_col = condition_col,
                                            donor_col = NULL,
                                            assay = NULL,
                                            layer = "counts",
                                            slot = NULL,
                                            max_fit_cells = NULL,
                                            sample_col = NULL,
                                            seed = 1,
                                            prefix = "fibrodynmix_",
                                            fit_args = list(),
                                            transfer_args = list(),
                                            prepare_args = list()) {
  require_seurat_object()
  species <- match.arg(species)
  marker_context <- match.arg(marker_context)
  if (is.null(marker_index)) {
    marker_index <- get_fibrodynmix_markers(species = species, context = marker_context)
  }
  assay <- resolve_seurat_assay(object, assay)
  metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  selected_cells <- rownames(metadata)
  if (!is.null(target_cell_type)) {
    if (is.null(cell_type_col) || !cell_type_col %in% colnames(metadata)) {
      stop("`cell_type_col` must name a metadata column when `target_cell_type` is supplied.", call. = FALSE)
    }
    selected_cells <- rownames(metadata)[as.character(metadata[[cell_type_col]]) %in% target_cell_type]
    if (length(selected_cells) == 0L) {
      stop("No cells match `target_cell_type`.", call. = FALSE)
    }
    object <- object[, selected_cells]
    metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  }

  counts <- extract_seurat_counts(object, assay = assay, layer = layer, slot = slot)
  metadata$.fibrodynmix_cell_id <- rownames(metadata)
  prepare_call <- c(
    list(
      counts = counts,
      cell_metadata = metadata,
      marker_index = marker_index,
      cell_id_col = ".fibrodynmix_cell_id",
      study_col = study_col,
      donor_col = donor_col
    ),
    prepare_args
  )
  prepared_all <- do.call(prepare_fibrodynmix_data, prepare_call)

  fit_cells <- colnames(prepared_all$counts)
  if (!is.null(max_fit_cells) && length(fit_cells) > max_fit_cells) {
    assert_positive_integer(max_fit_cells, "max_fit_cells")
    if (!is.null(seed)) {
      set.seed(seed)
    }
    sample_col <- resolve_workflow_sample_col(prepared_all$cell_metadata, sample_col, donor_col, condition_col)
    fit_cells <- sample_workflow_cells(prepared_all$cell_metadata, sample_col, max_fit_cells)
  }

  prepared_fit <- subset_fibrodynmix_prepared(prepared_all, fit_cells)
  fit <- do.call(fit_fibrodynmix_prepared, c(list(data = prepared_fit), fit_args))

  warm_start <- score_marker_baseline(
    counts = prepared_all$counts,
    marker_index = prepared_all$marker_index,
    library_size = prepared_all$library_size
  )$z_pred
  transfer_call <- utils::modifyList(
    list(
      counts = prepared_all$counts,
      fit = fit,
      library_size = prepared_all$library_size,
      z_init = warm_start[, rownames(fit$beta_hat), drop = FALSE],
      chunk_size = 500,
      maxit_z = 75
    ),
    transfer_args
  )
  transfer <- do.call(fit_fibrodynmix_transfer, transfer_call)

  object_out <- add_fibrodynmix_to_seurat(
    object = object,
    fit = transfer,
    cells = rownames(transfer$z_hat),
    assay = assay,
    prefix = prefix,
    add_reduction = TRUE,
    reduction_name = paste0(prefix, "transfer"),
    key = "FDMX_"
  )

  weights <- data.frame(
    cell_id = rownames(transfer$z_hat),
    transfer$z_hat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  eval_metadata <- as.data.frame(object_out[[]], stringsAsFactors = FALSE)
  eval_metadata$cell_id <- rownames(eval_metadata)
  eval_expression <- tryCatch(
    extract_seurat_counts(object_out, assay = assay, layer = "data", slot = "data"),
    error = function(e) extract_seurat_counts(object_out, assay = assay, layer = layer, slot = slot)
  )
  evaluation <- evaluate_fibrodynmix_annotation(
    weights,
    state_cols = colnames(transfer$z_hat),
    metadata = eval_metadata,
    metadata_cell_col = "cell_id",
    cluster_col = if ("seurat_clusters" %in% colnames(eval_metadata)) "seurat_clusters" else NULL,
    condition_col = condition_col,
    expression = eval_expression,
    marker_index = marker_index
  )

  out <- list(
    seurat = object_out,
    data = prepared_all,
    fit_data = prepared_fit,
    fit = fit,
    transfer = transfer,
    evaluation = evaluation,
    marker_index = marker_index,
    selected_cells = colnames(prepared_all$counts),
    fit_cells = fit_cells,
    assay = assay
  )
  class(out) <- c("FibroDynMixSeuratWorkflow", class(out))
  out
}

#' Attach FibroDynMix results to a Seurat object
#'
#' Adds state weights, dominant state, entropy, and FPI to Seurat metadata and
#' optionally stores the state-weight matrix as a Seurat dimensional reduction.
#'
#' @param object A Seurat object.
#' @param fit A fitted FibroDynMix model containing `z_hat`.
#' @param cells Optional cell identifiers corresponding to rows of `fit$z_hat`.
#' @param assay Assay name for the dimensional reduction. If `NULL`, the Seurat
#'   default assay is used.
#' @param prefix Prefix for added metadata columns.
#' @param add_reduction Whether to add a `DimReduc` containing state weights.
#' @param reduction_name Name of the Seurat reduction to create.
#' @param key Key prefix for reduction dimensions.
#'
#' @return The input Seurat object with FibroDynMix results attached.
#' @export
add_fibrodynmix_to_seurat <- function(object,
                                      fit,
                                      cells = NULL,
                                      assay = NULL,
                                      prefix = "fibrodynmix_",
                                      add_reduction = TRUE,
                                      reduction_name = "fibrodynmix",
                                      key = "FDM_") {
  require_seurat_object()
  if (is.null(fit$z_hat)) {
    stop("`fit` must contain a `z_hat` matrix.", call. = FALSE)
  }
  assay <- resolve_seurat_assay(object, assay)
  z <- as.matrix(fit$z_hat)
  if (is.null(colnames(z))) {
    colnames(z) <- sprintf("state_%d", seq_len(ncol(z)))
  }
  if (is.null(cells)) {
    cells <- rownames(z)
  }
  if (is.null(cells) || length(cells) != nrow(z)) {
    stop("`cells` must contain one identifier per row of `fit$z_hat`.", call. = FALSE)
  }
  if (anyDuplicated(cells)) {
    stop("`cells` must be unique.", call. = FALSE)
  }
  object_cells <- colnames(object)
  if (!all(cells %in% object_cells)) {
    missing_cells <- cells[!cells %in% object_cells]
    stop(
      sprintf("Seurat object is missing %d fitted cells, including: %s", length(missing_cells), paste(utils::head(missing_cells, 5), collapse = ", ")),
      call. = FALSE
    )
  }

  state_names <- colnames(z)
  state_columns <- paste0(prefix, "z_", sanitize_seurat_column_names(state_names))
  metadata <- data.frame(
    matrix(NA_real_, nrow = length(object_cells), ncol = length(state_names)),
    row.names = object_cells,
    check.names = FALSE
  )
  colnames(metadata) <- state_columns
  metadata[cells, state_columns] <- z

  dominant <- rep(NA_character_, length(object_cells))
  names(dominant) <- object_cells
  dominant[cells] <- state_names[max.col(z, ties.method = "first")]
  entropy <- rep(NA_real_, length(object_cells))
  names(entropy) <- object_cells
  entropy[cells] <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))
  max_weight <- rep(NA_real_, length(object_cells))
  names(max_weight) <- object_cells
  max_weight[cells] <- apply(z, 1L, max)
  normalized_entropy <- rep(NA_real_, length(object_cells))
  names(normalized_entropy) <- object_cells
  normalized_entropy[cells] <- entropy[cells] / log(ncol(z))
  confidence_class <- rep(NA_character_, length(object_cells))
  names(confidence_class) <- object_cells
  confidence_class[cells] <- as.character(classify_fibrodynmix_confidence(max_weight[cells]))

  fpi <- rep(NA_real_, length(object_cells))
  names(fpi) <- object_cells
  fpi_values <- compute_fpi(z)$fpi
  fpi[cells] <- fpi_values

  metadata[[paste0(prefix, "dominant_state")]] <- dominant
  metadata[[paste0(prefix, "entropy")]] <- entropy
  metadata[[paste0(prefix, "normalized_entropy")]] <- normalized_entropy
  metadata[[paste0(prefix, "max_weight")]] <- max_weight
  metadata[[paste0(prefix, "confidence_class")]] <- confidence_class
  metadata[[paste0(prefix, "fpi")]] <- fpi
  object <- SeuratObject::AddMetaData(object, metadata = metadata)

  if (isTRUE(add_reduction)) {
    embeddings <- matrix(NA_real_, nrow = length(object_cells), ncol = ncol(z))
    rownames(embeddings) <- object_cells
    colnames(embeddings) <- paste0(key, seq_len(ncol(z)))
    embeddings[cells, ] <- z
    reduction <- SeuratObject::CreateDimReducObject(
      embeddings = embeddings,
      assay = assay,
      key = key,
      misc = list(
        state_names = state_names,
        metadata_columns = state_columns,
        dominant_state_column = paste0(prefix, "dominant_state"),
        entropy_column = paste0(prefix, "entropy"),
        normalized_entropy_column = paste0(prefix, "normalized_entropy"),
        max_weight_column = paste0(prefix, "max_weight"),
        confidence_class_column = paste0(prefix, "confidence_class"),
        fpi_column = paste0(prefix, "fpi")
      )
    )
    object[[reduction_name]] <- reduction
  }

  object
}

require_seurat_object <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package `SeuratObject` is required for FibroDynMix Seurat integration.", call. = FALSE)
  }
  invisible(TRUE)
}

resolve_seurat_assay <- function(object, assay) {
  require_seurat_object()
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }
  if (length(assay) != 1L || is.na(assay) || !assay %in% names(object@assays)) {
    stop("`assay` must name an assay in the Seurat object.", call. = FALSE)
  }
  assay
}

extract_seurat_counts <- function(object, assay, layer = "counts", slot = NULL) {
  require_seurat_object()
  counts <- tryCatch(
    suppressWarnings(SeuratObject::LayerData(object, assay = assay, layer = layer)),
    error = function(e) NULL
  )
  if (!is.null(counts) && (nrow(counts) == 0L || ncol(counts) == 0L)) {
    counts <- NULL
  }
  if (is.null(counts)) {
    fallback_slot <- if (is.null(slot)) layer else slot
    counts <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, slot = fallback_slot),
      error = function(e) NULL
    )
  }
  if (is.null(counts)) {
    stop("Could not extract raw counts from the Seurat object.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Extracted Seurat counts must have gene and cell names.", call. = FALSE)
  }
  counts
}

sanitize_seurat_column_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "state")
}

resolve_workflow_sample_col <- function(metadata, sample_col, donor_col, condition_col) {
  candidates <- c(sample_col, donor_col, condition_col)
  if (length(candidates) == 0L) {
    return(NULL)
  }
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% colnames(metadata)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

sample_workflow_cells <- function(metadata, sample_col, max_fit_cells) {
  cell_id <- rownames(metadata)
  if (is.null(sample_col)) {
    return(sample(cell_id, max_fit_cells))
  }
  strata <- split(cell_id, metadata[[sample_col]])
  n_strata <- length(strata)
  base_n <- floor(max_fit_cells / n_strata)
  remainder <- max_fit_cells - base_n * n_strata
  sampled <- unlist(lapply(seq_along(strata), function(i) {
    target_n <- base_n + as.integer(i <= remainder)
    sample(strata[[i]], min(length(strata[[i]]), target_n))
  }), use.names = FALSE)
  if (length(sampled) > max_fit_cells) {
    sampled <- sample(sampled, max_fit_cells)
  }
  sampled
}

subset_fibrodynmix_prepared <- function(data, cells) {
  keep <- match(cells, colnames(data$counts))
  if (anyNA(keep)) {
    stop("Internal error: fit cells are not present in prepared data.", call. = FALSE)
  }
  out <- data
  out$counts <- data$counts[, keep, drop = FALSE]
  out$cell_metadata <- data$cell_metadata[keep, , drop = FALSE]
  out$library_size <- data$library_size[keep]
  out$study_id <- if (is.null(data$study_id)) NULL else data$study_id[keep]
  out$donor_id <- if (is.null(data$donor_id)) NULL else data$donor_id[keep]
  class(out) <- class(data)
  out
}
