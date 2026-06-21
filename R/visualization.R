#' Plot sample-level fibroblast state composition
#'
#' Draws a stacked composition bar plot from a long state-composition table.
#'
#' @param composition Data frame containing one row per group and state.
#' @param group_col Column used on the x-axis, such as `dataset_id`,
#'   `sample_id`, `donor_id`, or `condition`.
#' @param state_col Column containing fibroblast state names.
#' @param value_col Column containing state fractions or proportions.
#' @param facet_col Optional column used for faceting.
#'
#' @return A `ggplot` object.
#' @export
plot_state_composition <- function(composition,
                                   group_col = "dataset_id",
                                   state_col = "state",
                                   value_col = "composition",
                                   facet_col = NULL) {
  require_ggplot2()
  composition <- as.data.frame(composition, stringsAsFactors = FALSE)
  check_plot_columns(composition, c(group_col, state_col, value_col), "composition")
  if (!is.null(facet_col)) {
    check_plot_columns(composition, facet_col, "composition")
  }

  composition[[group_col]] <- factor(composition[[group_col]], levels = unique(composition[[group_col]]))
  composition[[state_col]] <- factor(composition[[state_col]], levels = unique(composition[[state_col]]))

  p <- ggplot2::ggplot(
    composition,
    aes_columns(x = group_col, y = value_col, fill = state_col)
  ) +
    ggplot2::geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
    ggplot2::scale_y_continuous(labels = ggplot_percent_labels()) +
    ggplot2::labs(x = NULL, y = "State composition", fill = "State") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_x")
  }
  p
}

#' Plot cell-level fibroblast state mixture weights
#'
#' Draws a cell-by-state heatmap from a table containing state weight columns.
#'
#' @param cell_weights Data frame containing one row per cell and state-weight
#'   columns.
#' @param state_cols Optional character vector naming the state-weight columns.
#'   If omitted, numeric non-metadata columns are used.
#' @param cell_col Column containing cell identifiers. If absent, row numbers
#'   are used.
#' @param order_by Cell ordering rule: dominant state, entropy, or original
#'   order.
#' @param max_cells Maximum number of cells to display.
#'
#' @return A `ggplot` object.
#' @export
plot_cell_state_heatmap <- function(cell_weights,
                                    state_cols = NULL,
                                    cell_col = "cell_id",
                                    order_by = c("dominant_state", "entropy", "none"),
                                    max_cells = 200) {
  require_ggplot2()
  order_by <- match.arg(order_by)
  cell_weights <- as.data.frame(cell_weights, stringsAsFactors = FALSE)
  state_cols <- infer_state_columns(cell_weights, state_cols)
  assert_plot_positive_integer(max_cells, "max_cells")

  z <- as.matrix(cell_weights[, state_cols, drop = FALSE])
  storage.mode(z) <- "double"
  if (anyNA(z) || any(z < 0)) {
    stop("State-weight columns must contain non-negative numeric values without NA.", call. = FALSE)
  }

  cell_id <- if (cell_col %in% colnames(cell_weights)) {
    as.character(cell_weights[[cell_col]])
  } else {
    sprintf("cell_%d", seq_len(nrow(cell_weights)))
  }
  dominant <- state_cols[max.col(z, ties.method = "first")]
  entropy <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))

  ord <- switch(
    order_by,
    dominant_state = order(dominant, -apply(z, 1, max), entropy, cell_id),
    entropy = order(entropy, cell_id),
    none = seq_len(nrow(cell_weights))
  )
  ord <- utils::head(ord, max_cells)

  heatmap_data <- long_state_weights(
    cell_id = cell_id[ord],
    z = z[ord, , drop = FALSE],
    state_cols = state_cols
  )
  heatmap_data$cell_id <- factor(heatmap_data$cell_id, levels = unique(heatmap_data$cell_id))
  heatmap_data$state <- factor(heatmap_data$state, levels = rev(state_cols))

  ggplot2::ggplot(
    heatmap_data,
    aes_columns(x = "cell_id", y = "state", fill = "state_weight")
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7fbff", "#9ecae1", "#3182bd", "#08519c"),
      limits = c(0, 1)
    ) +
    ggplot2::labs(x = "Cells", y = NULL, fill = "Weight") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

#' Plot cross-cohort transfer diagnostics
#'
#' Draws a horizontal diagnostic plot for leave-dataset-out or leave-donor-out
#' transfer analyses.
#'
#' @param transfer Data frame containing transfer diagnostics.
#' @param x_col Column naming held-out datasets or donors. If omitted, a common
#'   identifier column is inferred.
#' @param y_col Numeric diagnostic column to plot.
#' @param threshold Optional reference threshold.
#'
#' @return A `ggplot` object.
#' @export
plot_transfer_diagnostics <- function(transfer,
                                      x_col = NULL,
                                      y_col = "transfer_z_convergence_rate",
                                      threshold = 0.9) {
  require_ggplot2()
  transfer <- as.data.frame(transfer, stringsAsFactors = FALSE)
  if (is.null(x_col)) {
    x_col <- first_existing_column(
      transfer,
      c("heldout_dataset_id", "heldout_donor_id", "dataset_id", "donor_id", "sample_id")
    )
  }
  check_plot_columns(transfer, c(x_col, y_col), "transfer")
  transfer[[x_col]] <- stats::reorder(factor(transfer[[x_col]]), transfer[[y_col]])

  p <- ggplot2::ggplot(transfer, aes_columns(x = x_col, y = y_col)) +
    ggplot2::geom_col(width = 0.72, fill = "#4c78a8") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_col) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "#b2182b")
  }
  p
}

#' Plot fibroblast state transition flow
#'
#' Draws a source-state by target-state heatmap from a long transition-flow
#' table.
#'
#' @param flow Data frame containing transition-flow rows.
#' @param source_col Source-state column.
#' @param target_col Target-state column.
#' @param value_col Flow or probability column.
#'
#' @return A `ggplot` object.
#' @export
plot_transition_flow <- function(flow,
                                 source_col = "source_state",
                                 target_col = "target_state",
                                 value_col = "flow") {
  require_ggplot2()
  flow <- as.data.frame(flow, stringsAsFactors = FALSE)
  check_plot_columns(flow, c(source_col, target_col, value_col), "flow")
  flow[[source_col]] <- factor(flow[[source_col]], levels = unique(flow[[source_col]]))
  flow[[target_col]] <- factor(flow[[target_col]], levels = unique(flow[[target_col]]))

  ggplot2::ggplot(flow, aes_columns(x = target_col, y = source_col, fill = value_col)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradientn(colours = c("#f7fbff", "#9ecae1", "#3182bd", "#08519c")) +
    ggplot2::labs(x = "Target state", y = "Source state", fill = "Flow") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Plot Fibroblast Plasticity Index distributions
#'
#' Draws cell-level or sample-level FPI distributions by group.
#'
#' @param fpi_data Data frame containing FPI or entropy values.
#' @param group_col Grouping column, such as `condition`, `dataset_id`, or
#'   `donor_id`.
#' @param fpi_col Numeric FPI column. If absent and `entropy` exists, entropy is
#'   used.
#' @param point_alpha Alpha for overlaid points.
#'
#' @return A `ggplot` object.
#' @export
plot_fpi_distribution <- function(fpi_data,
                                  group_col = "condition",
                                  fpi_col = "fpi",
                                  point_alpha = 0.35) {
  require_ggplot2()
  fpi_data <- as.data.frame(fpi_data, stringsAsFactors = FALSE)
  if (!fpi_col %in% colnames(fpi_data) && "entropy" %in% colnames(fpi_data)) {
    fpi_col <- "entropy"
  }
  check_plot_columns(fpi_data, c(group_col, fpi_col), "fpi_data")
  fpi_data[[group_col]] <- factor(fpi_data[[group_col]], levels = unique(fpi_data[[group_col]]))

  ggplot2::ggplot(fpi_data, aes_columns(x = group_col, y = fpi_col, fill = group_col)) +
    ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.65) +
    ggplot2::geom_jitter(width = 0.14, height = 0, alpha = point_alpha, size = 0.8) +
    ggplot2::labs(x = NULL, y = fpi_col, fill = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Plot benchmark method rankings
#'
#' Draws benchmark metric comparisons across stress scenarios and methods.
#'
#' @param rankings Data frame containing benchmark rankings or summaries.
#' @param scenario_col Scenario column. If omitted, a common scenario column is
#'   inferred.
#' @param method_col Method column.
#' @param metric_col Numeric metric column.
#'
#' @return A `ggplot` object.
#' @export
plot_benchmark_rankings <- function(rankings,
                                    scenario_col = NULL,
                                    method_col = "method",
                                    metric_col = "rmse_mean") {
  require_ggplot2()
  rankings <- as.data.frame(rankings, stringsAsFactors = FALSE)
  if (is.null(scenario_col)) {
    scenario_col <- first_existing_column(rankings, c("stress_mode", "scenario", "benchmark", "setting"))
  }
  check_plot_columns(rankings, c(scenario_col, method_col, metric_col), "rankings")
  rankings[[scenario_col]] <- factor(rankings[[scenario_col]], levels = unique(rankings[[scenario_col]]))

  ggplot2::ggplot(
    rankings,
    aes_columns(x = scenario_col, y = metric_col, colour = method_col, group = method_col)
  ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::labs(x = NULL, y = metric_col, colour = "Method") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Plot fibroblast purity QC
#'
#' Draws dataset-level fibroblast purity margins from the bioinformatics
#' validation layer.
#'
#' @param purity Data frame containing purity QC metrics.
#' @param x_col Dataset or sample identifier column.
#' @param y_col Numeric purity-margin column.
#' @param colour_col Optional column used to colour bars, such as
#'   `low_purity_fraction`.
#'
#' @return A `ggplot` object.
#' @export
plot_purity_qc <- function(purity,
                           x_col = "dataset_id",
                           y_col = "purity_margin_mean",
                           colour_col = "low_purity_fraction") {
  require_ggplot2()
  purity <- as.data.frame(purity, stringsAsFactors = FALSE)
  check_plot_columns(purity, c(x_col, y_col), "purity")
  has_colour <- !is.null(colour_col) && colour_col %in% colnames(purity)
  purity[[x_col]] <- stats::reorder(factor(purity[[x_col]]), purity[[y_col]])

  mapping <- if (has_colour) {
    aes_columns(x = x_col, y = y_col, fill = colour_col)
  } else {
    aes_columns(x = x_col, y = y_col)
  }
  p <- ggplot2::ggplot(purity, mapping) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#6b6b6b") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_col, fill = colour_col) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (has_colour) {
    p <- p + ggplot2::scale_fill_gradient(low = "#d9f0d3", high = "#238b45")
  }
  p
}

#' Plot state-program pathway enrichment
#'
#' Draws top state-pathway enrichment results as a dot plot.
#'
#' @param enrichment Data frame containing pathway enrichment rows.
#' @param state_col State column.
#' @param pathway_col Pathway or gene-set column.
#' @param q_col Adjusted p-value or q-value column.
#' @param size_col Optional overlap-size column.
#' @param top_n Number of rows to display after sorting by `q_col`.
#'
#' @return A `ggplot` object.
#' @export
plot_pathway_enrichment <- function(enrichment,
                                    state_col = "state",
                                    pathway_col = "pathway",
                                    q_col = "q_value",
                                    size_col = "n_overlap",
                                    top_n = 20) {
  require_ggplot2()
  enrichment <- as.data.frame(enrichment, stringsAsFactors = FALSE)
  check_plot_columns(enrichment, c(state_col, pathway_col, q_col), "enrichment")
  assert_plot_positive_integer(top_n, "top_n")
  enrichment <- enrichment[order(enrichment[[q_col]], decreasing = FALSE), , drop = FALSE]
  enrichment <- utils::head(enrichment, top_n)
  enrichment$neg_log10_q <- -log10(pmax(as.numeric(enrichment[[q_col]]), .Machine$double.xmin))
  enrichment$state_pathway <- paste(enrichment[[state_col]], enrichment[[pathway_col]], sep = " | ")
  enrichment$state_pathway <- factor(enrichment$state_pathway, levels = rev(enrichment$state_pathway))

  mapping <- if (!is.null(size_col) && size_col %in% colnames(enrichment)) {
    aes_columns(x = "neg_log10_q", y = "state_pathway", colour = state_col, size = size_col)
  } else {
    aes_columns(x = "neg_log10_q", y = "state_pathway", colour = state_col)
  }

  ggplot2::ggplot(enrichment, mapping) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::labs(x = "-log10(q)", y = NULL, colour = "State", size = "Overlap") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot FibroDynMix marker support
#'
#' Draws state-level marker support returned by
#' `evaluate_fibrodynmix_annotation()`.
#'
#' @param marker_support Data frame containing marker support rows.
#' @param state_col State column.
#' @param support_col Numeric support column. Defaults to the log2 own-state
#'   versus other-state marker-score ratio.
#'
#' @return A `ggplot` object.
#' @export
plot_marker_support <- function(marker_support,
                                state_col = "state",
                                support_col = "log2_ratio_own_vs_other") {
  require_ggplot2()
  marker_support <- as.data.frame(marker_support, stringsAsFactors = FALSE)
  check_plot_columns(marker_support, c(state_col, support_col), "marker_support")
  marker_support <- marker_support[order(marker_support[[support_col]], decreasing = FALSE), , drop = FALSE]
  marker_support[[state_col]] <- factor(marker_support[[state_col]], levels = marker_support[[state_col]])

  ggplot2::ggplot(marker_support, aes_columns(x = support_col, y = state_col, fill = support_col)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#555555") +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0) +
    ggplot2::labs(x = support_col, y = NULL, fill = "Support") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot Seurat cluster and FibroDynMix state agreement
#'
#' Draws a heatmap from the `cluster_agreement` table returned by
#' `evaluate_fibrodynmix_annotation()`.
#'
#' @param cluster_agreement Data frame containing cluster-state rows.
#' @param cluster_col Cluster column.
#' @param state_col FibroDynMix state column.
#' @param value_col Numeric fraction column.
#'
#' @return A `ggplot` object.
#' @export
plot_cluster_state_agreement <- function(cluster_agreement,
                                         cluster_col = "cluster",
                                         state_col = "state",
                                         value_col = "fraction_within_cluster") {
  require_ggplot2()
  cluster_agreement <- as.data.frame(cluster_agreement, stringsAsFactors = FALSE)
  check_plot_columns(cluster_agreement, c(cluster_col, state_col, value_col), "cluster_agreement")
  cluster_agreement[[cluster_col]] <- factor(cluster_agreement[[cluster_col]], levels = unique(cluster_agreement[[cluster_col]]))
  cluster_agreement[[state_col]] <- factor(cluster_agreement[[state_col]], levels = unique(cluster_agreement[[state_col]]))

  ggplot2::ggplot(cluster_agreement, aes_columns(x = state_col, y = cluster_col, fill = value_col)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", labels = ggplot_percent_labels()) +
    ggplot2::labs(x = "FibroDynMix state", y = "Cluster", fill = "Fraction") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )
}

#' Plot fibroblast annotations on a Seurat embedding
#'
#' Draws a UMAP/tSNE/PCA-style annotation plot from a Seurat object. The default
#' annotation column is the FibroDynMix dominant-state column written by
#' `add_fibrodynmix_to_seurat()`.
#'
#' @param object A Seurat object.
#' @param reduction Reduction name. If `NULL`, common reductions are searched in
#'   order: `umap`, `tsne`, `pca`, `fibrodynmix`.
#' @param dims Two reduction dimensions to plot.
#' @param annotation_col Metadata column used for point colour.
#' @param facet_col Optional metadata column used for faceting.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot` object.
#' @export
plot_fibroblast_annotation <- function(object,
                                       reduction = NULL,
                                       dims = c(1, 2),
                                       annotation_col = "fibrodynmix_dominant_state",
                                       facet_col = NULL,
                                       point_size = 0.45,
                                       point_alpha = 0.8) {
  require_ggplot2()
  require_seurat_object()
  reduction <- resolve_seurat_reduction(object, reduction)
  dims <- validate_embedding_dims(dims)

  embeddings <- SeuratObject::Embeddings(object[[reduction]])
  if (ncol(embeddings) < max(dims)) {
    stop("Selected reduction does not contain the requested dimensions.", call. = FALSE)
  }
  metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  check_plot_columns(metadata, annotation_col, "Seurat metadata")
  if (!is.null(facet_col)) {
    check_plot_columns(metadata, facet_col, "Seurat metadata")
  }

  cells <- intersect(rownames(embeddings), rownames(metadata))
  if (length(cells) == 0L) {
    stop("Embedding cells do not match Seurat metadata row names.", call. = FALSE)
  }
  plot_data <- data.frame(
    cell_id = cells,
    dim_1 = embeddings[cells, dims[[1]]],
    dim_2 = embeddings[cells, dims[[2]]],
    annotation = metadata[cells, annotation_col],
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_col)) {
    plot_data$facet <- metadata[cells, facet_col]
  }
  plot_data$annotation <- factor(plot_data$annotation, levels = unique(plot_data$annotation))

  p <- ggplot2::ggplot(plot_data, aes_columns(x = "dim_1", y = "dim_2", colour = "annotation")) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::labs(
      x = paste0(reduction, "_", dims[[1]]),
      y = paste0(reduction, "_", dims[[2]]),
      colour = annotation_col
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula("~ facet"))
  }
  p
}

#' Plot fibroblast marker expression by annotation
#'
#' Draws a marker dot plot from a Seurat object, grouped by FibroDynMix state or
#' another fibroblast annotation column.
#'
#' @param object A Seurat object.
#' @param features Character vector of marker genes to plot.
#' @param group_col Metadata column used to group cells.
#' @param assay Assay name. If `NULL`, the Seurat default assay is used.
#' @param layer Assay layer to read. Defaults to `counts`.
#' @param slot Legacy Seurat slot fallback.
#' @param expression_transform Expression transform applied before averaging:
#'   `log1p` or `identity`.
#' @param min_pct_expression Expression threshold used for percentage expressed.
#'
#' @return A `ggplot` object.
#' @export
plot_fibroblast_marker_dot <- function(object,
                                       features,
                                       group_col = "fibrodynmix_dominant_state",
                                       assay = NULL,
                                       layer = "counts",
                                       slot = NULL,
                                       expression_transform = c("log1p", "identity"),
                                       min_pct_expression = 0) {
  require_ggplot2()
  require_seurat_object()
  expression_transform <- match.arg(expression_transform)
  if (!is.character(features) || length(features) == 0L || anyNA(features) || any(features == "")) {
    stop("`features` must be a non-empty character vector.", call. = FALSE)
  }
  if (length(min_pct_expression) != 1L || is.na(min_pct_expression) || min_pct_expression < 0) {
    stop("`min_pct_expression` must be a non-negative numeric scalar.", call. = FALSE)
  }

  assay <- resolve_seurat_assay(object, assay)
  expression <- extract_seurat_counts(object, assay = assay, layer = layer, slot = slot)
  metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  check_plot_columns(metadata, group_col, "Seurat metadata")
  cells <- intersect(colnames(expression), rownames(metadata))
  if (length(cells) == 0L) {
    stop("Expression cells do not match Seurat metadata row names.", call. = FALSE)
  }
  features <- unique(features)
  retained_features <- features[features %in% rownames(expression)]
  if (length(retained_features) == 0L) {
    stop("None of `features` are present in the selected Seurat assay.", call. = FALSE)
  }

  expression <- expression[retained_features, cells, drop = FALSE]
  if (expression_transform == "log1p") {
    expression <- log1p(expression)
  }
  groups <- as.character(metadata[cells, group_col])
  keep <- !is.na(groups) & groups != ""
  if (!any(keep)) {
    stop("`group_col` does not contain any non-missing group labels.", call. = FALSE)
  }
  expression <- expression[, keep, drop = FALSE]
  groups <- groups[keep]
  group_levels <- unique(groups)

  rows <- lapply(group_levels, function(group) {
    idx <- groups == group
    group_expression <- expression[, idx, drop = FALSE]
    data.frame(
      group = group,
      feature = retained_features,
      mean_expression = matrix_row_means(group_expression),
      pct_expressed = matrix_row_means(group_expression > min_pct_expression),
      stringsAsFactors = FALSE
    )
  })
  dot_data <- do.call(rbind, rows)
  dot_data$group <- factor(dot_data$group, levels = group_levels)
  dot_data$feature <- factor(dot_data$feature, levels = rev(retained_features))

  ggplot2::ggplot(
    dot_data,
    aes_columns(x = "group", y = "feature", colour = "mean_expression", size = "pct_expressed")
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_size_continuous(labels = ggplot_percent_labels(), range = c(1.2, 6)) +
    ggplot2::scale_colour_gradientn(colours = c("#f7fbff", "#9ecae1", "#3182bd", "#08519c")) +
    ggplot2::labs(x = NULL, y = NULL, colour = "Mean expression", size = "Expressed") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.15, colour = "#e5e5e5"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for FibroDynMix plotting functions.", call. = FALSE)
  }
  invisible(TRUE)
}

aes_columns <- function(...) {
  columns <- list(...)
  do.call(ggplot2::aes, lapply(columns, as.name))
}

check_plot_columns <- function(data, columns, data_name) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0L) {
    stop(
      sprintf("`%s` is missing required column(s): %s.", data_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

infer_state_columns <- function(data, state_cols = NULL) {
  if (!is.null(state_cols)) {
    check_plot_columns(data, state_cols, "cell_weights")
    return(state_cols)
  }
  metadata_cols <- c(
    "cell_id", "dataset_id", "sample_id", "donor_id", "condition", "disease",
    "dominant_state", "entropy", "fpi", "transition_potential"
  )
  numeric_cols <- colnames(data)[vapply(data, is.numeric, logical(1))]
  state_cols <- setdiff(numeric_cols, metadata_cols)
  if (length(state_cols) < 2L) {
    stop("Could not infer state-weight columns; provide `state_cols` explicitly.", call. = FALSE)
  }
  state_cols
}

first_existing_column <- function(data, candidates) {
  hit <- candidates[candidates %in% colnames(data)]
  if (length(hit) == 0L) {
    stop(
      sprintf("Could not infer an identifier column. Tried: %s.", paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  hit[[1]]
}

long_state_weights <- function(cell_id, z, state_cols) {
  rows <- lapply(seq_along(state_cols), function(k) {
    data.frame(
      cell_id = cell_id,
      state = state_cols[[k]],
      state_weight = z[, k],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ggplot_percent_labels <- function() {
  function(x) paste0(round(100 * x), "%")
}

assert_plot_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x <= 0 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

matrix_row_means <- function(x) {
  if (inherits(x, "sparseMatrix") || inherits(x, "Matrix")) {
    return(as.numeric(Matrix::rowMeans(x)))
  }
  x <- as.matrix(x)
  rowMeans(x)
}

resolve_seurat_reduction <- function(object, reduction) {
  reductions <- names(object@reductions)
  if (length(reductions) == 0L) {
    stop("Seurat object does not contain dimensional reductions.", call. = FALSE)
  }
  if (is.null(reduction)) {
    candidates <- c("umap", "tsne", "pca", "fibrodynmix")
    matched <- candidates[candidates %in% reductions]
    reduction <- if (length(matched) > 0L) matched[[1]] else reductions[[1]]
  }
  if (length(reduction) != 1L || is.na(reduction) || !reduction %in% reductions) {
    stop("`reduction` must name a dimensional reduction in the Seurat object.", call. = FALSE)
  }
  reduction
}

validate_embedding_dims <- function(dims) {
  if (length(dims) != 2L || anyNA(dims) || any(dims <= 0) || any(dims != as.integer(dims))) {
    stop("`dims` must contain two positive integer dimensions.", call. = FALSE)
  }
  as.integer(dims)
}
