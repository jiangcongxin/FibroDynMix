#' Fit an optional scVI latent baseline
#'
#' Fits a small scVI model through `reticulate`, then projects the learned latent
#' representation to a state simplex by soft cluster membership. Clusters are
#' aligned to fibroblast states using marker-score enrichment, so marker genes
#' orient labels but do not define the deep latent representation.
#'
#' This baseline is optional because it requires a Python environment with
#' `scvi-tools` and `anndata`.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of marker genes for each state.
#' @param library_size Optional vector of per-cell library sizes used only for
#'   marker-score alignment.
#' @param n_latent scVI latent dimensionality.
#' @param max_epochs Training epochs passed to `scvi.model.SCVI.train()`.
#' @param seed Optional seed passed to Python and R clustering.
#' @param projection_temperature Optional distance softmax temperature. If
#'   omitted, the median latent distance to cluster centers is used.
#'
#' @return A list with `latent`, `z_pred`, cluster assignments, cluster-state
#'   alignment, and projection metadata.
#' @export
fit_scvi_latent_baseline <- function(counts,
                                     marker_index,
                                     library_size = NULL,
                                     n_latent = length(marker_index),
                                     max_epochs = 40,
                                     seed = NULL,
                                     projection_temperature = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("`fit_scvi_latent_baseline()` requires the reticulate R package.", call. = FALSE)
  }
  if (!reticulate::py_module_available("scvi") || !reticulate::py_module_available("anndata")) {
    stop("`fit_scvi_latent_baseline()` requires Python modules scvi and anndata.", call. = FALSE)
  }
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("`counts` must contain finite non-negative numeric values.", call. = FALSE)
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
  assert_positive_integer(n_latent, "n_latent")
  assert_positive_integer(max_epochs, "max_epochs")

  if (!is.null(seed)) {
    reticulate::py_run_string(sprintf(
      paste(
        "import random, numpy as np",
        "random.seed(%d)",
        "np.random.seed(%d)",
        sep = "\n"
      ),
      as.integer(seed),
      as.integer(seed)
    ))
  }

  anndata <- reticulate::import("anndata", delay_load = FALSE)
  scvi <- reticulate::import("scvi", delay_load = FALSE)
  adata <- anndata$AnnData(X = t(counts))
  adata$var_names <- rownames(counts)
  adata$obs_names <- colnames(counts)
  scvi$model$SCVI$setup_anndata(adata)
  model <- scvi$model$SCVI(adata, n_latent = as.integer(n_latent), gene_likelihood = "nb")
  model$train(max_epochs = as.integer(max_epochs))
  latent <- as.matrix(model$get_latent_representation())
  rownames(latent) <- colnames(counts)
  colnames(latent) <- sprintf("scvi_latent_%d", seq_len(ncol(latent)))

  projection <- project_latent_to_state_simplex(
    latent = latent,
    counts = counts,
    marker_index = marker_index,
    library_size = library_size,
    seed = seed,
    projection_temperature = projection_temperature
  )

  c(
    list(
      latent = latent,
      max_epochs = max_epochs,
      n_latent = n_latent
    ),
    projection
  )
}

project_latent_to_state_simplex <- function(latent,
                                           counts,
                                           marker_index,
                                           library_size = NULL,
                                           seed = NULL,
                                           projection_temperature = NULL) {
  latent <- as.matrix(latent)
  n_states <- length(marker_index)
  assert_positive_integer(n_states, "n_states")
  if (nrow(latent) != ncol(counts)) {
    stop("`latent` must have one row per cell in `counts`.", call. = FALSE)
  }
  if (!is.null(seed)) {
    old_seed <- preserve_seed()
    on.exit(restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  km <- stats::kmeans(latent, centers = n_states, nstart = 5, iter.max = 50)
  distances <- matrix(
    NA_real_,
    nrow = nrow(latent),
    ncol = n_states,
    dimnames = list(rownames(latent), sprintf("cluster_%d", seq_len(n_states)))
  )
  for (cluster in seq_len(n_states)) {
    centered <- sweep(latent, 2, km$centers[cluster, ], "-")
    distances[, cluster] <- sqrt(rowSums(centered^2))
  }
  if (is.null(projection_temperature)) {
    projection_temperature <- stats::median(distances[is.finite(distances)])
    if (is.na(projection_temperature) || projection_temperature <= 0) {
      projection_temperature <- 1
    }
  }
  topic_prob <- softmax_rows(-distances / projection_temperature)
  colnames(topic_prob) <- colnames(distances)

  marker_scores <- score_marker_baseline(
    counts = counts,
    marker_index = marker_index,
    library_size = library_size
  )$scores
  cluster_state_scores <- matrix(
    NA_real_,
    nrow = n_states,
    ncol = length(marker_index),
    dimnames = list(colnames(topic_prob), names(marker_index))
  )
  for (cluster in seq_len(n_states)) {
    keep <- km$cluster == cluster
    cluster_state_scores[cluster, ] <- colMeans(marker_scores[keep, , drop = FALSE])
  }
  alignment <- greedy_align_clusters_to_states(cluster_state_scores)

  state_names <- names(marker_index)
  z_pred <- matrix(
    .Machine$double.eps,
    nrow = nrow(topic_prob),
    ncol = length(state_names),
    dimnames = list(rownames(topic_prob), state_names)
  )
  for (state in state_names) {
    cluster <- alignment$state_to_cluster[[state]]
    if (!is.na(cluster)) {
      z_pred[, state] <- topic_prob[, cluster]
    }
  }
  z_pred <- sweep(z_pred, 1, rowSums(z_pred), "/")

  list(
    z_pred = z_pred,
    cluster = km$cluster,
    centers = km$centers,
    cluster_state_scores = cluster_state_scores,
    cluster_to_state = alignment$cluster_to_state,
    state_to_cluster = alignment$state_to_cluster,
    projection_temperature = projection_temperature
  )
}

greedy_align_clusters_to_states <- function(cluster_state_scores) {
  cluster_names <- rownames(cluster_state_scores)
  state_names <- colnames(cluster_state_scores)
  available_clusters <- cluster_names
  available_states <- state_names
  assignments <- data.frame(cluster = character(), state = character(), stringsAsFactors = FALSE)
  while (length(available_clusters) > 0L && length(available_states) > 0L) {
    sub_scores <- cluster_state_scores[available_clusters, available_states, drop = FALSE]
    best <- which(sub_scores == max(sub_scores), arr.ind = TRUE)[1, , drop = FALSE]
    cluster <- rownames(sub_scores)[best[1, "row"]]
    state <- colnames(sub_scores)[best[1, "col"]]
    assignments <- rbind(assignments, data.frame(cluster = cluster, state = state, stringsAsFactors = FALSE))
    available_clusters <- setdiff(available_clusters, cluster)
    available_states <- setdiff(available_states, state)
  }
  cluster_to_state <- stats::setNames(rep(NA_character_, length(cluster_names)), cluster_names)
  state_to_cluster <- stats::setNames(rep(NA_character_, length(state_names)), state_names)
  for (i in seq_len(nrow(assignments))) {
    cluster_to_state[[assignments$cluster[i]]] <- assignments$state[i]
    state_to_cluster[[assignments$state[i]]] <- assignments$cluster[i]
  }
  list(cluster_to_state = cluster_to_state, state_to_cluster = state_to_cluster)
}
