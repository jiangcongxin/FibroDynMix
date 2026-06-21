#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(FibroDynMix)
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required to read real-data Seurat objects.", call. = FALSE)
}

out_dir <- file.path("analysis", "nb_transfer_bootstrap_realdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  GSE163973 = "projects/keloid_fibro_gse163973/results/fibro_seurat_with_fibrodynmix.rds",
  GSE243716 = "projects/scar_fibro_gse243716/results/gse243716_fibro_seurat_with_fibrodynmix.rds"
)
datasets <- datasets[file.exists(unlist(datasets))]
if (length(datasets) == 0L) {
  stop("No real-data Seurat result objects were found.", call. = FALSE)
}

n_boot <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_BOOT", "5"))
max_train_cells <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_TRAIN_CELLS", "140"))
max_eval_cells <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_EVAL_CELLS", "120"))
max_model_genes <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_MODEL_GENES", "450"))
maxit_beta <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_MAXIT_BETA", "15"))
maxit_fit_z <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_FIT_MAXIT_Z", "15"))
maxit_transfer_z <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_MAXIT_Z", "100"))
seed <- as.integer(Sys.getenv("FIBRODYNMIX_NB_TRANSFER_SEED", "20260606"))

read_counts <- function(object, assay = "RNA", layer = "counts") {
  assay <- if (assay %in% names(object@assays)) assay else SeuratObject::DefaultAssay(object)
  counts <- suppressWarnings(SeuratObject::LayerData(object, assay = assay, layer = layer))
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    counts <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
  }
  counts
}

metadata_with_ids <- function(object) {
  metadata <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  metadata$cell_id <- rownames(metadata)
  metadata
}

select_model_genes <- function(counts, marker_index, cells, max_model_genes) {
  marker_genes <- unique(unlist(marker_index, use.names = FALSE))
  detected <- Matrix::rowSums(counts[, cells, drop = FALSE])
  variable_genes <- names(sort(detected, decreasing = TRUE))
  genes <- unique(c(marker_genes, utils::head(variable_genes, max_model_genes)))
  intersect(genes, rownames(counts))
}

fit_transfer_bootstrap <- function(dataset, object_path) {
  object <- readRDS(object_path)
  counts <- read_counts(object)
  metadata <- metadata_with_ids(object)
  marker_index <- get_fibrodynmix_markers("human", "scar")
  marker_index <- lapply(marker_index, intersect, y = rownames(counts))
  marker_index <- marker_index[vapply(marker_index, length, integer(1)) > 0L]

  set.seed(seed + match(dataset, names(datasets)))
  all_cells <- intersect(colnames(counts), metadata$cell_id)
  eval_cells <- sample(all_cells, min(length(all_cells), max_eval_cells))
  candidate_train_cells <- setdiff(all_cells, eval_cells)
  if (length(candidate_train_cells) < 20L) {
    candidate_train_cells <- all_cells
  }
  model_genes <- select_model_genes(counts, marker_index, c(eval_cells, candidate_train_cells), max_model_genes)
  counts <- counts[model_genes, , drop = FALSE]
  marker_index <- lapply(marker_index, intersect, y = rownames(counts))
  marker_index <- marker_index[vapply(marker_index, length, integer(1)) > 0L]

  draws <- list()
  diagnostics <- list()
  for (b in seq_len(n_boot)) {
    set.seed(seed + b + 1000L * match(dataset, names(datasets)))
    train_cells <- sample(candidate_train_cells, min(length(candidate_train_cells), max_train_cells), replace = TRUE)
    train_counts <- counts[, train_cells, drop = FALSE]
    colnames(train_counts) <- paste0(train_cells, "__boot", b, "_", seq_along(train_cells))
    train_metadata <- metadata[match(train_cells, metadata$cell_id), , drop = FALSE]
    train_metadata$cell_id <- colnames(train_counts)

    prepared <- prepare_fibrodynmix_data(
      counts = train_counts,
      cell_metadata = train_metadata,
      marker_index = marker_index,
      cell_id_col = "cell_id",
      study_col = if ("condition" %in% colnames(train_metadata)) "condition" else NULL,
      donor_col = if ("sample_id" %in% colnames(train_metadata)) "sample_id" else NULL,
      min_cells_per_gene = 0,
      min_counts_per_gene = 0,
      require_all_states = TRUE
    )
    fit <- fit_fibrodynmix_prepared(
      prepared,
      fit_study_effect = "condition" %in% colnames(train_metadata),
      fit_donor_effect = FALSE,
      n_outer = 1,
      initializer_args = list(n_iter = 3),
      maxit_beta = maxit_beta,
      maxit_z = maxit_fit_z,
      early_stopping = TRUE
    )
    eval_counts <- counts[colnames(fit$beta_hat), eval_cells, drop = FALSE]
    eval_library <- as.numeric(Matrix::colSums(eval_counts))
    warm <- score_marker_baseline(eval_counts, marker_index = marker_index, library_size = eval_library)$z_pred
    warm <- warm[, rownames(fit$beta_hat), drop = FALSE]
    transfer <- fit_fibrodynmix_transfer(
      counts = eval_counts,
      fit = fit,
      library_size = eval_library,
      z_init = warm,
      chunk_size = 250,
      maxit_z = maxit_transfer_z,
      parallel = FALSE,
      return_cell_diagnostics = TRUE
    )

    z <- as.data.frame(transfer$z_hat, stringsAsFactors = FALSE)
    z$cell_id <- rownames(transfer$z_hat)
    z$replicate <- b
    z$dataset <- dataset
    draws[[b]] <- z[, c("dataset", "replicate", "cell_id", colnames(transfer$z_hat)), drop = FALSE]
    diagnostics[[b]] <- data.frame(
      dataset = dataset,
      replicate = b,
      n_train_cells = ncol(prepared$counts),
      n_eval_cells = ncol(eval_counts),
      n_genes = nrow(eval_counts),
      nb_objective = transfer$nb_objective,
      z_convergence_rate = transfer$z_convergence_rate,
      n_nonconverged_cells = transfer$n_nonconverged_cells,
      train_best_objective = fit$best_objective,
      train_stop_reason = fit$stop_reason,
      stringsAsFactors = FALSE
    )
  }
  list(draws = do.call(rbind, draws), diagnostics = do.call(rbind, diagnostics))
}

results <- lapply(names(datasets), function(dataset) fit_transfer_bootstrap(dataset, datasets[[dataset]]))
names(results) <- names(datasets)
draws <- do.call(rbind, lapply(results, `[[`, "draws"))
diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
state_cols <- setdiff(colnames(draws), c("dataset", "replicate", "cell_id"))

summary_rows <- list()
idx <- 1L
for (dataset in unique(draws$dataset)) {
  for (state in state_cols) {
    subset <- draws[draws$dataset == dataset, c("cell_id", "replicate", state), drop = FALSE]
    colnames(subset)[3L] <- "z"
    cells <- unique(subset$cell_id)
    cell_sd <- vapply(cells, function(cell) {
      stats::sd(subset$z[subset$cell_id == cell], na.rm = TRUE)
    }, numeric(1))
    summary_rows[[idx]] <- data.frame(
      dataset = dataset,
      state = state,
      bootstrap_method = "nb_fit_plus_transfer",
      bootstrap_n = n_boot,
      n_eval_cells = length(cells),
      median_cell_z_sd = stats::median(cell_sd, na.rm = TRUE),
      p90_cell_z_sd = as.numeric(stats::quantile(cell_sd, 0.9, na.rm = TRUE, names = FALSE)),
      max_cell_z_sd = max(cell_sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
summary <- do.call(rbind, summary_rows)

manifest <- data.frame(
  analysis = "nb_transfer_bootstrap_realdata",
  datasets = paste(names(datasets), collapse = ","),
  bootstrap_n = n_boot,
  max_train_cells = max_train_cells,
  max_eval_cells = max_eval_cells,
  max_model_genes = max_model_genes,
  maxit_beta = maxit_beta,
  maxit_fit_z = maxit_fit_z,
  maxit_transfer_z = maxit_transfer_z,
  median_transfer_convergence_rate = stats::median(diagnostics$z_convergence_rate, na.rm = TRUE),
  min_transfer_convergence_rate = min(diagnostics$z_convergence_rate, na.rm = TRUE),
  max_state_p90_cell_z_sd = max(summary$p90_cell_z_sd, na.rm = TRUE),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

utils::write.table(draws, file.path(out_dir, "nb_transfer_bootstrap_draws.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(diagnostics, file.path(out_dir, "nb_transfer_bootstrap_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(summary, file.path(out_dir, "nb_transfer_bootstrap_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(manifest, file.path(out_dir, "nb_transfer_bootstrap_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

print(summary)
print(manifest)
