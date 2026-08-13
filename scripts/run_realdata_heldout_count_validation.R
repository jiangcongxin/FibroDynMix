#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_realdata_heldout_count_validation.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

get_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0L) return(default)
  sub(key, "", hit[[length(hit)]], fixed = TRUE)
}

OUT <- normalizePath(get_arg("out", file.path(ROOT, "analysis", "realdata_heldout_count_validation")), mustWork = FALSE)
N_MARKER_SPLITS <- as.integer(get_arg("marker-splits", "4"))
MAX_GENES <- as.integer(get_arg("max-genes", "500"))
MAX_CELLS_PER_GROUP <- as.integer(get_arg("max-cells-per-group", "60"))
MAX_HOLDOUTS <- as.integer(get_arg("max-holdouts", "99"))
BASE_SEED <- as.integer(get_arg("seed", "20260812"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (file in c(
  "matrix_utils.R", "z_logistic_normal_prior.R", "simulate_fibrodynmix.R", "marker_sets.R",
  "baseline_marker_scoring.R", "fibrodynmix_initializer.R", "nb_likelihood.R", "fit_nb_model.R"
)) source(file.path(ROOT, "R", file))

write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
mean_abs <- function(x) mean(abs(x), na.rm = TRUE)

sample_columns_by_group <- function(counts, metadata, group, max_cells, seed) {
  set.seed(seed)
  keep <- unlist(lapply(split(seq_len(ncol(counts)), group), function(idx) sample(idx, min(length(idx), max_cells))), use.names = FALSE)
  keep <- sort(keep)
  list(counts = counts[, keep, drop = FALSE], metadata = metadata[keep, , drop = FALSE])
}

select_model_genes <- function(counts, markers, max_genes) {
  marker_genes <- intersect(unique(unlist(markers, use.names = FALSE)), rownames(counts))
  detected <- Matrix::rowSums(counts > 0)
  total <- Matrix::rowSums(counts)
  ranked <- order(detected, total, decreasing = TRUE)
  unique(c(marker_genes, rownames(counts)[head(ranked, max_genes)]))
}

load_gse167339 <- function() {
  manifest_path <- file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_human_fibroblast_dataset_manifest.tsv")
  manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    x <- readRDS(manifest$counts_path[i])
    if (is.null(dim(x)) || ncol(x) == 0L) return(NULL)
    colnames(x) <- paste(manifest$dataset_id[i], colnames(x), sep = "__")
    list(
      counts = x,
      metadata = data.frame(
        cell_id = colnames(x), holdout_group = manifest$donor_id[i],
        source_group = manifest$study_id[i], dataset_id = manifest$dataset_id[i],
        stringsAsFactors = FALSE
      )
    )
  })
  parts <- Filter(Negate(is.null), parts)
  common <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
  counts <- do.call(cbind, lapply(parts, function(x) x$counts[common, , drop = FALSE]))
  metadata <- do.call(rbind, lapply(parts, `[[`, "metadata"))
  markers <- get_fibrodynmix_markers("human", "scar")
  genes <- select_model_genes(counts, markers, MAX_GENES)
  counts <- counts[genes, , drop = FALSE]
  markers <- lapply(markers, intersect, y = rownames(counts))
  sampled <- sample_columns_by_group(counts, metadata, metadata$holdout_group, MAX_CELLS_PER_GROUP, BASE_SEED + 167339L)
  list(dataset = "GSE167339", counts = as.matrix(sampled$counts), metadata = sampled$metadata, markers = markers, holdout_unit = "donor")
}

load_gse246215 <- function() {
  base <- file.path(ROOT, "analysis", "gse246215_sensitivity", "runs", "raw_sample_seed246215")
  manifest <- read.delim(file.path(base, "dataset_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  selected_metadata <- read.delim(file.path(base, "selected_cell_metadata.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    x <- readRDS(manifest$counts_path[i])
    cells <- intersect(colnames(x), selected_metadata$CellName)
    x <- x[, cells, drop = FALSE]
    meta <- selected_metadata[match(cells, selected_metadata$CellName), , drop = FALSE]
    list(
      counts = x,
      metadata = data.frame(
        cell_id = cells, holdout_group = meta$CancerType_short,
        source_group = meta$CancerType_short, dataset_id = manifest$dataset_id[i],
        stringsAsFactors = FALSE
      )
    )
  })
  common <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
  counts <- do.call(cbind, lapply(parts, function(x) x$counts[common, , drop = FALSE]))
  metadata <- do.call(rbind, lapply(parts, `[[`, "metadata"))
  markers <- get_fibrodynmix_markers("human", "caf")
  genes <- select_model_genes(counts, markers, MAX_GENES)
  counts <- counts[genes, , drop = FALSE]
  markers <- lapply(markers, intersect, y = rownames(counts))
  sampled <- sample_columns_by_group(counts, metadata, metadata$holdout_group, MAX_CELLS_PER_GROUP, BASE_SEED + 246215L)
  list(dataset = "GSE246215", counts = as.matrix(sampled$counts), metadata = sampled$metadata, markers = markers, holdout_unit = "cancer_dataset")
}

fit_one <- function(dataset, holdout_group, marker_split) {
  counts <- dataset$counts
  metadata <- dataset$metadata
  markers <- dataset$markers
  set.seed(BASE_SEED + sum(utf8ToInt(dataset$dataset)) + marker_split * 10000L)
  input_markers <- lapply(markers, function(x) sample(x, floor(length(x) / 2L)))
  heldout_markers <- Map(setdiff, markers, input_markers)
  if (any(vapply(input_markers, length, integer(1)) < 2L) || any(vapply(heldout_markers, length, integer(1)) < 2L)) {
    stop("Every state requires at least two input and two held-out markers.", call. = FALSE)
  }
  train <- metadata$holdout_group != holdout_group
  test <- !train
  train_counts <- counts[, train, drop = FALSE]
  test_counts <- counts[, test, drop = FALSE]
  train_lib <- pmax(colSums(train_counts), 1)
  test_lib <- pmax(colSums(test_counts), 1)
  marker_train <- score_marker_baseline(train_counts, input_markers, library_size = train_lib)$z_pred
  marker_test <- score_marker_baseline(test_counts, input_markers, library_size = test_lib)$z_pred
  initializer <- fit_fibrodynmix_initializer(
    counts = train_counts, marker_index = input_markers,
    library_size = train_lib, n_iter = 3
  )
  fit <- fit_fibrodynmix_nb(
    counts = train_counts, marker_index = input_markers, library_size = train_lib,
    initial_state = list(z_hat = marker_train, beta_hat = initializer$beta_hat),
    n_outer = 3, estimate_phi = TRUE, beta_l2 = 0.01, marker_l2 = 0.05,
    beta_constraint = "sum_to_zero", fit_study_effect = FALSE,
    z_l2 = 0, z_anchor = "initializer_logit", z_anchor_sd = 0.1,
    maxit_beta = 25, maxit_z = 20, early_stopping = FALSE, rollback_to_best = TRUE
  )
  eval_genes <- unique(unlist(heldout_markers, use.names = FALSE))
  null_alpha <- initialize_nb_alpha(train_counts, train_lib)
  null_phi <- estimate_phi_moments(train_counts, train_lib)
  null_beta <- matrix(0, nrow = ncol(marker_test), ncol = length(eval_genes),
                      dimnames = list(colnames(marker_test), eval_genes))
  null_ll <- fibrodynmix_nb_loglik(
    counts = test_counts[eval_genes, , drop = FALSE], z = marker_test,
    beta = null_beta, alpha = null_alpha[eval_genes], phi = null_phi[eval_genes],
    library_size = test_lib
  ) / (length(eval_genes) * ncol(test_counts))
  anchored_ll <- fibrodynmix_nb_loglik(
    counts = test_counts[eval_genes, , drop = FALSE], z = marker_test,
    beta = fit$beta_hat[, eval_genes, drop = FALSE], alpha = fit$alpha_hat[eval_genes],
    phi = fit$phi_hat[eval_genes], library_size = test_lib
  ) / (length(eval_genes) * ncol(test_counts))
  thinned <- matrix(rbinom(length(test_counts), as.vector(test_counts), 0.5), nrow = nrow(test_counts), dimnames = dimnames(test_counts))
  thin_lib <- pmax(colSums(thinned), 1)
  marker_thin <- score_marker_baseline(thinned, input_markers, library_size = thin_lib)$z_pred
  data.frame(
    dataset = dataset$dataset, holdout_unit = dataset$holdout_unit,
    holdout_group = holdout_group, marker_split = marker_split,
    n_train_cells = sum(train), n_test_cells = sum(test), n_model_genes = nrow(counts),
    n_input_markers = length(unique(unlist(input_markers, use.names = FALSE))),
    n_heldout_state_genes = length(eval_genes),
    intercept_loglik_per_observation = null_ll,
    anchored_loglik_per_observation = anchored_ll,
    anchored_gain = anchored_ll - null_ll,
    train_drift_from_marker = mean_abs(fit$z_hat - marker_train),
    marker_downsample_mean_abs_delta = mean_abs(marker_thin - marker_test),
    stringsAsFactors = FALSE
  )
}

datasets <- list(load_gse167339(), load_gse246215())
rows <- list()
for (dataset in datasets) {
  holdouts <- head(sort(unique(dataset$metadata$holdout_group)), MAX_HOLDOUTS)
  for (holdout in holdouts) {
    for (marker_split in seq_len(N_MARKER_SPLITS)) {
      message(sprintf("%s: holdout=%s split=%d", dataset$dataset, holdout, marker_split))
      rows[[length(rows) + 1L]] <- fit_one(dataset, holdout, marker_split)
    }
  }
}
metrics <- do.call(rbind, rows)
write_tsv(metrics, file.path(OUT, "realdata_heldout_metrics.tsv"))

group_means <- aggregate(
  cbind(intercept_loglik_per_observation, anchored_loglik_per_observation, anchored_gain, train_drift_from_marker, marker_downsample_mean_abs_delta) ~ dataset + holdout_unit + holdout_group,
  data = metrics, FUN = mean
)
write_tsv(group_means, file.path(OUT, "realdata_holdout_group_means.tsv"))
summary <- aggregate(
  cbind(intercept_loglik_per_observation, anchored_loglik_per_observation, anchored_gain, train_drift_from_marker, marker_downsample_mean_abs_delta) ~ dataset,
  data = group_means, FUN = mean
)
summary$favorable_holdouts <- vapply(summary$dataset, function(x) sum(group_means$dataset == x & group_means$anchored_gain > 0), integer(1))
summary$total_holdouts <- vapply(summary$dataset, function(x) sum(group_means$dataset == x), integer(1))
write_tsv(summary, file.path(OUT, "realdata_heldout_summary.tsv"))

png(file.path(OUT, "realdata_heldout_validation.png"), width = 1400, height = 700, res = 180)
par(mfrow = c(1, 2), mar = c(8, 5, 3, 1))
boxplot(anchored_gain ~ dataset, data = group_means, las = 2, ylab = "Anchored NB gain in held-out log-likelihood", xlab = "", main = "Held-out group prediction")
abline(h = 0, lty = 2, col = "grey40")
boxplot(train_drift_from_marker ~ dataset, data = group_means, las = 2, ylab = "Mean absolute coordinate drift", xlab = "", main = "Marker-coordinate retention")
dev.off()

report <- c(
  "# Real-data held-out count validation",
  "",
  "## Design",
  "",
  sprintf("- %d random marker splits per held-out group.", N_MARKER_SPLITS),
  "- Half of each curated marker set defined the coordinate; complementary markers were excluded from initialization and used only as held-out state-associated genes.",
  "- GSE167339 used donor-level holdout. GSE246215 used cancer-dataset holdout.",
  "- Programs were trained without donor or study terms because the held-out group was unseen.",
  "",
  "## Results",
  "",
  unlist(lapply(seq_len(nrow(summary)), function(i) sprintf(
    "- %s: mean held-out gain %.4f; favorable holdouts %d/%d; mean coordinate drift %.4f.",
    summary$dataset[i], summary$anchored_gain[i], summary$favorable_holdouts[i], summary$total_holdouts[i], summary$train_drift_from_marker[i]
  ))),
  "",
  "## Claim boundary",
  "",
  "Positive held-out count likelihood shows that the marker-anchored program generalizes numerically to an unseen donor or dataset for curated state-associated genes. It does not establish that the fitted coordinates are biological truth, nor does it test disease effects or lineage."
)
writeLines(report, file.path(OUT, "realdata_heldout_report.md"))

stopifnot(
  nrow(metrics) > 0L,
  all(is.finite(metrics$anchored_gain)),
  all(metrics$train_drift_from_marker < 0.1),
  all(metrics$n_heldout_state_genes >= 12L)
)
message("Real-data held-out validation written to: ", OUT)
