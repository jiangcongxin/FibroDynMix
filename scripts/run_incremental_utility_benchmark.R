#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_incremental_utility_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

get_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0L) return(default)
  sub(key, "", hit[[length(hit)]], fixed = TRUE)
}

OUT <- normalizePath(get_arg("out", file.path(ROOT, "analysis", "incremental_utility_benchmark")), mustWork = FALSE)
N_REP <- as.integer(get_arg("replicates", "20"))
N_MARKER_SPLITS <- as.integer(get_arg("marker-splits", "4"))
BASE_SEED <- as.integer(get_arg("seed", "20260811"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (file in c(
  "matrix_utils.R", "z_logistic_normal_prior.R", "simulate_fibrodynmix.R",
  "benchmark_metrics.R", "baseline_marker_scoring.R", "fibrodynmix_initializer.R",
  "nb_likelihood.R", "fit_nb_model.R"
)) source(file.path(ROOT, "R", file))

mean_abs <- function(x) mean(abs(x), na.rm = TRUE)

stratified_split <- function(study_id, fraction = 0.25, seed = 1L) {
  set.seed(seed)
  holdout <- rep(FALSE, length(study_id))
  for (level in unique(study_id)) {
    idx <- which(study_id == level)
    holdout[sample(idx, max(1L, round(length(idx) * fraction)))] <- TRUE
  }
  holdout
}

center_group_effect <- function(x) sweep(x, 2L, colMeans(x), "-")

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) return(NA_real_)
  stats::cor(x[keep], y[keep])
}

heldout_nonmarker_loglik <- function(counts, z, fit, library_size, genes, study_id = NULL) {
  effect <- if (is.null(study_id)) NULL else fit$study_effect[, genes, drop = FALSE]
  sum(fibrodynmix_nb_loglik(
    counts = counts[genes, , drop = FALSE], z = z,
    beta = fit$beta_hat[, genes, drop = FALSE], alpha = fit$alpha_hat[genes],
    phi = fit$phi_hat[genes], library_size = library_size,
    study_effect = effect, study_id = study_id
  )) / (length(genes) * ncol(counts))
}

project_with_marker_anchor <- function(counts, marker_z, fit, library_size, study_id) {
  prior <- prepare_initializer_logit_anchor(marker_z, logit_sd = 0.1)
  update_z_nb(
    counts = counts, z = marker_z, beta = fit$beta_hat, alpha = fit$alpha_hat,
    phi = fit$phi_hat, library_size = library_size,
    study_effect = fit$study_effect, study_id = study_id,
    donor_effect = NULL, donor_id = NULL, z_l2 = 0,
    z_prior = prior, optimizer = "BFGS", optimizer_control = list(), maxit = 20
  )$z
}

metric_rows <- list()
fit_rows <- list()

for (replicate_id in seq_len(N_REP)) {
  seed <- BASE_SEED + replicate_id - 1L
  sim <- simulate_fibrodynmix(
    n_studies = 3, donors_per_study = 3, cells_per_donor = 10,
    n_genes = 216, marker_genes_per_state = 12,
    scenario = "batch_confounding", study_effect_sd = 0.30,
    donor_effect_sd = 0.08, seed = seed
  )
  holdout <- stratified_split(sim$cell_metadata$study_id, seed = seed + 1000L)
  train <- !holdout
  for (marker_split in seq_len(N_MARKER_SPLITS)) {
  set.seed(seed + marker_split * 10000L)
  input_marker_index <- lapply(sim$parameters$marker_index, function(x) sample(x, length(x) / 2L))
  heldout_marker_index <- Map(setdiff, sim$parameters$marker_index, input_marker_index)
  input_marker_genes <- rownames(sim$counts)[unique(unlist(input_marker_index, use.names = FALSE))]
  heldout_state_genes <- rownames(sim$counts)[unique(unlist(heldout_marker_index, use.names = FALSE))]
  noninput_genes <- setdiff(rownames(sim$counts), input_marker_genes)

  marker <- score_marker_baseline(
    sim$counts, input_marker_index,
    library_size = sim$cell_metadata$library_size
  )$z_pred
  init <- fit_fibrodynmix_initializer(
    counts = sim$counts[, train, drop = FALSE],
    marker_index = input_marker_index,
    library_size = sim$cell_metadata$library_size[train], n_iter = 3
  )
  initial_state <- list(z_hat = marker[train, , drop = FALSE], beta_hat = init$beta_hat)
  common <- list(
    counts = sim$counts[, train, drop = FALSE],
    marker_index = input_marker_index,
    library_size = sim$cell_metadata$library_size[train],
    initial_state = initial_state, n_outer = 3,
    estimate_phi = TRUE, beta_l2 = 0.01, marker_l2 = 0.05,
    beta_constraint = "sum_to_zero", z_l2 = 0,
    z_anchor = "initializer_logit", z_anchor_sd = 0.1,
    maxit_beta = 25, maxit_z = 20, early_stopping = FALSE,
    rollback_to_best = TRUE
  )
  fit_no_source <- do.call(fit_fibrodynmix_nb, c(common, list(fit_study_effect = FALSE)))
  fit_source <- do.call(fit_fibrodynmix_nb, c(common, list(
    study_id = sim$cell_metadata$study_id[train], fit_study_effect = TRUE, study_l2 = 0.1
  )))

  test_counts <- sim$counts[, holdout, drop = FALSE]
  test_lib <- sim$cell_metadata$library_size[holdout]
  test_study <- sim$cell_metadata$study_id[holdout]
  marker_test <- marker[holdout, , drop = FALSE]

  null_alpha <- initialize_nb_alpha(sim$counts[, train, drop = FALSE], sim$cell_metadata$library_size[train])
  null_phi <- estimate_phi_moments(sim$counts[, train, drop = FALSE], sim$cell_metadata$library_size[train])
  null_beta <- matrix(0, nrow = ncol(marker), ncol = length(heldout_state_genes),
                      dimnames = list(colnames(marker), heldout_state_genes))
  null_ll <- fibrodynmix_nb_loglik(
    counts = test_counts[heldout_state_genes, , drop = FALSE], z = marker_test,
    beta = null_beta, alpha = null_alpha[heldout_state_genes], phi = null_phi[heldout_state_genes],
    library_size = test_lib
  ) / (length(heldout_state_genes) * ncol(test_counts))
  no_source_ll <- heldout_nonmarker_loglik(test_counts, marker_test, fit_no_source, test_lib, heldout_state_genes)
  source_ll <- heldout_nonmarker_loglik(test_counts, marker_test, fit_source, test_lib, heldout_state_genes, test_study)

  truth_beta <- center_state_program_matrix(sim$parameters$beta_kg[, heldout_state_genes, drop = FALSE])
  beta_cor_no_source <- safe_cor(as.vector(truth_beta), as.vector(fit_no_source$beta_hat[, heldout_state_genes, drop = FALSE]))
  beta_cor_source <- safe_cor(as.vector(truth_beta), as.vector(fit_source$beta_hat[, heldout_state_genes, drop = FALSE]))
  true_study <- center_group_effect(sim$parameters$study_effect[, noninput_genes, drop = FALSE])
  estimated_study <- center_group_effect(fit_source$study_effect[, noninput_genes, drop = FALSE])
  study_cor <- safe_cor(as.vector(true_study), as.vector(estimated_study))

  thinned <- matrix(stats::rbinom(length(test_counts), size = as.vector(test_counts), prob = 0.5),
                    nrow = nrow(test_counts), dimnames = dimnames(test_counts))
  thin_lib <- pmax(colSums(thinned), 1)
  marker_thin <- score_marker_baseline(thinned, input_marker_index, library_size = thin_lib)$z_pred
  projection_full <- project_with_marker_anchor(test_counts, marker_test, fit_source, test_lib, test_study)
  projection_thin <- project_with_marker_anchor(thinned, marker_thin, fit_source, thin_lib, test_study)

  z_metrics <- list(
    marker_scoring = evaluate_state_weights(sim$z[train, , drop = FALSE], marker[train, , drop = FALSE]),
    anchored_no_source = evaluate_state_weights(sim$z[train, , drop = FALSE], fit_no_source$z_hat),
    anchored_source_aware = evaluate_state_weights(sim$z[train, , drop = FALSE], fit_source$z_hat)
  )
  for (method in names(z_metrics)) {
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      replicate = replicate_id, marker_split = marker_split, method = method,
      state_weight_rmse = z_metrics[[method]]$rmse,
      dominant_accuracy = z_metrics[[method]]$dominant_accuracy,
      drift_from_marker = if (method == "marker_scoring") 0 else mean_abs(
        (if (method == "anchored_no_source") fit_no_source$z_hat else fit_source$z_hat) - marker[train, , drop = FALSE]
      ), stringsAsFactors = FALSE
    )
  }
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(
    replicate = replicate_id, marker_split = marker_split,
    null_heldout_gene_loglik_per_observation = null_ll,
    anchored_no_source_heldout_gene_loglik_per_observation = no_source_ll,
    anchored_source_heldout_gene_loglik_per_observation = source_ll,
    source_gain_over_null = source_ll - null_ll,
    source_gain_over_no_source = source_ll - no_source_ll,
    beta_truth_correlation_no_source = beta_cor_no_source,
    beta_truth_correlation_source = beta_cor_source,
    study_effect_truth_correlation = study_cor,
    marker_downsample_mean_abs_delta = mean_abs(marker_thin - marker_test),
    anchored_projection_downsample_mean_abs_delta = mean_abs(projection_thin - projection_full),
    marker_holdout_rmse = evaluate_state_weights(sim$z[holdout, , drop = FALSE], marker_test)$rmse,
    anchored_projection_holdout_rmse = evaluate_state_weights(sim$z[holdout, , drop = FALSE], projection_full)$rmse,
    stringsAsFactors = FALSE
  )
  }
}

cell_metrics <- do.call(rbind, metric_rows)
fit_metrics <- do.call(rbind, fit_rows)
write.table(cell_metrics, file.path(OUT, "state_coordinate_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fit_metrics, file.path(OUT, "incremental_utility_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summarize_numeric <- function(data, fields) {
  do.call(rbind, lapply(fields, function(field) data.frame(
    metric = field, mean = mean(data[[field]], na.rm = TRUE),
    sd = stats::sd(data[[field]], na.rm = TRUE),
    median = stats::median(data[[field]], na.rm = TRUE), stringsAsFactors = FALSE
  )))
}
replicate_means <- aggregate(
  . ~ replicate,
  data = fit_metrics[, setdiff(colnames(fit_metrics), "marker_split"), drop = FALSE],
  FUN = mean
)
write.table(replicate_means, file.path(OUT, "simulation_replicate_means.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary_table <- summarize_numeric(replicate_means, setdiff(colnames(replicate_means), "replicate"))
write.table(summary_table, file.path(OUT, "incremental_utility_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
coordinate_summary <- aggregate(
  cbind(state_weight_rmse, dominant_accuracy, drift_from_marker) ~ method,
  data = cell_metrics, FUN = function(x) mean(x, na.rm = TRUE)
)
write.table(coordinate_summary, file.path(OUT, "state_coordinate_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
marker_split_variability <- aggregate(
  cbind(source_gain_over_null, source_gain_over_no_source, marker_downsample_mean_abs_delta, anchored_projection_downsample_mean_abs_delta) ~ replicate,
  data = fit_metrics, FUN = stats::sd
)
write.table(marker_split_variability, file.path(OUT, "marker_split_variability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
paired_checks <- data.frame(
  contrast = c("anchored_source_vs_intercept", "anchored_source_vs_no_source", "anchored_projection_vs_marker_downsampling"),
  mean_difference = c(
    mean(replicate_means$source_gain_over_null),
    mean(replicate_means$source_gain_over_no_source),
    mean(replicate_means$anchored_projection_downsample_mean_abs_delta - replicate_means$marker_downsample_mean_abs_delta)
  ),
  favorable_replicates = c(
    sum(replicate_means$source_gain_over_null > 0),
    sum(replicate_means$source_gain_over_no_source > 0),
    sum(replicate_means$anchored_projection_downsample_mean_abs_delta < replicate_means$marker_downsample_mean_abs_delta)
  ),
  total_replicates = N_REP,
  stringsAsFactors = FALSE
)
write.table(paired_checks, file.path(OUT, "paired_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(OUT, "incremental_utility_benchmark.png"), width = 1800, height = 800, res = 180)
par(mfrow = c(1, 3), mar = c(8, 5, 3, 1))
boxplot(
  fit_metrics[, c("null_heldout_gene_loglik_per_observation", "anchored_no_source_heldout_gene_loglik_per_observation", "anchored_source_heldout_gene_loglik_per_observation")],
  names = c("Intercept-only", "Anchored NB", "Anchored NB\n+ study"), las = 2, cex.axis = 0.85,
  ylab = "Held-out state-gene log-likelihood / observation", main = "Count-model utility"
)
cell_metrics$method_short <- factor(
  cell_metrics$method,
  levels = c("marker_scoring", "anchored_no_source", "anchored_source_aware"),
  labels = c("Marker score", "Anchored NB", "Anchored NB\n+ study")
)
boxplot(state_weight_rmse ~ method_short, data = cell_metrics, las = 2, cex.axis = 0.85,
        xlab = "", ylab = "Training state-weight RMSE", main = "Coordinate recovery")
boxplot(
  fit_metrics[, c("marker_downsample_mean_abs_delta", "anchored_projection_downsample_mean_abs_delta")],
  names = c("Marker score", "Anchored projection"), las = 2, cex.axis = 0.85,
  ylab = "Mean absolute change after 50% thinning", main = "Downsampling stability"
)
dev.off()

mean_value <- function(name) summary_table$mean[match(name, summary_table$metric)]
report <- c(
  "# Incremental utility beyond marker scoring",
  "",
  "## Design",
  "",
  sprintf("- %d independent batch-confounded simulations, each evaluated with %d random half-marker splits (%d fits in total); three studies, 90 cells, 216 genes per simulation.", N_REP, N_MARKER_SPLITS, N_REP * N_MARKER_SPLITS),
  "- All methods used the same marker-score simplex coordinate.",
  "- For every random split, half of each state's markers defined the coordinate; the complementary half was excluded from initialization and used only for held-out count prediction.",
  "- The source-aware model estimated study-by-gene effects; donor effects were intentionally omitted to keep the source decomposition identifiable.",
  "- Downsampling stability compared direct marker scoring with fixed-program projection under the same marker-logit anchor.",
  "",
  "## Main results",
  "",
  sprintf("- Intercept-only held-out-gene log-likelihood per observation: %.4f.", mean_value("null_heldout_gene_loglik_per_observation")),
  sprintf("- Anchored NB without source terms: %.4f (gain %.4f).", mean_value("anchored_no_source_heldout_gene_loglik_per_observation"), mean_value("anchored_no_source_heldout_gene_loglik_per_observation") - mean_value("null_heldout_gene_loglik_per_observation")),
  sprintf("- Anchored NB with study terms: %.4f (gain over no-source %.4f).", mean_value("anchored_source_heldout_gene_loglik_per_observation"), mean_value("source_gain_over_no_source")),
  sprintf("- Correlation between estimated and true study effects: %.3f.", mean_value("study_effect_truth_correlation")),
  sprintf("- Mean change after 50%% thinning: marker score %.4f; anchored projection %.4f.", mean_value("marker_downsample_mean_abs_delta"), mean_value("anchored_projection_downsample_mean_abs_delta")),
  sprintf("- Source-aware count prediction exceeded the intercept-only model in %d/%d replicates and the no-source model in %d/%d replicates.", paired_checks$favorable_replicates[1], N_REP, paired_checks$favorable_replicates[2], N_REP),
  sprintf("- Anchored projection was more downsampling-stable than direct marker scoring in %d/%d replicates.", paired_checks$favorable_replicates[3], N_REP),
  sprintf("- The mean within-simulation SD of the source-aware gain over the intercept model across marker splits was %.4f.", mean(marker_split_variability$source_gain_over_null)),
  "",
  "## Claim boundary",
  "",
  "This benchmark tests incremental numerical utility under a simulator with known count-generating structure. The evaluation genes are state-associated but excluded from coordinate initialization; they should therefore be described as held-out state-associated genes, not generic non-marker genes. The benchmark does not establish biological state truth in public datasets."
)
writeLines(report, file.path(OUT, "incremental_utility_report.md"))

stopifnot(
  nrow(fit_metrics) == N_REP * N_MARKER_SPLITS,
  all(is.finite(as.matrix(fit_metrics[, -c(1, 2), drop = FALSE]))),
  all(cell_metrics$drift_from_marker[cell_metrics$method != "marker_scoring"] < 0.10)
)
message("Incremental-utility benchmark written to: ", OUT)
