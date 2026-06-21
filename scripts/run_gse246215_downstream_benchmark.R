#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_gse246215_downstream_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

source_files <- c(
  "matrix_utils.R",
  "simulate_fibrodynmix.R",
  "marker_sets.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "topic_nmf_baseline.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "real_data_interface.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas"))
out_dir <- get_arg("out", file.path(ROOT, "analysis", "gse246215_downstream_benchmark"))
counts_path <- get_arg("counts", file.path(data_dir, "GSE246215_Fibroblast_counts.csv.gz"))
metadata_path <- get_arg("metadata", file.path(data_dir, "GSE246215_Fibroblast_metadata.csv.gz"))
seed <- as.integer(get_arg("seed", "246215"))
max_cells_per_cancer <- as.integer(get_arg("max-cells-per-cancer", "80"))
max_genes <- as.integer(get_arg("max-genes", "700"))
aggregation_col <- get_arg("aggregation-col", "PatientID")
n_outer <- as.integer(get_arg("n-outer", "2"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "12"))
maxit_z <- as.integer(get_arg("maxit-z", "12"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
topic_nmf_iter <- as.integer(get_arg("topic-nmf-iter", "40"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!file.exists(counts_path)) {
  stop(sprintf("Missing GSE246215 counts file: %s", counts_path), call. = FALSE)
}
if (!file.exists(metadata_path)) {
  stop(sprintf("Missing GSE246215 metadata file: %s", metadata_path), call. = FALSE)
}

metadata <- data.table::fread(metadata_path, data.table = FALSE)
metadata <- metadata[metadata$Phenotype == "Tumor", , drop = FALSE]
if (!aggregation_col %in% colnames(metadata)) {
  aggregation_col <- if ("SampleID" %in% colnames(metadata)) "SampleID" else "CancerType_short"
}
required_meta <- c("CellName", "CancerType_short", "CancerType", aggregation_col)
missing_meta <- setdiff(required_meta, colnames(metadata))
if (length(missing_meta) > 0L) {
  stop(sprintf("GSE246215 metadata is missing columns: %s", paste(missing_meta, collapse = ", ")), call. = FALSE)
}

groups <- split(metadata, metadata$CancerType_short)
groups <- groups[vapply(groups, nrow, integer(1)) > 0L]
selected <- do.call(rbind, lapply(sort(names(groups)), function(cancer) {
  group <- groups[[cancer]]
  group[sample(seq_len(nrow(group)), min(max_cells_per_cancer, nrow(group))), , drop = FALSE]
}))
selected_cells <- unique(selected$CellName)

message("Reading selected GSE246215 count matrix...")
counts_dt <- data.table::fread(
  counts_path,
  select = c("GeneName", selected_cells),
  data.table = FALSE,
  check.names = FALSE
)
gene_names <- make.unique(as.character(counts_dt$GeneName))
counts <- as.matrix(counts_dt[, -1L, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_names
colnames(counts) <- colnames(counts_dt)[-1L]

selected <- selected[match(colnames(counts), selected$CellName), , drop = FALSE]
markers <- get_fibrodynmix_markers("human", "caf")
marker_genes <- unique(unlist(markers, use.names = FALSE))
gene_rank <- order(rowSums(counts > 0), rowSums(counts), decreasing = TRUE)
keep_ranked <- rep(FALSE, nrow(counts))
keep_ranked[gene_rank[seq_len(min(max_genes, length(gene_rank)))]] <- TRUE
keep_genes <- keep_ranked | rownames(counts) %in% marker_genes
counts <- counts[keep_genes, , drop = FALSE]
markers <- lapply(markers, intersect, y = rownames(counts))
markers <- markers[vapply(markers, length, integer(1)) > 0L]

cell_metadata <- data.frame(
  cell_id = colnames(counts),
  dataset_id = selected$CancerType_short,
  study_id = selected$CancerType_short,
  donor_id = selected[[aggregation_col]],
  condition = selected$CancerType_short,
  cancer_type = selected$CancerType,
  cancer_type_short = selected$CancerType_short,
  aggregation_id = selected[[aggregation_col]],
  stringsAsFactors = FALSE
)

prepared <- prepare_fibrodynmix_data(
  counts = counts,
  cell_metadata = cell_metadata,
  marker_index = markers,
  cell_id_col = "cell_id",
  study_col = "study_id",
  donor_col = "donor_id",
  min_cells_per_gene = 0,
  min_counts_per_gene = 0,
  require_all_states = TRUE
)

message("Fitting FibroDynMix for GSE246215 downstream benchmark...")
fibro_fit <- fit_fibrodynmix_prepared(
  prepared,
  fit_study_effect = length(unique(prepared$study_id)) > 1L,
  fit_donor_effect = FALSE,
  n_outer = n_outer,
  initializer_args = list(n_iter = initializer_iter),
  study_l2 = study_l2,
  marker_l2 = marker_l2,
  maxit_beta = maxit_beta,
  maxit_z = maxit_z
)

message("Fitting marker-score and NMF baselines...")
marker_baseline <- score_marker_baseline(
  counts = prepared$counts,
  marker_index = prepared$marker_index,
  library_size = prepared$library_size
)
topic_backend <- if (requireNamespace("NMF", quietly = TRUE)) "nmf" else "multiplicative_update"
topic_fit <- fit_topic_nmf_baseline(
  counts = prepared$counts,
  marker_index = prepared$marker_index,
  n_topics = length(prepared$marker_index),
  n_iter = topic_nmf_iter,
  seed = seed,
  backend = topic_backend
)

aggregate_features <- function(cell_features, method, feature_type) {
  feature_df <- as.data.frame(cell_features, stringsAsFactors = FALSE)
  feature_df$aggregation_id <- prepared$cell_metadata$aggregation_id
  feature_df$cancer_type_short <- prepared$cell_metadata$cancer_type_short
  groups <- unique(feature_df[, c("aggregation_id", "cancer_type_short"), drop = FALSE])
  out <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    keep <- feature_df$aggregation_id == groups$aggregation_id[i] &
      feature_df$cancer_type_short == groups$cancer_type_short[i]
    values <- colMeans(as.matrix(feature_df[keep, setdiff(colnames(feature_df), c("aggregation_id", "cancer_type_short")), drop = FALSE]))
    data.frame(
      method = method,
      feature_type = feature_type,
      aggregation_col = aggregation_col,
      aggregation_id = groups$aggregation_id[i],
      cancer_type_short = groups$cancer_type_short[i],
      n_cells = sum(keep),
      t(values),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

feature_tables <- list(
  aggregate_features(fibro_fit$z_hat, "fibrodynmix_z", "state_composition"),
  aggregate_features(marker_baseline$z_pred, "marker_scoring_z", "state_composition"),
  aggregate_features(marker_baseline$scores, "marker_scoring_scores", "marker_scores"),
  aggregate_features(topic_fit$z_pred, "topic_nmf_z", "topic_composition")
)
features_all <- do.call(rbind, feature_tables)
write_tsv(features_all, file.path(out_dir, "gse246215_downstream_features.tsv"))

feature_cols_for_method <- function(df) {
  setdiff(colnames(df), c("method", "feature_type", "aggregation_col", "aggregation_id", "cancer_type_short", "n_cells"))
}

metrics <- do.call(rbind, lapply(split(features_all, features_all$method), function(df) {
  feature_cols <- feature_cols_for_method(df)
  clf <- evaluate_downstream_classification(
    features = as.matrix(df[, feature_cols, drop = FALSE]),
    labels = df$cancer_type_short,
    n_folds = min(5L, max(2L, min(table(df$cancer_type_short)))),
    seed = seed
  )
  data.frame(
    method = df$method[1L],
    feature_type = df$feature_type[1L],
    aggregation_col = aggregation_col,
    n_observations = clf$n_observations,
    n_evaluated = clf$n_evaluated,
    n_classes = clf$n_classes,
    n_folds = clf$n_folds,
    status = clf$status,
    accuracy = clf$accuracy,
    balanced_accuracy = clf$balanced_accuracy,
    macro_f1 = clf$macro_f1,
    macro_auroc = clf$macro_auroc,
    stringsAsFactors = FALSE
  )
}))
metrics <- metrics[order(-metrics$balanced_accuracy, -metrics$macro_f1, metrics$method), , drop = FALSE]
write_tsv(metrics, file.path(out_dir, "gse246215_downstream_classification_metrics.tsv"))

canonical_programs <- list(
  myofibroblast_gradient = intersect(c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1", "FAP", "PDPN"), rownames(prepared$counts)),
  ecm_gradient = intersect(c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "MMP2", "COL11A1", "MMP11"), rownames(prepared$counts)),
  inflammatory_gradient = intersect(c("IL6", "CXCL1", "CXCL8", "CXCL12", "CXCL14", "CCL2"), rownames(prepared$counts))
)
expected_gradient_state <- c(
  myofibroblast_gradient = "myofibroblast",
  ecm_gradient = "ECM-remodeling",
  inflammatory_gradient = "inflammatory"
)
normalized <- log_normalize_counts(prepared$counts, prepared$library_size, 10000)
gradient_rows <- list()
gradient_index <- 1L
for (program in names(canonical_programs)) {
  genes <- canonical_programs[[program]]
  if (length(genes) == 0L) {
    next
  }
  score <- colMeans(normalized[genes, , drop = FALSE])
  for (state in colnames(fibro_fit$z_hat)) {
    gradient_rows[[gradient_index]] <- data.frame(
      program = program,
      expected_state = expected_gradient_state[[program]],
      state = state,
      is_expected_state = identical(state, expected_gradient_state[[program]]),
      n_genes = length(genes),
      spearman_rho = suppressWarnings(stats::cor(score, fibro_fit$z_hat[, state], method = "spearman")),
      stringsAsFactors = FALSE
    )
    gradient_index <- gradient_index + 1L
  }
}
gradient_validation <- do.call(rbind, gradient_rows)
gradient_validation <- gradient_validation[order(gradient_validation$program, !gradient_validation$is_expected_state, -gradient_validation$spearman_rho), , drop = FALSE]
write_tsv(gradient_validation, file.path(out_dir, "gse246215_biological_gradient_validation.tsv"))

selection_summary <- aggregate(
  CellName ~ CancerType_short,
  data = selected,
  FUN = length
)
colnames(selection_summary) <- c("cancer_type_short", "n_selected_cells")
selection_summary$n_aggregation_units <- vapply(
  selection_summary$cancer_type_short,
  function(x) length(unique(selected[selected$CancerType_short == x, aggregation_col])),
  integer(1)
)
write_tsv(selection_summary, file.path(out_dir, "gse246215_downstream_selection_summary.tsv"))
write_tsv(prepared$marker_summary, file.path(out_dir, "gse246215_downstream_marker_coverage.tsv"))

manifest <- data.frame(
  analysis = "gse246215_downstream_benchmark",
  primary_claim = "GSE246215 cancer-type labels are predicted from patient/sample-level fibroblast-state representations, comparing FibroDynMix z, marker-scoring features, and NMF topics.",
  claim_boundary = "Processed public tumor-fibroblast count matrices from one GEO accession; classification evaluates representation utility, not diagnostic performance. Biological-gradient validation is marker-prior-oriented and should be interpreted as a case study rather than independent marker discovery.",
  source_accession = "GSE246215",
  aggregation_col = aggregation_col,
  n_cells = ncol(prepared$counts),
  n_genes = nrow(prepared$counts),
  n_states = length(prepared$marker_index),
  n_aggregation_units = length(unique(prepared$cell_metadata$aggregation_id)),
  n_cancer_types = length(unique(prepared$cell_metadata$cancer_type_short)),
  topic_backend = topic_backend,
  best_method_by_balanced_accuracy = metrics$method[which.max(metrics$balanced_accuracy)],
  best_balanced_accuracy = max(metrics$balanced_accuracy, na.rm = TRUE),
  fibrodynmix_balanced_accuracy = metrics$balanced_accuracy[metrics$method == "fibrodynmix_z"],
  marker_score_balanced_accuracy = metrics$balanced_accuracy[metrics$method == "marker_scoring_scores"],
  topic_nmf_balanced_accuracy = metrics$balanced_accuracy[metrics$method == "topic_nmf_z"],
  myofibroblast_gradient_expected_state_spearman = gradient_validation$spearman_rho[gradient_validation$program == "myofibroblast_gradient" & gradient_validation$is_expected_state][1L],
  ecm_gradient_expected_state_spearman = gradient_validation$spearman_rho[gradient_validation$program == "ecm_gradient" & gradient_validation$is_expected_state][1L],
  inflammatory_gradient_expected_state_spearman = gradient_validation$spearman_rho[gradient_validation$program == "inflammatory_gradient" & gradient_validation$is_expected_state][1L],
  seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "gse246215_downstream_manifest.tsv"))

message("GSE246215 downstream benchmark written to: ", normalizePath(out_dir))
