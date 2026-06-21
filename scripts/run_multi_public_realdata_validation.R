#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

has_flag <- function(name) {
  paste0("--", name) %in% args
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_multi_public_realdata_validation.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

source_package_files <- function(root) {
  r_files <- c(
    "simulate_fibrodynmix.R",
    "baseline_marker_scoring.R",
    "benchmark_metrics.R",
    "fibrodynmix_initializer.R",
    "nb_likelihood.R",
    "fit_nb_model.R",
    "transition_flow.R",
    "real_data_interface.R",
    "cross_cohort_transfer.R"
  )
  invisible(lapply(file.path(root, "R", r_files), source))
}

source_package_files(ROOT)

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_dryad_breast_fibroblast"))
dataset_manifest <- get_arg("dataset-manifest")
out_dir <- get_arg("out", file.path(ROOT, "analysis", "multi_public_realdata_validation"))
max_cells <- as.integer(get_arg("max-cells", "120"))
max_genes <- as.integer(get_arg("max-genes", "700"))
seed <- as.integer(get_arg("seed", "260606"))
n_outer <- as.integer(get_arg("n-outer", "2"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "15"))
maxit_z <- as.integer(get_arg("maxit-z", "15"))
transfer_maxit_z <- as.integer(get_arg("transfer-maxit-z", "80"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
donor_l2 <- as.numeric(get_arg("donor-l2", "0.1"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

default_public_registry <- function() {
  data.frame(
    dataset_id = c("dryad_mt3_caf", "dryad_normal_mammary_fibroblast"),
    study_id = c("Dryad_MT3_CAF", "Dryad_normal_mammary"),
    donor_id = c("MT3_CAF_pool", "normal_mammary_pool"),
    condition = c("disease", "normal"),
    organism = c("Mus musculus", "Mus musculus"),
    tissue = c("breast_tumor", "normal_mammary_gland"),
    counts_path = file.path(data_dir, c("MT3_CAFs_raw.txt", "Normal_mammary_fibroblasts_raw.txt")),
    public_record = "https://datadryad.org/dataset/doi:10.6071/M3238R",
    download_url = c(
      "https://zenodo.org/api/records/3977255/files/MT3_CAFs_raw.txt/content",
      "https://zenodo.org/api/records/3977255/files/Normal_mammary_fibroblasts_raw.txt/content"
    ),
    stringsAsFactors = FALSE
  )
}

read_registry <- function(path) {
  if (is.null(path)) {
    registry <- default_public_registry()
  } else {
    registry <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  required <- c("dataset_id", "study_id", "donor_id", "condition", "counts_path")
  missing <- setdiff(required, colnames(registry))
  if (length(missing) > 0L) {
    stop(sprintf("Dataset manifest is missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  registry$counts_path <- ifelse(
    grepl("^/", registry$counts_path),
    registry$counts_path,
    file.path(ROOT, registry$counts_path)
  )
  registry
}

read_table_matrix <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing public count matrix: %s", path), call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
    mat <- as.matrix(x)
    storage.mode(mat) <- "numeric"
    return(mat)
  }
  if (ext %in% c("txt", "tsv")) {
    tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext == "csv") {
    tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop(sprintf("Unsupported count format in multi-public validation: %s", path), call. = FALSE)
  }
  default_rownames <- identical(rownames(tab), as.character(seq_len(nrow(tab))))
  if (default_rownames) {
    if (ncol(tab) < 2L) {
      stop("Count table must contain one gene column and at least one cell column.", call. = FALSE)
    }
    genes <- tab[[1L]]
    mat <- as.matrix(tab[, -1L, drop = FALSE])
    rownames(mat) <- make.unique(as.character(genes))
  } else {
    mat <- as.matrix(tab)
    rownames(mat) <- make.unique(rownames(tab))
  }
  storage.mode(mat) <- "numeric"
  mat
}

default_fibroblast_markers <- function(gene_names) {
  human <- list(
    resident = c("DCN", "LUM", "COL14A1", "PDGFRA", "PI16"),
    inflammatory = c("IL6", "CXCL12", "CXCL14", "CCL2", "CXCL2"),
    myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1"),
    `ECM-remodeling` = c("COL1A1", "COL1A2", "FN1", "POSTN", "MMP2"),
    `antigen-presenting` = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"),
    `IFN-stress` = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1")
  )
  mouse <- list(
    resident = c("Dcn", "Lum", "Col14a1", "Pdgfra", "Pi16"),
    inflammatory = c("Il6", "Cxcl12", "Cxcl14", "Ccl2", "Cxcl2"),
    myofibroblast = c("Acta2", "Tagln", "Myl9", "Tpm2", "Cnn1"),
    `ECM-remodeling` = c("Col1a1", "Col1a2", "Fn1", "Postn", "Mmp2"),
    `antigen-presenting` = c("H2-Aa", "H2-Ab1", "Cd74", "H2-Eb1", "H2-DMa"),
    `IFN-stress` = c("Isg15", "Ifit1", "Ifit3", "Mx1", "Oas1")
  )
  human_hits <- sum(unique(unlist(human, use.names = FALSE)) %in% gene_names)
  mouse_hits <- sum(unique(unlist(mouse, use.names = FALSE)) %in% gene_names)
  if (mouse_hits > human_hits) mouse else human
}

subset_cells <- function(counts, metadata, max_cells) {
  if (ncol(counts) <= max_cells) {
    return(list(counts = counts, metadata = metadata))
  }
  selected <- sample(colnames(counts), max_cells)
  list(
    counts = counts[, selected, drop = FALSE],
    metadata = metadata[match(selected, metadata$cell_id), , drop = FALSE]
  )
}

simplex_entropy <- function(z) {
  z_safe <- pmax(z, .Machine$double.eps)
  -rowSums(z_safe * log(z_safe))
}

read_dataset <- function(row) {
  mat <- read_table_matrix(row$counts_path)
  colnames(mat) <- paste(row$study_id, colnames(mat), sep = "_")
  metadata <- data.frame(
    cell_id = colnames(mat),
    dataset_id = row$dataset_id,
    study_id = row$study_id,
    donor_id = row$donor_id,
    condition = row$condition,
    organism = if ("organism" %in% names(row)) row$organism else NA_character_,
    tissue = if ("tissue" %in% names(row)) row$tissue else NA_character_,
    stringsAsFactors = FALSE
  )
  subset <- subset_cells(mat, metadata, max_cells)
  list(counts = subset$counts, metadata = subset$metadata)
}

fit_prepared_fast <- function(prepared, fit_study_effect = TRUE, fit_donor_effect = FALSE) {
  fit_fibrodynmix_prepared(
    prepared,
    fit_study_effect = fit_study_effect,
    fit_donor_effect = fit_donor_effect,
    n_outer = n_outer,
    initializer_args = list(n_iter = initializer_iter),
    study_l2 = study_l2,
    donor_l2 = donor_l2,
    marker_l2 = marker_l2,
    maxit_beta = maxit_beta,
    maxit_z = maxit_z
  )
}

registry <- read_registry(dataset_manifest)
if (nrow(registry) < 2L) {
  stop("Multi-public validation requires at least two dataset rows.", call. = FALSE)
}
registry_source <- if (is.null(dataset_manifest)) {
  "default_dryad_breast_fibroblast_registry"
} else {
  dataset_manifest
}
registry_public_records <- if ("public_record" %in% colnames(registry)) {
  paste(unique(registry$public_record), collapse = ";")
} else {
  NA_character_
}
registry_organisms <- if ("organism" %in% colnames(registry)) {
  paste(unique(registry$organism), collapse = ";")
} else {
  NA_character_
}

parts <- lapply(seq_len(nrow(registry)), function(i) read_dataset(registry[i, , drop = FALSE]))
names(parts) <- registry$dataset_id

common_genes <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
if (length(common_genes) == 0L) {
  stop("No shared genes across public datasets.", call. = FALSE)
}

combined_counts <- do.call(cbind, lapply(parts, function(x) x$counts[common_genes, , drop = FALSE]))
combined_metadata <- do.call(rbind, lapply(parts, function(x) x$metadata))
markers <- default_fibroblast_markers(common_genes)
marker_genes <- unique(unlist(markers, use.names = FALSE))

gene_rank <- order(rowSums(combined_counts > 0), rowSums(combined_counts), decreasing = TRUE)
keep_ranked <- rep(FALSE, nrow(combined_counts))
keep_ranked[gene_rank[seq_len(min(max_genes, length(gene_rank)))]] <- TRUE
keep_genes <- keep_ranked | rownames(combined_counts) %in% marker_genes
combined_counts <- combined_counts[keep_genes, , drop = FALSE]
markers <- default_fibroblast_markers(rownames(combined_counts))

prepared <- prepare_fibrodynmix_data(
  counts = combined_counts,
  cell_metadata = combined_metadata,
  marker_index = markers,
  cell_id_col = "cell_id",
  study_col = "study_id",
  donor_col = "donor_id",
  min_cells_per_gene = 0,
  min_counts_per_gene = 0,
  require_all_states = TRUE
)

pooled_fit <- fit_prepared_fast(
  prepared,
  fit_study_effect = !has_flag("no-study-effect"),
  fit_donor_effect = has_flag("fit-donor-effect")
)

state_cols <- colnames(pooled_fit$z_hat)
cell_weights <- as.data.frame(pooled_fit$z_hat, stringsAsFactors = FALSE)
cell_weights$cell_id <- rownames(pooled_fit$z_hat)
cell_weights$dataset_id <- prepared$cell_metadata$dataset_id
cell_weights$study_id <- prepared$cell_metadata$study_id
cell_weights$condition <- prepared$cell_metadata$condition
cell_weights$dominant_state <- state_cols[max.col(pooled_fit$z_hat, ties.method = "first")]
cell_weights$entropy <- simplex_entropy(pooled_fit$z_hat)
cell_weights <- cell_weights[, c("dataset_id", "study_id", "condition", "cell_id", "dominant_state", "entropy", state_cols), drop = FALSE]
write_tsv(cell_weights, file.path(out_dir, "multi_public_cell_state_weights.tsv"))

dataset_groups <- unique(prepared$cell_metadata[, c("dataset_id", "study_id", "donor_id", "condition", "organism", "tissue"), drop = FALSE])
composition <- do.call(rbind, lapply(seq_len(nrow(dataset_groups)), function(i) {
  idx <- prepared$cell_metadata$dataset_id == dataset_groups$dataset_id[i]
  data.frame(
    dataset_id = dataset_groups$dataset_id[i],
    study_id = dataset_groups$study_id[i],
    donor_id = dataset_groups$donor_id[i],
    condition = dataset_groups$condition[i],
    organism = dataset_groups$organism[i],
    tissue = dataset_groups$tissue[i],
    state = state_cols,
    composition = colMeans(pooled_fit$z_hat[idx, , drop = FALSE]),
    mean_entropy = mean(simplex_entropy(pooled_fit$z_hat[idx, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
}))
write_tsv(composition, file.path(out_dir, "multi_public_dataset_state_composition.tsv"))

dataset_summary <- do.call(rbind, lapply(seq_len(nrow(dataset_groups)), function(i) {
  idx <- prepared$cell_metadata$dataset_id == dataset_groups$dataset_id[i]
  data.frame(
    dataset_id = dataset_groups$dataset_id[i],
    condition = dataset_groups$condition[i],
    n_cells = sum(idx),
    mean_library_size = mean(prepared$library_size[idx]),
    median_library_size = stats::median(prepared$library_size[idx]),
    mean_state_entropy = mean(simplex_entropy(pooled_fit$z_hat[idx, , drop = FALSE])),
    dominant_state_fraction = max(colMeans(pooled_fit$z_hat[idx, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
}))
write_tsv(dataset_summary, file.path(out_dir, "multi_public_dataset_summary.tsv"))

transition_summary <- data.frame(status = "not_run", reason = "Need normal and disease condition labels.", stringsAsFactors = FALSE)
if (all(c("normal", "disease") %in% prepared$cell_metadata$condition)) {
  normal <- colMeans(pooled_fit$z_hat[prepared$cell_metadata$condition == "normal", , drop = FALSE])
  disease <- colMeans(pooled_fit$z_hat[prepared$cell_metadata$condition == "disease", , drop = FALSE])
  cost <- compute_state_cost(pooled_fit$beta_hat)
  flow <- estimate_transition_flow(normal, disease, cost, lambda = 0.5)
  flow_long <- as.data.frame(as.table(flow$flow), stringsAsFactors = FALSE)
  colnames(flow_long) <- c("source_state", "target_state", "flow")
  write_tsv(flow_long, file.path(out_dir, "multi_public_transition_flow.tsv"))
  transition_summary <- data.frame(
    status = "ok",
    normal_label = "normal",
    disease_label = "disease",
    expected_cost = flow$expected_cost,
    entropy = flow$entropy,
    converged = flow$converged,
    iterations = flow$iterations,
    stringsAsFactors = FALSE
  )
}
write_tsv(transition_summary, file.path(out_dir, "multi_public_transition_summary.tsv"))

transfer_rows <- list()
transfer_composition_rows <- list()
for (heldout_id in registry$dataset_id) {
  train_ids <- setdiff(registry$dataset_id, heldout_id)
  train_cells <- prepared$cell_metadata$dataset_id %in% train_ids
  heldout_cells <- prepared$cell_metadata$dataset_id == heldout_id
  train_prepared <- prepare_fibrodynmix_data(
    counts = prepared$counts[, train_cells, drop = FALSE],
    cell_metadata = prepared$cell_metadata[train_cells, , drop = FALSE],
    marker_index = markers,
    cell_id_col = "cell_id",
    study_col = "study_id",
    donor_col = "donor_id",
    min_cells_per_gene = 0,
    min_counts_per_gene = 0,
    require_all_states = TRUE
  )
  train_fit <- fit_prepared_fast(train_prepared, fit_study_effect = length(unique(train_prepared$study_id)) > 1L, fit_donor_effect = FALSE)
  transfer <- fit_fibrodynmix_transfer(
    counts = prepared$counts[, heldout_cells, drop = FALSE],
    fit = train_fit,
    library_size = prepared$library_size[heldout_cells],
    maxit_z = transfer_maxit_z
  )
  heldout_meta <- unique(prepared$cell_metadata[heldout_cells, c("dataset_id", "condition", "study_id"), drop = FALSE])
  transfer_rows[[length(transfer_rows) + 1L]] <- data.frame(
    heldout_dataset_id = heldout_id,
    heldout_condition = heldout_meta$condition[1L],
    training_dataset_ids = paste(train_ids, collapse = ";"),
    n_train_cells = ncol(train_prepared$counts),
    n_heldout_cells = sum(heldout_cells),
    n_train_datasets = length(train_ids),
    n_shared_genes = length(transfer$shared_genes),
    train_initial_objective = train_fit$nb_objective_trace[1L],
    train_best_objective = train_fit$best_objective,
    train_stop_reason = train_fit$stop_reason,
    transfer_nb_objective = transfer$nb_objective,
    transfer_nb_loglik = transfer$nb_loglik,
    transfer_converged = transfer$converged,
    transfer_z_convergence_rate = transfer$z_convergence_rate,
    transfer_n_nonconverged_cells = transfer$n_nonconverged_cells,
    mean_transfer_entropy = mean(simplex_entropy(transfer$z_hat)),
    stringsAsFactors = FALSE
  )
  transfer_composition_rows[[length(transfer_composition_rows) + 1L]] <- data.frame(
    heldout_dataset_id = heldout_id,
    heldout_condition = heldout_meta$condition[1L],
    state = colnames(transfer$z_hat),
    transferred_composition = colMeans(transfer$z_hat),
    pooled_fit_composition = colMeans(pooled_fit$z_hat[heldout_cells, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}
transfer_diagnostics <- do.call(rbind, transfer_rows)
write_tsv(transfer_diagnostics, file.path(out_dir, "multi_public_transfer_diagnostics.tsv"))
write_tsv(do.call(rbind, transfer_composition_rows), file.path(out_dir, "multi_public_transfer_state_composition.tsv"))

fit_diagnostics <- data.frame(
  n_datasets = length(unique(prepared$cell_metadata$dataset_id)),
  n_conditions = length(unique(prepared$cell_metadata$condition)),
  n_cells = ncol(prepared$counts),
  n_genes = nrow(prepared$counts),
  n_states = length(prepared$marker_index),
  initial_objective = pooled_fit$nb_objective_trace[1L],
  best_objective = pooled_fit$best_objective,
  final_objective = utils::tail(pooled_fit$nb_objective_trace, 1),
  objective_improvement = pooled_fit$nb_objective_trace[1L] - pooled_fit$best_objective,
  best_iteration = pooled_fit$best_iteration,
  executed_iterations = pooled_fit$executed_iterations,
  stop_reason = pooled_fit$stop_reason,
  fit_study_effect = !is.null(pooled_fit$study_effect),
  fit_donor_effect = !is.null(pooled_fit$donor_effect),
  mean_transfer_z_convergence_rate = mean(transfer_diagnostics$transfer_z_convergence_rate),
  min_transfer_z_convergence_rate = min(transfer_diagnostics$transfer_z_convergence_rate),
  stringsAsFactors = FALSE
)
write_tsv(fit_diagnostics, file.path(out_dir, "multi_public_nb_fit_diagnostics.tsv"))
write_tsv(prepared$marker_summary, file.path(out_dir, "multi_public_marker_coverage.tsv"))
write_tsv(registry, file.path(out_dir, "multi_public_dataset_registry.tsv"))

manifest <- data.frame(
  analysis = "multi_public_realdata_validation",
  primary_claim = "FibroDynMix runs a pooled raw-count NB fit with study effects and leave-dataset-out transfer diagnostics across multiple public fibroblast count datasets.",
  claim_boundary = if (is.null(dataset_manifest)) {
    "Default registry uses two public mouse breast fibroblast count matrices from the Dryad record doi:10.6071/M3238R. This strengthens real-data validation beyond a single pooled smoke run but is not yet an independent multi-study human atlas validation."
  } else {
    "User-supplied public dataset registry validation. Interpret biological claims according to the registry source, organism, sample design, and whether matrices are independent studies or condition/cancer-type subsets."
  },
  registry_source = registry_source,
  public_records = registry_public_records,
  organisms = registry_organisms,
  n_public_datasets = length(unique(prepared$cell_metadata$dataset_id)),
  n_conditions = length(unique(prepared$cell_metadata$condition)),
  n_cells = ncol(prepared$counts),
  n_genes = nrow(prepared$counts),
  mean_transfer_z_convergence_rate = mean(transfer_diagnostics$transfer_z_convergence_rate),
  seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "multi_public_validation_manifest.tsv"))

message("Multi-public real-data validation written to: ", normalizePath(out_dir))
