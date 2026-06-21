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

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_gse167339_donor_robustness.R"
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

dataset_manifest <- get_arg("dataset-manifest", file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_human_fibroblast_dataset_manifest.tsv"))
data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast"))
validation_dir <- get_arg("validation-dir", file.path(ROOT, "analysis", "independent_geo_gse167339_validation"))
out_dir <- get_arg("out", file.path(ROOT, "analysis", "gse167339_donor_robustness"))
max_genes <- as.integer(get_arg("max-genes", "700"))
seed <- as.integer(get_arg("seed", "167339"))
n_outer <- as.integer(get_arg("n-outer", "2"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "15"))
maxit_z <- as.integer(get_arg("maxit-z", "15"))
transfer_maxit_z <- as.integer(get_arg("transfer-maxit-z", "120"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

simplex_entropy <- function(z) {
  z_safe <- pmax(z, .Machine$double.eps)
  -rowSums(z_safe * log(z_safe))
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

read_registry <- function(path) {
  registry <- read_tsv(path)
  required <- c("dataset_id", "study_id", "donor_id", "condition", "counts_path")
  missing <- setdiff(required, colnames(registry))
  if (length(missing) > 0L) {
    stop(sprintf("Dataset manifest is missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  registry$counts_path <- ifelse(grepl("^/", registry$counts_path), registry$counts_path, file.path(ROOT, registry$counts_path))
  registry
}

read_dataset <- function(row) {
  mat <- readRDS(row$counts_path)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  colnames(mat) <- paste(row$study_id, row$dataset_id, colnames(mat), sep = "_")
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
  list(counts = mat, metadata = metadata)
}

prepare_combined <- function(registry) {
  parts <- lapply(seq_len(nrow(registry)), function(i) read_dataset(registry[i, , drop = FALSE]))
  common_genes <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
  if (length(common_genes) == 0L) {
    stop("No shared genes across GSE167339 datasets.", call. = FALSE)
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
  prepare_fibrodynmix_data(
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
}

fit_prepared_fast <- function(prepared, fit_study_effect = TRUE) {
  fit_fibrodynmix_prepared(
    prepared,
    fit_study_effect = fit_study_effect,
    fit_donor_effect = FALSE,
    n_outer = n_outer,
    initializer_args = list(n_iter = initializer_iter),
    study_l2 = study_l2,
    marker_l2 = marker_l2,
    maxit_beta = maxit_beta,
    maxit_z = maxit_z
  )
}

registry <- read_registry(dataset_manifest)
prepared <- prepare_combined(registry)
state_cols <- names(prepared$marker_index)

cell_weights_path <- file.path(validation_dir, "multi_public_cell_state_weights.tsv")
if (!file.exists(cell_weights_path)) {
  stop(sprintf("Missing GSE167339 validation cell weights: %s", cell_weights_path), call. = FALSE)
}
cell_weights <- read_tsv(cell_weights_path)
cell_weights <- merge(
  cell_weights,
  unique(registry[, c("dataset_id", "donor_id", "condition"), drop = FALSE]),
  by = c("dataset_id", "condition"),
  all.x = TRUE,
  sort = FALSE
)

donor_condition_groups <- unique(cell_weights[, c("donor_id", "condition"), drop = FALSE])
donor_comp <- do.call(rbind, lapply(seq_len(nrow(donor_condition_groups)), function(i) {
  idx <- cell_weights$donor_id == donor_condition_groups$donor_id[i] & cell_weights$condition == donor_condition_groups$condition[i]
  z <- as.matrix(cell_weights[idx, state_cols, drop = FALSE])
  data.frame(
    donor_id = donor_condition_groups$donor_id[i],
    condition = donor_condition_groups$condition[i],
    n_cells = sum(idx),
    state = state_cols,
    composition = colMeans(z),
    mean_entropy = mean(cell_weights$entropy[idx]),
    plasticity_index = mean(cell_weights$entropy[idx]),
    stringsAsFactors = FALSE
  )
}))
write_tsv(donor_comp, file.path(out_dir, "gse167339_donor_state_composition.tsv"))

donor_summary <- do.call(rbind, lapply(split(donor_comp, donor_comp$donor_id), function(x) {
  x_wide <- reshape(
    x[, c("condition", "state", "composition")],
    idvar = "state",
    timevar = "condition",
    direction = "wide"
  )
  normal_col <- "composition.normal"
  disease_col <- "composition.disease"
  hash_col <- "composition.hash_unknown"
  data.frame(
    donor_id = unique(x$donor_id),
    n_conditions = length(unique(x$condition)),
    n_cells = sum(unique(x[, c("condition", "n_cells")])$n_cells),
    mean_entropy = mean(x$mean_entropy),
    has_normal = normal_col %in% colnames(x_wide),
    has_disease = disease_col %in% colnames(x_wide),
    has_hash_unknown = hash_col %in% colnames(x_wide),
    disease_minus_normal_l1 = if (all(c(normal_col, disease_col) %in% colnames(x_wide))) {
      sum(abs(x_wide[[disease_col]] - x_wide[[normal_col]]))
    } else {
      NA_real_
    },
    hash_minus_known_l1 = if (hash_col %in% colnames(x_wide) && any(c(normal_col, disease_col) %in% colnames(x_wide))) {
      known_cols <- intersect(c(normal_col, disease_col), colnames(x_wide))
      known_mean <- rowMeans(x_wide[, known_cols, drop = FALSE])
      sum(abs(x_wide[[hash_col]] - known_mean))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
write_tsv(donor_summary, file.path(out_dir, "gse167339_donor_robustness_summary.tsv"))

leave_donor_rows <- list()
leave_donor_comp_rows <- list()
for (heldout_donor in sort(unique(prepared$cell_metadata$donor_id))) {
  train_cells <- prepared$cell_metadata$donor_id != heldout_donor
  heldout_cells <- prepared$cell_metadata$donor_id == heldout_donor
  train_prepared <- prepare_fibrodynmix_data(
    counts = prepared$counts[, train_cells, drop = FALSE],
    cell_metadata = prepared$cell_metadata[train_cells, , drop = FALSE],
    marker_index = prepared$marker_index,
    cell_id_col = "cell_id",
    study_col = "study_id",
    donor_col = "donor_id",
    min_cells_per_gene = 0,
    min_counts_per_gene = 0,
    require_all_states = TRUE
  )
  train_fit <- fit_prepared_fast(train_prepared, fit_study_effect = length(unique(train_prepared$study_id)) > 1L)
  transfer <- fit_fibrodynmix_transfer(
    counts = prepared$counts[, heldout_cells, drop = FALSE],
    fit = train_fit,
    library_size = prepared$library_size[heldout_cells],
    maxit_z = transfer_maxit_z
  )
  heldout_meta <- prepared$cell_metadata[heldout_cells, , drop = FALSE]
  leave_donor_rows[[length(leave_donor_rows) + 1L]] <- data.frame(
    heldout_donor_id = heldout_donor,
    heldout_conditions = paste(sort(unique(heldout_meta$condition)), collapse = ";"),
    training_donor_ids = paste(sort(unique(train_prepared$donor_id)), collapse = ";"),
    n_train_cells = ncol(train_prepared$counts),
    n_heldout_cells = sum(heldout_cells),
    n_train_donors = length(unique(train_prepared$donor_id)),
    n_heldout_datasets = length(unique(heldout_meta$dataset_id)),
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
  leave_donor_comp_rows[[length(leave_donor_comp_rows) + 1L]] <- data.frame(
    heldout_donor_id = heldout_donor,
    state = colnames(transfer$z_hat),
    transferred_composition = colMeans(transfer$z_hat),
    stringsAsFactors = FALSE
  )
}
leave_donor_transfer <- do.call(rbind, leave_donor_rows)
write_tsv(leave_donor_transfer, file.path(out_dir, "gse167339_leave_donor_out_transfer.tsv"))
write_tsv(do.call(rbind, leave_donor_comp_rows), file.path(out_dir, "gse167339_leave_donor_out_state_composition.tsv"))

hash_matrix_path <- file.path(data_dir, "GSM5102538_H3-hash_matrix.mtx.gz")
hash_features_path <- file.path(data_dir, "GSM5102538_H3-hash_features.tsv.gz")
hash_sensitivity <- data.frame()
if (file.exists(hash_matrix_path) && file.exists(hash_features_path)) {
  hash_mat <- Matrix::readMM(hash_matrix_path)
  hash_features <- utils::read.delim(hash_features_path, header = FALSE, stringsAsFactors = FALSE)
  feature_type <- if (ncol(hash_features) >= 3L) hash_features[[3L]] else rep("Gene Expression", nrow(hash_features))
  antibody_rows <- feature_type == "Antibody Capture"
  antibody_mat <- as.matrix(hash_mat[antibody_rows, , drop = FALSE])
  antibody_names <- hash_features[[2L]][antibody_rows]
  rownames(antibody_mat) <- antibody_names
  top_hash_index <- max.col(t(antibody_mat), ties.method = "first")
  top_hash <- antibody_names[top_hash_index]
  top_count <- apply(antibody_mat, 2L, max)
  second_count <- apply(antibody_mat, 2L, function(x) sort(x, decreasing = TRUE)[2L])
  ratio <- top_count / pmax(second_count, 1)
  rows <- list()
  idx <- 1L
  for (min_count in c(5, 10, 20, 50)) {
    for (min_ratio in c(1.5, 2, 3)) {
      assigned <- top_count >= min_count & ratio >= min_ratio
      tab <- table(factor(top_hash[assigned], levels = antibody_names))
      rows[[idx]] <- data.frame(
        hash_min_count = min_count,
        hash_min_ratio = min_ratio,
        n_assigned_cells = sum(assigned),
        n_hash_groups_ge_40 = sum(tab >= 40),
        largest_hash_group_cells = if (length(tab)) max(tab) else 0L,
        smallest_retained_hash_group_cells = if (any(tab >= 40)) min(tab[tab >= 40]) else 0L,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  hash_sensitivity <- do.call(rbind, rows)
}
write_tsv(hash_sensitivity, file.path(out_dir, "gse167339_hash_threshold_sensitivity.tsv"))

manifest <- data.frame(
  analysis = "gse167339_donor_robustness",
  primary_claim = "FibroDynMix supports donor-level composition summaries, leave-donor-out transfer, and hash-demultiplexing sensitivity for the independent GSE167339 human fibroblast validation.",
  claim_boundary = "Human3 hash groups are included as hash_unknown donor-level robustness evidence because public files do not map HumanHashTag IDs to treatment labels; transition and treatment biology should be interpreted using Human1/Human2 labels only.",
  source_accession = "GSE167339",
  n_donors = length(unique(prepared$cell_metadata$donor_id)),
  n_datasets = length(unique(prepared$cell_metadata$dataset_id)),
  n_cells = ncol(prepared$counts),
  n_genes = nrow(prepared$counts),
  n_leave_donor_out_runs = nrow(leave_donor_transfer),
  mean_leave_donor_transfer_convergence = mean(leave_donor_transfer$transfer_z_convergence_rate),
  min_leave_donor_transfer_convergence = min(leave_donor_transfer$transfer_z_convergence_rate),
  n_hash_threshold_settings = nrow(hash_sensitivity),
  seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "gse167339_donor_robustness_manifest.tsv"))

message("GSE167339 donor robustness analysis written to: ", normalizePath(out_dir))
