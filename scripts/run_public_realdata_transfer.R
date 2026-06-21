#!/usr/bin/env Rscript

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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_public_realdata_transfer.R"
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
out_dir <- get_arg("out", file.path(ROOT, "analysis", "public_realdata_transfer"))
max_cells <- as.integer(get_arg("max-cells", "160"))
max_genes <- as.integer(get_arg("max-genes", "700"))
seed <- as.integer(get_arg("seed", "1207"))
n_outer <- as.integer(get_arg("n-outer", "2"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "15"))
maxit_z <- as.integer(get_arg("maxit-z", "15"))
transfer_maxit_z <- as.integer(get_arg("transfer-maxit-z", "100"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

dataset <- data.frame(
  file = c("MT3_CAFs_raw.txt", "Normal_mammary_fibroblasts_raw.txt"),
  condition = c("disease", "normal"),
  study_id = c("Dryad_MT3_CAF", "Dryad_normal_mammary"),
  donor_id = c("MT3_CAF_pool", "normal_mammary_pool"),
  stringsAsFactors = FALSE
)

resolve_data_path <- function(file) {
  primary <- file.path(data_dir, file)
  fallback <- file.path(ROOT, file)
  if (file.exists(primary)) {
    return(primary)
  }
  if (file.exists(fallback)) {
    return(fallback)
  }
  stop(sprintf("Missing public raw-count file: %s", file), call. = FALSE)
}

read_one <- function(path, condition, study_id, donor_id) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  default_rownames <- identical(rownames(tab), as.character(seq_len(nrow(tab))))
  if (default_rownames) {
    gene <- tab[[1L]]
    mat <- as.matrix(tab[, -1L, drop = FALSE])
    rownames(mat) <- make.unique(as.character(gene))
  } else {
    mat <- as.matrix(tab)
    rownames(mat) <- make.unique(rownames(tab))
  }
  storage.mode(mat) <- "numeric"
  colnames(mat) <- paste(study_id, colnames(mat), sep = "_")
  metadata <- data.frame(
    cell_id = colnames(mat),
    condition = condition,
    study_id = study_id,
    donor_id = donor_id,
    stringsAsFactors = FALSE
  )
  list(counts = mat, metadata = metadata, source_path = path)
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

subset_cells <- function(part, max_cells) {
  if (ncol(part$counts) <= max_cells) {
    return(part)
  }
  selected <- sample(colnames(part$counts), max_cells)
  part$counts <- part$counts[, selected, drop = FALSE]
  part$metadata <- part$metadata[match(selected, part$metadata$cell_id), , drop = FALSE]
  part
}

parts <- lapply(seq_len(nrow(dataset)), function(i) {
  path <- resolve_data_path(dataset$file[i])
  subset_cells(
    read_one(
      path = path,
      condition = dataset$condition[i],
      study_id = dataset$study_id[i],
      donor_id = dataset$donor_id[i]
    ),
    max_cells = max_cells
  )
})
names(parts) <- dataset$condition

common_genes <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
combined_counts <- do.call(cbind, lapply(parts, function(x) x$counts[common_genes, , drop = FALSE]))
markers <- default_fibroblast_markers(common_genes)
marker_genes <- unique(unlist(markers, use.names = FALSE))
gene_rank <- order(rowSums(combined_counts > 0), rowSums(combined_counts), decreasing = TRUE)
keep_ranked <- rep(FALSE, nrow(combined_counts))
keep_ranked[gene_rank[seq_len(min(max_genes, length(gene_rank)))]] <- TRUE
keep_genes <- keep_ranked | rownames(combined_counts) %in% marker_genes
common_genes <- rownames(combined_counts)[keep_genes]

parts <- lapply(parts, function(part) {
  part$counts <- part$counts[common_genes, , drop = FALSE]
  part
})

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

simplex_entropy <- function(z) {
  z_safe <- pmax(z, .Machine$double.eps)
  -rowSums(z_safe * log(z_safe))
}

fit_one_condition <- function(part) {
  prepared <- prepare_fibrodynmix_data(
    counts = part$counts,
    cell_metadata = part$metadata,
    marker_index = markers,
    cell_id_col = "cell_id",
    study_col = "study_id",
    donor_col = "donor_id",
    min_cells_per_gene = 0,
    min_counts_per_gene = 0,
    require_all_states = TRUE
  )
  fit <- fit_fibrodynmix_prepared(
    prepared,
    fit_study_effect = FALSE,
    fit_donor_effect = FALSE,
    n_outer = n_outer,
    initializer_args = list(n_iter = initializer_iter),
    marker_l2 = 0.05,
    maxit_beta = maxit_beta,
    maxit_z = maxit_z
  )
  list(prepared = prepared, fit = fit)
}

fits <- lapply(parts, fit_one_condition)

directions <- data.frame(
  train_condition = c("disease", "normal"),
  heldout_condition = c("normal", "disease"),
  stringsAsFactors = FALSE
)

direction_rows <- list()
composition_rows <- list()
cell_rows <- list()
marker_rows <- list()
row_id <- 1L
composition_id <- 1L
cell_id <- 1L
marker_id <- 1L

for (i in seq_len(nrow(directions))) {
  train_condition <- directions$train_condition[i]
  heldout_condition <- directions$heldout_condition[i]
  train <- fits[[train_condition]]
  heldout <- parts[[heldout_condition]]
  heldout_library <- colSums(heldout$counts)
  transfer <- fit_fibrodynmix_transfer(
    counts = heldout$counts,
    fit = train$fit,
    library_size = heldout_library,
    maxit_z = transfer_maxit_z
  )
  state_cols <- colnames(transfer$z_hat)
  entropy <- simplex_entropy(transfer$z_hat)
  dominant <- state_cols[max.col(transfer$z_hat, ties.method = "first")]

  direction_rows[[row_id]] <- data.frame(
    train_condition = train_condition,
    heldout_condition = heldout_condition,
    n_train_cells = ncol(train$prepared$counts),
    n_heldout_cells = ncol(heldout$counts),
    n_genes = nrow(train$prepared$counts),
    n_shared_genes = length(transfer$shared_genes),
    train_initial_objective = train$fit$nb_objective_trace[1L],
    train_best_objective = train$fit$best_objective,
    train_final_objective = utils::tail(train$fit$nb_objective_trace, 1),
    train_stop_reason = train$fit$stop_reason,
    transfer_nb_objective = transfer$nb_objective,
    transfer_nb_loglik = transfer$nb_loglik,
    transfer_converged = transfer$converged,
    transfer_z_convergence_rate = transfer$z_convergence_rate,
    transfer_n_nonconverged_cells = transfer$n_nonconverged_cells,
    mean_transfer_entropy = mean(entropy),
    stringsAsFactors = FALSE
  )
  row_id <- row_id + 1L

  composition_rows[[composition_id]] <- data.frame(
    train_condition = train_condition,
    heldout_condition = heldout_condition,
    source = "transferred_heldout",
    state = state_cols,
    composition = colMeans(transfer$z_hat),
    stringsAsFactors = FALSE
  )
  composition_id <- composition_id + 1L
  composition_rows[[composition_id]] <- data.frame(
    train_condition = train_condition,
    heldout_condition = heldout_condition,
    source = "training_fit",
    state = state_cols,
    composition = colMeans(train$fit$z_hat),
    stringsAsFactors = FALSE
  )
  composition_id <- composition_id + 1L

  cell_weight <- as.data.frame(transfer$z_hat, stringsAsFactors = FALSE)
  cell_weight$cell_id <- rownames(transfer$z_hat)
  cell_weight$train_condition <- train_condition
  cell_weight$heldout_condition <- heldout_condition
  cell_weight$dominant_state <- dominant
  cell_weight$entropy <- entropy
  cell_rows[[cell_id]] <- cell_weight[, c("train_condition", "heldout_condition", "cell_id", "dominant_state", "entropy", state_cols), drop = FALSE]
  cell_id <- cell_id + 1L

  marker_summary <- train$prepared$marker_summary
  marker_summary$train_condition <- train_condition
  marker_rows[[marker_id]] <- marker_summary[, c("train_condition", setdiff(colnames(marker_summary), "train_condition")), drop = FALSE]
  marker_id <- marker_id + 1L
}

dataset_paths <- data.frame(
  file = dataset$file,
  condition = dataset$condition,
  source_path = vapply(dataset$file, resolve_data_path, character(1)),
  stringsAsFactors = FALSE
)

write_tsv(dataset_paths, file.path(out_dir, "public_transfer_dataset_manifest.tsv"))
write_tsv(do.call(rbind, direction_rows), file.path(out_dir, "public_transfer_diagnostics.tsv"))
write_tsv(do.call(rbind, composition_rows), file.path(out_dir, "public_transfer_state_composition.tsv"))
write_tsv(do.call(rbind, cell_rows), file.path(out_dir, "public_transfer_cell_state_weights.tsv"))
write_tsv(do.call(rbind, marker_rows), file.path(out_dir, "public_transfer_marker_coverage.tsv"))

manifest <- data.frame(
  analysis = "public_realdata_transfer",
  primary_claim = "FibroDynMix can freeze a raw-count NB state program learned in one public fibroblast condition and infer held-out cell state mixtures in the other condition.",
  claim_boundary = "Real-data engineering smoke test on two public mouse breast fibroblast count files; not a full multi-donor human cross-cohort atlas validation or disease-mechanism claim.",
  n_directions = nrow(directions),
  max_cells_per_condition = max_cells,
  retained_genes = length(common_genes),
  seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "public_transfer_manifest.tsv"))

message("Public real-data transfer analysis written to: ", normalizePath(out_dir))
