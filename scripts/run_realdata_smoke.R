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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_realdata_smoke.R"
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
    "real_data_interface.R"
  )
  invisible(lapply(file.path(root, "R", r_files), source))
}

source_package_files(ROOT)

out_dir <- get_arg("out", file.path(ROOT, "analysis", "realdata_smoke"))
counts_path <- get_arg("counts")
metadata_path <- get_arg("metadata")
markers_path <- get_arg("markers")
genes_path <- get_arg("genes")
cells_path <- get_arg("cells")
condition_col <- get_arg("condition-col", "condition")
study_col <- get_arg("study-col", "study_id")
donor_col <- get_arg("donor-col", "donor_id")
cell_id_col <- get_arg("cell-id-col", "cell_id")
stratify_col <- get_arg("stratify-col", NULL)
normal_label <- get_arg("normal-label", "normal")
disease_label <- get_arg("disease-label", "disease")
max_cells <- as.integer(get_arg("max-cells", "300"))
max_genes <- as.integer(get_arg("max-genes", "800"))
min_cells_per_gene <- as.integer(get_arg("min-cells-per-gene", "3"))
min_counts_per_gene <- as.numeric(get_arg("min-counts-per-gene", "5"))
seed <- as.integer(get_arg("seed", "1001"))
n_outer <- as.integer(get_arg("n-outer", "2"))
initializer_iter <- as.integer(get_arg("initializer-iter", "4"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
donor_l2 <- as.numeric(get_arg("donor-l2", "0.1"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))

if (is.null(counts_path) || is.null(metadata_path)) {
  stop("Usage: Rscript scripts/run_realdata_smoke.R --counts=PATH --metadata=PATH [--markers=PATH] [--out=DIR]", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

read_table_matrix <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  default_rownames <- identical(rownames(tab), as.character(seq_len(nrow(tab))))
  if (!default_rownames) {
    mat <- as.matrix(tab)
    storage.mode(mat) <- "numeric"
    rownames(mat) <- make.unique(rownames(tab))
    return(mat)
  }
  if (ncol(tab) < 2L) {
    stop("Count table must contain one gene column and at least one cell column.", call. = FALSE)
  }
  gene_names <- tab[[1L]]
  mat <- as.matrix(tab[, -1L, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(as.character(gene_names))
  mat
}

read_counts <- function(path, genes_path = NULL, cells_path = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
    return(as.matrix(x))
  }
  if (ext %in% c("txt", "tsv")) {
    return(read_table_matrix(path))
  }
  if (ext == "csv") {
    tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    gene_names <- tab[[1L]]
    mat <- as.matrix(tab[, -1L, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- make.unique(as.character(gene_names))
    return(mat)
  }
  if (ext == "mtx") {
    if (is.null(genes_path) || is.null(cells_path)) {
      stop("MTX input requires --genes=PATH and --cells=PATH.", call. = FALSE)
    }
    mat <- Matrix::readMM(path)
    genes <- utils::read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)[[1L]]
    cells <- utils::read.delim(cells_path, header = FALSE, stringsAsFactors = FALSE)[[1L]]
    if (nrow(mat) == length(genes) && ncol(mat) == length(cells)) {
      rownames(mat) <- make.unique(as.character(genes))
      colnames(mat) <- make.unique(as.character(cells))
      return(as.matrix(mat))
    }
    if (nrow(mat) == length(cells) && ncol(mat) == length(genes)) {
      mat <- Matrix::t(mat)
      rownames(mat) <- make.unique(as.character(genes))
      colnames(mat) <- make.unique(as.character(cells))
      return(as.matrix(mat))
    }
    stop("MTX dimensions do not match supplied gene/cell files.", call. = FALSE)
  }
  stop(sprintf("Unsupported count format: %s", path), call. = FALSE)
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
  human_hits <- sum(unique(unlist(human)) %in% gene_names)
  mouse_hits <- sum(unique(unlist(mouse)) %in% gene_names)
  if (mouse_hits > human_hits) mouse else human
}

read_markers <- function(path, gene_names) {
  if (is.null(path)) {
    return(default_fibroblast_markers(gene_names))
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }
  tab <- utils::read.delim(path, stringsAsFactors = FALSE)
  if (!all(c("state", "gene") %in% colnames(tab))) {
    stop("Marker TSV must contain `state` and `gene` columns.", call. = FALSE)
  }
  split(tab$gene, tab$state)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

counts <- read_counts(counts_path, genes_path, cells_path)
metadata <- utils::read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!cell_id_col %in% colnames(metadata)) {
  if (!is.null(colnames(counts)) && nrow(metadata) == ncol(counts)) {
    metadata[[cell_id_col]] <- colnames(counts)
  } else {
    stop(sprintf("Metadata must contain `%s` or have one row per count-matrix cell.", cell_id_col), call. = FALSE)
  }
}

if (!study_col %in% colnames(metadata)) {
  metadata[[study_col]] <- "public_dataset"
}
if (!donor_col %in% colnames(metadata)) {
  metadata[[donor_col]] <- metadata[[study_col]]
}
if (!condition_col %in% colnames(metadata)) {
  metadata[[condition_col]] <- "condition_unknown"
}

metadata <- metadata[metadata[[cell_id_col]] %in% colnames(counts), , drop = FALSE]
counts <- counts[, colnames(counts) %in% metadata[[cell_id_col]], drop = FALSE]

if (ncol(counts) > max_cells) {
  if (!is.null(stratify_col)) {
    if (!stratify_col %in% colnames(metadata)) {
      stop("`--stratify-col` must name a metadata column.", call. = FALSE)
    }
    strata <- split(metadata[[cell_id_col]], metadata[[stratify_col]])
    n_strata <- length(strata)
    base_n <- floor(max_cells / n_strata)
    remainder <- max_cells - base_n * n_strata
    selected_cells <- unlist(lapply(seq_along(strata), function(i) {
      target_n <- base_n + as.integer(i <= remainder)
      sample(strata[[i]], min(length(strata[[i]]), target_n))
    }), use.names = FALSE)
  } else {
    selected_cells <- sample(colnames(counts), max_cells)
  }
  counts <- counts[, selected_cells, drop = FALSE]
  metadata <- metadata[metadata[[cell_id_col]] %in% selected_cells, , drop = FALSE]
}

gene_rank <- order(rowSums(counts > 0), rowSums(counts), decreasing = TRUE)
markers <- read_markers(markers_path, rownames(counts))
marker_genes <- unique(unlist(markers, use.names = FALSE))
keep_marker_genes <- rownames(counts) %in% marker_genes
if (nrow(counts) > max_genes) {
  keep_ranked <- rep(FALSE, nrow(counts))
  keep_ranked[gene_rank[seq_len(min(max_genes, length(gene_rank)))]] <- TRUE
  keep_genes <- keep_ranked | keep_marker_genes
  counts <- counts[keep_genes, , drop = FALSE]
}

prepared <- prepare_fibrodynmix_data(
  counts = counts,
  cell_metadata = metadata,
  marker_index = markers,
  cell_id_col = cell_id_col,
  study_col = study_col,
  donor_col = donor_col,
  min_cells_per_gene = min_cells_per_gene,
  min_counts_per_gene = min_counts_per_gene,
  require_all_states = TRUE
)

fit <- fit_fibrodynmix_prepared(
  prepared,
  fit_study_effect = !has_flag("no-study-effect"),
  fit_donor_effect = if (has_flag("fit-donor-effect")) TRUE else NULL,
  n_outer = n_outer,
  initializer_args = list(n_iter = initializer_iter),
  study_l2 = study_l2,
  donor_l2 = donor_l2,
  marker_l2 = marker_l2,
  maxit_beta = 20,
  maxit_z = 20
)

z <- as.data.frame(fit$z_hat, stringsAsFactors = FALSE)
z[[cell_id_col]] <- rownames(fit$z_hat)
z <- z[, c(cell_id_col, setdiff(colnames(z), cell_id_col)), drop = FALSE]
write_tsv(z, file.path(out_dir, "cell_state_weights.tsv"))

state_cols <- colnames(fit$z_hat)
group_df <- unique(prepared$cell_metadata[, c(study_col, donor_col, condition_col), drop = FALSE])
composition <- do.call(rbind, lapply(seq_len(nrow(group_df)), function(row_idx) {
  idx <- prepared$cell_metadata[[study_col]] == group_df[[study_col]][row_idx] &
    prepared$cell_metadata[[donor_col]] == group_df[[donor_col]][row_idx] &
    prepared$cell_metadata[[condition_col]] == group_df[[condition_col]][row_idx]
  data.frame(
    study_id = group_df[[study_col]][row_idx],
    donor_id = group_df[[donor_col]][row_idx],
    condition = group_df[[condition_col]][row_idx],
    state = state_cols,
    composition = colMeans(fit$z_hat[idx, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}))
write_tsv(composition, file.path(out_dir, "state_composition.tsv"))

condition_values <- unique(as.character(prepared$cell_metadata[[condition_col]]))
transition_summary <- data.frame(status = "not_run", reason = "Need at least two condition labels.", stringsAsFactors = FALSE)
fpi <- compute_fpi(fit$z_hat)
if (all(c(normal_label, disease_label) %in% condition_values)) {
  normal <- colMeans(fit$z_hat[prepared$cell_metadata[[condition_col]] == normal_label, , drop = FALSE])
  disease <- colMeans(fit$z_hat[prepared$cell_metadata[[condition_col]] == disease_label, , drop = FALSE])
  cost <- compute_state_cost(fit$beta_hat)
  flow <- estimate_transition_flow(normal, disease, cost, lambda = 0.5)
  fpi <- compute_fpi(fit$z_hat, flow = flow$flow)
  flow_long <- as.data.frame(as.table(flow$flow), stringsAsFactors = FALSE)
  colnames(flow_long) <- c("source_state", "target_state", "flow")
  write_tsv(flow_long, file.path(out_dir, "transition_flow.tsv"))
  transition_summary <- data.frame(
    status = "ok",
    normal_label = normal_label,
    disease_label = disease_label,
    expected_cost = flow$expected_cost,
    entropy = flow$entropy,
    converged = flow$converged,
    iterations = flow$iterations,
    stringsAsFactors = FALSE
  )
}
fpi[[cell_id_col]] <- rownames(fit$z_hat)
fpi <- merge(
  prepared$cell_metadata[, c(cell_id_col, study_col, donor_col, condition_col), drop = FALSE],
  fpi,
  by.x = cell_id_col,
  by.y = cell_id_col,
  all.x = TRUE,
  sort = FALSE
)
write_tsv(fpi, file.path(out_dir, "cell_fpi.tsv"))

diagnostics <- data.frame(
  n_cells = ncol(prepared$counts),
  n_genes = nrow(prepared$counts),
  n_states = length(prepared$marker_index),
  initial_objective = fit$nb_objective_trace[1L],
  best_objective = fit$best_objective,
  final_objective = utils::tail(fit$nb_objective_trace, 1),
  objective_improvement = fit$nb_objective_trace[1L] - fit$best_objective,
  best_iteration = fit$best_iteration,
  executed_iterations = fit$executed_iterations,
  stop_reason = fit$stop_reason,
  fit_study_effect = !is.null(fit$study_effect),
  fit_donor_effect = !is.null(fit$donor_effect),
  study_effect_l2_norm = if (is.null(fit$study_effect)) NA_real_ else sqrt(sum(fit$study_effect^2)),
  donor_effect_l2_norm = if (is.null(fit$donor_effect)) NA_real_ else sqrt(sum(fit$donor_effect^2)),
  stringsAsFactors = FALSE
)
write_tsv(diagnostics, file.path(out_dir, "nb_fit_diagnostics.tsv"))
write_tsv(prepared$marker_summary, file.path(out_dir, "marker_coverage.tsv"))
write_tsv(as.data.frame(prepared$filter_summary), file.path(out_dir, "filter_summary.tsv"))
write_tsv(transition_summary, file.path(out_dir, "transition_summary.tsv"))

metadata_summary <- data.frame(
  field = c("source_counts", "source_metadata", "condition_col", "study_col", "donor_col", "stratify_col", "seed"),
  value = c(counts_path, metadata_path, condition_col, study_col, donor_col, ifelse(is.null(stratify_col), "", stratify_col), seed),
  stringsAsFactors = FALSE
)
write_tsv(metadata_summary, file.path(out_dir, "run_manifest.tsv"))

sample_summary <- as.data.frame(table(metadata[[condition_col]]), stringsAsFactors = FALSE)
colnames(sample_summary) <- c("condition", "n_cells")
write_tsv(sample_summary, file.path(out_dir, "condition_cell_counts.tsv"))

message("Real-data smoke analysis written to: ", normalizePath(out_dir))
