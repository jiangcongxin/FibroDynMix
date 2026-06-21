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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_gse246215_sensitivity.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas"))
out_dir <- get_arg("out", file.path(ROOT, "analysis", "gse246215_sensitivity"))
counts_path <- get_arg("counts", file.path(data_dir, "GSE246215_Fibroblast_counts.csv.gz"))
metadata_path <- get_arg("metadata", file.path(data_dir, "GSE246215_Fibroblast_metadata.csv.gz"))
seeds <- as.integer(strsplit(get_arg("seeds", "246215,246216"), ",", fixed = TRUE)[[1]])
max_cells_per_group <- as.integer(get_arg("max-cells-per-group", "60"))
max_genes <- as.integer(get_arg("max-genes", "600"))
trim_quantile <- as.numeric(get_arg("trim-quantile", "0.95"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

if (!file.exists(counts_path)) {
  stop(sprintf("Missing GSE246215 counts file: %s", counts_path), call. = FALSE)
}
if (!file.exists(metadata_path)) {
  stop(sprintf("Missing GSE246215 metadata file: %s", metadata_path), call. = FALSE)
}
if (length(seeds) == 0L || anyNA(seeds)) {
  stop("`--seeds` must contain at least one integer seed.", call. = FALSE)
}
if (length(trim_quantile) != 1L || is.na(trim_quantile) || trim_quantile <= 0 || trim_quantile >= 1) {
  stop("`--trim-quantile` must be between 0 and 1.", call. = FALSE)
}

metadata <- data.table::fread(metadata_path, data.table = FALSE)
metadata <- metadata[metadata$Phenotype == "Tumor", , drop = FALSE]
groups <- split(metadata, metadata$CancerType_short)
groups <- groups[vapply(groups, nrow, integer(1)) > 0L]
selected_all <- unique(unlist(lapply(groups, function(group) group$CellName), use.names = FALSE))

message("Reading GSE246215 counts for library-size QC...")
counts_for_qc <- data.table::fread(
  counts_path,
  select = c("GeneName", selected_all),
  data.table = FALSE,
  check.names = FALSE
)
gene_names <- make.unique(as.character(counts_for_qc$GeneName))
qc_mat <- as.matrix(counts_for_qc[, -1L, drop = FALSE])
storage.mode(qc_mat) <- "numeric"
rownames(qc_mat) <- gene_names
library_size <- colSums(qc_mat)
library_df <- data.frame(
  CellName = colnames(qc_mat),
  library_size = as.numeric(library_size),
  stringsAsFactors = FALSE
)
metadata <- merge(metadata, library_df, by = "CellName", all.x = TRUE, sort = FALSE)
write_tsv(
  aggregate(
    library_size ~ CancerType_short,
    data = metadata,
    FUN = function(x) paste(c(n = length(x), mean = mean(x), median = stats::median(x), q95 = stats::quantile(x, 0.95)), collapse = ";")
  ),
  file.path(out_dir, "gse246215_library_size_summary.tsv")
)

select_cells <- function(seed, mode) {
  set.seed(seed)
  selected <- list()
  for (group_name in sort(names(groups))) {
    group <- metadata[metadata$CancerType_short == group_name, , drop = FALSE]
    if (mode == "library_trim_q95") {
      cutoff <- stats::quantile(group$library_size, trim_quantile, na.rm = TRUE)
      group <- group[group$library_size <= cutoff, , drop = FALSE]
    }
    if (nrow(group) == 0L) {
      stop(sprintf("No cells remain for group %s under mode %s.", group_name, mode), call. = FALSE)
    }
    idx <- sample(seq_len(nrow(group)), min(max_cells_per_group, nrow(group)))
    selected[[group_name]] <- group[idx, , drop = FALSE]
  }
  do.call(rbind, selected)
}

write_condition_inputs <- function(selected, run_dir, mode, seed) {
  selected_cells <- unique(selected$CellName)
  counts_dt <- counts_for_qc[, c("GeneName", selected_cells), drop = FALSE]
  genes <- make.unique(as.character(counts_dt$GeneName))
  dataset_rows <- list()
  selection_rows <- list()
  for (group_name in sort(unique(selected$CancerType_short))) {
    group_meta <- selected[selected$CancerType_short == group_name, , drop = FALSE]
    group_cells <- group_meta$CellName
    mat <- as.matrix(counts_dt[, group_cells, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- genes
    colnames(mat) <- group_cells
    counts_out <- file.path(run_dir, sprintf("GSE246215_%s_%s_seed%s_counts.rds", group_name, mode, seed))
    saveRDS(mat, counts_out)

    dataset_id <- sprintf("gse246215_%s_%s_seed%s", tolower(group_name), mode, seed)
    dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      study_id = sprintf("GSE246215_%s_%s", group_name, mode),
      donor_id = sprintf("GSE246215_%s_multi_patient", group_name),
      condition = "disease",
      organism = "Homo sapiens",
      tissue = group_meta$CancerType[1L],
      counts_path = counts_out,
      public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246215",
      download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE246215&format=file",
      stringsAsFactors = FALSE
    )
    selection_rows[[length(selection_rows) + 1L]] <- data.frame(
      sensitivity_mode = mode,
      seed = seed,
      dataset_id = dataset_id,
      cancer_type_short = group_name,
      n_selected_cells = nrow(group_meta),
      n_patients = length(unique(group_meta$PatientID)),
      n_samples = length(unique(group_meta$SampleID)),
      mean_library_size = mean(group_meta$library_size),
      median_library_size = stats::median(group_meta$library_size),
      max_library_size = max(group_meta$library_size),
      stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, dataset_rows)
  selection <- do.call(rbind, selection_rows)
  manifest_path <- file.path(run_dir, "dataset_manifest.tsv")
  write_tsv(manifest, manifest_path)
  write_tsv(selection, file.path(run_dir, "selection_summary.tsv"))
  write_tsv(selected, file.path(run_dir, "selected_cell_metadata.tsv"))
  manifest_path
}

modes <- c("raw_sample", "library_trim_q95")
run_rows <- list()
composition_rows <- list()
transfer_rows <- list()
selection_rows <- list()
run_index <- 1L

for (mode in modes) {
  for (seed in seeds) {
    run_name <- sprintf("%s_seed%s", mode, seed)
    run_dir <- file.path(out_dir, "runs", run_name)
    validation_dir <- file.path(run_dir, "validation")
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    selected <- select_cells(seed, mode)
    manifest_path <- write_condition_inputs(selected, run_dir, mode, seed)

    cmd <- c(
      file.path(ROOT, "scripts", "run_multi_public_realdata_validation.R"),
      paste0("--dataset-manifest=", manifest_path),
      paste0("--out=", validation_dir),
      paste0("--max-cells=", max_cells_per_group),
      paste0("--max-genes=", max_genes),
      "--n-outer=2",
      "--initializer-iter=3",
      "--maxit-beta=12",
      "--maxit-z=12",
      "--transfer-maxit-z=120"
    )
    status <- system2("Rscript", cmd)
    if (status != 0L) {
      stop(sprintf("GSE246215 sensitivity validation failed for %s.", run_name), call. = FALSE)
    }

    manifest <- read_tsv(file.path(validation_dir, "multi_public_validation_manifest.tsv"))
    fit <- read_tsv(file.path(validation_dir, "multi_public_nb_fit_diagnostics.tsv"))
    transfer <- read_tsv(file.path(validation_dir, "multi_public_transfer_diagnostics.tsv"))
    composition <- read_tsv(file.path(validation_dir, "multi_public_dataset_state_composition.tsv"))
    selection <- read_tsv(file.path(run_dir, "selection_summary.tsv"))

    run_rows[[run_index]] <- data.frame(
      sensitivity_mode = mode,
      seed = seed,
      n_public_datasets = manifest$n_public_datasets,
      n_cells = manifest$n_cells,
      n_genes = manifest$n_genes,
      initial_objective = fit$initial_objective,
      best_objective = fit$best_objective,
      objective_improvement = fit$objective_improvement,
      mean_transfer_z_convergence_rate = fit$mean_transfer_z_convergence_rate,
      min_transfer_z_convergence_rate = fit$min_transfer_z_convergence_rate,
      run_dir = validation_dir,
      stringsAsFactors = FALSE
    )
    composition$sensitivity_mode <- mode
    composition$seed <- seed
    transfer$sensitivity_mode <- mode
    transfer$seed <- seed
    composition_rows[[run_index]] <- composition
    transfer_rows[[run_index]] <- transfer
    selection_rows[[run_index]] <- selection
    run_index <- run_index + 1L
  }
}

run_summary <- do.call(rbind, run_rows)
composition_all <- do.call(rbind, composition_rows)
transfer_all <- do.call(rbind, transfer_rows)
selection_all <- do.call(rbind, selection_rows)

write_tsv(run_summary, file.path(out_dir, "gse246215_sensitivity_run_summary.tsv"))
write_tsv(composition_all, file.path(out_dir, "gse246215_sensitivity_state_composition.tsv"))
write_tsv(transfer_all, file.path(out_dir, "gse246215_sensitivity_transfer_diagnostics.tsv"))
write_tsv(selection_all, file.path(out_dir, "gse246215_sensitivity_selection_summary.tsv"))

composition_variability <- do.call(rbind, lapply(split(composition_all, paste(composition_all$sensitivity_mode, composition_all$tissue, composition_all$state, sep = "|")), function(df) {
  data.frame(
    sensitivity_mode = df$sensitivity_mode[1],
    tissue = df$tissue[1],
    state = df$state[1],
    composition_mean = mean(df$composition),
    composition_sd = stats::sd(df$composition),
    composition_min = min(df$composition),
    composition_max = max(df$composition),
    n_runs = nrow(df),
    stringsAsFactors = FALSE
  )
}))
rownames(composition_variability) <- NULL
write_tsv(composition_variability, file.path(out_dir, "gse246215_sensitivity_composition_variability.tsv"))

manifest <- data.frame(
  analysis = "gse246215_sensitivity",
  primary_claim = "GSE246215 human fibroblast validation was stress-tested across random downsampling seeds and library-size trimming.",
  claim_boundary = "Sensitivity analysis uses sampled processed public count matrices from one GEO study; it evaluates computational robustness and library-size QC sensitivity, not definitive cancer-type biology.",
  n_modes = length(modes),
  n_seeds = length(seeds),
  n_runs = nrow(run_summary),
  min_transfer_z_convergence_rate = min(run_summary$min_transfer_z_convergence_rate),
  mean_objective_improvement = mean(run_summary$objective_improvement),
  max_composition_sd = max(composition_variability$composition_sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "gse246215_sensitivity_manifest.tsv"))

message("GSE246215 sensitivity analysis written to: ", normalizePath(out_dir))
