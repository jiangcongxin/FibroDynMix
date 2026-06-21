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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/prepare_gse246215_fibroblast_inputs.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas"))
out_dir <- get_arg("out", data_dir)
max_cells_per_group <- as.integer(get_arg("max-cells-per-group", "80"))
seed <- as.integer(get_arg("seed", "246215"))

counts_path <- get_arg("counts", file.path(data_dir, "GSE246215_Fibroblast_counts.csv.gz"))
metadata_path <- get_arg("metadata", file.path(data_dir, "GSE246215_Fibroblast_metadata.csv.gz"))

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
required_cols <- c("CellName", "SampleID", "PatientID", "CancerType", "CancerType_short", "Phenotype")
missing <- setdiff(required_cols, colnames(metadata))
if (length(missing) > 0L) {
  stop(sprintf("GSE246215 metadata missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

metadata <- metadata[metadata$Phenotype == "Tumor", , drop = FALSE]
if (nrow(metadata) == 0L) {
  stop("No tumor fibroblast cells found in GSE246215 metadata.", call. = FALSE)
}

groups <- split(metadata, metadata$CancerType_short)
groups <- groups[vapply(groups, nrow, integer(1)) > 0L]
selected <- do.call(rbind, lapply(names(groups), function(group_name) {
  group <- groups[[group_name]]
  idx <- sample(seq_len(nrow(group)), min(max_cells_per_group, nrow(group)))
  group[idx, , drop = FALSE]
}))
selected <- selected[order(selected$CancerType_short, selected$SampleID, selected$CellName), , drop = FALSE]

selected_cells <- unique(selected$CellName)
counts_dt <- data.table::fread(
  counts_path,
  select = c("GeneName", selected_cells),
  data.table = FALSE,
  check.names = FALSE
)
genes <- make.unique(as.character(counts_dt$GeneName))

dataset_rows <- list()
summary_rows <- list()
for (group_name in sort(unique(selected$CancerType_short))) {
  group_meta <- selected[selected$CancerType_short == group_name, , drop = FALSE]
  group_cells <- group_meta$CellName
  mat <- as.matrix(counts_dt[, group_cells, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  colnames(mat) <- group_cells
  counts_out <- file.path(out_dir, sprintf("GSE246215_%s_fibroblast_counts.rds", group_name))
  saveRDS(mat, counts_out)

  meta_out <- group_meta[, c("CellName", "SampleID", "PatientID", "CancerType", "CancerType_short", "Phenotype", "Cluster"), drop = FALSE]
  meta_out$dataset_id <- sprintf("gse246215_%s_fibroblast", tolower(group_name))
  write_tsv(meta_out, file.path(out_dir, sprintf("GSE246215_%s_fibroblast_metadata.tsv", group_name)))

  dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
    dataset_id = sprintf("gse246215_%s_fibroblast", tolower(group_name)),
    study_id = sprintf("GSE246215_%s", group_name),
    donor_id = sprintf("GSE246215_%s_multi_patient", group_name),
    condition = "disease",
    organism = "Homo sapiens",
    tissue = group_meta$CancerType[1L],
    counts_path = counts_out,
    public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246215",
    download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE246215&format=file",
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    dataset_id = sprintf("gse246215_%s_fibroblast", tolower(group_name)),
    cancer_type_short = group_name,
    cancer_type = group_meta$CancerType[1L],
    n_selected_cells = nrow(group_meta),
    n_patients = length(unique(group_meta$PatientID)),
    n_samples = length(unique(group_meta$SampleID)),
    counts_path = counts_out,
    stringsAsFactors = FALSE
  )
}

dataset_manifest <- do.call(rbind, dataset_rows)
selection_summary <- do.call(rbind, summary_rows)
write_tsv(dataset_manifest, file.path(out_dir, "gse246215_fibroblast_dataset_manifest.tsv"))
write_tsv(selection_summary, file.path(out_dir, "gse246215_fibroblast_selection_summary.tsv"))
write_tsv(selected, file.path(out_dir, "gse246215_fibroblast_selected_cell_metadata.tsv"))

manifest <- data.frame(
  analysis = "prepare_gse246215_fibroblast_inputs",
  source_accession = "GSE246215",
  source_title = "Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation",
  n_groups = nrow(dataset_manifest),
  n_selected_cells = nrow(selected),
  n_genes = length(genes),
  max_cells_per_group = max_cells_per_group,
  seed = seed,
  claim_boundary = "Prepared sampled fibroblast count matrices from a processed GEO supplementary count matrix; raw FASTQ files are not provided for this GEO record.",
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "gse246215_prepare_manifest.tsv"))

message("Prepared GSE246215 FibroDynMix inputs under: ", normalizePath(out_dir))
