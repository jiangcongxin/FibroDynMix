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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/prepare_gse167339_human_fibroblast_inputs.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast"))
out_dir <- get_arg("out", data_dir)
max_cells_per_sample <- as.integer(get_arg("max-cells-per-sample", "80"))
seed <- as.integer(get_arg("seed", "167339"))
include_hash_pool <- identical(tolower(get_arg("include-hash-pool", "false")), "true")
min_hash_cells <- as.integer(get_arg("min-hash-cells", "40"))
hash_min_count <- as.numeric(get_arg("hash-min-count", "10"))
hash_min_ratio <- as.numeric(get_arg("hash-min-ratio", "2"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

sample_table <- data.frame(
  gsm = c(
    "GSM5102531", "GSM5102532", "GSM5102533",
    "GSM5102534", "GSM5102535", "GSM5102536", "GSM5102537",
    "GSM5102538"
  ),
  prefix = c(
    "GSM5102531_KC-C-HH25WDMXX",
    "GSM5102532_KC-S-HH25WDMXX",
    "GSM5102533_KC-T-HH25WDMXX",
    "GSM5102534_H2C-HCKLKBBXY",
    "GSM5102535_H2SB-HCKLKBBXY",
    "GSM5102536_H2SU-HCKLKBBXY",
    "GSM5102537_H2T-HCKLKBBXY",
    "GSM5102538_H3-hash"
  ),
  donor_id = c("Human1", "Human1", "Human1", "Human2", "Human2", "Human2", "Human2", "Human3"),
  condition = c("control", "strain", "strain_treatment", "control", "strain", "strain", "strain_treatment", "hash_pool"),
  treatment = c("No strain control", "Strain", "Strain plus treatment", "No strain control", "Strain v1", "Strain repeat", "Strain plus treatment", "Hash pooled mixed conditions"),
  use_by_default = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
sample_table <- sample_table[sample_table$use_by_default, , drop = FALSE]

read_10x_sample <- function(prefix) {
  matrix_path <- file.path(data_dir, paste0(prefix, "_matrix.mtx.gz"))
  features_path <- file.path(data_dir, paste0(prefix, "_features.tsv.gz"))
  barcodes_path <- file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz"))
  missing <- c(matrix_path, features_path, barcodes_path)[!file.exists(c(matrix_path, features_path, barcodes_path))]
  if (length(missing) > 0L) {
    stop(sprintf("Missing GSE167339 sample files: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  mat <- Matrix::readMM(matrix_path)
  features <- utils::read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- utils::read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  gene_name <- if (ncol(features) >= 2L) features[[2L]] else features[[1L]]
  rownames(mat) <- make.unique(as.character(gene_name))
  colnames(mat) <- as.character(barcodes[[1L]])
  mat
}

dataset_rows <- list()
summary_rows <- list()
cell_metadata_rows <- list()

for (i in seq_len(nrow(sample_table))) {
  row <- sample_table[i, , drop = FALSE]
  mat <- read_10x_sample(row$prefix)
  if (ncol(mat) > max_cells_per_sample) {
    selected_cells <- sample(colnames(mat), max_cells_per_sample)
    mat <- mat[, selected_cells, drop = FALSE]
  }
  cell_ids <- paste(row$gsm, colnames(mat), sep = "_")
  colnames(mat) <- cell_ids
  counts_out <- file.path(out_dir, sprintf("%s_human_fibroblast_counts.rds", row$gsm))
  saveRDS(mat, counts_out)

  cell_metadata <- data.frame(
    cell_id = cell_ids,
    original_barcode = sub(paste0("^", row$gsm, "_"), "", cell_ids),
    gsm = row$gsm,
    donor_id = row$donor_id,
    condition = row$condition,
    treatment = row$treatment,
    top_hash = NA_character_,
    top_hash_count = NA_real_,
    hash_ratio = NA_real_,
    stringsAsFactors = FALSE
  )
  cell_metadata_rows[[length(cell_metadata_rows) + 1L]] <- cell_metadata
  write_tsv(cell_metadata, file.path(out_dir, sprintf("%s_human_fibroblast_metadata.tsv", row$gsm)))

  dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
    dataset_id = sprintf("gse167339_%s_%s_%s", tolower(row$donor_id), row$condition, tolower(row$gsm)),
    study_id = sprintf("GSE167339_%s", row$condition),
    donor_id = row$donor_id,
    condition = ifelse(row$condition == "control", "normal", "disease"),
    organism = "Homo sapiens",
    tissue = "3D collagen scar-system human dermal fibroblast culture",
    counts_path = counts_out,
    public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167339",
    download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167339&format=file",
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    dataset_id = sprintf("gse167339_%s_%s_%s", tolower(row$donor_id), row$condition, tolower(row$gsm)),
    gsm = row$gsm,
    donor_id = row$donor_id,
    original_condition = row$condition,
    model_condition = ifelse(row$condition == "control", "normal", "disease"),
    treatment = row$treatment,
    n_selected_cells = ncol(mat),
    n_genes = nrow(mat),
    median_library_size = stats::median(Matrix::colSums(mat)),
    counts_path = counts_out,
    stringsAsFactors = FALSE
  )
}

if (include_hash_pool) {
  hash_prefix <- "GSM5102538_H3-hash"
  hash_mat <- read_10x_sample(hash_prefix)
  hash_features <- utils::read.delim(
    file.path(data_dir, paste0(hash_prefix, "_features.tsv.gz")),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  feature_type <- if (ncol(hash_features) >= 3L) hash_features[[3L]] else rep("Gene Expression", nrow(hash_features))
  gene_rows <- feature_type == "Gene Expression"
  antibody_rows <- feature_type == "Antibody Capture"
  if (!any(antibody_rows)) {
    stop("include-hash-pool=true requested, but GSM5102538 has no Antibody Capture rows.", call. = FALSE)
  }
  gene_mat <- hash_mat[gene_rows, , drop = FALSE]
  antibody_mat <- as.matrix(hash_mat[antibody_rows, , drop = FALSE])
  antibody_names <- rownames(hash_mat)[antibody_rows]
  rownames(antibody_mat) <- antibody_names

  top_hash_index <- max.col(t(antibody_mat), ties.method = "first")
  top_hash <- antibody_names[top_hash_index]
  top_count <- apply(antibody_mat, 2L, max)
  second_count <- apply(antibody_mat, 2L, function(x) sort(x, decreasing = TRUE)[2L])
  hash_ratio <- top_count / pmax(second_count, 1)
  assigned <- top_count >= hash_min_count & hash_ratio >= hash_min_ratio

  hash_summary <- data.frame(
    hash_tag = antibody_names,
    n_top_cells = as.integer(table(factor(top_hash, levels = antibody_names))),
    n_assigned_cells = as.integer(table(factor(top_hash[assigned], levels = antibody_names))),
    hash_min_count = hash_min_count,
    hash_min_ratio = hash_min_ratio,
    stringsAsFactors = FALSE
  )
  write_tsv(hash_summary, file.path(out_dir, "gse167339_human3_hash_demux_summary.tsv"))

  used_hashes <- antibody_names[hash_summary$n_assigned_cells >= min_hash_cells]
  for (hash_tag in used_hashes) {
    hash_cells <- colnames(gene_mat)[assigned & top_hash == hash_tag]
    if (length(hash_cells) > max_cells_per_sample) {
      hash_cells <- sample(hash_cells, max_cells_per_sample)
    }
    mat <- gene_mat[, hash_cells, drop = FALSE]
    cell_ids <- paste("GSM5102538", gsub("[^A-Za-z0-9]+", "_", hash_tag), colnames(mat), sep = "_")
    colnames(mat) <- cell_ids
    hash_id <- tolower(gsub("[^A-Za-z0-9]+", "_", sub("--.*$", "", hash_tag)))
    counts_out <- file.path(out_dir, sprintf("GSM5102538_%s_human_fibroblast_counts.rds", hash_id))
    saveRDS(mat, counts_out)

    cell_metadata <- data.frame(
      cell_id = cell_ids,
      original_barcode = hash_cells,
      gsm = "GSM5102538",
      donor_id = "Human3",
      condition = "hash_unknown",
      treatment = sprintf("Hash-demultiplexed %s; source condition not mapped in public workbook", hash_tag),
      top_hash = hash_tag,
      top_hash_count = top_count[hash_cells],
      hash_ratio = hash_ratio[hash_cells],
      stringsAsFactors = FALSE
    )
    cell_metadata_rows[[length(cell_metadata_rows) + 1L]] <- cell_metadata
    write_tsv(cell_metadata, file.path(out_dir, sprintf("GSM5102538_%s_human_fibroblast_metadata.tsv", hash_id)))

    dataset_id <- sprintf("gse167339_human3_%s_gsm5102538", hash_id)
    dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      study_id = "GSE167339_hash_unknown",
      donor_id = "Human3",
      condition = "hash_unknown",
      organism = "Homo sapiens",
      tissue = "3D collagen scar-system human dermal fibroblast culture",
      counts_path = counts_out,
      public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167339",
      download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167339&format=file",
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      gsm = "GSM5102538",
      donor_id = "Human3",
      original_condition = "hash_unknown",
      model_condition = "hash_unknown",
      treatment = sprintf("Hash-demultiplexed %s; source condition not mapped", hash_tag),
      n_selected_cells = ncol(mat),
      n_genes = nrow(mat),
      median_library_size = stats::median(Matrix::colSums(mat)),
      counts_path = counts_out,
      stringsAsFactors = FALSE
    )
  }
}

dataset_manifest <- do.call(rbind, dataset_rows)
selection_summary <- do.call(rbind, summary_rows)
selected_metadata <- do.call(rbind, cell_metadata_rows)

write_tsv(dataset_manifest, file.path(out_dir, "gse167339_human_fibroblast_dataset_manifest.tsv"))
write_tsv(selection_summary, file.path(out_dir, "gse167339_human_fibroblast_selection_summary.tsv"))
write_tsv(selected_metadata, file.path(out_dir, "gse167339_human_fibroblast_selected_cell_metadata.tsv"))

manifest <- data.frame(
  analysis = "prepare_gse167339_human_fibroblast_inputs",
  source_accession = "GSE167339",
  source_title = "Disrupting Mechanotransduction Promotes Regenerative Phenotypes in Human Cells",
  n_samples = nrow(dataset_manifest),
  n_donors = length(unique(dataset_manifest$donor_id)),
  n_selected_cells = nrow(selected_metadata),
  n_genes_min = min(selection_summary$n_genes),
  max_cells_per_sample = max_cells_per_sample,
  included_hash_pool = include_hash_pool,
  n_hash_demux_groups = if (include_hash_pool && exists("used_hashes")) length(used_hashes) else 0L,
  hash_min_count = hash_min_count,
  hash_min_ratio = hash_min_ratio,
  seed = seed,
  claim_boundary = if (include_hash_pool) {
    "Prepared sampled 10x-style processed MTX count matrices from GEO supplementary files. Human 3 is included as hash-demultiplexed pseudo-samples with condition hash_unknown because the public workbook does not map hash tags to treatment labels."
  } else {
    "Prepared sampled 10x-style processed MTX count matrices from GEO supplementary files. Default run excludes the hash-pooled Human 3 sample; use include-hash-pool=true to add hash-demultiplexed pseudo-samples without treatment-label interpretation."
  },
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "gse167339_prepare_manifest.tsv"))

message("Prepared GSE167339 FibroDynMix inputs under: ", normalizePath(out_dir))
