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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_public_realdata_smoke.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

data_dir <- get_arg("data-dir", file.path(ROOT, "data", "public_dryad_breast_fibroblast"))
out_dir <- get_arg("out", file.path(ROOT, "analysis", "public_realdata_smoke"))
max_cells <- get_arg("max-cells", "240")
max_genes <- get_arg("max-genes", "700")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset <- data.frame(
  file = c("MT3_CAFs_raw.txt", "Normal_mammary_fibroblasts_raw.txt"),
  condition = c("disease", "normal"),
  study_id = c("Dryad_MT3_CAF", "Dryad_normal_mammary"),
  donor_id = c("MT3_CAF_pool", "normal_mammary_pool"),
  download_url = c(
    "https://zenodo.org/api/records/3977255/files/MT3_CAFs_raw.txt/content",
    "https://zenodo.org/api/records/3977255/files/Normal_mammary_fibroblasts_raw.txt/content"
  ),
  stringsAsFactors = FALSE
)

utils::write.table(
  dataset,
  file.path(out_dir, "public_dataset_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

download_file <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 1000) {
    return(TRUE)
  }
  cmd <- sprintf(
    "curl -L --fail --retry 2 -A 'Mozilla/5.0' -o %s %s",
    shQuote(dest),
    shQuote(url)
  )
  status <- system(cmd)
  isTRUE(status == 0L) && file.exists(dest) && file.info(dest)$size > 1000
}

for (i in seq_len(nrow(dataset))) {
  dest <- file.path(data_dir, dataset$file[i])
  if (!download_file(dataset$download_url[i], dest)) {
    utils::write.table(
      data.frame(
        status = "download_blocked",
        doi = "10.6071/M3238R",
        zenodo_record = "https://zenodo.org/records/3977255",
        dryad_record = "https://datadryad.org/dataset/doi:10.6071/M3238R",
        data_dir = normalizePath(data_dir, mustWork = FALSE),
        required_files = paste(dataset$file, collapse = ", "),
        stringsAsFactors = FALSE
      ),
      file.path(out_dir, "public_download_status.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    stop(
      paste(
        "Public command-line download was blocked or failed.",
        "Download these public raw-count files in a browser from https://zenodo.org/records/3977255 or https://datadryad.org/dataset/doi:10.6071/M3238R",
        "and place them under:",
        normalizePath(data_dir, mustWork = FALSE),
        "Required files:",
        paste(dataset$file, collapse = ", "),
        sep = "\n"
      ),
      call. = FALSE
    )
  }
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
  list(counts = mat, metadata = metadata)
}

parts <- lapply(seq_len(nrow(dataset)), function(i) {
  read_one(
    file.path(data_dir, dataset$file[i]),
    condition = dataset$condition[i],
    study_id = dataset$study_id[i],
    donor_id = dataset$donor_id[i]
  )
})

common_genes <- Reduce(intersect, lapply(parts, function(x) rownames(x$counts)))
counts <- do.call(cbind, lapply(parts, function(x) x$counts[common_genes, , drop = FALSE]))
metadata <- do.call(rbind, lapply(parts, function(x) x$metadata))

counts_path <- file.path(data_dir, "dryad_breast_fibroblast_counts.rds")
metadata_path <- file.path(data_dir, "dryad_breast_fibroblast_metadata.tsv")
saveRDS(counts, counts_path)
utils::write.table(metadata, metadata_path, sep = "\t", quote = FALSE, row.names = FALSE)

cmd <- c(
  "Rscript",
  file.path(ROOT, "scripts", "run_realdata_smoke.R"),
  paste0("--counts=", counts_path),
  paste0("--metadata=", metadata_path),
  paste0("--out=", out_dir),
  "--condition-col=condition",
  "--study-col=study_id",
  "--donor-col=donor_id",
  "--cell-id-col=cell_id",
  "--stratify-col=condition",
  "--normal-label=normal",
  "--disease-label=disease",
  paste0("--max-cells=", max_cells),
  paste0("--max-genes=", max_genes),
  "--min-cells-per-gene=0",
  "--min-counts-per-gene=0",
  "--n-outer=2",
  "--initializer-iter=4"
)

status <- system2(cmd[[1]], cmd[-1])
if (status != 0L) {
  stop("Public real-data smoke analysis failed.", call. = FALSE)
}

message("Public real-data smoke analysis written to: ", normalizePath(out_dir))
