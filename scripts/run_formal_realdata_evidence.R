#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(FibroDynMix)
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required to read real-data evidence objects.", call. = FALSE)
}

out_dir <- file.path("analysis", "formal_realdata_evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_names <- names(get_fibrodynmix_markers("human", "scar"))
pathway_sets <- list(
  ECM_organization = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "POSTN", "SPARC", "THBS1", "MMP2", "MMP14", "TIMP1", "LOX"),
  TGF_beta_contractile = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "ACTA2", "TAGLN", "MYL9", "CNN1", "TPM2", "CTGF"),
  Interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "STAT1", "IRF7", "IFI6"),
  Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-A", "HLA-B", "HLA-C", "CD74", "B2M"),
  Cytokine_chemokine = c("IL6", "CXCL12", "CXCL14", "CCL2", "CXCL1", "CXCL2", "CXCL8", "CXCL10", "CXCL11"),
  Resident_matrix_homeostasis = c("DCN", "LUM", "COL14A1", "PI16", "PDGFRA", "DPT", "FBLN1", "FBLN2")
)

read_counts <- function(object, assay = "RNA", layer = "counts") {
  assay <- if (assay %in% names(object@assays)) assay else SeuratObject::DefaultAssay(object)
  counts <- suppressWarnings(SeuratObject::LayerData(object, assay = assay, layer = layer))
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    counts <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
  }
  counts
}

metadata_weights <- function(object) {
  md <- as.data.frame(object[[]], stringsAsFactors = FALSE)
  md$cell_id <- rownames(md)
  z_cols <- grep("^fibrodynmix_z_", colnames(md), value = TRUE)
  if (length(z_cols) == 0L) {
    z_cols <- grep("^fibrodynmix_(resident|inflammatory|myofibroblast|ECM|antigen|IFN)", colnames(md), value = TRUE)
  }
  weights <- md[, c("cell_id", z_cols), drop = FALSE]
  state_cols <- colnames(weights)[-1L]
  state_cols <- sub("^fibrodynmix_z_", "", state_cols)
  state_cols <- sub("^fibrodynmix_", "", state_cols)
  state_cols <- gsub("_", "-", state_cols, fixed = TRUE)
  state_cols <- gsub("\\.", "-", state_cols)
  colnames(weights) <- c("cell_id", state_cols)
  list(metadata = md, weights = weights, state_cols = setdiff(colnames(weights), "cell_id"))
}

top_state_genes <- function(counts, weights, state_cols, max_cells = 700L, top_n = 40L, seed = 20260606) {
  common_cells <- intersect(colnames(counts), weights$cell_id)
  set.seed(seed)
  if (length(common_cells) > max_cells) {
    common_cells <- sample(common_cells, max_cells)
  }
  counts <- counts[, common_cells, drop = FALSE]
  weights <- weights[match(common_cells, weights$cell_id), , drop = FALSE]
  lib <- matrix_col_sums(counts)
  lib[lib <= 0] <- 1
  norm <- as.matrix(log_normalize_counts(counts, lib))
  rows <- list()
  idx <- 1L
  for (state in state_cols) {
    z <- weights[[state]]
    cors <- apply(norm, 1L, function(g) suppressWarnings(stats::cor(g, z, method = "spearman")))
    cors[is.na(cors)] <- 0
    top <- sort(cors, decreasing = TRUE)
    top <- top[seq_len(min(top_n, length(top)))]
    rows[[idx]] <- data.frame(state = state, gene = names(top), rank = seq_along(top), spearman_rho = as.numeric(top), stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

pathway_enrichment <- function(programs, universe_genes) {
  rows <- list()
  idx <- 1L
  for (state in unique(programs$state)) {
    genes <- programs$gene[programs$state == state]
    for (pathway in names(pathway_sets)) {
      pathway_genes <- intersect(pathway_sets[[pathway]], universe_genes)
      overlap <- intersect(genes, pathway_genes)
      p <- if (length(pathway_genes) == 0L) NA_real_ else {
        stats::phyper(length(overlap) - 1L, length(pathway_genes), length(universe_genes) - length(pathway_genes), length(genes), lower.tail = FALSE)
      }
      rows[[idx]] <- data.frame(
        state = state,
        pathway = pathway,
        n_program_genes = length(genes),
        n_pathway_genes = length(pathway_genes),
        n_overlap = length(overlap),
        overlap_genes = paste(overlap, collapse = ";"),
        p_value = p,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$q_value <- stats::p.adjust(out$p_value, method = "BH")
  out
}

bootstrap_stability <- function(dataset, counts, metadata, marker_index, max_cells = 220L, max_genes = 700L, n_boot = 5L, seed = 20260606) {
  set.seed(seed)
  cells <- colnames(counts)
  if (length(cells) > max_cells) {
    cells <- sample(cells, max_cells)
  }
  marker_genes <- unique(unlist(marker_index, use.names = FALSE))
  variable_genes <- names(sort(matrix_row_sums(counts[, cells, drop = FALSE]), decreasing = TRUE))
  genes <- unique(c(marker_genes, utils::head(variable_genes, max_genes)))
  genes <- intersect(genes, rownames(counts))
  counts <- counts[genes, cells, drop = FALSE]
  metadata <- metadata[match(cells, metadata$cell_id), , drop = FALSE]
  boot <- bootstrap_fibrodynmix(
    counts = counts,
    marker_index = marker_index,
    library_size = matrix_col_sums(counts),
    cell_metadata = metadata,
    sample_col = if ("sample_id" %in% colnames(metadata)) "sample_id" else NULL,
    method = "initializer",
    n_boot = n_boot,
    seed = seed,
    fit_args = list(n_iter = 4),
    keep_fits = FALSE
  )
  rows <- lapply(boot$cell_summary$state |> unique(), function(state) {
    x <- boot$cell_summary[boot$cell_summary$state == state, , drop = FALSE]
    data.frame(
      dataset = dataset,
      state = state,
      bootstrap_method = "initializer",
      bootstrap_n = n_boot,
      bootstrap_cells = length(cells),
      bootstrap_genes = length(genes),
      median_cell_z_sd = stats::median(x$sd, na.rm = TRUE),
      p90_cell_z_sd = as.numeric(stats::quantile(x$sd, 0.9, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

analyze_dataset <- function(dataset, seurat_path, existing_state_summary = NULL, seed = 20260606) {
  object <- readRDS(seurat_path)
  counts <- read_counts(object)
  mw <- metadata_weights(object)
  marker_index <- get_fibrodynmix_markers("human", "scar")
  marker_index <- lapply(marker_index, intersect, y = rownames(counts))
  eval <- evaluate_fibrodynmix_annotation(
    mw$weights,
    state_cols = mw$state_cols,
    metadata = mw$metadata,
    metadata_cell_col = "cell_id",
    cluster_col = if ("seurat_clusters" %in% colnames(mw$metadata)) "seurat_clusters" else NULL,
    condition_col = if ("condition" %in% colnames(mw$metadata)) "condition" else NULL,
    expression = counts,
    marker_index = marker_index
  )
  programs <- top_state_genes(counts, mw$weights, mw$state_cols, seed = seed)
  programs$dataset <- dataset
  enrich <- pathway_enrichment(programs, rownames(counts))
  enrich$dataset <- dataset
  boot <- bootstrap_stability(dataset, counts, mw$metadata, marker_index, seed = seed)
  state <- eval$state_summary
  state$dataset <- dataset
  top_path <- do.call(rbind, lapply(split(enrich, enrich$state), function(x) x[order(x$q_value, -x$n_overlap), ][1L, , drop = FALSE]))
  top_path <- top_path[, c("state", "pathway", "n_overlap", "overlap_genes", "q_value"), drop = FALSE]
  merged <- merge(state, boot, by = c("state"), all.x = TRUE, sort = FALSE)
  merged <- merge(merged, top_path, by = "state", all.x = TRUE, sort = FALSE)
  merged$dataset <- dataset
  merged <- merged[, !colnames(merged) %in% c("dataset.x", "dataset.y"), drop = FALSE]
  list(state_evidence = merged, programs = programs, enrichment = enrich, bootstrap = boot, evaluation = eval)
}

datasets <- list(
  GSE163973 = "projects/keloid_fibro_gse163973/results/fibro_seurat_with_fibrodynmix.rds",
  GSE243716 = "projects/scar_fibro_gse243716/results/gse243716_fibro_seurat_with_fibrodynmix.rds"
)
available <- datasets[file.exists(unlist(datasets))]
if (length(available) == 0L) {
  stop("No real-data Seurat result objects were found.", call. = FALSE)
}

results <- lapply(names(available), function(dataset) analyze_dataset(dataset, available[[dataset]]))
names(results) <- names(available)

evidence <- do.call(rbind, lapply(results, `[[`, "state_evidence"))
programs <- do.call(rbind, lapply(results, `[[`, "programs"))
enrichment <- do.call(rbind, lapply(results, `[[`, "enrichment"))
bootstrap <- do.call(rbind, lapply(results, `[[`, "bootstrap"))

nb_transfer_path <- file.path("analysis", "nb_transfer_bootstrap_realdata", "nb_transfer_bootstrap_summary.tsv")
if (file.exists(nb_transfer_path)) {
  nb_transfer <- utils::read.delim(nb_transfer_path, stringsAsFactors = FALSE, check.names = FALSE)
  nb_transfer <- nb_transfer[, c("dataset", "state", "bootstrap_n", "n_eval_cells", "median_cell_z_sd", "p90_cell_z_sd", "max_cell_z_sd"), drop = FALSE]
  colnames(nb_transfer) <- c(
    "dataset",
    "state",
    "nb_transfer_bootstrap_n",
    "nb_transfer_eval_cells",
    "nb_transfer_median_cell_z_sd",
    "nb_transfer_p90_cell_z_sd",
    "nb_transfer_max_cell_z_sd"
  )
  evidence <- merge(evidence, nb_transfer, by = c("dataset", "state"), all.x = TRUE, sort = FALSE)
} else {
  evidence$nb_transfer_bootstrap_n <- NA_integer_
  evidence$nb_transfer_eval_cells <- NA_integer_
  evidence$nb_transfer_median_cell_z_sd <- NA_real_
  evidence$nb_transfer_p90_cell_z_sd <- NA_real_
  evidence$nb_transfer_max_cell_z_sd <- NA_real_
}

utils::write.table(evidence, file.path(out_dir, "formal_state_evidence_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(programs, file.path(out_dir, "state_associated_top_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(enrichment, file.path(out_dir, "state_pathway_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(bootstrap, file.path(out_dir, "bootstrap_stability_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

manifest <- data.frame(
  analysis = "formal_realdata_evidence",
  datasets = paste(names(available), collapse = ","),
  n_state_rows = nrow(evidence),
  n_supported_rows = sum(evidence$support_label == "supported", na.rm = TRUE),
  n_exploratory_rows = sum(evidence$support_label == "exploratory", na.rm = TRUE),
  max_median_bootstrap_sd = max(evidence$median_cell_z_sd, na.rm = TRUE),
  max_nb_transfer_p90_cell_z_sd = max(evidence$nb_transfer_p90_cell_z_sd, na.rm = TRUE),
  n_pathway_rows_q10 = sum(enrichment$q_value <= 0.10, na.rm = TRUE),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(out_dir, "formal_state_evidence_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(evidence[, c("dataset", "state", "n_cells", "support_label", "marker_log2_ratio_own_vs_other", "median_cell_z_sd", "pathway", "q_value"), drop = FALSE])
