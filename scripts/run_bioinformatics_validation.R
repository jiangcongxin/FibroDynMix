#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_bioinformatics_validation.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "bioinformatics_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

state_names <- c("resident", "inflammatory", "myofibroblast", "ECM-remodeling", "antigen-presenting", "IFN-stress")

marker_sets <- list(
  fibroblast_positive = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL14A1", "FBLN1"),
  endothelial_negative = c("PECAM1", "VWF", "KDR", "CDH5", "CLDN5"),
  epithelial_negative = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  immune_negative = c("PTPRC", "LYZ", "CD3D", "CD79A", "MS4A1", "NKG7"),
  pericyte_smc_negative = c("RGS5", "MCAM", "MYH11", "NOTCH3")
)

pathway_sets <- list(
  ECM_organization = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "POSTN", "SPARC", "THBS1", "MMP2", "MMP14", "TIMP1", "LOX"),
  TGF_beta_contractile = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "ACTA2", "TAGLN", "MYL9", "CNN1", "TPM2", "CTGF"),
  Interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "STAT1", "IRF7", "IFI6"),
  Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-A", "HLA-B", "HLA-C", "CD74", "B2M"),
  Cytokine_chemokine = c("IL6", "CXCL12", "CXCL14", "CCL2", "CXCL1", "CXCL2", "CXCL8", "CXCL10", "CXCL11"),
  Resident_matrix_homeostasis = c("DCN", "LUM", "COL14A1", "PI16", "PDGFRA", "DPT", "FBLN1", "FBLN2")
)

normalize_counts <- function(counts) {
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  log1p(t(t(counts) / lib * 10000))
}

read_dataset_counts <- function(registry_row) {
  counts <- readRDS(registry_row$counts_path)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  colnames(counts) <- paste(registry_row$study_id, colnames(counts), sep = "_")
  counts
}

score_gene_set <- function(norm_counts, genes) {
  present <- intersect(genes, rownames(norm_counts))
  if (length(present) == 0L) {
    return(rep(NA_real_, ncol(norm_counts)))
  }
  colMeans(norm_counts[present, , drop = FALSE])
}

dataset_qc <- function(dataset_label, registry_path) {
  registry <- read_tsv(registry_path)
  registry$counts_path <- ifelse(grepl("^/", registry$counts_path), registry$counts_path, file.path(ROOT, registry$counts_path))
  rows <- lapply(seq_len(nrow(registry)), function(i) {
    counts <- read_dataset_counts(registry[i, , drop = FALSE])
    norm <- normalize_counts(counts)
    scores <- as.data.frame(lapply(marker_sets, function(genes) score_gene_set(norm, genes)), stringsAsFactors = FALSE)
    positive <- scores$fibroblast_positive
    negative <- rowMeans(scores[, setdiff(colnames(scores), "fibroblast_positive"), drop = FALSE], na.rm = TRUE)
    purity_margin <- positive - negative
    data.frame(
      validation_dataset = dataset_label,
      dataset_id = registry$dataset_id[i],
      study_id = registry$study_id[i],
      donor_id = registry$donor_id[i],
      condition = registry$condition[i],
      organism = if ("organism" %in% colnames(registry)) registry$organism[i] else NA_character_,
      tissue = if ("tissue" %in% colnames(registry)) registry$tissue[i] else NA_character_,
      n_cells = ncol(counts),
      n_genes = nrow(counts),
      fibroblast_positive_score = mean(positive, na.rm = TRUE),
      endothelial_negative_score = mean(scores$endothelial_negative, na.rm = TRUE),
      epithelial_negative_score = mean(scores$epithelial_negative, na.rm = TRUE),
      immune_negative_score = mean(scores$immune_negative, na.rm = TRUE),
      pericyte_smc_negative_score = mean(scores$pericyte_smc_negative, na.rm = TRUE),
      purity_margin_mean = mean(purity_margin, na.rm = TRUE),
      purity_margin_median = stats::median(purity_margin, na.rm = TRUE),
      low_purity_fraction = mean(purity_margin < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

load_combined_counts_and_weights <- function(dataset_label, registry_path, weights_path) {
  registry <- read_tsv(registry_path)
  registry$counts_path <- ifelse(grepl("^/", registry$counts_path), registry$counts_path, file.path(ROOT, registry$counts_path))
  parts <- lapply(seq_len(nrow(registry)), function(i) read_dataset_counts(registry[i, , drop = FALSE]))
  common <- Reduce(intersect, lapply(parts, rownames))
  counts <- do.call(cbind, lapply(parts, function(x) x[common, , drop = FALSE]))
  weights <- read_tsv(weights_path)
  weights <- weights[weights$cell_id %in% colnames(counts), , drop = FALSE]
  counts <- counts[, weights$cell_id, drop = FALSE]
  list(dataset_label = dataset_label, counts = counts, weights = weights)
}

state_programs <- function(loaded, top_n = 30L) {
  norm <- normalize_counts(loaded$counts)
  rows <- list()
  idx <- 1L
  for (state in state_names) {
    z <- loaded$weights[[state]]
    cors <- apply(norm, 1L, function(g) suppressWarnings(stats::cor(g, z, method = "spearman")))
    cors[is.na(cors)] <- 0
    top <- sort(cors, decreasing = TRUE)
    top <- top[seq_len(min(top_n, length(top)))]
    rows[[idx]] <- data.frame(
      validation_dataset = loaded$dataset_label,
      state = state,
      gene = names(top),
      rank = seq_along(top),
      spearman_rho = as.numeric(top),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

enrich_programs <- function(program_table, universe_genes) {
  rows <- list()
  idx <- 1L
  for (dataset_label in unique(program_table$validation_dataset)) {
    for (state in unique(program_table$state)) {
      genes <- program_table$gene[program_table$validation_dataset == dataset_label & program_table$state == state]
      for (pathway in names(pathway_sets)) {
        pathway_genes <- intersect(pathway_sets[[pathway]], universe_genes)
        overlap <- intersect(genes, pathway_genes)
        p <- if (length(pathway_genes) == 0L) {
          NA_real_
        } else {
          stats::phyper(length(overlap) - 1L, length(pathway_genes), length(universe_genes) - length(pathway_genes), length(genes), lower.tail = FALSE)
        }
        rows[[idx]] <- data.frame(
          validation_dataset = dataset_label,
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
  }
  enrich <- do.call(rbind, rows)
  enrich$q_value <- stats::p.adjust(enrich$p_value, method = "BH")
  enrich
}

state_abundance <- function(dataset_label, weights_path, registry_path) {
  weights <- read_tsv(weights_path)
  registry <- read_tsv(registry_path)
  meta <- unique(registry[, c("dataset_id", "study_id", "donor_id", "condition", "tissue"), drop = FALSE])
  weights <- merge(weights, meta, by = c("dataset_id", "study_id", "condition"), all.x = TRUE, sort = FALSE)
  rows <- list()
  idx <- 1L
  groups <- unique(weights[, c("donor_id", "condition"), drop = FALSE])
  for (i in seq_len(nrow(groups))) {
    keep <- weights$donor_id == groups$donor_id[i] & weights$condition == groups$condition[i]
    z <- as.matrix(weights[keep, state_names, drop = FALSE])
    entropy <- weights$entropy[keep]
    for (state in state_names) {
      rows[[idx]] <- data.frame(
        validation_dataset = dataset_label,
        donor_id = groups$donor_id[i],
        condition = groups$condition[i],
        n_cells = sum(keep),
        state = state,
        composition = mean(z[, state]),
        mean_entropy = mean(entropy),
        mean_fpi = mean(entropy + (1 - apply(z, 1L, max))),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

donor_effects <- function(abundance) {
  rows <- list()
  idx <- 1L
  for (dataset_label in unique(abundance$validation_dataset)) {
    x <- abundance[abundance$validation_dataset == dataset_label, , drop = FALSE]
    for (state in unique(x$state)) {
      xs <- x[x$state == state, , drop = FALSE]
      donors <- unique(xs$donor_id)
      for (donor in donors) {
        xd <- xs[xs$donor_id == donor, , drop = FALSE]
        if (all(c("normal", "disease") %in% xd$condition)) {
          rows[[idx]] <- data.frame(
            validation_dataset = dataset_label,
            donor_id = donor,
            state = state,
            contrast = "disease_minus_normal",
            delta_composition = xd$composition[xd$condition == "disease"][1L] - xd$composition[xd$condition == "normal"][1L],
            delta_entropy = xd$mean_entropy[xd$condition == "disease"][1L] - xd$mean_entropy[xd$condition == "normal"][1L],
            delta_fpi = xd$mean_fpi[xd$condition == "disease"][1L] - xd$mean_fpi[xd$condition == "normal"][1L],
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

gse246215_registry <- file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas", "gse246215_fibroblast_dataset_manifest.tsv")
gse167339_registry <- file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_human_fibroblast_dataset_manifest.tsv")
gse246215_weights <- file.path(ROOT, "analysis", "independent_geo_gse246215_validation", "multi_public_cell_state_weights.tsv")
gse167339_weights <- file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_cell_state_weights.tsv")

purity <- rbind(
  dataset_qc("GSE246215", gse246215_registry),
  dataset_qc("GSE167339", gse167339_registry)
)
write_tsv(purity, file.path(OUT, "fibroblast_purity_qc.tsv"))

loaded <- list(
  load_combined_counts_and_weights("GSE246215", gse246215_registry, gse246215_weights),
  load_combined_counts_and_weights("GSE167339", gse167339_registry, gse167339_weights)
)
programs <- do.call(rbind, lapply(loaded, state_programs))
write_tsv(programs, file.path(OUT, "state_associated_top_genes.tsv"))

universe <- Reduce(union, lapply(loaded, function(x) rownames(x$counts)))
enrichment <- enrich_programs(programs, universe)
write_tsv(enrichment, file.path(OUT, "state_program_pathway_enrichment.tsv"))

abundance <- rbind(
  state_abundance("GSE246215", gse246215_weights, gse246215_registry),
  state_abundance("GSE167339", gse167339_weights, gse167339_registry)
)
write_tsv(abundance, file.path(OUT, "donor_state_abundance_fpi.tsv"))

effects <- donor_effects(abundance)
write_tsv(effects, file.path(OUT, "donor_aware_state_effects.tsv"))

manifest <- data.frame(
  analysis = "bioinformatics_validation",
  primary_claim = "FibroDynMix real-data outputs are supported by fibroblast purity QC, state-associated gene/pathway annotation, and donor-aware state abundance/FPI summaries.",
  claim_boundary = "Curated marker and pathway checks are lightweight bioinformatics validation, not exhaustive cell-type annotation, full GO enrichment, or causal disease biology.",
  n_validation_datasets = length(unique(purity$validation_dataset)),
  n_qc_dataset_rows = nrow(purity),
  min_purity_margin = min(purity$purity_margin_mean, na.rm = TRUE),
  max_low_purity_fraction = max(purity$low_purity_fraction, na.rm = TRUE),
  n_state_program_rows = nrow(programs),
  n_enrichment_rows = nrow(enrichment),
  n_enriched_state_pathways_q05 = sum(enrichment$q_value <= 0.05, na.rm = TRUE),
  n_abundance_rows = nrow(abundance),
  n_donor_effect_rows = nrow(effects),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "bioinformatics_validation_manifest.tsv"))

message("Bioinformatics validation written to: ", normalizePath(OUT))
