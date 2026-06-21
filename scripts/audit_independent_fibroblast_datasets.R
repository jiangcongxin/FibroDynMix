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
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/audit_independent_fibroblast_datasets.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "independent_fibroblast_dataset_screening"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(data.frame())
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

file_exists_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0L
}

rel <- function(path) {
  sub(paste0("^", ROOT, "/"), "", path)
}

present_paths <- function(paths) {
  paths <- paths[vapply(paths, file_exists_nonempty, logical(1))]
  if (length(paths) == 0L) {
    return(NA_character_)
  }
  paste(vapply(paths, rel, character(1)), collapse = ";")
}

first_value <- function(df, key_col, key, value_col, default = NA) {
  if (nrow(df) == 0L || !all(c(key_col, value_col) %in% colnames(df))) {
    return(default)
  }
  hit <- df[[value_col]][df[[key_col]] == key]
  hit <- hit[!is.na(hit) & nzchar(as.character(hit))]
  if (length(hit) == 0L) {
    return(default)
  }
  hit[[1]]
}

sum_value <- function(df, key_col, key, value_col, default = NA_real_) {
  if (nrow(df) == 0L || !all(c(key_col, value_col) %in% colnames(df))) {
    return(default)
  }
  hit <- suppressWarnings(as.numeric(df[[value_col]][df[[key_col]] == key]))
  hit <- hit[is.finite(hit)]
  if (length(hit) == 0L) {
    return(default)
  }
  sum(hit)
}

candidate_summary <- read_tsv(file.path(
  ROOT,
  "projects",
  "jtm_pathological_scar_fibroblast_states",
  "results",
  "fibroblast_extraction_final",
  "fibroblast_candidate_dataset_summary.tsv"
))
scar_transfer <- read_tsv(file.path(
  ROOT,
  "projects",
  "jtm_pathological_scar_fibroblast_states",
  "results",
  "state_modeling",
  "cross_cohort_transfer_final",
  "transfer_diagnostics.tsv"
))
scar_module_summary <- read_tsv(file.path(
  ROOT,
  "projects",
  "jtm_pathological_scar_fibroblast_states",
  "results",
  "disease_module",
  "module_score_validation",
  "scar_fibroblast_module_validation_summary.tsv"
))
dataset_intake <- read_tsv(file.path(
  ROOT,
  "projects",
  "jtm_pathological_scar_fibroblast_states",
  "results",
  "dataset_intake",
  "combined_dataset_manifest.tsv"
))
gse246215_manifest <- read_tsv(file.path(
  ROOT,
  "analysis",
  "independent_geo_gse246215_validation",
  "multi_public_validation_manifest.tsv"
))
gse167339_manifest <- read_tsv(file.path(
  ROOT,
  "analysis",
  "independent_geo_gse167339_validation",
  "multi_public_validation_manifest.tsv"
))
gse167339_donor_manifest <- read_tsv(file.path(
  ROOT,
  "analysis",
  "gse167339_donor_robustness",
  "gse167339_donor_robustness_manifest.tsv"
))
core_scvi_manifest <- read_tsv(file.path(
  ROOT,
  "analysis",
  "core_converged_method_benchmark_scvi_r5",
  "core_converged_benchmark_manifest.tsv"
))
selector_manifest <- read_tsv(file.path(
  ROOT,
  "analysis",
  "validation_aware_nb_selection",
  "validation_aware_nb_selection_manifest.tsv"
))

local_status <- function(paths) {
  if (any(vapply(paths, file_exists_nonempty, logical(1)))) {
    return("local_evidence_available")
  }
  "not_found_locally"
}

manifest_n_cells <- function(manifest) {
  if (nrow(manifest) == 0L || !"n_cells" %in% colnames(manifest)) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(manifest$n_cells[[1]]))
  if (is.finite(out)) {
    return(out)
  }
  NA_real_
}

screened_n_cells <- function(dataset_id) {
  candidate_cells <- sum_value(candidate_summary, "dataset_id", dataset_id, "n_cells")
  if (is.finite(candidate_cells)) {
    return(candidate_cells)
  }
  if (dataset_id == "gse246215") {
    return(manifest_n_cells(gse246215_manifest))
  }
  if (dataset_id == "gse167339") {
    return(manifest_n_cells(gse167339_manifest))
  }
  NA_real_
}

row <- function(priority_rank, dataset_id, accession, disease_context, organism, modality,
                submission_role, independence_class, current_status, biological_question,
                endpoint, limitations, recommended_next_action, figure_candidate,
                evidence_paths, readiness_score, iScience_value) {
  data.frame(
    priority_rank = priority_rank,
    dataset_id = dataset_id,
    accession = accession,
    disease_context = disease_context,
    organism = organism,
    modality = modality,
    submission_role = submission_role,
    independence_class = independence_class,
    current_status = current_status,
    local_status = local_status(evidence_paths),
    n_cells = screened_n_cells(dataset_id),
    n_fibroblast_candidates = sum_value(candidate_summary, "dataset_id", dataset_id, "n_fibroblast_candidates"),
    transfer_convergence_rate = first_value(scar_transfer, "dataset_id", dataset_id, "z_convergence_rate"),
    module_delta_positive_minus_negative = first_value(
      scar_module_summary,
      "dataset_id",
      dataset_id,
      "delta_positive_minus_negative"
    ),
    module_auc_positive_high = first_value(scar_module_summary, "dataset_id", dataset_id, "auc_positive_high"),
    biological_question = biological_question,
    endpoint = endpoint,
    limitations = limitations,
    recommended_next_action = recommended_next_action,
    figure_candidate = figure_candidate,
    evidence_paths = present_paths(evidence_paths),
    readiness_score = readiness_score,
    iScience_value = iScience_value,
    stringsAsFactors = FALSE
  )
}

rows <- list(
  row(
    0,
    "gse246215",
    "GSE246215",
    "multi-cancer CAF atlas",
    "Homo sapiens",
    "scRNA-seq processed count matrix",
    "current cancer real-data case study",
    "same GEO record; multiple cancer-type subsets, not fully independent multi-study validation",
    "already analyzed with leave-dataset-out transfer, sensitivity/QC, and cancer-type downstream benchmark",
    "Do fibroblast-state compositions encode cancer-context differences while preserving marker gradients?",
    "Cancer-type classification and marker-gradient validation.",
    "Useful but insufficient as the only biological validation because all subsets come from one GEO study and HCC remains QC-sensitive.",
    "Keep as main cancer case study, but explicitly pair with at least one non-cancer external fibroblast dataset.",
    "Figure 4",
    c(
      file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas", "gse246215_fibroblast_dataset_manifest.tsv"),
      file.path(ROOT, "analysis", "independent_geo_gse246215_validation", "multi_public_validation_manifest.tsv"),
      file.path(ROOT, "analysis", "gse246215_downstream_benchmark", "gse246215_downstream_classification_metrics.tsv")
    ),
    4,
    "medium"
  ),
  row(
    1,
    "gse167339",
    "GSE167339",
    "human dermal fibroblast collagen scar/strain perturbation",
    "Homo sapiens",
    "10x-style scRNA-seq MTX",
    "promote now as independent non-cancer validation",
    "independent GEO study; separate disease/perturbation context from GSE246215",
    "already prepared and validated with leave-dataset-out plus leave-donor-out transfer",
    "Are fibroblast-state gradients transferable to a non-cancer scar/strain context?",
    "Donor-level state composition, leave-donor-out transfer, and treatment/strain gradient where labels are available.",
    "Human3 hash groups lack public treatment mapping, so Human3 should support donor robustness only.",
    "Promote to external validation panel and rerun selected-NB/transfer summary if final selector changes the default.",
    "Figure 5",
    c(
      file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_human_fibroblast_dataset_manifest.tsv"),
      file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_validation_manifest.tsv"),
      file.path(ROOT, "analysis", "gse167339_donor_robustness", "gse167339_donor_robustness_manifest.tsv")
    ),
    5,
    "high"
  ),
  row(
    2,
    "gse156326",
    "GSE156326",
    "hypertrophic scar versus normal skin",
    "Homo sapiens",
    "scRNA-seq",
    "highest-priority new external scar validation",
    "independent GEO study; non-cancer disease context",
    "fibroblast candidates and final transfer outputs available locally",
    "Does FibroDynMix recover ECM-remodeling/inflammatory/resident-state shifts in hypertrophic scar?",
    "Disease versus control pseudobulk composition and module score validation.",
    "Current module result is positive-direction but small donor/sample count; treat as external validation, not clinical diagnostic evidence.",
    "Integrate into FibroDynMix manuscript as the first independent non-cancer dataset and regenerate concise state-composition/gradient panels.",
    "Figure 5",
    c(
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "fibroblast_extraction_final", "gse156326_fibroblast_candidates_final.rds"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "state_modeling", "cross_cohort_transfer_final", "gse156326_transfer_state_composition.tsv"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "disease_module", "module_score_validation", "scar_fibroblast_module_validation_summary.tsv")
    ),
    5,
    "high"
  ),
  row(
    3,
    "gse181316",
    "GSE181316",
    "keloid, normal scar, and normal skin",
    "Homo sapiens",
    "scRNA-seq",
    "second high-priority external scar validation",
    "independent GEO study; non-cancer disease context",
    "fibroblast candidates and final transfer outputs available locally",
    "Does the inferred fibroblast plasticity axis separate keloid from normal scar/skin?",
    "Disease versus control pseudobulk composition, module score, and external transfer diagnostics.",
    "Strong local signal but still observational; claims should remain composition/gradient, not mechanism or therapy.",
    "Use as the stronger scar validation after GSE156326, with disease-control statistics and state-gradient figure source data.",
    "Figure 5 or Figure 6",
    c(
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "fibroblast_extraction_final", "gse181316_fibroblast_candidates_final.rds"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "state_modeling", "cross_cohort_transfer_final", "gse181316_transfer_state_composition.tsv"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "disease_module", "module_score_validation", "scar_fibroblast_module_validation_summary.tsv")
    ),
    5,
    "high"
  ),
  row(
    4,
    "gse163973",
    "GSE163973",
    "keloid versus normal scar/skin reference source",
    "Homo sapiens",
    "scRNA-seq",
    "possible discovery/reference source",
    "independent public dataset, but not independent if used to derive scar programs",
    "local raw/meta and FibroDynMix project outputs available",
    "Can a reference FibroDynMix scar-state basis be learned and transferred to other scar datasets?",
    "Reference fit, state definitions, and transfer basis.",
    "If used to define programs, it cannot also serve as independent validation for those programs.",
    "Use either as discovery/reference or exclude from validation claims; document the choice before manuscript drafting.",
    "Methods/Figure 5 reference panel",
    c(
      file.path(ROOT, "projects", "keloid_fibro_gse163973", "data", "GSE163973_RAW.tar"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "state_modeling", "cross_cohort_transfer_final", "gse163973_reference_fit_final_capped.rds")
    ),
    3,
    "medium"
  ),
  row(
    5,
    "gse175866",
    "GSE175866",
    "sorted keloid fibroblast subsets",
    "Homo sapiens",
    "bulk RNA-seq",
    "external endpoint support",
    "independent bulk fibroblast dataset; not single-cell composition validation",
    "bulk signature/module validation available locally",
    "Do FibroDynMix-linked programs have external expression support in sorted fibroblast subsets?",
    "Subset-level module score contrast.",
    "Only two sorted profiles in the current summary; descriptive support only, not inferential validation.",
    "Use as supplementary external endpoint support, not as the minimum independent scRNA-seq validation.",
    "Supplementary Figure",
    c(
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "bulk_validation", "gse175866", "gse175866_signature_validation_manifest.tsv"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "bulk_validation", "gse175866", "gse175866_signature_scores.tsv")
    ),
    3,
    "medium"
  ),
  row(
    6,
    "gse181297",
    "GSE181297",
    "keloid spatial/scar context",
    "Homo sapiens",
    "single-cell/spatial-context public data",
    "supporting sensitivity dataset",
    "independent GEO study; limited control structure in current local summary",
    "fibroblast candidates and final transfer outputs available locally",
    "Is the inferred scar-state axis recoverable in an additional keloid-context dataset?",
    "Transfer stability and state composition.",
    "Current summary has disease samples without matched control rows, so endpoint evidence is incomplete.",
    "Keep as supplementary transfer stability evidence unless control labels are curated.",
    "Supplementary Figure",
    c(
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "fibroblast_extraction_final", "gse181297_fibroblast_candidates_final.rds"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "state_modeling", "cross_cohort_transfer_final", "gse181297_transfer_state_composition.tsv")
    ),
    3,
    "medium"
  ),
  row(
    7,
    "gse243716",
    "GSE243716",
    "keloid versus hypertrophic scar",
    "Homo sapiens",
    "scRNA-seq",
    "supporting sensitivity dataset",
    "independent GEO study; two-condition scar comparison",
    "fibroblast candidates, project outputs, and separate FibroDynMix workflow available locally",
    "Can the state axis distinguish keloid from hypertrophic scar in a compact independent dataset?",
    "Transfer stability and keloid-versus-hypertrophic composition contrast.",
    "Small two-sample structure limits statistical claims.",
    "Use as sensitivity/supporting panel after GSE156326 and GSE181316.",
    "Supplementary Figure",
    c(
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "fibroblast_extraction_final", "gse243716_fibroblast_candidates_final.rds"),
      file.path(ROOT, "projects", "jtm_pathological_scar_fibroblast_states", "results", "state_modeling", "cross_cohort_transfer_final", "gse243716_transfer_state_composition.tsv"),
      file.path(ROOT, "projects", "scar_fibro_gse243716", "results", "gse243716_fibrodynmix_workflow.rds")
    ),
    3,
    "medium"
  ),
  row(
    8,
    "dryad_breast_fibroblast",
    "Dryad breast fibroblast count matrices",
    "mouse normal mammary fibroblast versus CAF",
    "Mus musculus",
    "public count matrix",
    "low-priority supporting public execution",
    "independent public record but non-human and limited condition design",
    "already analyzed in multi-public real-data validation",
    "Can the raw-count workflow execute on an additional public fibroblast count source?",
    "Technical public-data transfer and state composition.",
    "Non-human and narrow source; weak support for an iScience human disease-atlas claim.",
    "Retain as technical reproducibility evidence, not a central biological validation.",
    "Supplementary Methods",
    c(
      file.path(ROOT, "data", "public_dryad_breast_fibroblast", "dryad_breast_fibroblast_counts.rds"),
      file.path(ROOT, "analysis", "multi_public_realdata_validation", "multi_public_validation_manifest.tsv")
    ),
    2,
    "low"
  )
)

readiness <- do.call(rbind, rows)
readiness <- readiness[order(readiness$priority_rank), , drop = FALSE]
write_tsv(readiness, file.path(OUT, "independent_fibroblast_dataset_readiness.tsv"))

next_actions <- readiness[readiness$readiness_score >= 4 | readiness$priority_rank %in% c(1, 2, 3), , drop = FALSE]
next_actions <- next_actions[, c(
  "priority_rank",
  "dataset_id",
  "accession",
  "submission_role",
  "recommended_next_action",
  "figure_candidate",
  "limitations",
  "evidence_paths"
), drop = FALSE]
write_tsv(next_actions, file.path(OUT, "independent_fibroblast_dataset_next_actions.tsv"))

manifest <- data.frame(
  analysis = "independent_fibroblast_dataset_screening",
  primary_claim = "Independent real-data expansion should promote GSE167339 and scar/fibrosis-context human fibroblast datasets, while keeping GSE246215 as a cancer case study and acknowledging NMF/marker competitiveness.",
  claim_boundary = "Readiness audit only. It does not rerun heavy real-data fitting and does not convert scar-project outputs into FibroDynMix manuscript claims without dedicated figure/source-data integration.",
  n_screened_datasets = nrow(readiness),
  n_high_value_ready = sum(readiness$iScience_value == "high" & readiness$local_status == "local_evidence_available"),
  highest_priority_external_scrna = paste(readiness$dataset_id[readiness$priority_rank %in% c(1, 2, 3)], collapse = ";"),
  validation_aware_selector_available = file_exists_nonempty(file.path(ROOT, "R", "nb_model_selection.R")) &&
    file_exists_nonempty(file.path(ROOT, "analysis", "validation_aware_nb_selection", "validation_aware_nb_selection_manifest.tsv")),
  selector_selected_n_outer_values = if (nrow(selector_manifest) > 0L && "selected_n_outer_values" %in% colnames(selector_manifest)) selector_manifest$selected_n_outer_values[[1]] else NA_character_,
  scvi_core_benchmark_available = if (nrow(core_scvi_manifest) > 0L && "include_scvi" %in% colnames(core_scvi_manifest)) as.character(core_scvi_manifest$include_scvi[[1]]) else "FALSE",
  gse246215_public_datasets = if (nrow(gse246215_manifest) > 0L && "n_public_datasets" %in% colnames(gse246215_manifest)) gse246215_manifest$n_public_datasets[[1]] else NA,
  gse167339_public_datasets = if (nrow(gse167339_manifest) > 0L && "n_public_datasets" %in% colnames(gse167339_manifest)) gse167339_manifest$n_public_datasets[[1]] else NA,
  gse167339_leave_donor_runs = if (nrow(gse167339_donor_manifest) > 0L && "n_leave_donor_out_runs" %in% colnames(gse167339_donor_manifest)) gse167339_donor_manifest$n_leave_donor_out_runs[[1]] else NA,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "independent_fibroblast_dataset_screening_manifest.tsv"))

message("Independent fibroblast dataset screening written to: ", normalizePath(OUT))
