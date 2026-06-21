#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_project_maturity_audit.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "project_maturity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

rbind_fill <- function(rows) {
  all_cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_cols, colnames(row))
    for (col in missing) {
      row[[col]] <- NA
    }
    row[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

file_exists_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0
}

namespace_lines <- readLines(file.path(ROOT, "NAMESPACE"), warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace_lines, value = TRUE))
man_files <- list.files(file.path(ROOT, "man"), pattern = "\\.Rd$", full.names = TRUE)
r_files <- list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE)
test_files <- list.files(file.path(ROOT, "tests", "testthat"), pattern = "\\.R$", full.names = TRUE)
script_files <- list.files(file.path(ROOT, "scripts"), pattern = "\\.(R|py)$", full.names = TRUE)
doc_files <- list.files(file.path(ROOT, "docs"), pattern = "\\.md$", full.names = TRUE)

function_source <- data.frame(
  function_name = exports,
  source_file = NA_character_,
  man_page = file.path("man", paste0(exports, ".Rd")),
  has_man_page = file.exists(file.path(ROOT, "man", paste0(exports, ".Rd"))),
  stringsAsFactors = FALSE
)
for (i in seq_along(exports)) {
  pattern <- paste0("(^|[^A-Za-z0-9_.])", exports[i], "\\s*<-\\s*function")
  hits <- r_files[vapply(r_files, function(path) any(grepl(pattern, readLines(path, warn = FALSE))), logical(1))]
  if (length(hits) > 0L) {
    function_source$source_file[i] <- file.path("R", basename(hits[1]))
  }
}
write_tsv(function_source, file.path(OUT, "exported_function_index.tsv"))

figure_manifest_files <- list.files(file.path(ROOT, "figures"), pattern = "figure_manifest.tsv$", recursive = TRUE, full.names = TRUE)
figure_rows <- lapply(figure_manifest_files, function(path) {
  manifest <- read_tsv(path)
  manifest$figure_package <- basename(dirname(path))
  manifest$manifest_path <- sub(paste0("^", ROOT, "/"), "", path)
  manifest
})
figure_matrix <- rbind_fill(figure_rows)
write_tsv(figure_matrix, file.path(OUT, "figure_claim_matrix.tsv"))

evidence_matrix <- data.frame(
  claim_area = c(
    "Raw-count generative model",
    "Latent simplex optimization",
    "Study and donor hierarchy",
    "Posterior uncertainty",
    "VI calibration",
    "Extended method benchmark",
    "Validation-aware NB selection",
    "scVI-inclusive core benchmark",
    "Marker stress benchmark",
    "Transition flow and FPI",
    "Cross-cohort transfer",
    "Multi-public real-data validation",
    "Independent human GEO validation",
    "Independent human donor perturbation validation",
    "GSE167339 donor robustness",
    "Independent fibroblast dataset screening",
    "Figure 5 external biological validation",
    "Selected-NB final summaries",
    "Bioinformatics validation",
    "GSE246215 sensitivity and QC",
    "Public raw-count execution",
    "Reproducibility audit",
    "Runtime and dependency lock",
    "Package maturity",
    "Manuscript figure reproducibility"
  ),
  implemented_evidence = c(
    "fibrodynmix_nb_loglik(), fibrodynmix_nb_objective(), fit_fibrodynmix_nb()",
    "fit_fibrodynmix_initializer(), fit_fibrodynmix_nb(), fit_fibrodynmix_vi()",
    "fit_fibrodynmix_nb(study_id, donor_id), study/donor effect diagnostics",
    "bootstrap_fibrodynmix(), fit_fibrodynmix_vi(), posterior interval summaries",
    "evaluate_posterior_intervals(), calibrate_posterior_interval_scale(), analysis/vi_benchmark",
    "fit_topic_nmf_baseline(), run_extended_method_benchmark.R, marker scoring, NMF/topic baseline, FibroDynMix NB/VI comparisons",
    "select_fibrodynmix_nb_model(), run_validation_aware_nb_selection.R, held-out likelihood, z stability, marker-gradient preservation, and downstream validation",
    "run_core_converged_method_benchmark.R with include_scvi, project-local scvi-tools/anndata/torch environment, and scVI latent projection baseline",
    "run_marker_stress_benchmark.R corrupted, missing, shared, swapped, and hidden-program marker-prior stress modes",
    "compute_state_cost(), estimate_transition_flow(), compute_fpi(), Figure 6",
    "fit_fibrodynmix_transfer(), run_cross_cohort_transfer_benchmark(), public real-data transfer smoke",
    "run_multi_public_realdata_validation.R pooled public fit, study-effect adjustment, transition flow, and leave-dataset-out transfer",
    "prepare_gse246215_fibroblast_inputs.R plus run_multi_public_realdata_validation.R on GSE246215 human fibroblast count matrices",
    "prepare_gse167339_human_fibroblast_inputs.R plus run_multi_public_realdata_validation.R on GSE167339 human fibroblast MTX samples",
    "run_gse167339_donor_robustness.R donor composition, leave-donor-out transfer, and hash-threshold sensitivity",
    "audit_independent_fibroblast_datasets.R ranks GSE167339, GSE156326, GSE181316, GSE246215, scar sensitivity datasets, and bulk endpoint support",
    "make_figure5.R integrates GSE167339 donor/transfer evidence with GSE156326/GSE181316 scar state composition and module validation",
    "run_selected_nb_final_summaries.R reruns VI, transfer, and GSE246215 downstream summaries using validation-aware selected NB defaults",
    "run_bioinformatics_validation.R fibroblast purity QC, state-associated genes/pathways, donor-aware abundance, and FPI",
    "run_gse246215_sensitivity.R random downsampling, library-size q95 trimming, transfer convergence, and composition variability",
    "run_public_realdata_smoke.R, run_public_realdata_transfer.R, Figure 2, public smoke figure",
    "run_reproducibility_audit.R analysis catalog, script inventory, documentation inventory, and runbook",
    "run_runtime_lock.R package dependency lock, R session info, external data lock, and heavy-analysis I/O contract",
    "NEWS.md, inst/CITATION, cran-comments.md, _pkgdown.yml, function index, R CMD check OK",
    "figures/figure1-6, public_realdata_smoke, panel source-data manifests, export QC"
  ),
  primary_outputs = c(
    "R/nb_likelihood.R; R/fit_nb_model.R",
    "R/fibrodynmix_initializer.R; R/fit_nb_model.R; R/vi_posterior.R",
    "analysis/study_effect_sensitivity; figures/figure4",
    "analysis/bootstrap_uncertainty; analysis/vi_posterior; figures/supplementary_uncertainty",
    "analysis/vi_benchmark",
    "analysis/extended_method_benchmark; docs/extended-method-benchmark.md",
    "R/nb_model_selection.R; analysis/validation_aware_nb_selection; docs/core-converged-method-benchmark.md",
    "R/scvi_baseline.R; analysis/core_converged_method_benchmark_scvi_r5; figures/core_converged_tradeoff_scvi_r5",
    "analysis/marker_stress_benchmark; docs/marker-stress-benchmark.md",
    "analysis/transition_flow; figures/figure6",
    "analysis/cross_cohort_transfer; analysis/public_realdata_transfer",
    "analysis/multi_public_realdata_validation; docs/multi-public-realdata-validation.md",
    "data/public_geo_gse246215_fibroblast_atlas; analysis/independent_geo_gse246215_validation; docs/independent-geo-gse246215-validation.md",
    "data/public_geo_gse167339_human_fibroblast; analysis/independent_geo_gse167339_validation; docs/independent-geo-gse167339-validation.md",
    "analysis/gse167339_donor_robustness; docs/gse167339-donor-robustness.md",
    "analysis/independent_fibroblast_dataset_screening; docs/iscience-method-biological-execution.md",
    "figures/figure5; docs/figure5-external-validation.md",
    "analysis/selected_nb_final_summaries; docs/selected-nb-final-summaries.md",
    "analysis/bioinformatics_validation; docs/bioinformatics-validation.md",
    "analysis/gse246215_sensitivity; docs/gse246215-sensitivity.md",
    "analysis/public_realdata_smoke; analysis/public_realdata_figure2; figures/figure2",
    "analysis/reproducibility_audit; docs/reproducibility-runbook.md",
    "analysis/runtime_lock; docs/package-runtime-lock.md",
    "NEWS.md; inst/CITATION; docs/package-function-index.md; inst/doc/fibrodynmix-workflow.md",
    "figures/*/source_data; figures/*/qc/export_image_qc.tsv"
  ),
  claim_boundary = c(
    "Current optimizer is a pragmatic NB implementation, not full joint Bayesian inference.",
    "Latent z is inferred, but full amortized inference remains future work.",
    "Ridge hierarchy effects are implemented; full posterior over effects remains future work.",
    "Bootstrap and lightweight VI posterior are implemented; full all-parameter posterior remains future work.",
    "Calibration uses simulated truth and should not be presented as real-data calibration.",
    "Simulation method comparison; marker scoring remains competitive in marker-aligned low-noise settings and should not be dismissed.",
    "Selector is simulation-calibrated; truth metrics are audit columns only and final real-data reruns still need selected defaults.",
    "scVI is included as a bounded latent baseline; it is not a full scANVI or cNMF replacement and should not be overinterpreted from small simulated runs.",
    "Stress benchmark supports bounded failure-mode claims, not universal dominance over marker scoring.",
    "Flow is cross-sectional optimal transport, not lineage tracing.",
    "Transfer mechanics are implemented; completed human cross-atlas generalization remains future work.",
    "Default registry uses two public condition-specific count matrices from one Dryad record; a completed human disease-atlas claim remains future work.",
    "GSE246215 uses public processed count matrices and sampled cancer-type subsets; raw FASTQ-level reprocessing and definitive cancer biology claims remain future work.",
    "GSE167339 uses bounded sampled 10x-style processed MTX matrices from three donors; Human 3 hash groups are included as hash_unknown without treatment-label interpretation.",
    "GSE167339 donor robustness supports transfer and threshold stability, not Human3 treatment-specific biology.",
    "Readiness audit does not rerun heavy fitting and does not convert scar-project outputs into FibroDynMix claims without figure/source-data integration.",
    "External Figure 5 supports reproducibility and biological utility; GSE156326/GSE181316 module contrasts are descriptive small-sample effects and do not pass FDR significance.",
    "Selected-NB final summaries are simulation-calibrated and preserve baseline-competitive language for GSE246215.",
    "Bioinformatics validation uses curated marker/pathway checks and donor summaries, not exhaustive annotation or causal biology.",
    "QC sensitivity supports computational robustness but HCC remains a library-size outlier and should not be overinterpreted biologically.",
    "Public analyses are technical smoke tests, not disease-mechanism claims.",
    "Audit checks file presence and metadata coverage; it does not rerun every heavy analysis.",
    "Runtime lock is a package-level snapshot and analysis contract, not a full renv.lock or container image.",
    "Package is development version 0.0.0.9000, not yet submitted to CRAN/Bioconductor.",
    "Figures are manuscript-facing computational outputs; external biological panels remain bounded public count-matrix validations, not standalone disease-atlas proof."
  ),
  target_relevance = c(
    "Core distinction from gene-set scoring for iScience/Communications Biology.",
    "Supports mixed-state plasticity claim.",
    "Addresses cross-cohort confounding reviewers will expect.",
    "Supports uncertainty reporting required for mature probabilistic methods.",
    "Quantifies posterior calibration rather than asserting uncertainty qualitatively.",
    "Adds a mature-package NMF/topic comparator and clarifies what FibroDynMix improves beyond ordinary topic models.",
    "Addresses objective-only NB stopping as a reviewer risk and gives a defensible default-selection rule.",
    "Closes the deep generative baseline gap for the core simulation benchmark.",
    "Clarifies when marker scoring fails and when FibroDynMix provides competitive or improved recovery plus richer model outputs.",
    "Connects state mixtures to fibroblast plasticity hypothesis.",
    "Supports cross-dataset migration claim with bounded evidence.",
    "Moves the real-data layer beyond one pooled smoke run and gives reviewers explicit leave-dataset-out diagnostics.",
    "Adds an independent human multi-patient fibroblast atlas source beyond the mouse Dryad validation.",
    "Adds an independent human donor-level perturbation dataset beyond the cancer atlas validation.",
    "Converts the GSE167339 real-data layer into explicit donor-level reproducibility evidence.",
    "Identifies the external-data path for a method-first submission with bounded biological-utility evidence.",
    "Promotes independent public human fibroblast data into a main-figure-level biological validation panel.",
    "Anchors final summary claims to explicit model-selection defaults instead of fixed smoke settings.",
    "Shows input fibroblast purity, biological interpretability of state programs, and donor-aware state/FPI shifts.",
    "Addresses real-data robustness and library-size confounding concerns without relying on new figures.",
    "Shows raw-count execution outside simulation.",
    "Supports reproducible research handoff and manuscript evidence traceability.",
    "Supports dependency/runtime review and public-data handoff for a mature R package.",
    "Supports mature R package review and reproducibility expectations.",
    "Supports source-backed figure review."
  ),
  stringsAsFactors = FALSE
)
write_tsv(evidence_matrix, file.path(OUT, "manuscript_evidence_matrix.tsv"))

readiness <- data.frame(
  item = c(
    "R CMD check",
    "Unit tests",
    "Integrity gate",
    "Function documentation",
    "Citation metadata",
    "Release notes",
    "Figure source data",
    "Public data claim boundary",
    "Multi-public real-data validation",
    "Independent human GEO validation",
    "Independent human donor perturbation validation",
    "GSE167339 donor robustness",
    "Bioinformatics validation",
    "GSE246215 sensitivity and QC",
    "Extended method benchmark",
    "Validation-aware NB selection",
    "scVI-inclusive core benchmark",
    "Marker stress benchmark",
    "Independent fibroblast dataset screening",
    "Figure 5 external biological validation",
    "Selected-NB final summaries",
    "Reproducibility audit",
    "Runtime and dependency lock",
    "VI calibration claim boundary",
    "Remaining atlas-level validation",
    "Manuscript skeleton",
    "Reviewer risk matrix",
    "Manuscript draft"
  ),
  status = c(
    "pass_when_last_run",
    "pass_when_last_run",
    "pass_when_last_run",
    ifelse(all(function_source$has_man_page), "complete", "incomplete"),
    ifelse(file_exists_nonempty(file.path(ROOT, "inst", "CITATION")), "complete", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "NEWS.md")), "complete", "missing"),
    ifelse(all(file.exists(file.path(dirname(figure_manifest_files), "panel_source_data_manifest.tsv"))), "complete", "incomplete"),
    "bounded",
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "multi_public_realdata_validation", "multi_public_validation_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "independent_geo_gse246215_validation", "multi_public_validation_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_validation_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "gse167339_donor_robustness", "gse167339_donor_robustness_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "bioinformatics_validation", "bioinformatics_validation_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "gse246215_sensitivity", "gse246215_sensitivity_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "extended_method_benchmark", "extended_method_benchmark_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "validation_aware_nb_selection", "validation_aware_nb_selection_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "core_converged_method_benchmark_scvi_r5", "core_converged_benchmark_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "marker_stress_benchmark", "marker_stress_benchmark_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "independent_fibroblast_dataset_screening", "independent_fibroblast_dataset_screening_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "figures", "figure5", "figure_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "selected_nb_final_summaries", "selected_nb_final_summaries_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "reproducibility_audit", "reproducibility_audit_manifest.tsv")), "complete_bounded", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "analysis", "runtime_lock", "runtime_lock_manifest.tsv")), "complete_bounded", "missing"),
    "bounded",
    "remaining",
    ifelse(file_exists_nonempty(file.path(ROOT, "docs", "manuscript-skeleton-isci-commbio.md")), "complete", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "docs", "reviewer-risk-matrix.md")), "complete", "missing"),
    ifelse(file_exists_nonempty(file.path(ROOT, "docs", "manuscript-draft-isci-commbio.md")), "complete", "missing")
  ),
  evidence = c(
    "R CMD check --no-manual reported Status: OK in the latest validation run.",
    "testthat::test_local() reported all tests passing in the latest validation run.",
    "scripts/check_project_integrity.R reported project integrity passed.",
    "NAMESPACE exports and man/*.Rd files are indexed in exported_function_index.tsv.",
    "inst/CITATION",
    "NEWS.md",
    "figures/*/panel_source_data_manifest.tsv and figures/*/source_data",
    "figure manifests and public real-data docs state smoke-test boundaries.",
    "analysis/multi_public_realdata_validation and docs/multi-public-realdata-validation.md",
    "analysis/independent_geo_gse246215_validation and docs/independent-geo-gse246215-validation.md",
    "analysis/independent_geo_gse167339_validation and docs/independent-geo-gse167339-validation.md",
    "analysis/gse167339_donor_robustness and docs/gse167339-donor-robustness.md",
    "analysis/bioinformatics_validation and docs/bioinformatics-validation.md",
    "analysis/gse246215_sensitivity and docs/gse246215-sensitivity.md",
    "analysis/extended_method_benchmark and docs/extended-method-benchmark.md",
    "analysis/validation_aware_nb_selection and R/nb_model_selection.R",
    "analysis/core_converged_method_benchmark_scvi_r5 and R/scvi_baseline.R",
    "analysis/marker_stress_benchmark and docs/marker-stress-benchmark.md",
    "analysis/independent_fibroblast_dataset_screening and docs/iscience-method-biological-execution.md",
    "figures/figure5 and docs/figure5-external-validation.md",
    "analysis/selected_nb_final_summaries and docs/selected-nb-final-summaries.md",
    "analysis/reproducibility_audit and docs/reproducibility-runbook.md",
    "analysis/runtime_lock and docs/package-runtime-lock.md",
    "docs/vi-posterior.md and analysis/vi_benchmark state simulation-calibration boundary.",
    "Completed human disease-atlas and cross-atlas generalization remain outside current evidence.",
    "docs/manuscript-skeleton-isci-commbio.md",
    "docs/reviewer-risk-matrix.md",
    "docs/manuscript-draft-isci-commbio.md"
  ),
  stringsAsFactors = FALSE
)
write_tsv(readiness, file.path(OUT, "submission_readiness_checklist.tsv"))

description <- read.dcf(file.path(ROOT, "DESCRIPTION"))[1, ]
manifest <- data.frame(
  package = unname(description[["Package"]]),
  version = unname(description[["Version"]]),
  n_exports = length(exports),
  n_man_pages = length(man_files),
  n_r_files = length(r_files),
  n_test_files = length(test_files),
  n_scripts = length(script_files),
  n_docs = length(doc_files),
  n_figure_packages = length(figure_manifest_files),
  n_analysis_dirs = length(list.dirs(file.path(ROOT, "analysis"), recursive = FALSE)),
  has_news = file_exists_nonempty(file.path(ROOT, "NEWS.md")),
  has_citation = file_exists_nonempty(file.path(ROOT, "inst", "CITATION")),
  has_pkgdown_config = file_exists_nonempty(file.path(ROOT, "_pkgdown.yml")),
  all_exports_documented = all(function_source$has_man_page),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "project_maturity_manifest.tsv"))

message("Project maturity audit written to: ", normalizePath(OUT))
