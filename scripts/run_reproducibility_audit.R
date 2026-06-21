#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_reproducibility_audit.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "reproducibility_audit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

file_exists_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

analysis_specs <- data.frame(
  analysis_id = c(
    "public_realdata_smoke",
    "public_realdata_transfer",
    "multi_public_realdata_validation",
    "independent_geo_gse246215_validation",
    "independent_geo_gse167339_validation",
    "gse167339_donor_robustness",
    "bioinformatics_validation",
    "gse246215_sensitivity",
    "gse246215_downstream_benchmark",
    "gse246215_downstream_benchmark_sample",
    "extended_method_benchmark",
    "nb_convergence_benchmark",
    "core_converged_method_benchmark",
    "marker_stress_benchmark",
    "vi_benchmark",
    "vi_posterior",
    "cross_cohort_transfer",
    "study_effect_sensitivity",
    "bootstrap_uncertainty",
    "transition_flow",
    "project_maturity"
  ),
  run_script = c(
    "scripts/run_public_realdata_smoke.R",
    "scripts/run_public_realdata_transfer.R",
    "scripts/run_multi_public_realdata_validation.R",
    "scripts/run_multi_public_realdata_validation.R --dataset-manifest=data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv --out=analysis/independent_geo_gse246215_validation",
    "scripts/prepare_gse167339_human_fibroblast_inputs.R --max-cells-per-sample=80 --include-hash-pool=true; scripts/run_multi_public_realdata_validation.R --dataset-manifest=data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv --out=analysis/independent_geo_gse167339_validation --max-cells=80 --max-genes=700 --transfer-maxit-z=120 --seed=167339",
    "scripts/run_gse167339_donor_robustness.R",
    "scripts/run_bioinformatics_validation.R",
    "scripts/run_gse246215_sensitivity.R",
    "scripts/run_gse246215_downstream_benchmark.R",
    "scripts/run_gse246215_downstream_benchmark.R --aggregation-col=SampleID --out=analysis/gse246215_downstream_benchmark_sample",
    "scripts/run_extended_method_benchmark.R",
    "scripts/run_nb_convergence_benchmark.R",
    "scripts/run_core_converged_method_benchmark.R",
    "scripts/run_marker_stress_benchmark.R",
    "scripts/run_vi_benchmark.R",
    "scripts/run_vi_posterior.R",
    "scripts/run_cross_cohort_transfer.R",
    "scripts/run_study_sensitivity.R",
    "scripts/run_bootstrap_uncertainty.R",
    "scripts/run_transition_flow.R",
    "scripts/run_project_maturity_audit.R"
  ),
  output_dir = file.path("analysis", c(
    "public_realdata_smoke",
    "public_realdata_transfer",
    "multi_public_realdata_validation",
    "independent_geo_gse246215_validation",
    "independent_geo_gse167339_validation",
    "gse167339_donor_robustness",
    "bioinformatics_validation",
    "gse246215_sensitivity",
    "gse246215_downstream_benchmark",
    "gse246215_downstream_benchmark_sample",
    "extended_method_benchmark",
    "nb_convergence_benchmark",
    "core_converged_method_benchmark",
    "marker_stress_benchmark",
    "vi_benchmark",
    "vi_posterior",
    "cross_cohort_transfer",
    "study_effect_sensitivity",
    "bootstrap_uncertainty",
    "transition_flow",
    "project_maturity"
  )),
  primary_manifest = c(
    "run_manifest.tsv",
    "public_transfer_manifest.tsv",
    "multi_public_validation_manifest.tsv",
    "multi_public_validation_manifest.tsv",
    "multi_public_validation_manifest.tsv",
    "gse167339_donor_robustness_manifest.tsv",
    "bioinformatics_validation_manifest.tsv",
    "gse246215_sensitivity_manifest.tsv",
    "gse246215_downstream_manifest.tsv",
    "gse246215_downstream_manifest.tsv",
    "extended_method_benchmark_manifest.tsv",
    "nb_convergence_manifest.tsv",
    "core_converged_benchmark_manifest.tsv",
    "marker_stress_benchmark_manifest.tsv",
    "vi_benchmark_manifest.tsv",
    "vi_manifest.tsv",
    "transfer_manifest.tsv",
    "study_effect_sensitivity_summary.tsv",
    "sample_composition_uncertainty.tsv",
    "transition_flow_summary.tsv",
    "project_maturity_manifest.tsv"
  ),
  research_role = c(
    "public raw-count execution smoke test",
    "bidirectional public transfer smoke test",
    "multi-public Dryad count-matrix validation",
    "independent human GSE246215 validation",
    "independent human GSE167339 donor perturbation validation",
    "GSE167339 donor-level robustness and hash-threshold sensitivity",
    "fibroblast purity, state program annotation, and donor-aware state/FPI validation",
    "GSE246215 downsampling and library-size QC sensitivity",
    "GSE246215 patient-level cancer-type downstream representation benchmark",
    "GSE246215 sample-level cancer-type downstream representation benchmark",
    "simulation comparison against marker scoring and NMF/topic baseline",
    "NB outer-iteration convergence adequacy benchmark",
    "core method rerun under n_outer 10/20",
    "simulation stress test for marker-prior failure modes",
    "simulation posterior interval calibration",
    "posterior state uncertainty demonstration",
    "leave-study-out simulation transfer benchmark",
    "study-effect penalty sensitivity",
    "bootstrap uncertainty analysis",
    "cross-sectional transition-flow and FPI analysis",
    "package and manuscript evidence audit"
  ),
  claim_boundary = c(
    "technical smoke test, not disease mechanism",
    "transfer mechanics, not full atlas generalization",
    "two condition-specific matrices from one Dryad record",
    "processed public count matrix and within-study cancer-type subsets",
    "processed 10x-style MTX matrices from three donors; Human 3 hash groups included as hash_unknown without treatment-label interpretation",
    "donor-level transfer robustness; Human3 remains hash_unknown without treatment-label interpretation",
    "curated marker/pathway checks; not exhaustive annotation or causal biology",
    "computational robustness and QC sensitivity, not cancer biology",
    "processed public count matrix from one GEO accession; representation utility, not diagnostic performance",
    "processed public count matrix from one GEO accession; sample-level representation utility, not diagnostic performance",
    "simulation truth only; marker scoring remains competitive in some settings",
    "bounded convergence sensitivity; n_outer=2 is not convergence-level evidence",
    "higher-outer rerun exposes objective/RMSE mismatch; not final optimizer selection",
    "bounded stress simulation; does not prove universal superiority over marker scoring",
    "simulation truth only; lightweight VI, not full posterior",
    "posterior skeleton around NB mode",
    "simulation transfer truth only",
    "simulation-calibrated penalty choice",
    "bootstrap uncertainty, not full Bayesian posterior",
    "cross-sectional OT flow, not lineage tracing",
    "audit summary, not independent validation"
  ),
  stringsAsFactors = FALSE
)

analysis_rows <- lapply(seq_len(nrow(analysis_specs)), function(i) {
  output_dir <- file.path(ROOT, analysis_specs$output_dir[i])
  manifest_path <- file.path(output_dir, analysis_specs$primary_manifest[i])
  files <- if (dir.exists(output_dir)) list.files(output_dir, recursive = TRUE, full.names = TRUE) else character()
  data.frame(
    analysis_id = analysis_specs$analysis_id[i],
    run_script = analysis_specs$run_script[i],
    output_dir = analysis_specs$output_dir[i],
    output_dir_exists = dir.exists(output_dir),
    n_output_files = length(files),
    primary_manifest = file.path(analysis_specs$output_dir[i], analysis_specs$primary_manifest[i]),
    primary_manifest_exists = file_exists_nonempty(manifest_path),
    research_role = analysis_specs$research_role[i],
    claim_boundary = analysis_specs$claim_boundary[i],
    stringsAsFactors = FALSE
  )
})
analysis_catalog <- do.call(rbind, analysis_rows)
write_tsv(analysis_catalog, file.path(OUT, "analysis_catalog.tsv"))

script_files <- list.files(file.path(ROOT, "scripts"), pattern = "\\.(R|py)$", full.names = TRUE)
script_inventory <- do.call(rbind, lapply(script_files, function(path) {
  text <- readLines(path, warn = FALSE)
  data.frame(
    script = sub(paste0("^", ROOT, "/"), "", path),
    language = tools::file_ext(path),
    n_lines = length(text),
    has_usage_args = any(grepl("get_arg|argparse|commandArgs", text)),
    sources_package_r = any(grepl("source\\(|source_files|source_package_files", text)),
    writes_analysis_output = any(grepl("write\\.table|write_tsv|saveRDS|ggsave", text)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(script_inventory, file.path(OUT, "script_inventory.tsv"))

doc_files <- list.files(file.path(ROOT, "docs"), pattern = "\\.md$", full.names = TRUE)
doc_inventory <- do.call(rbind, lapply(doc_files, function(path) {
  text <- readLines(path, warn = FALSE)
  data.frame(
    document = sub(paste0("^", ROOT, "/"), "", path),
    n_lines = length(text),
    has_claim_boundary = any(grepl("Claim Boundary|claim boundary|Claim boundary", text)),
    mentions_outputs = any(grepl("analysis/|figures/|scripts/", text)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(doc_inventory, file.path(OUT, "doc_inventory.tsv"))

manifest <- data.frame(
  analysis = "reproducibility_audit",
  primary_claim = "FibroDynMix has a source-backed reproducibility map linking scripts, analysis outputs, documents, and claim boundaries.",
  claim_boundary = "Audit checks file presence and metadata coverage; it does not rerun every heavy analysis.",
  n_registered_analyses = nrow(analysis_catalog),
  n_existing_analysis_dirs = sum(analysis_catalog$output_dir_exists),
  n_primary_manifests = sum(analysis_catalog$primary_manifest_exists),
  n_scripts = nrow(script_inventory),
  n_docs = nrow(doc_inventory),
  all_registered_outputs_present = all(analysis_catalog$output_dir_exists & analysis_catalog$primary_manifest_exists),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "reproducibility_audit_manifest.tsv"))

message("Reproducibility audit written to: ", normalizePath(OUT))
