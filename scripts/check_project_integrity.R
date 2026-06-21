#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/check_project_integrity.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))

required_figures <- c(
  "figure1",
  "figure2",
  "figure3",
  "figure4",
  "figure5",
  "figure6",
  "public_realdata_smoke"
)

failures <- character()
notes <- character()

add_failure <- function(...) {
  failures <<- c(failures, sprintf(...))
}

add_note <- function(...) {
  notes <<- c(notes, sprintf(...))
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

is_nonempty <- function(x) {
  length(x) > 0L && !all(is.na(x)) && all(nzchar(as.character(x)))
}

check_figure_package <- function(fig) {
  fig_dir <- file.path(ROOT, "figures", fig)
  exports_dir <- file.path(fig_dir, "exports")
  qc_path <- file.path(fig_dir, "qc", "export_image_qc.tsv")
  manifest_path <- file.path(fig_dir, "figure_manifest.tsv")
  panel_manifest_path <- file.path(fig_dir, "panel_source_data_manifest.tsv")
  legend_path <- file.path(fig_dir, "main_figure_legends.md")

  if (!dir.exists(fig_dir)) {
    add_failure("Missing figure directory: %s", fig_dir)
    return(invisible(FALSE))
  }

  stem <- if (fig == "public_realdata_smoke") "public_realdata_smoke" else fig
  required_exports <- file.path(exports_dir, paste0(stem, c(".pdf", ".svg", ".png", ".tiff")))
  required_exports <- c(required_exports, file.path(exports_dir, "contact_sheet.png"))
  missing_exports <- required_exports[!file.exists(required_exports)]
  if (length(missing_exports) > 0L) {
    add_failure("%s missing exports: %s", fig, paste(basename(missing_exports), collapse = ", "))
  }

  if (!file.exists(qc_path)) {
    add_failure("%s missing QC file.", fig)
  } else {
    qc <- read_tsv(qc_path)
    if (!all(qc$exists %in% TRUE)) {
      add_failure("%s QC reports missing export files.", fig)
    }
    png_qc <- qc[qc$file %in% c(paste0(stem, ".png"), "contact_sheet.png"), , drop = FALSE]
    if (nrow(png_qc) == 0L || any(is.na(png_qc$width)) || any(png_qc$width < 2500) || any(png_qc$height < 2000)) {
      add_failure("%s PNG/contact-sheet dimensions are below review threshold or missing.", fig)
    }
  }

  if (!file.exists(manifest_path)) {
    add_failure("%s missing figure_manifest.tsv.", fig)
  } else {
    manifest <- read_tsv(manifest_path)
    if (!"claim_boundary" %in% colnames(manifest) || !is_nonempty(manifest$claim_boundary)) {
      add_failure("%s manifest lacks a non-empty claim_boundary.", fig)
    }
    if (!"primary_claim" %in% colnames(manifest) || !is_nonempty(manifest$primary_claim)) {
      add_failure("%s manifest lacks a non-empty primary_claim.", fig)
    }
  }

  if (!file.exists(panel_manifest_path)) {
    add_failure("%s missing panel_source_data_manifest.tsv.", fig)
  } else {
    panel_manifest <- read_tsv(panel_manifest_path)
    if (!all(c("panel", "source_data", "claim") %in% colnames(panel_manifest))) {
      add_failure("%s panel manifest lacks required columns.", fig)
    } else {
      for (i in seq_len(nrow(panel_manifest))) {
        sources <- trimws(unlist(strsplit(panel_manifest$source_data[i], ";", fixed = TRUE)))
        sources <- sources[nzchar(sources)]
        for (source in sources) {
          source_path <- file.path(fig_dir, source)
          if (!file.exists(source_path)) {
            add_failure("%s panel %s missing source data: %s", fig, panel_manifest$panel[i], source)
          }
        }
      }
      if (!is_nonempty(panel_manifest$claim)) {
        add_failure("%s panel manifest has empty claims.", fig)
      }
    }
  }

  if (!file.exists(legend_path) || file.info(legend_path)$size == 0L) {
    add_failure("%s missing non-empty main_figure_legends.md.", fig)
  }
}

for (fig in required_figures) {
  check_figure_package(fig)
}

ds_store <- list.files(ROOT, pattern = "^\\.DS_Store$", recursive = TRUE, all.files = TRUE, full.names = TRUE)
if (length(ds_store) > 0L) {
  add_failure("Repository contains .DS_Store files: %s", paste(ds_store, collapse = ", "))
}

rbuildignore <- file.path(ROOT, ".Rbuildignore")
if (!file.exists(rbuildignore)) {
  add_failure("Missing .Rbuildignore.")
} else {
  ignored <- readLines(rbuildignore, warn = FALSE)
  for (pattern in c("^data$", "^analysis$", "^figures$", "^cran-comments\\.md$", "^_pkgdown\\.yml$", "^MT3_CAFs_raw\\.txt$", "^Normal_mammary_fibroblasts_raw\\.txt$")) {
    if (!pattern %in% ignored) {
      add_failure(".Rbuildignore missing required pattern: %s", pattern)
    }
  }
}

required_release_files <- file.path(
  ROOT,
  c(
    "NEWS.md",
    "cran-comments.md",
    "inst/CITATION",
    "_pkgdown.yml",
    "docs/package-function-index.md",
    "docs/submission-readiness-isci-commbio.md",
    "docs/manuscript-skeleton-isci-commbio.md",
    "docs/manuscript-draft-isci-commbio.md",
    "docs/reviewer-risk-matrix.md",
    "docs/multi-public-realdata-validation.md",
    "docs/independent-geo-gse246215-validation.md",
    "docs/independent-geo-gse167339-validation.md",
    "docs/gse167339-donor-robustness.md",
    "docs/bioinformatics-validation.md",
    "docs/gse246215-sensitivity.md",
    "docs/extended-method-benchmark.md",
    "docs/marker-stress-benchmark.md",
    "docs/reproducibility-runbook.md",
    "docs/package-runtime-lock.md",
    "inst/doc/fibrodynmix-workflow.md"
  )
)
missing_release_files <- required_release_files[!file.exists(required_release_files)]
if (length(missing_release_files) > 0L) {
  add_failure("Missing package maturity files: %s", paste(missing_release_files, collapse = ", "))
}

for (release_file in required_release_files[file.exists(required_release_files)]) {
  if (file.info(release_file)$size == 0L) {
    add_failure("Package maturity file is empty: %s", release_file)
  }
}

news_path <- file.path(ROOT, "NEWS.md")
if (file.exists(news_path)) {
  news <- readLines(news_path, warn = FALSE)
  if (!any(grepl("Claim Boundary", news, fixed = TRUE))) {
    add_failure("NEWS.md must include a Claim Boundary section.")
  }
}

citation_path <- file.path(ROOT, "inst", "CITATION")
if (file.exists(citation_path)) {
  citation_text <- paste(readLines(citation_path, warn = FALSE), collapse = "\n")
  if (!grepl("bibentry", citation_text, fixed = TRUE)) {
    add_failure("inst/CITATION must contain a bibentry.")
  }
}

required_analysis <- file.path(
  ROOT,
  c(
    "analysis/public_realdata_figure2/nb_fit_diagnostics.tsv",
    "analysis/public_realdata_figure2/condition_cell_counts.tsv",
    "analysis/public_realdata_smoke/nb_fit_diagnostics.tsv",
    "analysis/cross_cohort_transfer/transfer_benchmark.tsv",
    "analysis/public_realdata_transfer/public_transfer_diagnostics.tsv",
    "analysis/multi_public_realdata_validation/multi_public_validation_manifest.tsv",
    "analysis/multi_public_realdata_validation/multi_public_nb_fit_diagnostics.tsv",
    "analysis/multi_public_realdata_validation/multi_public_transfer_diagnostics.tsv",
    "analysis/multi_public_realdata_validation/multi_public_dataset_state_composition.tsv",
    "data/public_geo_gse246215_fibroblast_atlas/gse246215_prepare_manifest.tsv",
    "data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv",
    "data/public_geo_gse167339_human_fibroblast/gse167339_prepare_manifest.tsv",
    "data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv",
    "analysis/independent_geo_gse246215_validation/multi_public_validation_manifest.tsv",
    "analysis/independent_geo_gse246215_validation/multi_public_nb_fit_diagnostics.tsv",
    "analysis/independent_geo_gse246215_validation/multi_public_transfer_diagnostics.tsv",
    "analysis/independent_geo_gse246215_validation/multi_public_marker_coverage.tsv",
    "analysis/independent_geo_gse167339_validation/multi_public_validation_manifest.tsv",
    "analysis/independent_geo_gse167339_validation/multi_public_nb_fit_diagnostics.tsv",
    "analysis/independent_geo_gse167339_validation/multi_public_transfer_diagnostics.tsv",
    "analysis/independent_geo_gse167339_validation/multi_public_transition_summary.tsv",
    "analysis/independent_geo_gse167339_validation/multi_public_dataset_state_composition.tsv",
    "analysis/gse167339_donor_robustness/gse167339_donor_robustness_manifest.tsv",
    "analysis/gse167339_donor_robustness/gse167339_donor_state_composition.tsv",
    "analysis/gse167339_donor_robustness/gse167339_donor_robustness_summary.tsv",
    "analysis/gse167339_donor_robustness/gse167339_leave_donor_out_transfer.tsv",
    "analysis/gse167339_donor_robustness/gse167339_leave_donor_out_state_composition.tsv",
    "analysis/gse167339_donor_robustness/gse167339_hash_threshold_sensitivity.tsv",
    "analysis/bioinformatics_validation/bioinformatics_validation_manifest.tsv",
    "analysis/bioinformatics_validation/fibroblast_purity_qc.tsv",
    "analysis/bioinformatics_validation/state_associated_top_genes.tsv",
    "analysis/bioinformatics_validation/state_program_pathway_enrichment.tsv",
    "analysis/bioinformatics_validation/donor_state_abundance_fpi.tsv",
    "analysis/bioinformatics_validation/donor_aware_state_effects.tsv",
    "analysis/gse246215_sensitivity/gse246215_sensitivity_manifest.tsv",
    "analysis/gse246215_sensitivity/gse246215_sensitivity_run_summary.tsv",
    "analysis/gse246215_sensitivity/gse246215_sensitivity_composition_variability.tsv",
    "analysis/gse246215_sensitivity/gse246215_sensitivity_transfer_diagnostics.tsv",
    "analysis/extended_method_benchmark/extended_method_benchmark_manifest.tsv",
    "analysis/extended_method_benchmark/extended_method_benchmark_metrics.tsv",
    "analysis/extended_method_benchmark/extended_method_benchmark_summary.tsv",
    "analysis/extended_method_benchmark/extended_method_rankings.tsv",
    "analysis/marker_stress_benchmark/marker_stress_benchmark_manifest.tsv",
    "analysis/marker_stress_benchmark/marker_stress_benchmark_metrics.tsv",
    "analysis/marker_stress_benchmark/marker_stress_benchmark_summary.tsv",
    "analysis/marker_stress_benchmark/marker_stress_rankings.tsv",
    "analysis/marker_stress_benchmark/marker_stress_contrast_vs_marker.tsv",
    "analysis/marker_stress_benchmark/marker_stress_prior_audit.tsv",
    "analysis/reproducibility_audit/reproducibility_audit_manifest.tsv",
    "analysis/reproducibility_audit/analysis_catalog.tsv",
    "analysis/reproducibility_audit/script_inventory.tsv",
    "analysis/reproducibility_audit/doc_inventory.tsv",
    "analysis/runtime_lock/runtime_lock_manifest.tsv",
    "analysis/runtime_lock/r_session_info.tsv",
    "analysis/runtime_lock/package_dependency_lock.tsv",
    "analysis/runtime_lock/external_data_lock.tsv",
    "analysis/runtime_lock/heavy_analysis_io.tsv",
    "analysis/vi_posterior/vi_manifest.tsv",
    "analysis/vi_posterior/vi_elbo_trace.tsv",
    "analysis/vi_benchmark/vi_benchmark_manifest.tsv",
    "analysis/vi_benchmark/vi_benchmark_metrics.tsv",
    "analysis/project_maturity/project_maturity_manifest.tsv",
    "analysis/project_maturity/manuscript_evidence_matrix.tsv",
    "analysis/project_maturity/submission_readiness_checklist.tsv"
  )
)
missing_analysis <- required_analysis[!file.exists(required_analysis)]
if (length(missing_analysis) > 0L) {
  add_failure("Missing required analysis outputs: %s", paste(missing_analysis, collapse = ", "))
}

multi_public_manifest_path <- file.path(ROOT, "analysis", "multi_public_realdata_validation", "multi_public_validation_manifest.tsv")
if (file.exists(multi_public_manifest_path)) {
  multi_public_manifest <- read_tsv(multi_public_manifest_path)
  if (!"claim_boundary" %in% colnames(multi_public_manifest) || !is_nonempty(multi_public_manifest$claim_boundary)) {
    add_failure("Multi-public validation manifest lacks a non-empty claim_boundary.")
  }
  if (!"n_public_datasets" %in% colnames(multi_public_manifest) || any(multi_public_manifest$n_public_datasets < 2L)) {
    add_failure("Multi-public validation must report at least two public datasets.")
  }
}

multi_public_transfer_path <- file.path(ROOT, "analysis", "multi_public_realdata_validation", "multi_public_transfer_diagnostics.tsv")
if (file.exists(multi_public_transfer_path)) {
  multi_public_transfer <- read_tsv(multi_public_transfer_path)
  if (!all(c("heldout_dataset_id", "transfer_z_convergence_rate", "transfer_nb_objective") %in% colnames(multi_public_transfer))) {
    add_failure("Multi-public transfer diagnostics lack required transfer columns.")
  }
  if (nrow(multi_public_transfer) < 2L) {
    add_failure("Multi-public transfer diagnostics must include at least two leave-dataset-out rows.")
  }
}

geo_manifest_path <- file.path(ROOT, "analysis", "independent_geo_gse246215_validation", "multi_public_validation_manifest.tsv")
if (file.exists(geo_manifest_path)) {
  geo_manifest <- read_tsv(geo_manifest_path)
  if (!"organisms" %in% colnames(geo_manifest) || !any(grepl("Homo sapiens", geo_manifest$organisms, fixed = TRUE))) {
    add_failure("Independent GSE246215 validation manifest must report Homo sapiens.")
  }
  if (!"n_public_datasets" %in% colnames(geo_manifest) || any(geo_manifest$n_public_datasets < 4L)) {
    add_failure("Independent GSE246215 validation must report at least four cancer-type datasets.")
  }
  if (!"claim_boundary" %in% colnames(geo_manifest) || !is_nonempty(geo_manifest$claim_boundary)) {
    add_failure("Independent GSE246215 validation manifest lacks a non-empty claim_boundary.")
  }
}

geo_transfer_path <- file.path(ROOT, "analysis", "independent_geo_gse246215_validation", "multi_public_transfer_diagnostics.tsv")
if (file.exists(geo_transfer_path)) {
  geo_transfer <- read_tsv(geo_transfer_path)
  if (nrow(geo_transfer) < 4L) {
    add_failure("Independent GSE246215 transfer diagnostics must include four leave-dataset-out rows.")
  }
  if (!all(c("heldout_dataset_id", "transfer_z_convergence_rate", "transfer_nb_objective") %in% colnames(geo_transfer))) {
    add_failure("Independent GSE246215 transfer diagnostics lack required transfer columns.")
  }
}

gse167339_prepare_path <- file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_prepare_manifest.tsv")
if (file.exists(gse167339_prepare_path)) {
  gse167339_prepare <- read_tsv(gse167339_prepare_path)
  if (!"source_accession" %in% colnames(gse167339_prepare) || !"GSE167339" %in% gse167339_prepare$source_accession) {
    add_failure("GSE167339 prepare manifest must report source accession GSE167339.")
  }
  if (!"n_donors" %in% colnames(gse167339_prepare) || any(gse167339_prepare$n_donors < 3L)) {
    add_failure("GSE167339 prepare manifest must report at least three donors after hash demultiplexing.")
  }
  if (!"n_samples" %in% colnames(gse167339_prepare) || any(gse167339_prepare$n_samples < 11L)) {
    add_failure("GSE167339 prepare manifest must report at least eleven sample-level datasets.")
  }
  if (!"n_hash_demux_groups" %in% colnames(gse167339_prepare) || any(gse167339_prepare$n_hash_demux_groups < 4L)) {
    add_failure("GSE167339 prepare manifest must report at least four Human3 hash-demux groups.")
  }
  if (!"claim_boundary" %in% colnames(gse167339_prepare) || !is_nonempty(gse167339_prepare$claim_boundary)) {
    add_failure("GSE167339 prepare manifest lacks a non-empty claim_boundary.")
  }
}

gse167339_manifest_path <- file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_validation_manifest.tsv")
if (file.exists(gse167339_manifest_path)) {
  gse167339_manifest <- read_tsv(gse167339_manifest_path)
  if (!"organisms" %in% colnames(gse167339_manifest) || !any(grepl("Homo sapiens", gse167339_manifest$organisms, fixed = TRUE))) {
    add_failure("Independent GSE167339 validation manifest must report Homo sapiens.")
  }
  if (!"n_public_datasets" %in% colnames(gse167339_manifest) || any(gse167339_manifest$n_public_datasets < 11L)) {
    add_failure("Independent GSE167339 validation must report at least eleven sample-level datasets.")
  }
  if (!"mean_transfer_z_convergence_rate" %in% colnames(gse167339_manifest) || any(gse167339_manifest$mean_transfer_z_convergence_rate < 0.9)) {
    add_failure("Independent GSE167339 validation mean transfer convergence is below 0.9.")
  }
  if (!"claim_boundary" %in% colnames(gse167339_manifest) || !is_nonempty(gse167339_manifest$claim_boundary)) {
    add_failure("Independent GSE167339 validation manifest lacks a non-empty claim_boundary.")
  }
}

gse167339_transfer_path <- file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_transfer_diagnostics.tsv")
if (file.exists(gse167339_transfer_path)) {
  gse167339_transfer <- read_tsv(gse167339_transfer_path)
  if (nrow(gse167339_transfer) < 11L) {
    add_failure("Independent GSE167339 transfer diagnostics must include eleven leave-dataset-out rows.")
  }
  if (!all(c("heldout_dataset_id", "transfer_z_convergence_rate", "transfer_nb_objective") %in% colnames(gse167339_transfer))) {
    add_failure("Independent GSE167339 transfer diagnostics lack required transfer columns.")
  }
  if ("transfer_z_convergence_rate" %in% colnames(gse167339_transfer) && any(gse167339_transfer$transfer_z_convergence_rate < 0.9)) {
    add_failure("Independent GSE167339 transfer diagnostics contain convergence rates below 0.9.")
  }
}

gse167339_transition_path <- file.path(ROOT, "analysis", "independent_geo_gse167339_validation", "multi_public_transition_summary.tsv")
if (file.exists(gse167339_transition_path)) {
  gse167339_transition <- read_tsv(gse167339_transition_path)
  if (!"status" %in% colnames(gse167339_transition) || !"ok" %in% gse167339_transition$status) {
    add_failure("Independent GSE167339 transition summary must report status ok.")
  }
  if (!"converged" %in% colnames(gse167339_transition) || !all(gse167339_transition$converged %in% TRUE)) {
    add_failure("Independent GSE167339 transition summary must report convergence.")
  }
}

gse167339_donor_manifest_path <- file.path(ROOT, "analysis", "gse167339_donor_robustness", "gse167339_donor_robustness_manifest.tsv")
if (file.exists(gse167339_donor_manifest_path)) {
  gse167339_donor_manifest <- read_tsv(gse167339_donor_manifest_path)
  if (!"source_accession" %in% colnames(gse167339_donor_manifest) || !"GSE167339" %in% gse167339_donor_manifest$source_accession) {
    add_failure("GSE167339 donor robustness manifest must report source accession GSE167339.")
  }
  if (!"n_donors" %in% colnames(gse167339_donor_manifest) || any(gse167339_donor_manifest$n_donors < 3L)) {
    add_failure("GSE167339 donor robustness manifest must report at least three donors.")
  }
  if (!"n_leave_donor_out_runs" %in% colnames(gse167339_donor_manifest) || any(gse167339_donor_manifest$n_leave_donor_out_runs < 3L)) {
    add_failure("GSE167339 donor robustness manifest must report three leave-donor-out runs.")
  }
  if (!"min_leave_donor_transfer_convergence" %in% colnames(gse167339_donor_manifest) || any(gse167339_donor_manifest$min_leave_donor_transfer_convergence < 0.9)) {
    add_failure("GSE167339 donor robustness min leave-donor transfer convergence is below 0.9.")
  }
  if (!"claim_boundary" %in% colnames(gse167339_donor_manifest) || !is_nonempty(gse167339_donor_manifest$claim_boundary)) {
    add_failure("GSE167339 donor robustness manifest lacks a non-empty claim_boundary.")
  }
}

gse167339_leave_donor_path <- file.path(ROOT, "analysis", "gse167339_donor_robustness", "gse167339_leave_donor_out_transfer.tsv")
if (file.exists(gse167339_leave_donor_path)) {
  gse167339_leave_donor <- read_tsv(gse167339_leave_donor_path)
  if (nrow(gse167339_leave_donor) < 3L) {
    add_failure("GSE167339 leave-donor-out transfer must include three donor rows.")
  }
  if (!all(c("heldout_donor_id", "transfer_z_convergence_rate", "transfer_nb_objective") %in% colnames(gse167339_leave_donor))) {
    add_failure("GSE167339 leave-donor-out transfer lacks required columns.")
  }
  if ("transfer_z_convergence_rate" %in% colnames(gse167339_leave_donor) && any(gse167339_leave_donor$transfer_z_convergence_rate < 0.9)) {
    add_failure("GSE167339 leave-donor-out transfer contains convergence rates below 0.9.")
  }
}

gse167339_hash_sensitivity_path <- file.path(ROOT, "analysis", "gse167339_donor_robustness", "gse167339_hash_threshold_sensitivity.tsv")
if (file.exists(gse167339_hash_sensitivity_path)) {
  gse167339_hash_sensitivity <- read_tsv(gse167339_hash_sensitivity_path)
  if (!all(c("hash_min_count", "hash_min_ratio", "n_assigned_cells", "n_hash_groups_ge_40") %in% colnames(gse167339_hash_sensitivity))) {
    add_failure("GSE167339 hash-threshold sensitivity lacks required columns.")
  }
  if (nrow(gse167339_hash_sensitivity) < 12L) {
    add_failure("GSE167339 hash-threshold sensitivity must include at least twelve threshold settings.")
  }
  if ("n_hash_groups_ge_40" %in% colnames(gse167339_hash_sensitivity) && any(gse167339_hash_sensitivity$n_hash_groups_ge_40 < 4L)) {
    add_failure("GSE167339 hash-threshold sensitivity must retain four hash groups across settings.")
  }
}

bioinfo_manifest_path <- file.path(ROOT, "analysis", "bioinformatics_validation", "bioinformatics_validation_manifest.tsv")
if (file.exists(bioinfo_manifest_path)) {
  bioinfo_manifest <- read_tsv(bioinfo_manifest_path)
  if (!"claim_boundary" %in% colnames(bioinfo_manifest) || !is_nonempty(bioinfo_manifest$claim_boundary)) {
    add_failure("Bioinformatics validation manifest lacks a non-empty claim_boundary.")
  }
  if (!"n_validation_datasets" %in% colnames(bioinfo_manifest) || any(bioinfo_manifest$n_validation_datasets < 2L)) {
    add_failure("Bioinformatics validation must include at least two validation datasets.")
  }
  if (!"min_purity_margin" %in% colnames(bioinfo_manifest) || any(bioinfo_manifest$min_purity_margin <= 0)) {
    add_failure("Bioinformatics validation must report positive minimum fibroblast purity margin.")
  }
  if (!"n_enriched_state_pathways_q05" %in% colnames(bioinfo_manifest) || any(bioinfo_manifest$n_enriched_state_pathways_q05 < 2L)) {
    add_failure("Bioinformatics validation must report at least two enriched state/pathway pairs at q <= 0.05.")
  }
  if (!"n_donor_effect_rows" %in% colnames(bioinfo_manifest) || any(bioinfo_manifest$n_donor_effect_rows < 12L)) {
    add_failure("Bioinformatics validation must report donor-aware state effects for GSE167339.")
  }
}

bioinfo_purity_path <- file.path(ROOT, "analysis", "bioinformatics_validation", "fibroblast_purity_qc.tsv")
if (file.exists(bioinfo_purity_path)) {
  bioinfo_purity <- read_tsv(bioinfo_purity_path)
  if (!all(c("validation_dataset", "dataset_id", "fibroblast_positive_score", "purity_margin_mean", "low_purity_fraction") %in% colnames(bioinfo_purity))) {
    add_failure("Fibroblast purity QC lacks required columns.")
  }
  if ("purity_margin_mean" %in% colnames(bioinfo_purity) && any(bioinfo_purity$purity_margin_mean <= 0)) {
    add_failure("Fibroblast purity QC contains non-positive purity margins.")
  }
}

bioinfo_enrichment_path <- file.path(ROOT, "analysis", "bioinformatics_validation", "state_program_pathway_enrichment.tsv")
if (file.exists(bioinfo_enrichment_path)) {
  bioinfo_enrichment <- read_tsv(bioinfo_enrichment_path)
  if (!all(c("validation_dataset", "state", "pathway", "n_overlap", "q_value") %in% colnames(bioinfo_enrichment))) {
    add_failure("State program pathway enrichment lacks required columns.")
  }
  if ("q_value" %in% colnames(bioinfo_enrichment) && sum(bioinfo_enrichment$q_value <= 0.05, na.rm = TRUE) < 2L) {
    add_failure("State program pathway enrichment must include at least two q <= 0.05 rows.")
  }
}

bioinfo_effects_path <- file.path(ROOT, "analysis", "bioinformatics_validation", "donor_aware_state_effects.tsv")
if (file.exists(bioinfo_effects_path)) {
  bioinfo_effects <- read_tsv(bioinfo_effects_path)
  if (!all(c("validation_dataset", "donor_id", "state", "delta_composition", "delta_fpi") %in% colnames(bioinfo_effects))) {
    add_failure("Donor-aware state effects lack required columns.")
  }
  if (nrow(bioinfo_effects) < 12L) {
    add_failure("Donor-aware state effects must include at least 12 rows.")
  }
}

extended_benchmark_manifest_path <- file.path(ROOT, "analysis", "extended_method_benchmark", "extended_method_benchmark_manifest.tsv")
if (file.exists(extended_benchmark_manifest_path)) {
  extended_manifest <- read_tsv(extended_benchmark_manifest_path)
  if (!"topic_backend" %in% colnames(extended_manifest) || !is_nonempty(extended_manifest$topic_backend)) {
    add_failure("Extended method benchmark manifest must report topic_backend.")
  }
  if (!"claim_boundary" %in% colnames(extended_manifest) || !is_nonempty(extended_manifest$claim_boundary)) {
    add_failure("Extended method benchmark manifest lacks a non-empty claim_boundary.")
  }
}

extended_benchmark_metrics_path <- file.path(ROOT, "analysis", "extended_method_benchmark", "extended_method_benchmark_metrics.tsv")
if (file.exists(extended_benchmark_metrics_path)) {
  extended_metrics <- read_tsv(extended_benchmark_metrics_path)
  if (!"method" %in% colnames(extended_metrics) || !"topic_nmf" %in% extended_metrics$method) {
    add_failure("Extended method benchmark metrics must include the topic_nmf baseline.")
  }
}

marker_stress_manifest_path <- file.path(ROOT, "analysis", "marker_stress_benchmark", "marker_stress_benchmark_manifest.tsv")
if (file.exists(marker_stress_manifest_path)) {
  marker_stress_manifest <- read_tsv(marker_stress_manifest_path)
  if (!"claim_boundary" %in% colnames(marker_stress_manifest) || !is_nonempty(marker_stress_manifest$claim_boundary)) {
    add_failure("Marker stress benchmark manifest lacks a non-empty claim_boundary.")
  }
  if (!"n_stress_modes" %in% colnames(marker_stress_manifest) || any(marker_stress_manifest$n_stress_modes < 7L)) {
    add_failure("Marker stress benchmark must report at least seven stress modes.")
  }
  if (!"n_methods" %in% colnames(marker_stress_manifest) || any(marker_stress_manifest$n_methods < 5L)) {
    add_failure("Marker stress benchmark must report at least five methods.")
  }
  if (!"n_modes_where_any_fibrodynmix_beats_marker_rmse" %in% colnames(marker_stress_manifest) || any(marker_stress_manifest$n_modes_where_any_fibrodynmix_beats_marker_rmse < 2L)) {
    add_failure("Marker stress benchmark must show at least two stress modes where a FibroDynMix variant beats marker scoring by RMSE.")
  }
}

marker_stress_metrics_path <- file.path(ROOT, "analysis", "marker_stress_benchmark", "marker_stress_benchmark_metrics.tsv")
if (file.exists(marker_stress_metrics_path)) {
  marker_stress_metrics <- read_tsv(marker_stress_metrics_path)
  required_stress_methods <- c("marker_scoring", "topic_nmf", "fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_vi")
  if (!"method" %in% colnames(marker_stress_metrics) || length(setdiff(required_stress_methods, marker_stress_metrics$method)) > 0L) {
    add_failure("Marker stress benchmark metrics must include marker_scoring, topic_nmf, initializer, NB, and VI methods.")
  }
  required_stress_modes <- c("clean_prior", "missing_markers", "corrupted_prior", "shared_markers", "swapped_markers", "hidden_program_corrupted_prior", "hidden_program_missing_markers")
  if (!"stress_mode" %in% colnames(marker_stress_metrics) || length(setdiff(required_stress_modes, marker_stress_metrics$stress_mode)) > 0L) {
    add_failure("Marker stress benchmark metrics are missing required stress modes.")
  }
}

marker_stress_contrast_path <- file.path(ROOT, "analysis", "marker_stress_benchmark", "marker_stress_contrast_vs_marker.tsv")
if (file.exists(marker_stress_contrast_path)) {
  marker_stress_contrast <- read_tsv(marker_stress_contrast_path)
  if (!all(c("method", "rmse_delta_vs_marker", "dominant_accuracy_delta_vs_marker") %in% colnames(marker_stress_contrast))) {
    add_failure("Marker stress contrast table lacks required delta columns.")
  } else {
    fibro_contrast <- marker_stress_contrast[grepl("^fibrodynmix", marker_stress_contrast$method), , drop = FALSE]
    if (!any(fibro_contrast$rmse_delta_vs_marker < 0, na.rm = TRUE)) {
      add_failure("Marker stress contrast must include at least one FibroDynMix RMSE improvement over marker scoring.")
    }
    if (!any(fibro_contrast$dominant_accuracy_delta_vs_marker > 0, na.rm = TRUE)) {
      add_failure("Marker stress contrast must include at least one FibroDynMix dominant-accuracy improvement over marker scoring.")
    }
  }
}

gse_sensitivity_manifest_path <- file.path(ROOT, "analysis", "gse246215_sensitivity", "gse246215_sensitivity_manifest.tsv")
if (file.exists(gse_sensitivity_manifest_path)) {
  gse_sensitivity_manifest <- read_tsv(gse_sensitivity_manifest_path)
  if (!"claim_boundary" %in% colnames(gse_sensitivity_manifest) || !is_nonempty(gse_sensitivity_manifest$claim_boundary)) {
    add_failure("GSE246215 sensitivity manifest lacks a non-empty claim_boundary.")
  }
  if (!"n_runs" %in% colnames(gse_sensitivity_manifest) || any(gse_sensitivity_manifest$n_runs < 4L)) {
    add_failure("GSE246215 sensitivity analysis must report at least four validation runs.")
  }
}

gse_sensitivity_transfer_path <- file.path(ROOT, "analysis", "gse246215_sensitivity", "gse246215_sensitivity_transfer_diagnostics.tsv")
if (file.exists(gse_sensitivity_transfer_path)) {
  gse_sensitivity_transfer <- read_tsv(gse_sensitivity_transfer_path)
  if (!all(c("sensitivity_mode", "seed", "transfer_z_convergence_rate") %in% colnames(gse_sensitivity_transfer))) {
    add_failure("GSE246215 sensitivity transfer diagnostics lack required sensitivity columns.")
  }
}

repro_manifest_path <- file.path(ROOT, "analysis", "reproducibility_audit", "reproducibility_audit_manifest.tsv")
if (file.exists(repro_manifest_path)) {
  repro_manifest <- read_tsv(repro_manifest_path)
  if (!"all_registered_outputs_present" %in% colnames(repro_manifest) || !all(repro_manifest$all_registered_outputs_present %in% TRUE)) {
    add_failure("Reproducibility audit must report all registered outputs present.")
  }
  if (!"claim_boundary" %in% colnames(repro_manifest) || !is_nonempty(repro_manifest$claim_boundary)) {
    add_failure("Reproducibility audit manifest lacks a non-empty claim_boundary.")
  }
}

analysis_catalog_path <- file.path(ROOT, "analysis", "reproducibility_audit", "analysis_catalog.tsv")
if (file.exists(analysis_catalog_path)) {
  analysis_catalog <- read_tsv(analysis_catalog_path)
  if (!all(c("analysis_id", "run_script", "primary_manifest_exists", "claim_boundary") %in% colnames(analysis_catalog))) {
    add_failure("Reproducibility analysis catalog lacks required columns.")
  }
  if (!all(analysis_catalog$primary_manifest_exists %in% TRUE)) {
    add_failure("Reproducibility analysis catalog reports missing primary manifests.")
  }
}

runtime_manifest_path <- file.path(ROOT, "analysis", "runtime_lock", "runtime_lock_manifest.tsv")
if (file.exists(runtime_manifest_path)) {
  runtime_manifest <- read_tsv(runtime_manifest_path)
  if (!"claim_boundary" %in% colnames(runtime_manifest) || !is_nonempty(runtime_manifest$claim_boundary)) {
    add_failure("Runtime lock manifest lacks a non-empty claim_boundary.")
  }
  if (!"n_missing_suggests" %in% colnames(runtime_manifest) || any(runtime_manifest$n_missing_suggests > 0L)) {
    add_failure("Runtime lock reports missing Suggests packages.")
  }
  if (!"all_external_files_present" %in% colnames(runtime_manifest) || !all(runtime_manifest$all_external_files_present %in% TRUE)) {
    add_failure("Runtime lock reports missing external data files.")
  }
}

package_lock_path <- file.path(ROOT, "analysis", "runtime_lock", "package_dependency_lock.tsv")
if (file.exists(package_lock_path)) {
  package_lock <- read_tsv(package_lock_path)
  if (!all(c("package", "dependency_field", "installed", "version") %in% colnames(package_lock))) {
    add_failure("Package dependency lock lacks required columns.")
  } else {
    required_suggests <- c("data.table", "GEOquery", "ggplot2", "irlba", "lsa", "Matrix", "NMF", "scater", "scran", "scuttle", "SeuratObject", "SingleCellExperiment", "testthat")
    missing_suggest_rows <- setdiff(required_suggests, package_lock$package)
    if (length(missing_suggest_rows) > 0L) {
      add_failure("Package dependency lock is missing Suggests rows: %s", paste(missing_suggest_rows, collapse = ", "))
    }
    suggest_rows <- package_lock[package_lock$package %in% required_suggests, , drop = FALSE]
    if (nrow(suggest_rows) > 0L && !all(suggest_rows$installed %in% TRUE)) {
      add_failure("Package dependency lock reports uninstalled Suggests packages.")
    }
    if (nrow(suggest_rows) > 0L && !is_nonempty(suggest_rows$version)) {
      add_failure("Package dependency lock has empty Suggests versions.")
    }
  }
}

external_lock_path <- file.path(ROOT, "analysis", "runtime_lock", "external_data_lock.tsv")
if (file.exists(external_lock_path)) {
  external_lock <- read_tsv(external_lock_path)
  if (!all(c("source_id", "accession_or_doi", "public_record", "local_files_present", "claim_boundary") %in% colnames(external_lock))) {
    add_failure("External data lock lacks required columns.")
  } else {
    if (!"10.6071/M3238R" %in% external_lock$accession_or_doi) {
      add_failure("External data lock must include Dryad DOI 10.6071/M3238R.")
    }
    if (!"GSE246215" %in% external_lock$accession_or_doi) {
      add_failure("External data lock must include GEO accession GSE246215.")
    }
    if (!"GSE167339" %in% external_lock$accession_or_doi) {
      add_failure("External data lock must include GEO accession GSE167339.")
    }
    if (!all(external_lock$local_files_present %in% TRUE)) {
      add_failure("External data lock reports missing local files.")
    }
    if (!is_nonempty(external_lock$claim_boundary)) {
      add_failure("External data lock has empty claim boundaries.")
    }
  }
}

heavy_io_path <- file.path(ROOT, "analysis", "runtime_lock", "heavy_analysis_io.tsv")
if (file.exists(heavy_io_path)) {
  heavy_io <- read_tsv(heavy_io_path)
  if (!all(c("analysis_id", "expected_inputs", "expected_outputs", "runtime_class", "rerun_trigger") %in% colnames(heavy_io))) {
    add_failure("Heavy analysis I/O lock lacks required columns.")
  } else {
    if (nrow(heavy_io) < 13L) {
      add_failure("Heavy analysis I/O lock must cover at least 13 analyses.")
    }
    if (!is_nonempty(heavy_io$expected_inputs) || !is_nonempty(heavy_io$expected_outputs)) {
      add_failure("Heavy analysis I/O lock has empty expected input/output fields.")
    }
  }
}

if (length(failures) > 0L) {
  cat("Project integrity check failed:\n")
  cat(paste0("- ", failures, "\n"), sep = "")
  quit(status = 1)
}

cat("Project integrity check passed.\n")
cat(sprintf("Checked %d figure packages and required analysis outputs.\n", length(required_figures)))
if (length(notes) > 0L) {
  cat(paste0("- ", notes, "\n"), sep = "")
}
