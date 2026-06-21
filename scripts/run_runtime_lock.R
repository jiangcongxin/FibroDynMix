#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_runtime_lock.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "runtime_lock")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

file_exists_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0
}

collapse_paths <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  paste(paths, collapse = ";")
}

description <- read.dcf(file.path(ROOT, "DESCRIPTION"))[1, ]
suggests_raw <- unname(description[["Suggests"]])
suggest_pkgs <- trimws(unlist(strsplit(gsub("\n", " ", suggests_raw), ",")))
suggest_pkgs <- sub("\\s*\\(.*\\)$", "", suggest_pkgs)
suggest_pkgs <- suggest_pkgs[nzchar(suggest_pkgs)]
runtime_pkgs <- unique(c("R", suggest_pkgs))

package_lock <- data.frame(
  package = runtime_pkgs,
  dependency_field = c("runtime", rep("Suggests", length(suggest_pkgs))),
  installed = c(TRUE, vapply(suggest_pkgs, requireNamespace, logical(1), quietly = TRUE)),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(suggest_pkgs, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    }, character(1))
  ),
  source = c("base R runtime", rep("DESCRIPTION", length(suggest_pkgs))),
  required_for = c(
    "package tests, analysis scripts, and R CMD check",
    rep("optional benchmark, public data, single-cell container, and test workflows", length(suggest_pkgs))
  ),
  stringsAsFactors = FALSE
)
write_tsv(package_lock, file.path(OUT, "package_dependency_lock.tsv"))

runtime_session <- data.frame(
  key = c(
    "r_version_string",
    "platform",
    "os_type",
    "blas",
    "lapack",
    "locale",
    "working_directory",
    "package_version",
    "generated_at"
  ),
  value = c(
    R.version.string,
    R.version$platform,
    .Platform$OS.type,
    sessionInfo()$BLAS,
    sessionInfo()$LAPACK,
    paste(Sys.getlocale(), collapse = ";"),
    ROOT,
    unname(description[["Version"]]),
    as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)
write_tsv(runtime_session, file.path(OUT, "r_session_info.tsv"))

external_rows <- list(
  data.frame(
    source_id = "Dryad DOI 10.6071/M3238R",
    source_type = "DOI",
    accession_or_doi = "10.6071/M3238R",
    public_record = "https://datadryad.org/dataset/doi:10.6071/M3238R",
    download_url = "https://zenodo.org/records/3977255",
    local_files = collapse_paths(file.path(ROOT, c(
      "MT3_CAFs_raw.txt",
      "Normal_mammary_fibroblasts_raw.txt",
      "data/public_dryad_breast_fibroblast/MT3_CAFs_raw.txt",
      "data/public_dryad_breast_fibroblast/Normal_mammary_fibroblasts_raw.txt"
    ))),
    local_files_present = all(file.exists(file.path(ROOT, c("MT3_CAFs_raw.txt", "Normal_mammary_fibroblasts_raw.txt")))),
    role = "mouse breast CAF and normal mammary fibroblast raw count matrices for public smoke, transfer, and multi-public validation",
    claim_boundary = "Two condition-specific count matrices from one Dryad record; not an independent multi-study human atlas.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    source_id = "GEO GSE246215",
    source_type = "GEO",
    accession_or_doi = "GSE246215",
    public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246215",
    download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE246215&format=file",
    local_files = collapse_paths(file.path(ROOT, c(
      "data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_counts.csv.gz",
      "data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_metadata.csv.gz",
      "data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv"
    ))),
    local_files_present = all(file.exists(file.path(ROOT, c(
      "data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_counts.csv.gz",
      "data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_metadata.csv.gz",
      "data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv"
    )))),
    role = "independent human fibroblast count-matrix validation and GSE246215 sensitivity/QC analysis",
    claim_boundary = "Uses processed GEO supplementary count and metadata files; no FASTQ-level reprocessing.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    source_id = "GEO GSE167339",
    source_type = "GEO",
    accession_or_doi = "GSE167339",
    public_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167339",
    download_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167339&format=file",
    local_files = collapse_paths(file.path(ROOT, c(
      "data/public_geo_gse167339_human_fibroblast/GSE167339_RAW.tar",
      "data/public_geo_gse167339_human_fibroblast/gse167339_prepare_manifest.tsv",
      "data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv"
    ))),
    local_files_present = all(file.exists(file.path(ROOT, c(
      "data/public_geo_gse167339_human_fibroblast/GSE167339_RAW.tar",
      "data/public_geo_gse167339_human_fibroblast/gse167339_prepare_manifest.tsv",
      "data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv"
    )))),
    role = "independent human donor-level fibroblast perturbation validation from processed 10x-style MTX files",
    claim_boundary = "Uses bounded sampled GEO supplementary MTX matrices from three donors; Human 3 hash groups are included as hash_unknown because treatment-label mapping is not present in the public workbook.",
    stringsAsFactors = FALSE
  )
)

dryad_registry <- file.path(ROOT, "analysis", "multi_public_realdata_validation", "multi_public_dataset_registry.tsv")
if (file_exists_nonempty(dryad_registry)) {
  registry <- read_tsv(dryad_registry)
  if (all(c("counts_path", "public_record", "download_url") %in% colnames(registry))) {
    external_rows[[1]]$local_files <- collapse_paths(c(external_rows[[1]]$local_files, registry$counts_path))
    external_rows[[1]]$local_files_present <- all(file.exists(registry$counts_path))
  }
}

geo_registry <- file.path(ROOT, "data", "public_geo_gse246215_fibroblast_atlas", "gse246215_fibroblast_dataset_manifest.tsv")
if (file_exists_nonempty(geo_registry)) {
  registry <- read_tsv(geo_registry)
  if (all(c("counts_path", "public_record", "download_url") %in% colnames(registry))) {
    external_rows[[2]]$local_files <- collapse_paths(c(external_rows[[2]]$local_files, registry$counts_path))
    external_rows[[2]]$local_files_present <- all(file.exists(registry$counts_path))
  }
}
gse167339_registry <- file.path(ROOT, "data", "public_geo_gse167339_human_fibroblast", "gse167339_human_fibroblast_dataset_manifest.tsv")
if (file_exists_nonempty(gse167339_registry)) {
  registry <- read_tsv(gse167339_registry)
  if (all(c("counts_path", "public_record", "download_url") %in% colnames(registry))) {
    external_rows[[3]]$local_files <- collapse_paths(c(external_rows[[3]]$local_files, registry$counts_path))
    external_rows[[3]]$local_files_present <- all(file.exists(registry$counts_path))
  }
}
external_data_lock <- do.call(rbind, external_rows)
write_tsv(external_data_lock, file.path(OUT, "external_data_lock.tsv"))

analysis_catalog_path <- file.path(ROOT, "analysis", "reproducibility_audit", "analysis_catalog.tsv")
if (file_exists_nonempty(analysis_catalog_path)) {
  analysis_catalog <- read_tsv(analysis_catalog_path)
} else {
  analysis_catalog <- data.frame(
    analysis_id = character(),
    run_script = character(),
    output_dir = character(),
    primary_manifest = character(),
    research_role = character(),
    claim_boundary = character(),
    stringsAsFactors = FALSE
  )
}

io_notes <- data.frame(
  analysis_id = c(
    "public_realdata_smoke",
    "public_realdata_transfer",
    "multi_public_realdata_validation",
    "independent_geo_gse246215_validation",
    "independent_geo_gse167339_validation",
    "gse167339_donor_robustness",
    "bioinformatics_validation",
    "gse246215_sensitivity",
    "extended_method_benchmark",
    "marker_stress_benchmark",
    "vi_benchmark",
    "vi_posterior",
    "cross_cohort_transfer",
    "study_effect_sensitivity",
    "bootstrap_uncertainty",
    "transition_flow",
    "project_maturity"
  ),
  expected_inputs = c(
    "MT3_CAFs_raw.txt; Normal_mammary_fibroblasts_raw.txt; marker prior list in script",
    "MT3_CAFs_raw.txt; Normal_mammary_fibroblasts_raw.txt; fitted source/target subsets",
    "data/public_dryad_breast_fibroblast dataset registry; public raw count matrices",
    "data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv; sampled RDS count matrices",
    "data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv; sampled RDS count matrices prepared from 10x-style MTX files",
    "GSE167339 sample-level validation outputs; GSE167339 dataset manifest; Human3 hash-demux count matrices",
    "GSE246215 and GSE167339 count matrices; FibroDynMix cell state weights; curated marker and pathway sets",
    "GSE246215 processed counts and metadata; sampling seeds; library-size trimming parameters",
    "simulated FibroDynMix count data across discrete, continuous, batch-confounded, and rare-transition scenarios",
    "simulated batch-confounded count data with corrupted, missing, shared, swapped, and hidden-program marker-prior stress modes",
    "simulated truth z and calibrated VI interval scale grid",
    "simulated count matrix and NB mode fit",
    "simulated leave-study-out count matrices and metadata",
    "simulated batch-confounded count matrices with study identifiers",
    "simulated count matrix and bootstrap seeds",
    "state composition summaries and state marker/program estimates",
    "DESCRIPTION; NAMESPACE; man; tests; scripts; docs; figure manifests"
  ),
  expected_outputs = c(
    "nb_fit_diagnostics.tsv; state_composition.tsv; run_manifest.tsv",
    "public_transfer_diagnostics.tsv; public_transfer_manifest.tsv; source/target state composition tables",
    "multi_public_validation_manifest.tsv; NB diagnostics; transfer diagnostics; transition flow tables",
    "multi_public_validation_manifest.tsv; NB diagnostics; transfer diagnostics; marker coverage",
    "multi_public_validation_manifest.tsv; NB diagnostics; transfer diagnostics; transition summary; donor/sample composition",
    "donor robustness manifest; donor state composition; leave-donor-out transfer; hash-threshold sensitivity",
    "bioinformatics validation manifest; purity QC; state top genes; pathway enrichment; donor abundance/FPI effects",
    "sensitivity manifest; composition variability; transfer diagnostics; run summaries",
    "extended benchmark manifest; method metrics; summary; rankings",
    "marker stress benchmark manifest; stress metrics; prior audit; rankings; contrast versus marker scoring",
    "VI benchmark manifest and interval calibration metrics",
    "VI manifest; ELBO trace; posterior z summaries",
    "transfer benchmark metrics and manifest",
    "study-effect sensitivity summary tables",
    "sample composition uncertainty and bootstrap intervals",
    "transition flow summary and FPI tables",
    "project maturity manifest; evidence matrix; submission readiness checklist"
  ),
  runtime_class = c(
    "medium real-data smoke",
    "medium real-data transfer",
    "heavy public count-matrix validation",
    "heavy independent human public validation",
    "heavy independent human donor perturbation validation",
    "heavy donor-level robustness analysis",
    "medium bioinformatics validation analysis",
    "heavy repeated public sensitivity analysis",
    "heavy simulation benchmark",
    "heavy marker-prior stress simulation benchmark",
    "medium simulation calibration",
    "medium posterior demonstration",
    "medium simulation transfer benchmark",
    "medium simulation sensitivity",
    "medium bootstrap analysis",
    "medium transition-flow analysis",
    "lightweight audit"
  ),
  rerun_trigger = c(
    "NB optimizer, public data parsing, or marker prior changes",
    "transfer optimizer or public data parsing changes",
    "public registry, NB optimizer, transfer optimizer, or transition-flow changes",
    "GSE246215 preparation, NB optimizer, transfer optimizer, or marker coverage changes",
    "GSE167339 preparation, hash demultiplexing, MTX parsing, NB optimizer, transfer optimizer, or transition-flow changes",
    "GSE167339 validation outputs, donor metadata, hash demultiplex thresholds, or transfer optimizer changes",
    "public count matrices, state weights, curated marker/pathway sets, or FPI definitions change",
    "GSE246215 preparation, QC filtering, transfer optimizer, or real-data claim changes",
    "benchmark method set, NMF/topic backend, optimizer, or simulation changes",
    "marker prior perturbation design, hidden-program injection, benchmark method set, or optimizer changes",
    "VI posterior, interval calibration, or simulation changes",
    "VI posterior implementation or NB mode changes",
    "transfer implementation or simulation changes",
    "study-effect penalty or optimizer changes",
    "bootstrap uncertainty or composition summary changes",
    "transition cost, OT flow, or FPI changes",
    "package files, docs, figures, or analysis evidence changes"
  ),
  stringsAsFactors = FALSE
)

heavy_analysis_io <- merge(
  analysis_catalog,
  io_notes,
  by = "analysis_id",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(heavy_analysis_io) == 0L) {
  heavy_analysis_io <- io_notes
}
write_tsv(heavy_analysis_io, file.path(OUT, "heavy_analysis_io.tsv"))

manifest <- data.frame(
  analysis = "runtime_lock",
  primary_claim = "FibroDynMix has a package-level runtime, dependency, public-data, and heavy-analysis I/O lock suitable for handoff and review.",
  claim_boundary = "This is a package-level runtime snapshot and analysis contract, not a full renv.lock or container image.",
  r_version = R.version.string,
  n_locked_packages = nrow(package_lock),
  n_missing_suggests = sum(!package_lock$installed),
  n_external_sources = nrow(external_data_lock),
  all_external_files_present = all(external_data_lock$local_files_present),
  n_heavy_analysis_contracts = nrow(heavy_analysis_io),
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "runtime_lock_manifest.tsv"))

message("Runtime lock written to: ", normalizePath(OUT))
