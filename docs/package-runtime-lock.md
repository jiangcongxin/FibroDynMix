# Package Runtime Lock

This document records the package-level runtime contract for FibroDynMix. It is
intended for manuscript review, analyst handoff, and release-readiness checks.
The generated tables live under `analysis/runtime_lock/`.

## Scope

The runtime lock captures:

- the R runtime used for the current package and analysis snapshot;
- installed versions of key `Suggests` dependencies;
- public data sources, DOI/GEO accessions, download locations, and local files;
- expected inputs and outputs for heavy analyses that are intentionally kept
  outside unit tests.

This is not a full `renv.lock`, Docker image, or Bioconductor release snapshot.
It is a package-level runtime record and analysis I/O contract. A future release
can add `renv` or a container once the public-data workflow is frozen.

## Generated Files

Regenerate the lock with:

```bash
Rscript scripts/run_runtime_lock.R
```

The script writes:

- `analysis/runtime_lock/runtime_lock_manifest.tsv`
- `analysis/runtime_lock/r_session_info.tsv`
- `analysis/runtime_lock/package_dependency_lock.tsv`
- `analysis/runtime_lock/external_data_lock.tsv`
- `analysis/runtime_lock/heavy_analysis_io.tsv`

## R Runtime

The current validated snapshot was generated under:

```text
R version 4.6.0 (2026-04-24)
```

Exact platform, BLAS/LAPACK, locale, package version, and generation time are
stored in `analysis/runtime_lock/r_session_info.tsv`.

## Key Suggests Packages

The runtime lock records the installed version and availability of each package
listed in `DESCRIPTION` `Suggests`:

- `data.table`
- `GEOquery`
- `ggplot2`
- `irlba`
- `lsa`
- `Matrix`
- `NMF`
- `scater`
- `scran`
- `scuttle`
- `SeuratObject`
- `SingleCellExperiment`
- `testthat`

These packages are optional at package installation time but important for the
full research workflow. In particular, `NMF` supports the mature topic baseline,
`GEOquery` supports GEO metadata access, `ggplot2` supports package-level
visualization helpers, and the Bioconductor single-cell packages support public
count-matrix handling.

## External Public Data

The real-data layer currently uses two public sources:

1. Dryad DOI `10.6071/M3238R`

   Public record: `https://datadryad.org/dataset/doi:10.6071/M3238R`

   Local raw count files:

   - `MT3_CAFs_raw.txt`
   - `Normal_mammary_fibroblasts_raw.txt`
   - `data/public_dryad_breast_fibroblast/MT3_CAFs_raw.txt`
   - `data/public_dryad_breast_fibroblast/Normal_mammary_fibroblasts_raw.txt`

   Claim boundary: two condition-specific mouse count matrices from one Dryad
   record; this supports public raw-count execution and transfer smoke tests,
   not a broad human atlas claim.

2. GEO accession `GSE246215`

   Public record:
   `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246215`

   Local processed files:

   - `data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_counts.csv.gz`
   - `data/public_geo_gse246215_fibroblast_atlas/GSE246215_Fibroblast_metadata.csv.gz`
   - sampled cancer-type RDS count matrices listed in
     `data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv`

   Claim boundary: processed GEO supplementary count and metadata files are
   used; raw FASTQ-level reprocessing is outside the current package workflow.

3. GEO accession `GSE167339`

   Public record:
   `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167339`

   Local processed files:

   - `data/public_geo_gse167339_human_fibroblast/GSE167339_RAW.tar`
   - `data/public_geo_gse167339_human_fibroblast/gse167339_prepare_manifest.tsv`
   - `data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv`
   - sampled sample-level RDS count matrices listed in the dataset manifest

   Claim boundary: processed 10x-style GEO supplementary MTX files are used.
   Human 3 is included as hash-demultiplexed `hash_unknown` pseudo-samples
   because the public workbook does not map HumanHashTag IDs to treatment
   labels.

## Heavy Analysis I/O

Heavy analyses are tracked in `analysis/runtime_lock/heavy_analysis_io.tsv`.
Each row records:

- `analysis_id`
- script command from the reproducibility catalog;
- expected input objects or files;
- expected output tables/manifests;
- approximate runtime class;
- rerun trigger;
- research role and claim boundary.

The GSE167339 donor robustness contract records the donor-level composition,
leave-donor-out transfer, and hash-threshold sensitivity layer used to support
three-donor public validation.

The marker stress benchmark contract records deliberately imperfect marker-prior
simulation scenarios used to compare FibroDynMix, marker scoring, and NMF/topic
baselines under marker failure modes.

The bioinformatics validation contract records fibroblast purity QC,
state-associated top genes, curated pathway enrichment, donor-level state
abundance, and FPI summaries for the current public real-data layer.

This makes the project reviewable without hiding substantial analyses inside
unit tests. Unit tests check package behavior; heavy analyses provide manuscript
evidence and should be rerun before a final manuscript freeze.

## Quality Gate

After changing dependencies, public-data scripts, or heavy-analysis contracts,
run:

```bash
Rscript scripts/run_runtime_lock.R
Rscript scripts/run_reproducibility_audit.R
Rscript scripts/run_project_maturity_audit.R
Rscript scripts/check_project_integrity.R
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual FibroDynMix_0.0.0.9000.tar.gz
```

Then remove generated check artifacts:

```bash
rm -rf FibroDynMix.Rcheck FibroDynMix_0.0.0.9000.tar.gz Rplots.pdf .DS_Store tests/testthat/_snaps
```
