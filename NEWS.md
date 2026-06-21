# FibroDynMix 0.0.0.9000

## Major Changes

- Added raw-count negative-binomial likelihood, deviance, and optimization
  objective functions.
- Added alternating NB optimizer with latent simplex updates, marker-orientation
  penalty, rollback diagnostics, and early stopping.
- Added ridge-penalized study-level and donor-level gene effects.
- Added bootstrap uncertainty summaries for cell state weights, entropy,
  sample-level composition, and marker-program stability.
- Added logistic-normal VI posterior skeleton with posterior draws, credible
  intervals, sample composition intervals, and ELBO-like trace.
- Added posterior interval calibration utilities for simulation studies.
- Added mature-package NMF/topic baseline support through
  `fit_topic_nmf_baseline()` and extended method benchmark outputs.
- Added cross-cohort transfer API and leave-study-out simulation benchmark.
- Added public raw-count smoke analyses and bidirectional public real-data
  transfer smoke analysis.
- Added registry-driven multi-public real-data validation with pooled public
  raw-count fitting, study-effect adjustment, transition flow, and
  leave-dataset-out transfer diagnostics.
- Added independent GEO GSE246215 human fibroblast atlas preparation and
  validation outputs with four cancer-type count matrices and leave-dataset-out
  transfer diagnostics.
- Added manuscript-facing Figure 1-6 and public real-data smoke figure packages
  with source data, manifests, legends, and export QC.

## Package Quality

- Added project integrity gate for figure packages, required analysis outputs,
  and release metadata.
- Added CRAN-style check notes, citation metadata, pkgdown configuration, and
  workflow documentation.
- Added project maturity audit, submission-readiness documentation, manuscript
  skeleton, manuscript draft, and reviewer-risk matrix.
- Current validation target: `testthat::test_local()`,
  `scripts/check_project_integrity.R`, and `R CMD check --no-manual`.

## Claim Boundary

The current package supports a mature method-development workflow for a
raw-count FibroDynMix prototype. It is not yet a full amortized inference engine
or a fully hierarchical variational posterior over every model parameter, and
public real-data analyses are bounded count-matrix validations and case studies
rather than diagnostic, causal, therapeutic, or completed human disease-atlas
evidence. The multi-public validation layer currently uses two public
condition-specific count matrices from one Dryad record and should not be
described as an independent human multi-study atlas. The GSE246215 layer uses
public processed count matrices rather than FASTQ-level reprocessing, so
cancer-type composition results require additional sample-level QC before
biological interpretation. GSE167339, GSE156326, and GSE181316 provide bounded
external biological-utility evidence, not definitive cross-atlas generalization.
