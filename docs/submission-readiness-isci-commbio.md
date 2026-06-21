# Submission Readiness for iScience / Communications Biology

## Current Positioning

FibroDynMix is now positioned as a computational methods manuscript rather than
a marker-score application. The defensible central claim is:

```text
FibroDynMix models fibroblast state plasticity as latent compositional programs
under a raw-count negative-binomial generative framework, with hierarchy-aware
effects, posterior uncertainty, transition-flow summaries, and cross-cohort
transfer diagnostics.
```

The final manuscript framing should remain **method-first**. Public mouse,
GSE246215, GSE167339, GSE156326, and GSE181316 analyses can be used as bounded
count-matrix validation and biological-utility evidence, but not as a completed
human disease atlas or a cross-atlas generalization claim.

## Evidence Already Present

- Raw-count likelihood and optimizer:
  - `fibrodynmix_nb_loglik()`
  - `fibrodynmix_nb_objective()`
  - `fit_fibrodynmix_nb()`
- Latent state simplex and weak marker priors:
  - `fit_fibrodynmix_initializer()`
  - `fit_fibrodynmix_nb()`
- Study and donor effects:
  - `fit_fibrodynmix_nb(..., fit_study_effect = TRUE, fit_donor_effect = TRUE)`
  - Figure 4 study-effect calibration
- Posterior uncertainty:
  - `bootstrap_fibrodynmix()`
  - `fit_fibrodynmix_vi()`
  - `evaluate_posterior_intervals()`
  - `calibrate_posterior_interval_scale()`
- Transition and plasticity:
  - `estimate_transition_flow()`
  - `compute_fpi()`
  - Figure 6
- Cross-cohort transfer:
  - `fit_fibrodynmix_transfer()`
  - `run_cross_cohort_transfer_benchmark()`
  - public real-data bidirectional transfer smoke
- Multi-public real-data validation:
  - `scripts/run_multi_public_realdata_validation.R`
  - pooled public raw-count NB fit with study effects
  - leave-dataset-out public transfer diagnostics
  - dataset-level state composition and transition-flow summaries
- Independent human GEO validation:
  - `scripts/prepare_gse246215_fibroblast_inputs.R`
  - `analysis/independent_geo_gse246215_validation`
  - four GSE246215 human cancer-type fibroblast count matrices
  - leave-dataset-out transfer diagnostics across GC, HCC, NSCLC, and TNBC
- GSE246215 sensitivity and QC:
  - `scripts/run_gse246215_sensitivity.R`
  - random downsampling and library-size q95 trimming
  - transfer convergence and state-composition variability summaries
  - explicit HCC library-size outlier boundary
- Extended method benchmark:
  - `fit_topic_nmf_baseline()`
  - `scripts/run_extended_method_benchmark.R`
  - marker scoring, NMF/topic baseline, FibroDynMix NB, study-effect NB, and VI
  - explicit boundary that marker scoring remains competitive in some
    marker-aligned simulations
- Core converged method benchmark:
  - `scripts/run_core_converged_method_benchmark.R`
  - `n_outer = 10` and `20` rerun for the core NB-family methods
  - explicit evidence that lower NB objective can worsen latent-state RMSE and
    downstream balanced accuracy
  - supports validation-aware stopping/model selection rather than
    objective-only selection
- Validation-aware NB selection:
  - `select_fibrodynmix_nb_model()`
  - `scripts/run_validation_aware_nb_selection.R`
  - candidate `n_outer = 2, 5, 10, 20`
  - held-out NB objective, z stability, marker-gradient preservation, and
    downstream validation are used for selection
  - training objective is reported but not used as a standalone selector
- scVI-inclusive core benchmark:
  - project-local `.venv-scvi` environment with `scvi-tools`, `anndata`, and
    `torch`
  - `analysis/core_converged_method_benchmark_scvi_r5`
  - marker scoring, NMF, FibroDynMix initializer, FibroDynMix NB, FibroDynMix
    NB plus study effects, and scVI latent projection are compared
- GSE246215 downstream representation benchmark:
  - `scripts/run_gse246215_downstream_benchmark.R`
  - patient- and sample-level cancer-type classification from inferred
    representations
  - biological-gradient validation for myofibroblast, ECM-remodeling, and
    inflammatory programs
  - explicit boundary that NMF topics are strongest for cancer-type
    classification in the current GSE246215 task
- Independent fibroblast dataset screening:
  - `scripts/audit_independent_fibroblast_datasets.R`
  - `analysis/independent_fibroblast_dataset_screening`
  - prioritizes GSE167339, GSE156326, and GSE181316 for independent non-cancer
    fibroblast validation
  - keeps GSE246215 as the cancer-context case study rather than the sole
    biological validation
- Figure 5 external biological validation:
  - `scripts/make_figure5.R`
  - `figures/figure5`
  - integrates GSE167339 donor/transfer robustness with GSE156326/GSE181316
    scar-context state composition and module validation
- Selected-NB final summaries:
  - `scripts/run_selected_nb_final_summaries.R`
  - `analysis/selected_nb_final_summaries`
  - final VI, transfer, and downstream summaries use validation-aware selected
    `n_outer` defaults
- Reproducibility audit:
  - `scripts/run_reproducibility_audit.R`
  - `analysis/reproducibility_audit/analysis_catalog.tsv`
  - `analysis/reproducibility_audit/script_inventory.tsv`
  - `docs/reproducibility-runbook.md`
- Reproducible figure package:
  - Figures 1-6
  - public real-data smoke figure
  - source data, manifests, legends, QC exports
- R package maturity:
  - `NEWS.md`
  - `inst/CITATION`
  - `cran-comments.md`
  - `_pkgdown.yml`
  - `inst/doc/fibrodynmix-workflow.md`
  - `R CMD check --no-manual` status OK in the latest validation

## Evidence Matrix

Run:

```bash
Rscript scripts/run_project_maturity_audit.R
```

Expected outputs:

- `analysis/project_maturity/project_maturity_manifest.tsv`
- `analysis/project_maturity/exported_function_index.tsv`
- `analysis/project_maturity/figure_claim_matrix.tsv`
- `analysis/project_maturity/manuscript_evidence_matrix.tsv`
- `analysis/project_maturity/submission_readiness_checklist.tsv`

These files provide the audit trail for package maturity and manuscript claim
discipline.

## Manuscript Claim Boundaries

The current manuscript should not claim:

- full amortized inference;
- a joint posterior over every model parameter;
- lineage tracing or observed temporal state conversion;
- diagnostic, causal, therapeutic, or mechanistic disease conclusions;
- a completed human disease atlas or full cross-atlas/multi-cohort
  generalization;
- disease mechanism discovery from the public mouse pooled-count smoke data.

The current manuscript can claim:

- raw-count generative modeling rather than normalized gene-set scoring;
- inferred latent state mixtures rather than fixed signature scores;
- study/donor effect modeling under ridge-penalized hierarchy terms;
- posterior state uncertainty under bootstrap and lightweight logistic-normal VI;
- simulation-calibrated posterior interval diagnostics;
- simulation comparison against marker scoring and mature-package NMF/topic
  baselines;
- downstream disease-label prediction in simulation, where FibroDynMix-family
  features outperform the marker/NMF baseline average under the current
  bounded benchmark;
- high-outer NB optimizer sensitivity showing that objective improvement,
  latent-state RMSE, and downstream utility can diverge;
- transition-flow summaries as cross-sectional optimal transport;
- reproducible public raw-count smoke execution;
- registry-driven multi-public count-matrix validation with bounded
  leave-dataset-out transfer diagnostics;
- independent public human fibroblast count-matrix validation from GSE246215,
  bounded by processed-matrix and sampling caveats;
- GSE246215 robustness under random downsampling and simple high-library-size
  trimming, with HCC treated as QC-sensitive;
- GSE246215 downstream-feature case study, including cancer-type prediction
  from patient/sample-level features and expected fibroblast marker-gradient
  recovery;
- bounded external biological utility from GSE167339 donor/transfer robustness
  and GSE156326/GSE181316 scar-context module trends, without diagnostic,
  causal, or therapeutic interpretation;
- source-backed reproducibility audit linking scripts, analysis outputs,
  documents, and claim boundaries;
- mature R package engineering gates.

## Remaining Work Before Submission

- Final Figure 5 ordering is resolved: `figures/figure5` is the external
  biological validation package, and the previous bootstrap uncertainty Figure 5
  is now `figures/supplementary_uncertainty`.
- Rerun study-effect summaries around validation-aware selected NB defaults if
  study-effect calibration remains a main figure claim.
- Decide whether GSE163973 is a discovery/reference dataset or excluded from
  independent validation claims; it should not serve both roles.
- Keep scVI, NMF, and marker scoring in the comparison tables and state
  explicitly where they are competitive or stronger.
- Polish the manuscript draft into journal style with Methods, Results,
  limitations, and data/code availability matched to the generated evidence
  matrices.
- Keep the final target as a method-first paper; use external human datasets as
  bounded biological-utility evidence rather than as an atlas-level claim.

## Drafting Files

- Manuscript skeleton: `docs/manuscript-skeleton-isci-commbio.md`
- Manuscript draft: `docs/manuscript-draft-isci-commbio.md`
- Reviewer risk matrix: `docs/reviewer-risk-matrix.md`
- Multi-public validation: `docs/multi-public-realdata-validation.md`
- Independent GSE246215 validation:
  `docs/independent-geo-gse246215-validation.md`
- GSE246215 sensitivity: `docs/gse246215-sensitivity.md`
- GSE246215 downstream benchmark:
  `docs/gse246215-downstream-benchmark.md`
- Extended method benchmark: `docs/extended-method-benchmark.md`
- Core converged benchmark: `docs/core-converged-method-benchmark.md`
- iScience execution log:
  `docs/iscience-method-biological-execution.md`
- Figure 5 external validation:
  `docs/figure5-external-validation.md`
- Selected-NB final summaries:
  `docs/selected-nb-final-summaries.md`
- NB convergence benchmark: `docs/nb-convergence-benchmark.md`
- Reproducibility runbook: `docs/reproducibility-runbook.md`
