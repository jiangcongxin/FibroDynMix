# iScience Method Plus Biological Insight Execution

This note records the execution status for moving FibroDynMix from a
method-first package manuscript toward an iScience-style
method-plus-biological-insight manuscript.

## Completed on 2026-06-14

### Validation-Aware NB Selection

Implemented a validation-aware NB outer-iteration selector:

- `R/nb_model_selection.R`
- `scripts/run_validation_aware_nb_selection.R`
- `tests/testthat/test-nb_model_selection.R`

Primary output:

- `analysis/validation_aware_nb_selection/validation_aware_nb_selection_manifest.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_selection_candidates.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_selected_summary.tsv`

The selector uses held-out NB objective, z stability, marker-gradient
preservation, and downstream validation. Training objective is reported but is
not used as a standalone selection criterion. Truth RMSE and dominant-state
accuracy are audit columns only when simulated truth is available.

Current light calibration selected `n_outer` values of `2`, `5`, and `10`
across scenarios and NB variants. This supports the manuscript position that
FibroDynMix should not chase count reconstruction objective alone when latent
biological composition is the target.

### scVI Baseline

Created a project-local scVI environment and ran the scVI-inclusive
core-converged benchmark:

- `.venv-scvi`
- `analysis/core_converged_method_benchmark_scvi_r5`
- `R/scvi_baseline.R`
- `scripts/run_core_converged_method_benchmark.R`

The benchmark includes marker scoring, topic NMF, FibroDynMix initializer,
FibroDynMix NB, FibroDynMix NB plus study effects, and scVI latent projection.
It should be reported conservatively: marker scoring and NMF remain competitive
or stronger for several metrics in the current bounded simulations.

### Objective-RMSE Tradeoff Figure

Generated the requested objective-versus-RMSE panel:

- `figures/core_converged_tradeoff_scvi_r5/exports/core_converged_objective_rmse_tradeoff.png`
- `figures/core_converged_tradeoff_scvi_r5/exports/core_converged_objective_rmse_tradeoff.pdf`
- `figures/core_converged_tradeoff_scvi_r5/source_data/objective_rmse_tradeoff_points.tsv`
- `figures/core_converged_tradeoff_scvi_r5/source_data/objective_rmse_tradeoff_pairs.tsv`

The 5-replicate run found 40/40 lower-objective/worse-RMSE tradeoff pairs
among tracked NB comparisons. Paired tests should be interpreted by effect
direction and magnitude rather than strong significance because finite Wilcoxon
p values were near, but not below, 0.05.

### Scale-Up Smoke

Completed one scale-up smoke run:

- `analysis/core_converged_scaleup_c50`

This used `cells_per_donor = 50` for one `batch_confounding` scenario and one
replicate. It is scale-readiness evidence, not performance evidence.

### Independent Dataset Screening

Added a reproducible dataset readiness audit:

- `scripts/audit_independent_fibroblast_datasets.R`
- `analysis/independent_fibroblast_dataset_screening/independent_fibroblast_dataset_readiness.tsv`
- `analysis/independent_fibroblast_dataset_screening/independent_fibroblast_dataset_next_actions.tsv`
- `analysis/independent_fibroblast_dataset_screening/independent_fibroblast_dataset_screening_manifest.tsv`

Highest-priority external scRNA-seq datasets for the iScience revision:

- `GSE167339`: promote now as independent non-cancer validation.
- `GSE156326`: integrate as high-priority hypertrophic scar versus normal skin
  validation.
- `GSE181316`: integrate as high-priority keloid, normal scar, and normal skin
  validation.

GSE246215 should remain the cancer-context case study, but not the only
biological validation. GSE163973 should be used either as a discovery/reference
source or excluded from independent validation claims. GSE175866 is useful as
bulk endpoint support, not as a substitute for independent scRNA-seq validation.

### Package Validation

The build-package check completed with `Status: OK`:

```bash
R CMD build .
R CMD check --no-manual --no-build-vignettes --no-multiarch FibroDynMix_0.0.0.9000.tar.gz
```

The project-local scVI virtual environment and local analysis/data folders are
excluded from the package tarball through `.Rbuildignore`.

## Follow-Up Completed

### Figure 5 External Validation

Generated the main Figure 5 from independent human fibroblast datasets:

- `scripts/make_figure5.R`
- `figures/figure5/exports/figure5.png`
- `figures/figure5/source_data/`
- `docs/figure5-external-validation.md`

The figure integrates GSE167339 donor-level validation with GSE156326 and
GSE181316 scar-context transfer/module validation. Current readout: GSE167339
leave-donor transfer convergence is 1.000, 1.000, and 0.978 across Human1,
Human2, and Human3; scar transfer convergence is 0.954 for GSE156326 and 0.960
for GSE181316; module deltas are positive in both scar cohorts but do not pass
FDR correction, so they support descriptive directionality rather than
statistically significant disease separation.

### Selected-NB Final Summaries

Reran final VI, transfer, and downstream summaries with validation-aware
selected defaults:

- `scripts/run_selected_nb_final_summaries.R`
- `analysis/selected_nb_final_summaries/selected_nb_defaults_used.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_final_summaries_manifest.tsv`
- `docs/selected-nb-final-summaries.md`

Defaults used: VI continuous `n_outer=2`, VI batch-confounding `n_outer=10`,
VI rare-transition `n_outer=2`, transfer `n_outer=10`, and real-data downstream
NB-study `n_outer=2`.

## Manuscript Position

Recommended central biological-method claim:

```text
FibroDynMix combines raw-count compositional inference with validation-aware
model selection to map conserved and disease-context-specific fibroblast state
plasticity across simulated benchmarks and public human disease datasets.
```

Recommended claim boundary:

```text
FibroDynMix improves interpretability and workflow transparency for
marker-oriented fibroblast composition, but it should not be claimed as
uniformly superior to marker scoring, NMF, or scVI embeddings across all tasks.
```

## Figure-Level Implications

- Figure 2: selected-NB rule plus simulation benchmark.
- Figure 3: marker scoring, NMF, scVI, initializer, selected NB, and VI
  comparison with explicit competitive-baseline language.
- Figure 4: GSE246215 downstream classification and marker-gradient biology.
- Figure 5: GSE167339 plus GSE156326/GSE181316 independent non-cancer
  fibroblast validation.
- Figure 6: cross-cohort transfer, uncertainty, and transition-flow summaries.

## Immediate Next Work

1. Promote `analysis/independent_fibroblast_dataset_screening` into the
   manuscript evidence map and figure source index.
2. Rerun study-effect summaries around selected defaults if the final figure
   order keeps study-effect calibration as a main claim.
3. Rewrite the Results narrative around biological plasticity axes and
   downstream utility, while proactively stating where marker scoring and NMF
   are competitive.
