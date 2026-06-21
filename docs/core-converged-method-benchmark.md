# Core Converged Method Benchmark

This benchmark reruns the core simulation comparison under higher NB
outer-iteration budgets, separating fast smoke settings from candidate
submission settings. It follows the convergence sensitivity analysis in
`docs/nb-convergence-benchmark.md`.

## Command

```bash
Rscript scripts/run_core_converged_method_benchmark.R
```

Default design:

- scenarios: `continuous`, `discrete`, `batch_confounding`, `rare_transition`
- outer grid: `n_outer = 10, 20`
- replicates: 2
- methods: `marker_scoring`, `topic_nmf`, `fibrodynmix_initializer`,
  `fibrodynmix_nb`, `fibrodynmix_nb_study`
- VI is excluded by default because the purpose is to isolate NB optimizer
  behavior before rerunning posterior summaries around a selected NB mode.

## Outputs

Files are written to `analysis/core_converged_method_benchmark/`:

- `core_converged_benchmark_metrics.tsv`
- `core_converged_benchmark_summary.tsv`
- `core_converged_optimizer_diagnostics.tsv`
- `core_converged_benchmark_rankings.tsv`
- `core_converged_objective_rmse_tradeoff.tsv`
- `core_converged_paired_tests.tsv`
- `core_converged_benchmark_manifest.tsv`

## Extended Runs Executed

Additional core-converged runs were executed on 2026-06-14:

```bash
RETICULATE_PYTHON=/Volumes/bioSSD/RCode/FibroDynMix/.venv-scvi/bin/python \
  Rscript scripts/run_core_converged_method_benchmark.R \
  --out=analysis/core_converged_method_benchmark_scvi_r5 \
  --n-replicates=5 \
  --include-scvi \
  --scvi-max-epochs=20
```

The scVI environment is project-local under `.venv-scvi` and contains
`scvi-tools 1.4.3`, `anndata 0.12.16`, and `torch 2.12.0`.

The 5-replicate scVI-inclusive run wrote:

- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_benchmark_metrics.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_benchmark_summary.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_optimizer_diagnostics.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_benchmark_rankings.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_objective_rmse_tradeoff.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_paired_tests.tsv`
- `analysis/core_converged_method_benchmark_scvi_r5/core_converged_benchmark_manifest.tsv`

Manifest-level summary for this run:

- 240 method rows
- methods: marker scoring, topic NMF, FibroDynMix initializer, FibroDynMix NB,
  FibroDynMix NB + study, and scVI latent
- mean FibroDynMix NB-family RMSE: 0.2567
- mean marker/NMF baseline RMSE: 0.2457
- mean FibroDynMix NB-family downstream balanced accuracy: 0.6555
- mean marker/NMF baseline downstream balanced accuracy: 0.7055
- objective/RMSE tradeoff pairs: 40 / 40 showed lower objective and worse RMSE
- paired-test rows: 480

The paired Wilcoxon tests should be interpreted conservatively at
`n_replicates = 5`; the minimum finite Wilcoxon p value across tracked metrics
was approximately 0.054 to 0.059. The paired-test table is still useful for
effect directions and paired deltas, but not as a standalone significance claim.

An Objective-RMSE tradeoff figure was generated from the 5-replicate run:

- `figures/core_converged_tradeoff_scvi_r5/exports/core_converged_objective_rmse_tradeoff.png`
- `figures/core_converged_tradeoff_scvi_r5/exports/core_converged_objective_rmse_tradeoff.pdf`
- `figures/core_converged_tradeoff_scvi_r5/source_data/objective_rmse_tradeoff_points.tsv`
- `figures/core_converged_tradeoff_scvi_r5/source_data/objective_rmse_tradeoff_pairs.tsv`

The scale-up smoke run used one `batch_confounding` replicate with
`cells_per_donor = 50` and `n_outer = 10`:

```bash
Rscript scripts/run_core_converged_method_benchmark.R \
  --out=analysis/core_converged_scaleup_c50 \
  --scenarios=batch_confounding \
  --n-replicates=1 \
  --n-outer-grid=10 \
  --cells-per-donor=50 \
  --n-genes=90
```

It completed successfully with 5 method rows. This supports feasibility of a
larger-cell simulation smoke, but because it is one scenario and one replicate,
it should be used as scale-readiness evidence rather than performance evidence.

A validation-aware NB selector was then implemented and run:

```bash
Rscript scripts/run_validation_aware_nb_selection.R \
  --out=analysis/validation_aware_nb_selection \
  --scenarios=continuous,discrete,batch_confounding,rare_transition \
  --n-replicates=1 \
  --n-outer-grid=2,5,10,20 \
  --variants=nb,nb_study
```

Primary outputs:

- `analysis/validation_aware_nb_selection/validation_aware_nb_selection_candidates.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_selected_summary.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_marker_gradient_components.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_splits.tsv`
- `analysis/validation_aware_nb_selection/validation_aware_nb_selection_manifest.tsv`

The selector combines held-out NB objective, z stability, marker-gradient
preservation, and downstream validation. Training objective is reported but is
not used alone for selection. Simulated truth metrics are included only as
audit columns.

## Current Result

The higher-outer run does not rescue the current NB optimizer as a primary
state-recovery winner. Across `n_outer = 10` and `20`, marker scoring and the
initializer remain stronger than full NB in RMSE for most scenarios.

Manifest-level summary:

- 80 method rows
- mean FibroDynMix NB-family RMSE: 0.2622
- mean marker/NMF baseline RMSE: 0.2514
- mean FibroDynMix NB-family downstream balanced accuracy: 0.6338
- mean marker/NMF baseline downstream balanced accuracy: 0.6836
- objective/RMSE tradeoff pairs: 16 / 16 showed lower objective and worse RMSE

Mean NB-family behavior:

| method | n_outer | mean RMSE | mean downstream balanced accuracy |
|---|---:|---:|---:|
| fibrodynmix_nb | 10 | 0.2541 | 0.6758 |
| fibrodynmix_nb | 20 | 0.2663 | 0.6328 |
| fibrodynmix_nb_study | 10 | 0.2582 | 0.6289 |
| fibrodynmix_nb_study | 20 | 0.2701 | 0.5977 |

Thus `n_outer = 10` is a better candidate than `n_outer = 20` in the current
bounded simulation, but it should be selected through an explicit
model-selection rule rather than by training objective alone.

## Why Lower NB Objective Can Worsen z RMSE

The NB objective is a penalized training count likelihood. It is not a supervised
loss on the true simulated latent state weights. Therefore, improving the NB
objective can move the fit toward better reconstruction of sampled counts while
moving the inferred `z` away from the data-generating `z`.

Mechanistically, several factors contribute:

- The model is partly non-identifiable: `alpha`, `beta`, study effects, and
  cell-level `z` can trade off while producing similar fitted means.
- Marker priors are weak orientation terms, not hard constraints; later
  iterations can fit high-count or noisy genes at the expense of marker-aligned
  latent truth.
- The benchmark evaluates latent truth RMSE, but the optimizer sees only counts,
  penalties, and library offsets.
- The fitted model is misspecified relative to finite simulated data because
  dispersion is estimated, genes are filtered/bounded, and local optimization
  uses approximate alternating updates.
- Current early stopping is objective-based, so it cannot detect that latent
  `z` quality or downstream utility is degrading while training objective
  improves.

The practical implication is that fixed high `n_outer` should not become the
default merely because it lowers the training objective. The implemented
validation-aware selector now provides the default rule for final reruns.

## Claim Boundary

This run supports a stricter manuscript position: full NB optimization is
implemented and auditable, but objective-only stopping is not a
submission-grade state-recovery rule. Existing NB, VI, transfer, and
study-effect summaries should be interpreted as bounded until final reruns are
anchored to the validation-aware selected defaults.
