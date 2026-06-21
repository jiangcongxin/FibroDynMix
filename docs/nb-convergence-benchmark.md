# NB Convergence Benchmark

This benchmark tests whether the current bounded `n_outer = 2` NB optimizer
setting is sufficient for manuscript-level benchmark claims.

## Command

```bash
Rscript scripts/run_nb_convergence_benchmark.R
```

The default bounded run evaluates:

- scenarios: `continuous`, `discrete`, `batch_confounding`, `rare_transition`
- methods: `fibrodynmix_nb`, `fibrodynmix_nb_study`
- outer grid: `2, 5, 10, 20`
- replicates: 2

Each fit records state-recovery metrics, downstream disease-label metrics,
optimizer diagnostics, objective traces, and `z` changes relative to the
`n_outer = 20` reference fit on the same simulated data.

## Outputs

Files are written to `analysis/nb_convergence_benchmark/`:

- `nb_convergence_metrics.tsv`
- `nb_convergence_summary.tsv`
- `nb_convergence_objective_traces.tsv`
- `nb_convergence_z_delta_vs_reference.tsv`
- `nb_convergence_delta_vs_min_n_outer.tsv`
- `nb_convergence_manifest.tsv`

## Current Result

The current bounded run does **not** support treating `n_outer = 2` as a
convergence-adequate main benchmark setting. It also shows that the NB issue is
not solved by simply increasing the outer-iteration budget.

Manifest-level diagnostics:

- 64 method/setting rows and 656 objective-trace rows
- maximum best-objective gap before the `n_outer = 20` reference: 0.0345
- maximum mean absolute `z` difference at `n_outer = 2` vs `n_outer = 20`:
  0.1072
- minimum dominant-state agreement at `n_outer = 2` vs `n_outer = 20`: 0.5000
- maximum absolute RMSE change across the grid relative to `n_outer = 2`:
  0.0552
- maximum absolute downstream balanced-accuracy change relative to
  `n_outer = 2`: 0.2344
- no `n_outer = 20` run stopped by the current early-stopping criterion

The objective continues improving as outer iterations increase. However, lower
NB objective does not consistently improve latent-state RMSE or downstream
classification. In several scenarios, RMSE worsens as the optimizer moves to
lower objective values. This means the issue is not only "run more iterations";
it also exposes a model-selection/regularization problem for the current NB
objective.

## Core High-Outer Follow-Up

The core method benchmark has now been rerun under higher NB outer budgets:
`docs/core-converged-method-benchmark.md`.

That follow-up uses `n_outer = 10` and `20` for the core NB-family methods and
writes results to `analysis/core_converged_method_benchmark/`.

Manifest-level summary from the follow-up:

- 80 method rows
- mean FibroDynMix NB-family RMSE: 0.2622
- mean marker/NMF baseline RMSE: 0.2514
- mean FibroDynMix NB-family downstream balanced accuracy: 0.6338
- mean marker/NMF baseline downstream balanced accuracy: 0.6836
- objective/RMSE tradeoff pairs: 16 / 16 showed lower objective and worse RMSE

This upgrades the conclusion: the convergence benchmark first showed that
`n_outer = 2` is inadequate as submission evidence; the high-outer core rerun
confirms a deeper objective/latent-quality tradeoff. The current NB training
objective can improve while the target manuscript metrics degrade.

## Interpretation

`n_outer = 2` remains appropriate for smoke tests and fast execution checks. It
is not sufficient for primary iScience/Communications Biology claims about NB
optimizer performance, VI behavior around the NB mode, transfer initialized from
NB fits, or study-effect calibration.

The manuscript should not present existing `n_outer = 2` NB results as
converged optimizer results. It also should not choose `n_outer = 20` merely
because it has the lower training objective. The defensible options are:

- define an explicit early-stopping/model-selection rule and report objective
  traces, `z` stability, marker-gradient preservation, and downstream behavior;
- rerun primary NB-related benchmarks using the selected rule or a justified
  higher-outer configuration such as the current `n_outer = 10` candidate;
- treat current NB/VI/transfer results as bounded prototype evidence until the
  convergence and regularization behavior is resolved.

## Claim Boundary

This is a bounded convergence sensitivity analysis with small simulated
datasets. Together with the core high-outer follow-up, it establishes that
`n_outer = 2` should not be used as the sole evidence base for NB-related
submission claims and that objective-only selection is not adequate for the
current optimizer. The next submission-grade step is validation-aware NB
selection followed by a rerun of NB, VI, transfer, and study-effect summaries
around the selected mode.
