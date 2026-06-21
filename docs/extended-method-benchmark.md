# Extended Method Benchmark

This non-visual benchmark adds mature-package topic/NMF comparison, optional
scVI latent comparison, and a downstream disease-label task to the simulation
evidence layer. The purpose is to separate FibroDynMix from direct marker
scoring, ordinary unsupervised topic/matrix factorization, and deep latent
representation baselines when those dependencies are available.

## Command

```bash
Rscript scripts/run_extended_method_benchmark.R
```

## Methods

The benchmark runs four simulated scenarios:

- continuous mixtures
- discrete states
- batch confounding
- rare transition-like cells

It compares:

- `marker_scoring`
- `topic_nmf`
- `scvi_latent` when `reticulate`, Python `scvi`, and Python `anndata` are
  available
- `fibrodynmix_initializer`
- `fibrodynmix_nb`
- `fibrodynmix_nb_study`
- `fibrodynmix_vi`

The topic/NMF comparator uses the `NMF` package when available and falls back to
an internal KL-NMF implementation otherwise. Marker genes are used only after
factorization to align learned topics to fibroblast-state labels. The optional
scVI comparator fits a small scVI model, clusters the latent space, projects
cluster distances to a simplex, and aligns clusters to states by marker-score
enrichment. This projection is a representation readout rather than a native
scVI mixture model.

Each method row also includes donor-grouped downstream classification of the
simulated `disease` label from the inferred feature matrix:

- `downstream_balanced_accuracy`
- `downstream_macro_f1`
- `downstream_macro_auroc`

## Outputs

The analysis writes to `analysis/extended_method_benchmark/`:

- `extended_method_benchmark_metrics.tsv`
- `extended_method_benchmark_summary.tsv`
- `extended_method_optimizer_diagnostics.tsv`
- `extended_method_rankings.tsv`
- `extended_method_benchmark_manifest.tsv`

## Current Result

The current run contains 48 benchmark rows across four scenarios and six
methods. The topic baseline used the `NMF` backend. scVI was skipped because
Python `scvi` and `anndata` were not available in the current runtime.
NB-family methods in this run use the bounded runtime setting `n_outer = 2`;
see `docs/nb-convergence-benchmark.md` before using these rows as
convergence-level optimizer evidence.

Mean RMSE across FibroDynMix-family methods is lower than the marker/NMF/scVI
baseline average:

- mean FibroDynMix-family RMSE: 0.2296
- mean baseline RMSE: 0.2568
- mean FibroDynMix-family dominant-state accuracy: 0.3815
- mean baseline dominant-state accuracy: 0.3086
- mean FibroDynMix-family downstream balanced accuracy: 0.7383
- mean baseline downstream balanced accuracy: 0.7090
- mean FibroDynMix-family downstream macro-F1: 0.7371
- mean baseline downstream macro-F1: 0.7049

The NMF/topic baseline performs worse than FibroDynMix-family methods across all
four scenarios in RMSE. Marker scoring remains competitive in several low-noise
simulation settings, so the manuscript should not claim uniform superiority over
marker scoring. The stronger defensible claim is that FibroDynMix provides
generative raw-count modeling, hierarchy-aware effects, uncertainty, transfer
diagnostics, and competitive downstream disease-label prediction that marker
scoring and ordinary NMF/topic models do not provide as a unified framework.

## Claim Boundary

This is a simulation benchmark with known latent truth. It supports
method-comparison claims against marker scoring and topic/NMF baselines, but it
does not replace independent biological validation. It also shows that marker
scoring can remain a strong baseline when the simulation exactly matches the
marker-prior structure, which should be stated directly in the Results and
Discussion.

The scVI baseline is optional and environment-dependent. When included, its
simplex projection should be described as a post hoc readout from a deep latent
space, not as a native scVI state-composition estimate.

NB optimizer rows in this benchmark are bounded runtime results. The separate
NB convergence benchmark shows that `n_outer = 2` is not sufficient as the sole
evidence base for primary NB optimizer claims.
