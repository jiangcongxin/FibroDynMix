# Variational Posterior Layer

## Purpose

`fit_fibrodynmix_vi()` adds an auditable posterior layer to the implemented
raw-count FibroDynMix optimizer:

```text
raw counts
-> fit NB FibroDynMix mode
-> fix alpha, beta, phi and hierarchy effects
-> define q(eta_i) as diagonal logistic-normal
-> z_i = softmax(eta_i)
-> return posterior draws, credible intervals, and ELBO-like trace
```

## Current Scope

The current implementation is a lightweight mean-field posterior over cell
state logits. It supports:

- cell-level posterior state draws;
- state-weight credible intervals;
- entropy credible intervals;
- optional sample-level composition intervals;
- an ELBO-like Monte Carlo trace for posterior scale refinement.
- simulation interval coverage and post-hoc interval-scale calibration.

## Reproducible Smoke Run

```bash
Rscript scripts/run_vi_posterior.R
```

Expected outputs:

- `analysis/vi_posterior/vi_elbo_trace.tsv`
- `analysis/vi_posterior/vi_cell_state_intervals.tsv`
- `analysis/vi_posterior/vi_cell_entropy_intervals.tsv`
- `analysis/vi_posterior/vi_sample_composition_intervals.tsv`
- `analysis/vi_posterior/vi_state_recovery_metrics.tsv`
- `analysis/vi_posterior/vi_manifest.tsv`

For simulation calibration against known latent truth:

```bash
Rscript scripts/run_vi_benchmark.R
```

Expected benchmark outputs:

- `analysis/vi_benchmark/vi_benchmark_metrics.tsv`
- `analysis/vi_benchmark/vi_benchmark_summary.tsv`
- `analysis/vi_benchmark/vi_optimizer_diagnostics.tsv`
- `analysis/vi_benchmark/vi_benchmark_manifest.tsv`

## Claim Boundary

This layer supports the claim that FibroDynMix now exposes posterior state
uncertainty under a raw-count NB generative program. It is not yet a full
amortized inference engine, and it does not yet place a joint variational
posterior over all model parameters, study effects, donor effects, and
dispersion parameters. The interval-scale calibration uses simulated truth for
method diagnostics and should not be presented as real-data posterior
calibration.
