# Marker Stress Benchmark

This benchmark redesigns the simulation comparison around marker-scoring
failure modes. The goal is not to show universal superiority, but to test when
direct marker scoring becomes brittle and whether FibroDynMix retains useful
state recovery, likelihood diagnostics, uncertainty, and transfer-compatible
state programs.

## Run

```bash
Rscript scripts/run_marker_stress_benchmark.R
```

Outputs:

- `analysis/marker_stress_benchmark/marker_stress_benchmark_manifest.tsv`
- `analysis/marker_stress_benchmark/marker_stress_benchmark_metrics.tsv`
- `analysis/marker_stress_benchmark/marker_stress_benchmark_summary.tsv`
- `analysis/marker_stress_benchmark/marker_stress_rankings.tsv`
- `analysis/marker_stress_benchmark/marker_stress_contrast_vs_marker.tsv`
- `analysis/marker_stress_benchmark/marker_stress_prior_audit.tsv`

## Stress Modes

The benchmark uses simulated truth with deliberately imperfect marker priors:

- `clean_prior`: original marker priors.
- `missing_markers`: only one third of true markers retained.
- `corrupted_prior`: half of marker priors replaced by non-marker genes.
- `shared_markers`: non-marker decoys added to every state.
- `swapped_markers`: half of markers swapped across neighboring states.
- `hidden_program_corrupted_prior`: corrupted priors plus diffuse hidden
  state signal in non-marker genes.
- `hidden_program_missing_markers`: missing priors plus diffuse hidden state
  signal in non-marker genes.

## Current Result

Current bounded benchmark:

- 7 stress modes
- 5 methods
- 2 replicates per stress mode
- 70 method-level rows
- topic baseline backend: `NMF`

Summary:

- NMF/topic remains the weakest comparator by RMSE in all stress modes.
- Marker scoring remains strong under clean priors and some missing-marker
  settings.
- FibroDynMix initializer has lower RMSE than marker scoring in
  `corrupted_prior`, `shared_markers`, and `hidden_program_corrupted_prior`.
- FibroDynMix NB/VI improve dominant-state accuracy over marker scoring in
  `swapped_markers` and match or improve it in `hidden_program_corrupted_prior`.

Manifest-level summary:

- modes where any FibroDynMix variant beats marker scoring by RMSE: 3 / 7
- modes where any FibroDynMix variant beats marker scoring by dominant
  accuracy: 2 / 7

## Claim Boundary

This benchmark does not support the claim that FibroDynMix universally
outperforms marker scoring on all state-weight RMSE metrics. It supports a more
precise claim: FibroDynMix is competitive with marker scoring, clearly stronger
than NMF/topic in these settings, and provides additional generative likelihood,
uncertainty, study/donor, transfer, and transition-flow outputs that marker
scoring does not provide.
