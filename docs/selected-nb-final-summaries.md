# Selected-NB Final Summaries

This analysis reruns final VI, transfer, and downstream summaries using
validation-aware selected NB defaults rather than fixed smoke settings.

## Command

```bash
Rscript scripts/run_selected_nb_final_summaries.R
```

## Outputs

- `analysis/selected_nb_final_summaries/selected_nb_defaults_used.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_vi_benchmark_metrics.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_vi_benchmark_summary.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_transfer_benchmark.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_transfer_benchmark_summary.tsv`
- `analysis/selected_nb_final_summaries/selected_nb_gse246215_downstream_classification_metrics.tsv`
- `analysis/selected_nb_final_summaries/gse246215_downstream_selected_nb/`
- `analysis/selected_nb_final_summaries/selected_nb_final_summaries_manifest.tsv`

## Defaults Used

The defaults are read from
`analysis/validation_aware_nb_selection/validation_aware_nb_selected_summary.tsv`.

| task | scenario | variant | selected n_outer |
|---|---|---|---:|
| VI continuous | continuous | nb | 2 |
| VI batch confounding | batch_confounding | nb | 10 |
| VI rare transition | rare_transition | nb | 2 |
| Cross-cohort transfer | batch_confounding | nb | 10 |
| Real-data downstream | batch_confounding | nb_study | 2 |

## Current Readout

- Mean selected-NB VI RMSE: 0.202.
- Mean selected-NB VI downstream balanced accuracy: 0.766.
- Mean selected-NB transfer RMSE: 0.215.
- Mean selected-NB transfer z convergence rate: 1.000.
- GSE246215 downstream benchmark remains baseline-competitive: topic NMF is the
  best cancer-type representation in this case study
  (balanced accuracy = 0.810), while FibroDynMix z reaches 0.631.

## Claim Boundary

These summaries support the move away from objective-only or smoke-only NB
settings. They do not establish universal superiority over marker scoring, NMF,
or scVI. The GSE246215 downstream rerun remains a representation-utility case
study from one GEO accession, not diagnostic validation.
