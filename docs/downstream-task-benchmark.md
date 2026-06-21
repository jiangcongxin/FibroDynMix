# Downstream Task Benchmark

This evidence layer asks whether inferred state features are useful beyond
direct state-weight RMSE. It adds label-prediction tasks to the simulation
benchmark and to the public GSE246215 fibroblast atlas.

## Simulation Command

```bash
Rscript scripts/run_extended_method_benchmark.R
```

The updated benchmark now records downstream donor-grouped prediction of the
simulated `disease` label from each method's inferred feature matrix. The
classifier is a lightweight nearest-centroid model with grouped
cross-validation by donor.

Additional output columns in
`analysis/extended_method_benchmark/extended_method_benchmark_metrics.tsv`:

- `downstream_status`
- `downstream_balanced_accuracy`
- `downstream_macro_f1`
- `downstream_macro_auroc`

Current manifest-level summary:

- mean FibroDynMix-family downstream balanced accuracy: 0.7383
- mean marker/NMF/scVI-baseline downstream balanced accuracy: 0.7090
- mean FibroDynMix-family downstream macro-F1: 0.7371
- mean marker/NMF/scVI-baseline downstream macro-F1: 0.7049
- scVI status: `not_available_skipped` because Python `scvi` and `anndata`
  are not available in the current runtime

This supports the claim that FibroDynMix-family features are useful for a
simulated downstream disease task on average, while marker scoring remains
strong in some clean-prior scenarios.

## GSE246215 Command

A dedicated record for the GSE246215 downstream layer is maintained in
`docs/gse246215-downstream-benchmark.md`.

Patient-level aggregation:

```bash
Rscript scripts/run_gse246215_downstream_benchmark.R
```

Sample-level aggregation:

```bash
Rscript scripts/run_gse246215_downstream_benchmark.R \
  --aggregation-col=SampleID \
  --out=analysis/gse246215_downstream_benchmark_sample
```

The script compares:

- `fibrodynmix_z`: FibroDynMix state composition
- `marker_scoring_z`: softmax-normalized marker-score composition
- `marker_scoring_scores`: raw marker-score features
- `topic_nmf_z`: state-aligned NMF topic composition

Primary outputs:

- `gse246215_downstream_classification_metrics.tsv`
- `gse246215_downstream_features.tsv`
- `gse246215_biological_gradient_validation.tsv`
- `gse246215_downstream_marker_coverage.tsv`
- `gse246215_downstream_manifest.tsv`

Current PatientID-aggregated cancer-type classification:

- `topic_nmf_z`: balanced accuracy 0.8100, macro-F1 0.8021
- `marker_scoring_z`: balanced accuracy 0.7074, macro-F1 0.6131
- `fibrodynmix_z`: balanced accuracy 0.6314, macro-F1 0.7585
- `marker_scoring_scores`: balanced accuracy 0.3919, macro-F1 0.5218

Current SampleID-aggregated cancer-type classification is similar:

- `topic_nmf_z`: balanced accuracy 0.8084, macro-F1 0.8036
- `marker_scoring_z`: balanced accuracy 0.6762, macro-F1 0.5682
- `fibrodynmix_z`: balanced accuracy 0.6106, macro-F1 0.7641
- `marker_scoring_scores`: balanced accuracy 0.4172, macro-F1 0.5758

These real-data results do not support a blanket claim that FibroDynMix `z`
beats all baselines for cancer-type classification. They do show that
FibroDynMix `z` is competitive with marker-score representations and has higher
macro-F1/AUROC than raw marker scores in this bounded run, while NMF topics are
the strongest cancer-type separator.

## Biological Gradient Case Study

The GSE246215 script also tests whether canonical fibroblast programs vary
continuously along expected FibroDynMix state weights:

- myofibroblast marker score vs myofibroblast `z`: Spearman rho 0.4236
- ECM marker score vs ECM-remodeling `z`: Spearman rho 0.6714
- inflammatory marker score vs inflammatory `z`: Spearman rho 0.3475

This supports a within-dataset biological gradient sanity check. Because the
same marker families orient the model, it should be presented as a case study
and not as independent marker discovery.

## Claim Boundary

The downstream layer strengthens the manuscript by answering "what is `z` good
for?" with explicit prediction tasks and feature comparisons. The defensible
claim is representation utility and biological-gradient recovery under bounded
settings, not universal dominance over marker scoring or NMF in every
downstream task.
