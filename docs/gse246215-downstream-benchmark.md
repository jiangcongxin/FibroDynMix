# GSE246215 Downstream Benchmark

This document records the dedicated downstream representation benchmark for
the public GSE246215 fibroblast atlas. It complements the broader downstream
benchmark in `docs/downstream-task-benchmark.md` and the independent GEO
validation in `docs/independent-geo-gse246215-validation.md`.

## Commands

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

## Outputs

Patient-level outputs are written to
`analysis/gse246215_downstream_benchmark/`.

Sample-level outputs are written to
`analysis/gse246215_downstream_benchmark_sample/`.

Primary files:

- `gse246215_downstream_classification_metrics.tsv`
- `gse246215_downstream_features.tsv`
- `gse246215_biological_gradient_validation.tsv`
- `gse246215_downstream_marker_coverage.tsv`
- `gse246215_downstream_manifest.tsv`

## Feature Sets

The benchmark compares four representations:

- `fibrodynmix_z`: FibroDynMix state composition.
- `marker_scoring_z`: softmax-normalized marker-score composition.
- `marker_scoring_scores`: raw marker-score features.
- `topic_nmf_z`: state-aligned NMF topic composition.

Cancer type is predicted with grouped cross-validation after aggregating
fibroblast-level features by patient or sample.

## Classification Metrics

Patient-level cancer-type classification:

| method | balanced accuracy | macro-F1 | macro-AUROC |
|---|---:|---:|---:|
| `topic_nmf_z` | 0.8100 | 0.8021 | 0.8624 |
| `marker_scoring_z` | 0.7074 | 0.6131 | 0.7803 |
| `fibrodynmix_z` | 0.6314 | 0.7585 | 0.8141 |
| `marker_scoring_scores` | 0.3919 | 0.5218 | 0.6875 |

Sample-level cancer-type classification:

| method | balanced accuracy | macro-F1 | macro-AUROC |
|---|---:|---:|---:|
| `topic_nmf_z` | 0.8084 | 0.8036 | 0.8390 |
| `marker_scoring_z` | 0.6762 | 0.5682 | 0.7704 |
| `fibrodynmix_z` | 0.6106 | 0.7641 | 0.8154 |
| `marker_scoring_scores` | 0.4172 | 0.5758 | 0.6757 |

These results should not be described as a universal downstream win for
FibroDynMix. NMF topics are the strongest cancer-type separator in this
bounded task. FibroDynMix `z` is nevertheless useful as a representation:
it exceeds raw marker-score features on macro-F1 and macro-AUROC and remains
competitive with normalized marker-score composition.

## Biological Gradient Validation

The same script evaluates whether canonical fibroblast marker programs vary
along the expected FibroDynMix state weights.

Expected-state correlations:

| marker gradient | expected state | genes | Spearman rho |
|---|---|---:|---:|
| `ecm_gradient` | `ECM-remodeling` | 8 | 0.6714 |
| `myofibroblast_gradient` | `myofibroblast` | 7 | 0.4236 |
| `inflammatory_gradient` | `inflammatory` | 6 | 0.3475 |

This supports a biological-gradient sanity check: FibroDynMix `z` recovers
continuous variation in expected fibroblast programs. Because related marker
families orient the latent states, this should be framed as gradient recovery
and representation validation, not as independent marker discovery.

## Manuscript Use

Recommended claim:

```text
In GSE246215, FibroDynMix state compositions provide biologically interpretable
downstream features: they recover expected myofibroblast, ECM-remodeling, and
inflammatory marker gradients and improve over raw marker-score features on
macro-F1 and macro-AUROC, although NMF topics remain the strongest cancer-type
classifier in this bounded task.
```

Avoid claiming:

- FibroDynMix universally outperforms marker scoring in real-data downstream
  classification.
- Cancer-type classification proves disease mechanism discovery.
- GSE246215 is an independent multi-study atlas validation.

## Claim Boundary

This is a processed public count-matrix validation from one GEO study. It is
useful as a downstream-feature and biological-gradient case study, but stronger
biology-first claims require independent human datasets from separate studies
or platforms.
