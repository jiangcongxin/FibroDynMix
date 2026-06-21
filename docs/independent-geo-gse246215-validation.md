# Independent GEO GSE246215 Human Fibroblast Validation

This analysis adds an independent human fibroblast public dataset source beyond
the Dryad mouse breast fibroblast count matrices.

## Source

- GEO accession: `GSE246215`
- Public record: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246215>
- Title: Cross-tissue human fibroblast atlas reveals myofibroblast subtypes
  with distinct roles in immune modulation
- Available supplementary files:
  - `GSE246215_Fibroblast_counts.csv.gz`
  - `GSE246215_Fibroblast_metadata.csv.gz`

The GEO record provides processed fibroblast count and metadata matrices. The
GEO page states that raw files are not provided for this record, so this
analysis should be described as public count-matrix validation rather than
FASTQ-level reprocessing.

## Preparation

```bash
Rscript scripts/prepare_gse246215_fibroblast_inputs.R \
  --max-cells-per-group=80
```

The preparation script samples tumor fibroblast cells by `CancerType_short` and
writes one FibroDynMix count matrix per cancer type:

- GC: 80 cells from 7 patients
- HCC: 80 cells from 2 patients
- NSCLC: 80 cells from 13 patients
- TNBC: 80 cells from 6 patients and 8 samples

Outputs are stored under
`data/public_geo_gse246215_fibroblast_atlas/`, including:

- `gse246215_fibroblast_dataset_manifest.tsv`
- `gse246215_fibroblast_selection_summary.tsv`
- `gse246215_fibroblast_selected_cell_metadata.tsv`
- `gse246215_prepare_manifest.tsv`

## Validation Command

```bash
Rscript scripts/run_multi_public_realdata_validation.R \
  --dataset-manifest=data/public_geo_gse246215_fibroblast_atlas/gse246215_fibroblast_dataset_manifest.tsv \
  --out=analysis/independent_geo_gse246215_validation \
  --max-cells=80 \
  --max-genes=700 \
  --n-outer=2 \
  --initializer-iter=3 \
  --maxit-beta=15 \
  --maxit-z=15 \
  --transfer-maxit-z=160
```

## Current Result

The current run fits 320 human fibroblast cells across four cancer-type
datasets and 721 retained genes. All weak-prior marker states retain five
markers. The pooled NB objective improves from 3.7101 to 3.4894.

Leave-dataset-out transfer converges for all held-out cancer-type datasets:

- GC held out: `transfer_z_convergence_rate = 1`
- HCC held out: `transfer_z_convergence_rate = 1`
- NSCLC held out: `transfer_z_convergence_rate = 1`
- TNBC held out: `transfer_z_convergence_rate = 1`

The analysis produces:

- `multi_public_validation_manifest.tsv`
- `multi_public_nb_fit_diagnostics.tsv`
- `multi_public_dataset_summary.tsv`
- `multi_public_dataset_state_composition.tsv`
- `multi_public_cell_state_weights.tsv`
- `multi_public_marker_coverage.tsv`
- `multi_public_transfer_diagnostics.tsv`
- `multi_public_transfer_state_composition.tsv`

## Downstream Classification and Gradient Case Study

The validation now includes a downstream representation benchmark:
`docs/gse246215-downstream-benchmark.md`.

```bash
Rscript scripts/run_gse246215_downstream_benchmark.R
```

This evaluates cancer-type classification from patient-level fibroblast
features and compares:

- FibroDynMix state composition: `fibrodynmix_z`
- marker-scoring softmax composition: `marker_scoring_z`
- raw marker scores: `marker_scoring_scores`
- state-aligned NMF topics: `topic_nmf_z`

Patient-level current result:

- `topic_nmf_z`: balanced accuracy 0.8100, macro-F1 0.8021
- `marker_scoring_z`: balanced accuracy 0.7074, macro-F1 0.6131
- `fibrodynmix_z`: balanced accuracy 0.6314, macro-F1 0.7585
- `marker_scoring_scores`: balanced accuracy 0.3919, macro-F1 0.5218

The same script also writes
`gse246215_biological_gradient_validation.tsv`. In the current bounded run,
canonical marker programs correlate with their expected FibroDynMix states:

- myofibroblast marker gradient vs myofibroblast `z`: Spearman rho 0.4236
- ECM marker gradient vs ECM-remodeling `z`: Spearman rho 0.6714
- inflammatory marker gradient vs inflammatory `z`: Spearman rho 0.3475

These results support downstream-feature evaluation and a biological-gradient
sanity check. They do not show that FibroDynMix `z` is the best cancer-type
classifier in GSE246215; NMF topics are strongest for that specific task.

## QC Notes

The HCC subset has much larger library sizes than the other cancer-type subsets
in the sampled public count matrix. The validation therefore supports
cross-dataset execution and transfer diagnostics, but HCC-specific state
composition should not be interpreted biologically without deeper source-data
normalization and sample-level QC.

## Claim Boundary

This is an independent human public count-matrix validation using the same
FibroDynMix registry-driven workflow as the Dryad public validation. It
strengthens the real-data evidence from one public mouse source to multiple
public sources and a human multi-patient fibroblast atlas. It remains a sampled
computational validation, not a fully reprocessed FASTQ-level atlas or a
definitive biological claim about cancer-type-specific fibroblast states.
