# Independent GEO GSE167339 Validation

This document records the second independent human real-data validation layer
for FibroDynMix.

## Dataset

Source: GEO `GSE167339`

Public record:
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167339`

Title: *Disrupting Mechanotransduction Promotes Regenerative Phenotypes in
Human Cells*

The GEO record describes human fibroblast scRNA-seq from three human patients
in a 3D collagen scar-system perturbation experiment. The current bounded
FibroDynMix run uses the seven explicit donor/condition MTX samples plus four
Human 3 hash-demultiplexed pseudo-samples. Human 3 hash groups are labeled
`hash_unknown` because the public workbook does not map HumanHashTag IDs to
treatment labels.

## Preparation

```bash
Rscript scripts/prepare_gse167339_human_fibroblast_inputs.R \
  --max-cells-per-sample=80 \
  --include-hash-pool=true
```

Generated inputs:

- `data/public_geo_gse167339_human_fibroblast/gse167339_prepare_manifest.tsv`
- `data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv`
- `data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_selection_summary.tsv`
- sampled sample-level count matrices as `*_human_fibroblast_counts.rds`
- `data/public_geo_gse167339_human_fibroblast/gse167339_human3_hash_demux_summary.tsv`

Preparation summary from the current bounded run:

- 11 sample-level datasets
- 3 donors
- 880 selected cells
- 33,538 genes before validation gene filtering
- 4 Human 3 hash-demultiplexed pseudo-samples, labeled `hash_unknown`

## Validation Run

```bash
Rscript scripts/run_multi_public_realdata_validation.R \
  --dataset-manifest=data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv \
  --out=analysis/independent_geo_gse167339_validation \
  --max-cells=80 \
  --max-genes=700 \
  --transfer-maxit-z=120 \
  --seed=167339
```

Primary outputs:

- `analysis/independent_geo_gse167339_validation/multi_public_validation_manifest.tsv`
- `analysis/independent_geo_gse167339_validation/multi_public_nb_fit_diagnostics.tsv`
- `analysis/independent_geo_gse167339_validation/multi_public_transfer_diagnostics.tsv`
- `analysis/independent_geo_gse167339_validation/multi_public_transition_summary.tsv`
- `analysis/independent_geo_gse167339_validation/multi_public_dataset_state_composition.tsv`

Current validation summary:

- 11 public sample-level datasets
- 3 model conditions (`normal`, `disease`, `hash_unknown`)
- 880 cells
- 720 retained genes after marker-preserving gene filtering
- pooled NB objective improved from 2.6149 to 2.3037
- mean leave-dataset-out transfer convergence rate: 0.9955
- minimum leave-dataset-out transfer convergence rate: 0.9750
- cross-sectional normal-to-disease transition flow converged

## Claim Boundary

This analysis strengthens the real-data layer by adding an independent
three-donor human fibroblast perturbation dataset outside GSE246215. It supports
computational transfer, donor/sample robustness, hash-demultiplexed Human 3
inclusion, and cross-sectional plasticity-flow feasibility.

It should not be presented as a definitive in vivo disease atlas. The run uses
bounded downsampling, maps strain/treatment conditions to a broad `disease`
label for model-level normal-to-perturbed comparison, and includes Human 3 hash
groups only as `hash_unknown` donor-level robustness evidence because the public
workbook does not map hash tags to treatment labels.
