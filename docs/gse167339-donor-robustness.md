# GSE167339 Donor Robustness

This analysis turns the GSE167339 validation from sample-level transfer into
donor-level robustness evidence.

## Inputs

Primary input:

- `data/public_geo_gse167339_human_fibroblast/gse167339_human_fibroblast_dataset_manifest.tsv`

The current manifest contains 11 sample-level datasets from three donors:

- Human1: explicit normal and perturbed samples
- Human2: explicit normal and perturbed samples
- Human3: four hash-demultiplexed `hash_unknown` pseudo-samples

Human3 is included only for donor-level robustness because the public workbook
does not map HumanHashTag IDs to treatment labels.

## Run

```bash
Rscript scripts/run_gse167339_donor_robustness.R
```

Outputs:

- `analysis/gse167339_donor_robustness/gse167339_donor_robustness_manifest.tsv`
- `analysis/gse167339_donor_robustness/gse167339_donor_state_composition.tsv`
- `analysis/gse167339_donor_robustness/gse167339_donor_robustness_summary.tsv`
- `analysis/gse167339_donor_robustness/gse167339_leave_donor_out_transfer.tsv`
- `analysis/gse167339_donor_robustness/gse167339_leave_donor_out_state_composition.tsv`
- `analysis/gse167339_donor_robustness/gse167339_hash_threshold_sensitivity.tsv`

## Current Result

The current bounded run used 880 cells and 720 retained genes.

Leave-donor-out transfer:

- Human1 held out: convergence rate 1.000
- Human2 held out: convergence rate 1.000
- Human3 held out: convergence rate 0.978
- mean convergence rate: 0.993
- minimum convergence rate: 0.978

Hash-demultiplexing sensitivity:

- 12 threshold settings were evaluated.
- All settings retained four hash groups with at least 40 assigned cells.
- The default threshold (`top count >= 10`, `top/second >= 2`) assigned 2,358
  Human3 cells before bounded sampling.

## Claim Boundary

This analysis supports donor-level computational robustness and transfer across
GSE167339 donors. It does not assign Human3 hash groups to treatment labels and
should not be used for Human3 treatment-specific biological claims.
