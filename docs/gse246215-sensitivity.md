# GSE246215 Sensitivity and Library-Size QC

This analysis stress-tests the independent GSE246215 human fibroblast validation
without adding new figures. It evaluates whether the real-data conclusions are
stable across random downsampling and simple library-size trimming.

## Command

```bash
Rscript scripts/run_gse246215_sensitivity.R \
  --seeds=246215,246216 \
  --max-cells-per-group=60 \
  --max-genes=600
```

## Design

The script starts from the public GSE246215 processed fibroblast count and
metadata files. It runs two sensitivity modes:

- `raw_sample`: random per-cancer-type downsampling
- `library_trim_q95`: random per-cancer-type downsampling after removing cells
  above the within-cancer-type 95th library-size percentile

For each mode and seed, the script generates a dataset manifest and calls
`scripts/run_multi_public_realdata_validation.R`. Outputs are written to
`analysis/gse246215_sensitivity/`.

## Outputs

- `gse246215_library_size_summary.tsv`
- `gse246215_sensitivity_run_summary.tsv`
- `gse246215_sensitivity_selection_summary.tsv`
- `gse246215_sensitivity_state_composition.tsv`
- `gse246215_sensitivity_composition_variability.tsv`
- `gse246215_sensitivity_transfer_diagnostics.tsv`
- `gse246215_sensitivity_manifest.tsv`

## Current Result

The current run includes four validations: two random seeds and two QC modes.
Each validation uses four cancer-type datasets and 240 cells. The minimum
leave-dataset-out transfer convergence rate is 0.9833, and the mean objective
improvement is 0.2521. The maximum state-composition standard deviation across
seed repeats is 0.0420.

HCC remains a library-size outlier even after q95 trimming, with mean selected
library sizes around 1.0 million counts compared with thousands to tens of
thousands in the other sampled cancer-type subsets. Its high myofibroblast
composition is therefore a real output of the processed public matrix but should
be treated as QC-sensitive until deeper sample-level normalization or
independent HCC validation is added.

## Claim Boundary

This sensitivity analysis supports computational robustness of the GSE246215
workflow under random downsampling and simple high-library-size trimming. It
does not establish cancer-type-specific fibroblast biology. HCC-specific state
composition should be described as a QC flag and hypothesis-generating result,
not as a definitive biological finding.
