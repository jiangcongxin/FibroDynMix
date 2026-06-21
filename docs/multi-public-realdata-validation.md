# Multi-Public Real-Data Validation

This analysis upgrades the public real-data layer from a single pooled smoke
test to a registry-driven validation workflow across multiple public
fibroblast count datasets.

## Command

```bash
Rscript scripts/run_multi_public_realdata_validation.R \
  --max-cells=160 \
  --max-genes=700 \
  --n-outer=2 \
  --initializer-iter=3 \
  --maxit-beta=15 \
  --maxit-z=15 \
  --transfer-maxit-z=160
```

The default registry uses two public mouse breast fibroblast raw-count matrices
from Dryad doi:10.6071/M3238R:

- `MT3_CAFs_raw.txt`
- `Normal_mammary_fibroblasts_raw.txt`

The script can also read a user-supplied manifest through
`--dataset-manifest=PATH`. Required columns are `dataset_id`, `study_id`,
`donor_id`, `condition`, and `counts_path`. Optional columns such as
`organism`, `tissue`, `public_record`, and `download_url` are preserved in the
registry output.

## Outputs

The analysis writes to `analysis/multi_public_realdata_validation/`:

- `multi_public_dataset_registry.tsv`
- `multi_public_nb_fit_diagnostics.tsv`
- `multi_public_dataset_summary.tsv`
- `multi_public_dataset_state_composition.tsv`
- `multi_public_cell_state_weights.tsv`
- `multi_public_marker_coverage.tsv`
- `multi_public_transition_flow.tsv`
- `multi_public_transition_summary.tsv`
- `multi_public_transfer_diagnostics.tsv`
- `multi_public_transfer_state_composition.tsv`
- `multi_public_validation_manifest.tsv`

## Current Result

The current run fits a pooled raw-count negative-binomial model with study
effects across 320 cells and 719 retained genes. The pooled objective improves
from 2.2291 to 2.1045. Leave-dataset-out transfer converges in both directions:

- held-out MT3 CAF cells: `transfer_z_convergence_rate = 1`
- held-out normal mammary fibroblast cells: `transfer_z_convergence_rate = 1`

The dataset-level state composition separates the disease and normal public
fibroblast count matrices in the expected direction: the MT3 CAF dataset has
higher ECM-remodeling/myofibroblast weight, whereas the normal mammary
fibroblast dataset has higher resident/inflammatory mixture weight under the
current weak-prior state basis.

## Claim Boundary

This is real raw-count validation across multiple public count matrices, with
pooled fitting, study-effect adjustment, state composition summaries, transition
flow, and leave-dataset-out transfer diagnostics. The current default registry
still comes from one public Dryad record and two condition-specific matrices.
It should therefore be presented as a strengthened public real-data validation
and transfer stress test, not as an independent multi-study human atlas
validation.
