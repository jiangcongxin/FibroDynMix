# Public Real-Data Transfer Smoke Analysis

## Purpose

This analysis checks that FibroDynMix can train a raw-count negative-binomial
state program in one public fibroblast condition and transfer that fixed program
to held-out cells from the other condition:

```text
condition A raw counts
-> fit NB state program
-> freeze alpha, beta, phi
-> optimize held-out z under raw-count NB likelihood in condition B
```

The analysis is bidirectional: disease to normal and normal to disease.

## Command

```bash
Rscript scripts/run_public_realdata_transfer.R \
  --data-dir=data/public_dryad_breast_fibroblast \
  --out=analysis/public_realdata_transfer
```

The script also falls back to root-level `MT3_CAFs_raw.txt` and
`Normal_mammary_fibroblasts_raw.txt` when the files are not present under
`data/public_dryad_breast_fibroblast/`.

## Expected Outputs

- `public_transfer_dataset_manifest.tsv`
- `public_transfer_diagnostics.tsv`
- `public_transfer_state_composition.tsv`
- `public_transfer_cell_state_weights.tsv`
- `public_transfer_marker_coverage.tsv`
- `public_transfer_manifest.tsv`

`public_transfer_diagnostics.tsv` reports the training objective, transfer
negative-binomial objective, held-out z convergence status, convergence rate,
non-converged cell count, and transferred-state entropy for each direction.

## Claim Boundary

This is an engineering smoke test on two public mouse breast fibroblast raw-count
files. It supports the claim that the implemented transfer layer runs on real
raw counts with auditable diagnostics. It does not establish a disease mechanism
and does not replace a full multi-donor human cross-cohort atlas validation.
