# Cross-Cohort Transfer

## Purpose

The transfer layer turns cross-cohort migration from a manuscript claim into an
executable benchmark:

```text
fit FibroDynMix on training studies
-> freeze alpha, beta, phi
-> optimize held-out z under raw-count NB likelihood
-> evaluate z recovery against simulated truth
```

## Current API

```r
transfer <- fit_fibrodynmix_transfer(
  counts = heldout_counts,
  fit = train_fit,
  library_size = heldout_library_size
)

benchmark <- run_cross_cohort_transfer_benchmark()
```

For public raw-count smoke testing, run:

```bash
Rscript scripts/run_public_realdata_transfer.R --out=analysis/public_realdata_transfer
```

## Claim Boundary

This is a leave-study-out simulation transfer benchmark. It evaluates whether
learned state programs can be applied to held-out count data. The public
real-data transfer smoke extends that check to two public mouse breast
fibroblast count files. Neither analysis is a full real-human cross-cohort atlas
validation.
