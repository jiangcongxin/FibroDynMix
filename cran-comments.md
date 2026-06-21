## Test Environments

- Local macOS Tahoe 26.4.1, R 4.6.0, aarch64-apple-darwin23

## R CMD Check Results

Current local result:

```text
Status: OK
```

The standard local release gate is:

```bash
Rscript -e 'testthat::test_local()'
Rscript scripts/check_project_integrity.R
R CMD build .
R CMD check --no-manual FibroDynMix_0.0.0.9000.tar.gz
```

## Package Scope

FibroDynMix is currently a method-development package for raw-count
negative-binomial fibroblast state mixture modeling. It includes simulation,
benchmarking, NB optimization, bootstrap uncertainty, a lightweight
logistic-normal VI posterior layer, transition-flow utilities, and public
real-data smoke analyses.

## Large Generated Files

The following directories and local public raw-count files are intentionally
excluded from R package builds:

- `data/`
- `analysis/`
- `figures/`
- `MT3_CAFs_raw.txt`
- `Normal_mammary_fibroblasts_raw.txt`

These files support manuscript reproducibility and local analysis but should not
be bundled into the source package tarball.
