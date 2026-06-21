# FibroDynMix Workflow

This workflow is a lightweight package-facing guide. It avoids executing large
analyses during package checks; reproducible manuscript analyses are run through
the scripts under `scripts/`.

## 1. Simulate Raw Counts

```r
sim <- simulate_fibrodynmix(
  n_studies = 2,
  donors_per_study = 3,
  cells_per_donor = 20,
  n_genes = 200,
  marker_genes_per_state = 8,
  scenario = "batch_confounding",
  seed = 1
)
```

## 2. Fit the NB Model

```r
fit <- fit_fibrodynmix_nb(
  counts = sim$counts,
  marker_index = sim$parameters$marker_index,
  library_size = sim$cell_metadata$library_size,
  study_id = sim$cell_metadata$study_id,
  donor_id = sim$cell_metadata$donor_id,
  fit_study_effect = TRUE,
  fit_donor_effect = TRUE,
  n_outer = 2
)
```

## 3. Estimate Posterior State Uncertainty

```r
vi <- fit_fibrodynmix_vi(
  counts = sim$counts,
  marker_index = sim$parameters$marker_index,
  library_size = sim$cell_metadata$library_size,
  nb_args = list(n_outer = 2),
  n_draws = 50
)
```

## 4. Benchmark Against Known Truth

```r
benchmark <- run_simulation_benchmark(
  scenarios = c("continuous", "batch_confounding"),
  methods = c("marker_scoring", "fibrodynmix_nb", "fibrodynmix_vi"),
  n_replicates = 2
)
```

## 5. Run Manuscript Smoke Analyses

```bash
Rscript scripts/run_vi_benchmark.R
Rscript scripts/run_cross_cohort_transfer.R
Rscript scripts/run_public_realdata_transfer.R
Rscript scripts/check_project_integrity.R
```

## Claim Boundary

The package currently supports raw-count method development, simulation
calibration, public real-data smoke testing, and manuscript figure
reproducibility. It is not yet a complete full-parameter Bayesian inference
engine or a full human multi-cohort atlas validation.
