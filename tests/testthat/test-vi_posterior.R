test_that("fit_fibrodynmix_vi returns posterior state summaries", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 6,
    n_genes = 50,
    marker_genes_per_state = 3,
    seed = 91
  )

  vi <- fit_fibrodynmix_vi(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    nb_args = list(n_outer = 1, initializer_args = list(n_iter = 1), maxit_beta = 8, maxit_z = 8),
    n_draws = 8,
    n_elbo_draws = 3,
    n_vi_iter = 2,
    seed = 92
  )

  expect_s3_class(vi, "FibroDynMixVI")
  expect_equal(dim(vi$z_mean), dim(sim$z))
  expect_true(all(abs(rowSums(vi$z_mean) - 1) < 1e-8))
  expect_true(all(is.finite(vi$elbo_trace$elbo)))
  expect_true(is.finite(vi$best_elbo))
  expect_true(all(c("cell_id", "state", "mean", "lower", "upper") %in% colnames(vi$cell_summary$z)))
  expect_true(all(vi$cell_summary$z$lower <= vi$cell_summary$z$upper))
  expect_equal(length(unique(vi$cell_draws$draw)), 8)
})

test_that("fit_fibrodynmix_vi returns sample-level credible intervals", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 4,
    n_genes = 48,
    marker_genes_per_state = 3,
    seed = 93
  )
  metadata <- sim$cell_metadata
  rownames(metadata) <- colnames(sim$counts)

  vi <- fit_fibrodynmix_vi(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = metadata$library_size,
    nb_args = list(n_outer = 1, initializer_args = list(n_iter = 1), maxit_beta = 6, maxit_z = 6),
    n_draws = 6,
    n_elbo_draws = 2,
    n_vi_iter = 1,
    cell_metadata = metadata,
    sample_col = "donor_id",
    seed = 94,
    keep_draws = FALSE
  )

  expect_null(vi$cell_draws)
  expect_true(all(c("sample_id", "state", "mean", "lower", "upper") %in% colnames(vi$sample_summary)))
  expect_equal(sort(unique(vi$sample_summary$sample_id)), sort(unique(metadata$donor_id)))
})

test_that("fit_fibrodynmix_nb retains effect identifiers for posterior evaluation", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 3,
    n_genes = 45,
    marker_genes_per_state = 3,
    seed = 95
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    study_id = sim$cell_metadata$study_id,
    donor_id = sim$cell_metadata$donor_id,
    fit_study_effect = TRUE,
    fit_donor_effect = TRUE,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 5,
    maxit_z = 5
  )

  expect_equal(fit$study_id, sim$cell_metadata$study_id)
  expect_equal(fit$donor_id, sim$cell_metadata$donor_id)
})

test_that("evaluate_posterior_intervals measures interval coverage", {
  z_true <- matrix(
    c(0.7, 0.3, 0.2, 0.8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("s1", "s2"))
  )
  intervals <- data.frame(
    cell_id = c("c1", "c1", "c2", "c2"),
    state = c("s1", "s2", "s1", "s2"),
    mean = c(0.65, 0.35, 0.25, 0.75),
    lower = c(0.6, 0.2, 0.1, 0.7),
    upper = c(0.8, 0.4, 0.3, 0.9),
    stringsAsFactors = FALSE
  )

  metrics <- evaluate_posterior_intervals(z_true, intervals)

  expect_equal(metrics$interval_coverage, 1)
  expect_true(is.finite(metrics$mean_interval_width))
  expect_equal(metrics$n_interval_rows, 4)
})

test_that("calibrate_posterior_interval_scale expands intervals toward target coverage", {
  z_true <- matrix(
    c(0.7, 0.3, 0.2, 0.8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("s1", "s2"))
  )
  intervals <- data.frame(
    cell_id = c("c1", "c1", "c2", "c2"),
    state = c("s1", "s2", "s1", "s2"),
    mean = c(0.55, 0.45, 0.35, 0.65),
    lower = c(0.5, 0.4, 0.3, 0.6),
    upper = c(0.6, 0.5, 0.4, 0.7),
    stringsAsFactors = FALSE
  )

  raw <- evaluate_posterior_intervals(z_true, intervals)
  calibrated <- calibrate_posterior_interval_scale(
    z_true,
    intervals,
    scale_grid = c(1, 2, 4),
    target_coverage = 0.75
  )

  expect_true(calibrated$interval_coverage >= raw$interval_coverage)
  expect_true(calibrated$selected_scale >= 1)
  expect_true(nrow(calibrated$calibration_table) == 3)
})
