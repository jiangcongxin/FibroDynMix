test_that("fit_fibrodynmix_nb returns coherent NB model estimates", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 15
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 8,
    maxit_z = 8
  )

  expect_equal(dim(fit$z_hat), dim(sim$z))
  expect_equal(dim(fit$beta_hat), dim(sim$parameters$beta_kg))
  expect_equal(length(fit$alpha_hat), nrow(sim$counts))
  expect_equal(length(fit$phi_hat), nrow(sim$counts))
  expect_true(all(rowSums(fit$z_hat) > 0.999))
  expect_true(all(rowSums(fit$z_hat) < 1.001))
  expect_true(all(fit$phi_hat > 0))
  expect_equal(length(fit$nb_objective_trace), 2)
  expect_true(all(is.finite(fit$nb_objective_trace)))
  expect_true(all(diff(fit$nb_objective_trace) <= 1e-8))
  expect_true("rolled_back" %in% names(fit$convergence))
  expect_equal(fit$executed_iterations, 1)
  expect_true(fit$best_iteration >= 0)
})

test_that("fit_fibrodynmix_nb produces a finite NB objective score", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 5,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 16
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    estimate_phi = FALSE,
    phi_init = rep(10, nrow(sim$counts)),
    maxit_beta = 8,
    maxit_z = 8
  )

  objective <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = fit$z_hat,
    beta = fit$beta_hat,
    alpha = fit$alpha_hat,
    phi = fit$phi_hat,
    library_size = sim$cell_metadata$library_size,
    average = TRUE
  )

  expect_true(is.finite(objective))
  expect_gt(objective, 0)
})

test_that("fit_fibrodynmix_nb validates phi_init", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 17
  )

  expect_error(
    fit_fibrodynmix_nb(
      counts = sim$counts,
      marker_index = sim$parameters$marker_index,
      library_size = sim$cell_metadata$library_size,
      n_outer = 1,
      phi_init = rep(-1, nrow(sim$counts))
    ),
    "positive"
  )
})

test_that("fit_fibrodynmix_nb supports early stopping and rollback diagnostics", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 18
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    n_outer = 3,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6,
    early_stopping = TRUE,
    patience = 1,
    rollback_to_best = TRUE
  )

  expect_lte(fit$executed_iterations, 3)
  expect_equal(length(fit$nb_objective_trace), fit$executed_iterations + 1)
  expect_true(all(diff(fit$nb_objective_trace) <= 1e-8))
  expect_true(fit$stop_reason %in% c("max_iterations", "early_stopping", "non_finite_objective"))
  expect_equal(dim(fit$marker_target), dim(fit$beta_hat))
})

test_that("fit_fibrodynmix_nb can fit study effects", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    scenario = "batch_confounding",
    seed = 22
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    study_id = sim$cell_metadata$study_id,
    fit_study_effect = TRUE,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6
  )

  expect_equal(dim(fit$study_effect), c(2, nrow(sim$counts)))
  expect_equal(rownames(fit$study_effect), sort(unique(sim$cell_metadata$study_id)))
  expect_true(all(is.finite(fit$study_effect)))
  expect_true(is.finite(fit$best_objective))
})

test_that("fit_fibrodynmix_nb can fit donor effects", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 3,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 23
  )

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    donor_id = sim$cell_metadata$donor_id,
    fit_donor_effect = TRUE,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6
  )

  expect_equal(dim(fit$donor_effect), c(3, nrow(sim$counts)))
  expect_equal(rownames(fit$donor_effect), sort(unique(sim$cell_metadata$donor_id)))
  expect_true(all(is.finite(fit$donor_effect)))
  expect_equal(length(fit$convergence$donor), fit$executed_iterations)
  expect_true(is.finite(fit$best_objective))
})

test_that("fit_fibrodynmix_nb validates donor effect identifiers", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 3,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 24
  )

  expect_error(
    fit_fibrodynmix_nb(
      counts = sim$counts,
      marker_index = sim$parameters$marker_index,
      library_size = sim$cell_metadata$library_size,
      fit_donor_effect = TRUE,
      n_outer = 1
    ),
    "donor_id"
  )
})
