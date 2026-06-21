test_that("fit_fibrodynmix_transfer estimates held-out simplex weights", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 5,
    n_genes = 48,
    marker_genes_per_state = 3,
    seed = 41
  )
  studies <- sort(unique(sim$cell_metadata$study_id))
  train_cells <- sim$cell_metadata$study_id == studies[1]
  test_cells <- sim$cell_metadata$study_id == studies[2]

  fit <- fit_fibrodynmix_nb(
    counts = sim$counts[, train_cells, drop = FALSE],
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size[train_cells],
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6
  )
  transfer <- fit_fibrodynmix_transfer(
    counts = sim$counts[, test_cells, drop = FALSE],
    fit = fit,
    library_size = sim$cell_metadata$library_size[test_cells],
    maxit_z = 6
  )

  expect_equal(dim(transfer$z_hat), c(sum(test_cells), ncol(sim$z)))
  expect_true(all(rowSums(transfer$z_hat) > 0.999))
  expect_true(all(rowSums(transfer$z_hat) < 1.001))
  expect_true(is.finite(transfer$nb_objective))
  expect_true(is.finite(transfer$nb_loglik))
  expect_equal(transfer$n_genes, ncol(fit$beta_hat))
  expect_true(is.finite(transfer$z_convergence_rate))
  expect_true(transfer$z_convergence_rate >= 0)
  expect_true(transfer$z_convergence_rate <= 1)
  expect_equal(length(transfer$cell_converged), nrow(transfer$z_hat))
  expect_true(is.data.frame(transfer$cell_diagnostics))
  expect_true(all(c("cell_id", "convergence_code", "objective", "max_state_weight", "entropy") %in% colnames(transfer$cell_diagnostics)))
})

test_that("fit_fibrodynmix_transfer supports warm starts, chunks, and parallel workers", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 411
  )
  train_cells <- sim$cell_metadata$study_id == sort(unique(sim$cell_metadata$study_id))[1]
  test_cells <- !train_cells
  fit <- fit_fibrodynmix_nb(
    counts = sim$counts[, train_cells, drop = FALSE],
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size[train_cells],
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 5,
    maxit_z = 5
  )
  z_init <- sim$z[test_cells, colnames(fit$z_hat), drop = FALSE]
  rownames(z_init) <- colnames(sim$counts)[test_cells]

  transfer <- fit_fibrodynmix_transfer(
    counts = sim$counts[, test_cells, drop = FALSE],
    fit = fit,
    library_size = sim$cell_metadata$library_size[test_cells],
    z_init = z_init,
    chunk_size = 2,
    parallel = TRUE,
    n_workers = 2,
    maxit_z = 4
  )

  expect_equal(dim(transfer$z_hat), dim(z_init))
  expect_equal(transfer$chunk_size, 2)
  expect_equal(transfer$n_chunks, 2)
  expect_equal(nrow(transfer$cell_diagnostics), nrow(z_init))
  expect_true(all(is.finite(transfer$cell_diagnostics$objective)))
  expect_true(all(transfer$cell_diagnostics$max_state_weight <= 1))
})

test_that("fit_fibrodynmix_transfer requires shared genes", {
  fit <- list(
    beta_hat = matrix(0, nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), c("g1", "g2"))),
    alpha_hat = c(g1 = 0, g2 = 0),
    phi_hat = c(g1 = 1, g2 = 1)
  )
  counts <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("x1", "x2"), c("c1", "c2")))

  expect_error(
    fit_fibrodynmix_transfer(counts, fit),
    "shared"
  )
})

test_that("run_cross_cohort_transfer_benchmark returns held-out metrics", {
  result <- run_cross_cohort_transfer_benchmark(
    n_replicates = 1,
    seed = 42,
    simulation_args = list(
      n_studies = 2,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 48,
      marker_genes_per_state = 3
    ),
    train_fit_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 6,
      maxit_z = 6,
      fit_donor_effect = FALSE
    ),
    transfer_args = list(maxit_z = 6)
  )

  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$transfer_rmse)))
  expect_true(all(result$transfer_dominant_accuracy >= 0))
  expect_true(all(result$transfer_dominant_accuracy <= 1))
  expect_true(all(is.finite(result$transfer_nb_objective)))
  expect_true(all(is.finite(result$transfer_z_convergence_rate)))
  expect_true(all(result$n_shared_genes > 0))
})
