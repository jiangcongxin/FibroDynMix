test_that("fit_topic_nmf_baseline returns state-aligned simplex weights", {
  sim <- simulate_fibrodynmix(
    scenario = "continuous",
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 5,
    n_genes = 48,
    marker_genes_per_state = 3,
    seed = 101
  )

  fit <- fit_topic_nmf_baseline(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    n_iter = 8,
    seed = 102,
    backend = "multiplicative_update"
  )

  expect_equal(dim(fit$z_pred), dim(sim$z))
  expect_equal(colnames(fit$z_pred), colnames(sim$z))
  expect_true(all(rowSums(fit$z_pred) > 0.999))
  expect_true(all(rowSums(fit$z_pred) < 1.001))
  expect_true(is.finite(fit$final_objective))
  expect_equal(fit$backend, "multiplicative_update")
  expect_equal(length(fit$topic_to_state), ncol(sim$z))
})

test_that("fit_topic_nmf_baseline can use the NMF package backend when available", {
  skip_if_not_installed("NMF")
  sim <- simulate_fibrodynmix(
    scenario = "discrete",
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 8,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 103
  )

  fit <- fit_topic_nmf_baseline(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    seed = 104,
    backend = "nmf"
  )

  expect_equal(dim(fit$z_pred), dim(sim$z))
  expect_equal(fit$backend, "nmf")
  expect_true(is.finite(fit$final_objective))
})
