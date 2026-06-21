test_that("fit_fibrodynmix_initializer returns coherent estimates", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 6,
    n_genes = 90,
    marker_genes_per_state = 4,
    seed = 6
  )

  fit <- fit_fibrodynmix_initializer(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    n_iter = 3
  )

  expect_equal(dim(fit$z_hat), dim(sim$z))
  expect_equal(colnames(fit$z_hat), colnames(sim$z))
  expect_true(all(rowSums(fit$z_hat) > 0.999))
  expect_true(all(rowSums(fit$z_hat) < 1.001))
  expect_equal(dim(fit$beta_hat), dim(sim$parameters$beta_kg))
  expect_equal(length(fit$objective), 4)
  expect_true(all(is.finite(fit$objective)))
  expect_equal(dim(fit$prior_scores), dim(sim$z))
})

test_that("fit_fibrodynmix_initializer accepts gene-name markers", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 5,
    n_genes = 90,
    marker_genes_per_state = 4,
    seed = 7
  )
  marker_names <- lapply(
    sim$parameters$marker_index,
    function(idx) rownames(sim$counts)[idx]
  )

  fit <- fit_fibrodynmix_initializer(
    counts = sim$counts,
    marker_index = marker_names,
    n_iter = 2
  )

  expect_equal(dim(fit$z_hat), dim(sim$z))
  expect_false(anyNA(fit$z_hat))
})

test_that("fit_fibrodynmix_initializer rejects empty marker states", {
  counts <- matrix(1, nrow = 5, ncol = 3, dimnames = list(paste0("g", 1:5), NULL))
  marker_index <- list(state_a = "not_present")

  expect_error(
    fit_fibrodynmix_initializer(counts, marker_index, n_iter = 1),
    "no markers present"
  )
})
