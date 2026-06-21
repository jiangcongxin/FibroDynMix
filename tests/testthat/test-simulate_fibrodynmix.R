test_that("simulate_fibrodynmix returns coherent count data and truth", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 5,
    n_genes = 120,
    marker_genes_per_state = 5,
    seed = 1
  )

  expect_equal(dim(sim$counts), c(120, 20))
  expect_equal(nrow(sim$cell_metadata), 20)
  expect_true("is_transition" %in% colnames(sim$cell_metadata))
  expect_false(any(sim$cell_metadata$is_transition))
  expect_equal(dim(sim$z), c(20, 6))
  expect_true(all(rowSums(sim$z) > 0.999))
  expect_true(all(rowSums(sim$z) < 1.001))
  expect_true(all(sim$counts >= 0))
  expect_true(all(sim$counts == round(sim$counts)))
  expect_equal(nrow(sim$gene_metadata), 120)
  expect_equal(dim(sim$parameters$beta_kg), c(6, 120))
})

test_that("discrete scenario produces one-hot state memberships", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 8,
    n_genes = 80,
    marker_genes_per_state = 4,
    scenario = "discrete",
    seed = 2
  )

  expect_true(all(sim$z %in% c(0, 1)))
  expect_true(all(rowSums(sim$z) == 1))
})

test_that("rare_transition scenario labels transition cells", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 25,
    n_genes = 80,
    marker_genes_per_state = 4,
    scenario = "rare_transition",
    seed = 9
  )

  expect_true(any(sim$cell_metadata$is_transition))
  expect_true(!all(sim$cell_metadata$is_transition))
})
