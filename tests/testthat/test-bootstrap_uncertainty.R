test_that("bootstrap_fibrodynmix returns uncertainty summaries", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 4,
    n_genes = 48,
    marker_genes_per_state = 3,
    seed = 50
  )

  boot <- bootstrap_fibrodynmix(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    cell_metadata = sim$cell_metadata,
    sample_col = "donor_id",
    method = "initializer",
    n_boot = 2,
    seed = 51,
    fit_args = list(n_iter = 1)
  )

  expect_equal(boot$n_boot, 2)
  expect_true(is.data.frame(boot$cell_draws))
  expect_true(is.data.frame(boot$sample_draws))
  expect_true(is.data.frame(boot$marker_draws))
  expect_true(all(c("mean", "lower", "upper") %in% colnames(boot$sample_summary)))
  expect_true(all(boot$sample_summary$lower <= boot$sample_summary$upper))
})

test_that("bootstrap_fibrodynmix can use NB study method", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 3,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 52
  )

  boot <- bootstrap_fibrodynmix(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    cell_metadata = sim$cell_metadata,
    sample_col = "study_id",
    method = "nb_study",
    n_boot = 1,
    seed = 53,
    fit_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 6,
      maxit_z = 6
    )
  )

  expect_equal(boot$method, "nb_study")
  expect_true(!is.null(boot$base_fit$study_effect))
})

test_that("summarize_bootstrap_uncertainty validates probabilities", {
  draws <- data.frame(replicate = 1, cell_id = "c1", state = "s1", z = 0.5, entropy = 1)
  sample <- data.frame(replicate = 1, sample_id = "a", state = "s1", composition = 0.5)
  marker <- data.frame(replicate = 1, state = "s1", gene = "g1", beta_abs = 1)

  expect_error(
    summarize_bootstrap_uncertainty(draws, sample, marker, probs = c(0.9, 0.1)),
    "probs"
  )
})
