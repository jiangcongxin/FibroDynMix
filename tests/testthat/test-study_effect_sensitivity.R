test_that("run_study_effect_sensitivity returns baseline and grid rows", {
  result <- run_study_effect_sensitivity(
    study_l2_grid = c(0.1, 1),
    marker_l2_grid = c(0.05),
    n_replicates = 1,
    seed = 30,
    simulation_args = list(
      n_studies = 2,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 36,
      marker_genes_per_state = 3
    ),
    nb_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 6,
      maxit_z = 6
    )
  )

  expect_equal(nrow(result), 3)
  expect_equal(sum(result$method == "fibrodynmix_nb"), 1)
  expect_equal(sum(result$method == "fibrodynmix_nb_study"), 2)
  expect_true(all(is.finite(result$rmse)))
  expect_true(all(is.finite(result$nb_best_objective)))
  expect_true(all(is.na(result$study_effect_l2_norm[result$method == "fibrodynmix_nb"])))
  expect_true(all(is.finite(result$study_effect_l2_norm[result$method == "fibrodynmix_nb_study"])))
})

test_that("run_study_effect_sensitivity validates grids", {
  expect_error(run_study_effect_sensitivity(study_l2_grid = -1), "study_l2_grid")
  expect_error(run_study_effect_sensitivity(marker_l2_grid = numeric()), "marker_l2_grid")
})
