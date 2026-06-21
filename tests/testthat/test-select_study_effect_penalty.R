test_that("select_study_effect_penalty recommends a grid row", {
  result <- run_study_effect_sensitivity(
    study_l2_grid = c(0.1, 1),
    marker_l2_grid = c(0.05, 0.1),
    n_replicates = 1,
    seed = 40,
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

  selected <- select_study_effect_penalty(result)

  expect_true(selected$recommended_study_l2 %in% c(0.1, 1))
  expect_true(selected$recommended_marker_l2 %in% c(0.05, 0.1))
  expect_true(is.character(selected$selection_reason))
  expect_true(nrow(selected$tradeoff_table) >= 1)
  expect_true("selection_score" %in% colnames(selected$tradeoff_table))
})

test_that("select_study_effect_penalty validates inputs", {
  bad <- data.frame(method = "fibrodynmix_nb")
  expect_error(select_study_effect_penalty(bad), "missing columns")

  no_study <- data.frame(
    method = "fibrodynmix_nb",
    study_l2 = NA_real_,
    marker_l2 = 0.05,
    rmse = 0.1,
    nb_best_objective = 1,
    study_effect_l2_norm = NA_real_,
    nb_any_rollback = FALSE
  )
  expect_error(select_study_effect_penalty(no_study), "fibrodynmix_nb_study")
})
