test_that("select_fibrodynmix_nb_model returns a validation-aware choice", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 101
  )

  result <- select_fibrodynmix_nb_model(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    candidate_n_outer = c(1, 2),
    labels = sim$cell_metadata$disease,
    groups = sim$cell_metadata$donor_id,
    seed = 101,
    nb_args = list(
      initializer_args = list(n_iter = 1),
      maxit_beta = 5,
      maxit_z = 5
    ),
    transfer_args = list(maxit_z = 5),
    truth_z = sim$z
  )

  expect_true(result$selected_n_outer %in% c(1, 2))
  expect_equal(nrow(result$candidate_scores), 2)
  expect_true(all(c(
    "heldout_nb_objective",
    "z_stability_delta_vs_previous",
    "marker_gradient_mean_spearman",
    "downstream_balanced_accuracy",
    "selection_score",
    "truth_rmse"
  ) %in% colnames(result$candidate_scores)))
  expect_true(all(is.finite(result$candidate_scores$heldout_nb_objective)))
  expect_true(all(is.finite(result$candidate_scores$selection_score)))
  expect_equal(nrow(result$split), ncol(sim$counts))
  expect_true(any(result$split$split == "holdout"))
  expect_true(any(result$split$split == "train"))
})
