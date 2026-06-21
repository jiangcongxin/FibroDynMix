test_that("run_simulation_benchmark returns one row per scenario replicate method", {
  result <- run_simulation_benchmark(
    scenarios = c("continuous", "discrete"),
    n_replicates = 2,
    seed = 10,
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 2,
      cells_per_donor = 5,
      n_genes = 80,
      marker_genes_per_state = 4
    )
  )

  expect_equal(nrow(result), 4)
  expect_equal(sort(unique(result$scenario)), c("continuous", "discrete"))
  expect_equal(unique(result$method), "marker_scoring")
  expect_true(all(is.finite(result$rmse)))
  expect_true(all(result$dominant_accuracy >= 0))
  expect_true(all(result$dominant_accuracy <= 1))
  expect_true(all(is.na(result$nb_best_objective)))
  expect_true("downstream_balanced_accuracy" %in% colnames(result))
  expect_true(all(result$downstream_status %in% c("ok", "insufficient_class_replicates")))
})

test_that("run_simulation_benchmark can run the FibroDynMix initializer", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 11,
    methods = c("marker_scoring", "fibrodynmix_initializer"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 2,
      cells_per_donor = 5,
      n_genes = 80,
      marker_genes_per_state = 4
    ),
    initializer_args = list(n_iter = 2)
  )

  expect_equal(nrow(result), 2)
  expect_equal(
    sort(result$method),
    c("fibrodynmix_initializer", "marker_scoring")
  )
  expect_true(all(is.finite(result$rmse)))
  expect_true(all(is.na(result$nb_best_objective)))
})

test_that("run_simulation_benchmark can run the topic NMF baseline", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 111,
    methods = c("marker_scoring", "topic_nmf"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 42,
      marker_genes_per_state = 3
    ),
    topic_nmf_args = list(n_iter = 6, backend = "multiplicative_update")
  )

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$method), c("marker_scoring", "topic_nmf"))
  expect_true(all(is.finite(result$rmse)))
  expect_true(all(is.na(result$nb_best_objective)))
  expect_true("downstream_macro_f1" %in% colnames(result))
})

test_that("run_simulation_benchmark can run FibroDynMix NB", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 12,
    methods = c("marker_scoring", "fibrodynmix_initializer", "fibrodynmix_nb"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 36,
      marker_genes_per_state = 3
    ),
    initializer_args = list(n_iter = 1),
    nb_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 8,
      maxit_z = 8
    )
  )

  expect_equal(nrow(result), 3)
  expect_equal(
    sort(result$method),
    c("fibrodynmix_initializer", "fibrodynmix_nb", "marker_scoring")
  )
  expect_true(all(is.finite(result$rmse)))
  nb_row <- result[result$method == "fibrodynmix_nb", , drop = FALSE]
  expect_true(is.finite(nb_row$nb_initial_objective))
  expect_true(is.finite(nb_row$nb_final_objective))
  expect_true(is.finite(nb_row$nb_best_objective))
  expect_true(nb_row$nb_executed_iterations >= 1)
  expect_true(nb_row$nb_stop_reason %in% c("max_iterations", "early_stopping", "non_finite_objective"))
  expect_false(is.na(nb_row$nb_any_rollback))
})

test_that("run_simulation_benchmark can run study-effect FibroDynMix NB", {
  result <- run_simulation_benchmark(
    scenarios = "batch_confounding",
    n_replicates = 1,
    seed = 22,
    methods = c("fibrodynmix_nb", "fibrodynmix_nb_study"),
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

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$method), c("fibrodynmix_nb", "fibrodynmix_nb_study"))
  expect_true(all(is.finite(result$nb_best_objective)))
})

test_that("run_simulation_benchmark can run donor-effect FibroDynMix NB", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 23,
    methods = c("fibrodynmix_nb_donor", "fibrodynmix_nb_study_donor"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 2,
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

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$method), c("fibrodynmix_nb_donor", "fibrodynmix_nb_study_donor"))
  expect_true(all(is.finite(result$nb_best_objective)))
  expect_true(all(is.finite(result$nb_donor_effect_l2_norm)))
})

test_that("run_simulation_benchmark can run FibroDynMix VI with interval diagnostics", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 24,
    methods = c("fibrodynmix_nb", "fibrodynmix_vi"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 42,
      marker_genes_per_state = 3
    ),
    nb_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 6,
      maxit_z = 6
    ),
    vi_args = list(
      n_draws = 6,
      n_elbo_draws = 2,
      n_vi_iter = 1,
      seed = 25,
      keep_draws = FALSE
    )
  )

  expect_equal(nrow(result), 2)
  vi_row <- result[result$method == "fibrodynmix_vi", , drop = FALSE]
  nb_row <- result[result$method == "fibrodynmix_nb", , drop = FALSE]
  expect_true(is.finite(vi_row$vi_best_elbo))
  expect_true(is.finite(vi_row$vi_interval_coverage))
  expect_true(vi_row$vi_interval_coverage >= 0)
  expect_true(vi_row$vi_interval_coverage <= 1)
  expect_true(is.finite(vi_row$vi_mean_interval_width))
  expect_true(is.finite(vi_row$vi_calibrated_interval_coverage))
  expect_true(vi_row$vi_calibrated_interval_coverage >= vi_row$vi_interval_coverage)
  expect_true(is.finite(vi_row$vi_calibrated_interval_scale))
  expect_true(is.na(nb_row$vi_best_elbo))
})

test_that("run_simulation_benchmark can retain full results", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 10,
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 80,
      marker_genes_per_state = 4
    ),
    keep_results = TRUE
  )

  full_results <- attr(result, "results")
  expect_equal(length(full_results), 1)
  expect_true(is.list(full_results[[1]]$simulation))
  expect_true(is.list(full_results[[1]]$baseline))
})

test_that("summarize_benchmark_results aggregates metrics", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 2,
    seed = 20,
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 80,
      marker_genes_per_state = 4
    )
  )

  summary <- summarize_benchmark_results(result)

  expect_equal(nrow(summary), 1)
  expect_equal(summary$n_replicates, 2)
  expect_true("rmse_mean" %in% colnames(summary))
  expect_true("dominant_accuracy_sd" %in% colnames(summary))
})

test_that("summarize_optimizer_diagnostics aggregates NB diagnostics", {
  result <- run_simulation_benchmark(
    scenarios = "continuous",
    n_replicates = 1,
    seed = 21,
    methods = c("marker_scoring", "fibrodynmix_nb"),
    simulation_args = list(
      n_studies = 1,
      donors_per_study = 1,
      cells_per_donor = 4,
      n_genes = 36,
      marker_genes_per_state = 3
    ),
    nb_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 8,
      maxit_z = 8
    )
  )

  summary <- summarize_optimizer_diagnostics(result)

  expect_equal(nrow(summary), 2)
  expect_true("objective_improvement_mean" %in% colnames(summary))
  expect_true("study_effect_l2_norm_mean" %in% colnames(summary))
  expect_true("donor_effect_l2_norm_mean" %in% colnames(summary))
  expect_equal(summary$n_nb_runs[summary$method == "marker_scoring"], 0)
  expect_equal(summary$n_nb_runs[summary$method == "fibrodynmix_nb"], 1)
  expect_true(is.finite(summary$objective_improvement_mean[summary$method == "fibrodynmix_nb"]))
})

test_that("run_simulation_benchmark rejects unknown methods", {
  expect_error(
    run_simulation_benchmark(methods = "not_a_method"),
    "'arg' should be"
  )
})
