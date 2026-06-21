test_that("score_marker_baseline returns simplex predictions", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 5,
    n_genes = 80,
    marker_genes_per_state = 4,
    seed = 3
  )

  baseline <- score_marker_baseline(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size
  )

  expect_equal(dim(baseline$z_pred), dim(sim$z))
  expect_equal(colnames(baseline$z_pred), colnames(sim$z))
  expect_true(all(rowSums(baseline$z_pred) > 0.999))
  expect_true(all(rowSums(baseline$z_pred) < 1.001))
  expect_equal(dim(baseline$scores), dim(sim$z))
})

test_that("score_marker_baseline accepts gene-name markers", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 80,
    marker_genes_per_state = 4,
    seed = 4
  )
  marker_names <- lapply(
    sim$parameters$marker_index,
    function(idx) rownames(sim$counts)[idx]
  )

  baseline <- score_marker_baseline(sim$counts, marker_names)

  expect_equal(dim(baseline$z_pred), dim(sim$z))
  expect_false(anyNA(baseline$z_pred))
})

test_that("run_marker_scoring_benchmark returns metrics", {
  result <- run_marker_scoring_benchmark(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 5,
    n_genes = 80,
    marker_genes_per_state = 4,
    seed = 5
  )

  expect_true(is.list(result$simulation))
  expect_true(is.list(result$baseline))
  expect_true(is.list(result$metrics))
  expect_true(is.numeric(result$metrics$rmse))
  expect_true(result$metrics$dominant_accuracy >= 0)
  expect_true(result$metrics$dominant_accuracy <= 1)
})

test_that("score_marker_baseline rejects missing state markers", {
  counts <- matrix(1, nrow = 5, ncol = 3, dimnames = list(paste0("g", 1:5), NULL))
  marker_index <- list(state_a = "not_present")

  expect_error(score_marker_baseline(counts, marker_index), "no markers present")
})
