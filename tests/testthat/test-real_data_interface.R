test_that("prepare_fibrodynmix_data aligns metadata and filters raw data", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 48,
    marker_genes_per_state = 3,
    seed = 31
  )
  counts <- sim$counts
  counts[1, ] <- 0
  counts[, 1] <- 0
  metadata <- sim$cell_metadata
  metadata <- metadata[nrow(metadata):1, , drop = FALSE]

  prepared <- prepare_fibrodynmix_data(
    counts = counts,
    cell_metadata = metadata,
    marker_index = sim$parameters$marker_index,
    cell_id_col = "cell_id",
    study_col = "study_id",
    donor_col = "donor_id",
    min_cells_per_gene = 1,
    min_counts_per_gene = 1
  )

  expect_s3_class(prepared, "FibroDynMixData")
  expect_equal(ncol(prepared$counts), ncol(counts) - 1)
  expect_lt(nrow(prepared$counts), nrow(counts))
  expect_true(all(rowSums(prepared$counts > 0) >= 1))
  expect_true(all(rowSums(prepared$counts) >= 1))
  expect_equal(rownames(prepared$cell_metadata), colnames(prepared$counts))
  expect_equal(prepared$cell_metadata$cell_id, colnames(prepared$counts))
  expect_equal(length(prepared$library_size), ncol(prepared$counts))
  expect_true(all(prepared$library_size > 0))
  expect_equal(length(prepared$study_id), ncol(prepared$counts))
  expect_equal(length(prepared$donor_id), ncol(prepared$counts))
  expect_true(all(prepared$marker_summary$retained_markers > 0))
  expect_equal(prepared$filter_summary$dropped_zero_library_cells, 1)
})

test_that("prepare_fibrodynmix_data rejects missing marker states after filtering", {
  counts <- matrix(
    c(10, 0, 0, 5, 1, 1),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), paste0("c", 1:2))
  )
  metadata <- data.frame(cell_id = paste0("c", 1:2), study = "s1")
  markers <- list(state_a = "g1", state_b = "g2")

  expect_error(
    prepare_fibrodynmix_data(
      counts = counts,
      cell_metadata = metadata,
      marker_index = markers,
      cell_id_col = "cell_id",
      min_cells_per_gene = 2,
      min_counts_per_gene = 1
    ),
    "no retained markers"
  )
})

test_that("fit_fibrodynmix_prepared dispatches to NB optimizer", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 42,
    marker_genes_per_state = 3,
    scenario = "batch_confounding",
    seed = 32
  )
  prepared <- prepare_fibrodynmix_data(
    counts = sim$counts,
    cell_metadata = sim$cell_metadata,
    marker_index = sim$parameters$marker_index,
    cell_id_col = "cell_id",
    study_col = "study_id"
  )

  fit <- fit_fibrodynmix_prepared(
    prepared,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6
  )

  expect_equal(dim(fit$z_hat), c(ncol(prepared$counts), length(prepared$marker_index)))
  expect_false(is.null(fit$study_effect))
  expect_true(is.null(fit$donor_effect))
  expect_true(is.finite(fit$best_objective))
})

test_that("fit_fibrodynmix_prepared infers donor effects for finer donor hierarchy", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 4,
    n_genes = 42,
    marker_genes_per_state = 3,
    seed = 33
  )
  prepared <- prepare_fibrodynmix_data(
    counts = sim$counts,
    cell_metadata = sim$cell_metadata,
    marker_index = sim$parameters$marker_index,
    cell_id_col = "cell_id",
    study_col = "study_id",
    donor_col = "donor_id"
  )

  fit <- fit_fibrodynmix_prepared(
    prepared,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 6,
    maxit_z = 6
  )

  expect_false(is.null(fit$study_effect))
  expect_false(is.null(fit$donor_effect))
  expect_equal(nrow(fit$donor_effect), length(unique(sim$cell_metadata$donor_id)))
})

test_that("fit_fibrodynmix_prepared requires prepared data", {
  expect_error(
    fit_fibrodynmix_prepared(list()),
    "prepare_fibrodynmix_data"
  )
})
