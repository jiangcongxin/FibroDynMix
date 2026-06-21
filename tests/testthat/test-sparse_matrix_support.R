test_that("NB likelihood matches dense and sparse counts", {
  skip_if_not_installed("Matrix")
  counts_dense <- matrix(
    c(0, 2, 1, 0, 4, 0, 3, 0, 0, 1, 0, 2),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:3))
  )
  counts_sparse <- Matrix::Matrix(counts_dense, sparse = TRUE)
  z <- matrix(
    c(0.8, 0.2, 0.3, 0.7, 0.5, 0.5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(counts_dense), c("s1", "s2"))
  )
  beta <- matrix(
    c(0.1, -0.2, 0.3, 0.0, -0.1, 0.2, 0.05, 0.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("s1", "s2"), rownames(counts_dense))
  )
  alpha <- setNames(rep(-4, 4), rownames(counts_dense))
  phi <- setNames(rep(5, 4), rownames(counts_dense))
  library_size <- colSums(counts_dense) + 100

  dense_loglik <- fibrodynmix_nb_loglik(counts_dense, z, beta, alpha, phi, library_size)
  sparse_loglik <- fibrodynmix_nb_loglik(counts_sparse, z, beta, alpha, phi, library_size)
  expect_equal(sparse_loglik, dense_loglik)

  dense_objective <- fibrodynmix_nb_objective(counts_dense, z, beta, alpha, phi, library_size, average = TRUE)
  sparse_objective <- fibrodynmix_nb_objective(counts_sparse, z, beta, alpha, phi, library_size, average = TRUE)
  expect_equal(sparse_objective, dense_objective)
})

test_that("prepare data, marker scoring, and transfer accept sparse counts", {
  skip_if_not_installed("Matrix")
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 36,
    marker_genes_per_state = 3,
    seed = 601
  )
  sparse_counts <- Matrix::Matrix(sim$counts, sparse = TRUE)
  prepared <- prepare_fibrodynmix_data(
    sparse_counts,
    sim$cell_metadata,
    sim$parameters$marker_index,
    cell_id_col = "cell_id",
    study_col = "study_id",
    min_cells_per_gene = 0,
    min_counts_per_gene = 0
  )
  expect_s3_class(prepared, "FibroDynMixData")
  expect_s4_class(prepared$counts, "dgCMatrix")

  baseline <- score_marker_baseline(prepared$counts, prepared$marker_index, prepared$library_size)
  expect_equal(nrow(baseline$z_pred), ncol(prepared$counts))

  train <- prepared$study_id == unique(prepared$study_id)[1]
  fit <- fit_fibrodynmix_prepared(
    subset_fibrodynmix_prepared(prepared, colnames(prepared$counts)[train]),
    fit_study_effect = FALSE,
    n_outer = 1,
    initializer_args = list(n_iter = 1),
    maxit_beta = 4,
    maxit_z = 4
  )
  transfer <- fit_fibrodynmix_transfer(
    counts = prepared$counts[, !train, drop = FALSE],
    fit = fit,
    library_size = prepared$library_size[!train],
    maxit_z = 4
  )
  expect_equal(nrow(transfer$z_hat), sum(!train))
  expect_true(is.data.frame(transfer$cell_diagnostics))
})
