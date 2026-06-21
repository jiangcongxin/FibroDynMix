test_that("prepare_fibrodynmix_seurat extracts counts and metadata", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    c(
      5, 0, 1, 0, 3, 0,
      0, 4, 0, 2, 0, 1,
      3, 1, 4, 1, 2, 1,
      1, 3, 1, 4, 1, 3,
      2, 2, 2, 2, 2, 2,
      1, 1, 1, 1, 1, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:6), paste0("cell", 1:6))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  metadata <- data.frame(
    study = rep(c("s1", "s2"), each = 3),
    donor = rep(c("d1", "d2", "d3"), each = 2),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = metadata)
  markers <- list(resident = c("gene1", "gene3"), inflammatory = c("gene2", "gene4"))

  prepared <- prepare_fibrodynmix_seurat(
    object = object,
    marker_index = markers,
    study_col = "study",
    donor_col = "donor",
    min_cells_per_gene = 1,
    min_counts_per_gene = 1
  )

  expect_s3_class(prepared, "FibroDynMixData")
  expect_equal(dim(prepared$counts), dim(counts))
  expect_equal(colnames(prepared$counts), colnames(object))
  expect_equal(prepared$cell_metadata$study, metadata$study)
  expect_equal(prepared$study_id, metadata$study)
  expect_equal(prepared$donor_id, metadata$donor)
})

test_that("add_fibrodynmix_to_seurat writes metadata and reduction", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    rpois(40, lambda = 3),
    nrow = 8,
    dimnames = list(paste0("gene", 1:8), paste0("cell", 1:5))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  z <- matrix(
    c(0.8, 0.2, 0.6, 0.4, 0.1, 0.9),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("cell", 1:3), c("resident", "ECM-remodeling"))
  )
  fit <- list(z_hat = z)

  object <- add_fibrodynmix_to_seurat(
    object = object,
    fit = fit,
    prefix = "fdm_",
    reduction_name = "fdm",
    key = "FDM_"
  )

  meta <- object[[]]
  expect_true(all(c(
    "fdm_z_resident", "fdm_z_ECM_remodeling", "fdm_dominant_state",
    "fdm_entropy", "fdm_normalized_entropy", "fdm_max_weight",
    "fdm_confidence_class", "fdm_fpi"
  ) %in% colnames(meta)))
  expect_equal(meta["cell1", "fdm_dominant_state"], "resident")
  expect_equal(meta["cell3", "fdm_dominant_state"], "ECM-remodeling")
  expect_equal(meta["cell1", "fdm_confidence_class"], "high")
  expect_true(meta["cell1", "fdm_max_weight"] >= 0.7)
  expect_true(meta["cell1", "fdm_normalized_entropy"] >= 0)
  expect_true(meta["cell1", "fdm_normalized_entropy"] <= 1)
  expect_true(is.na(meta["cell4", "fdm_z_resident"]))
  expect_true("fdm" %in% names(object@reductions))
  expect_equal(ncol(object[["fdm"]]@cell.embeddings), 2)
  expect_equal(object[["fdm"]]@misc$state_names, colnames(z))
})

test_that("fit_fibrodynmix_seurat returns fitted Seurat workflow object", {
  skip_if_not_installed("SeuratObject")

  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 1,
    cells_per_donor = 4,
    n_genes = 42,
    marker_genes_per_state = 3,
    seed = 101
  )
  metadata <- sim$cell_metadata
  rownames(metadata) <- metadata$cell_id
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(sim$counts, sparse = TRUE), meta.data = metadata)

  result <- fit_fibrodynmix_seurat(
    object = object,
    marker_index = sim$parameters$marker_index,
    study_col = "study_id",
    attach = TRUE,
    fit_args = list(
      n_outer = 1,
      initializer_args = list(n_iter = 1),
      maxit_beta = 5,
      maxit_z = 5
    )
  )

  expect_s3_class(result, "FibroDynMixSeuratFit")
  expect_s3_class(result$data, "FibroDynMixData")
  expect_equal(nrow(result$fit$z_hat), ncol(sim$counts))
  expect_true("fibrodynmix_dominant_state" %in% colnames(result$seurat[[]]))
  expect_true("fibrodynmix" %in% names(result$seurat@reductions))
})

test_that("run_fibrodynmix_seurat_workflow fits, transfers, and evaluates", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    c(
      5, 4, 0, 0, 4, 5, 0, 0,
      0, 0, 5, 4, 0, 0, 4, 5,
      3, 3, 1, 1, 3, 3, 1, 1,
      1, 1, 3, 3, 1, 1, 3, 3,
      2, 2, 2, 2, 2, 2, 2, 2,
      1, 1, 1, 1, 1, 1, 1, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:6), paste0("cell", 1:8))
  )
  metadata <- data.frame(
    cell_type = rep("Fibroblast", 8),
    condition = rep(c("normal", "scar"), each = 4),
    sample_id = rep(c("s1", "s2"), each = 4),
    seurat_clusters = rep(c("0", "1"), each = 4),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE), meta.data = metadata)
  markers <- list(state_a = c("gene1", "gene3"), state_b = c("gene2", "gene4"))

  workflow <- run_fibrodynmix_seurat_workflow(
    object = object,
    marker_index = markers,
    cell_type_col = "cell_type",
    target_cell_type = "Fibroblast",
    condition_col = "condition",
    donor_col = "sample_id",
    max_fit_cells = 6,
    fit_args = list(n_outer = 1, initializer_args = list(n_iter = 1), maxit_beta = 4, maxit_z = 4),
    transfer_args = list(maxit_z = 4, chunk_size = 3, parallel = FALSE),
    prepare_args = list(min_cells_per_gene = 0, min_counts_per_gene = 0)
  )

  expect_s3_class(workflow, "FibroDynMixSeuratWorkflow")
  expect_equal(nrow(workflow$transfer$z_hat), ncol(object))
  expect_true("fibrodynmix_dominant_state" %in% colnames(workflow$seurat[[]]))
  expect_true(is.data.frame(workflow$evaluation$state_summary))
  expect_true(is.data.frame(workflow$transfer$cell_diagnostics))
})
