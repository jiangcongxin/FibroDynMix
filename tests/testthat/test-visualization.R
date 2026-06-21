test_that("visualization functions return ggplot objects", {
  skip_if_not_installed("ggplot2")

  composition <- data.frame(
    dataset_id = rep(c("D1", "D2"), each = 2),
    state = rep(c("resident", "inflammatory"), 2),
    composition = c(0.7, 0.3, 0.4, 0.6),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_state_composition(composition), "ggplot")

  cell_weights <- data.frame(
    cell_id = paste0("cell", 1:4),
    dominant_state = c("resident", "resident", "inflammatory", "inflammatory"),
    entropy = c(0.2, 0.4, 0.5, 0.3),
    resident = c(0.9, 0.7, 0.2, 0.1),
    inflammatory = c(0.1, 0.3, 0.8, 0.9),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_cell_state_heatmap(cell_weights), "ggplot")

  transfer <- data.frame(
    heldout_dataset_id = c("D1", "D2"),
    transfer_z_convergence_rate = c(0.94, 0.97),
    transfer_nb_objective = c(10.2, 9.8),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_transfer_diagnostics(transfer), "ggplot")

  flow <- data.frame(
    source_state = rep(c("resident", "inflammatory"), each = 2),
    target_state = rep(c("resident", "inflammatory"), 2),
    flow = c(0.5, 0.2, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_transition_flow(flow), "ggplot")

  fpi <- data.frame(
    condition = rep(c("normal", "disease"), each = 4),
    fpi = c(0.2, 0.3, 0.25, 0.35, 0.5, 0.6, 0.55, 0.65),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_fpi_distribution(fpi), "ggplot")

  rankings <- data.frame(
    stress_mode = rep(c("marker_dropout", "marker_contamination"), each = 2),
    method = rep(c("marker_scoring", "fibrodynmix_nb"), 2),
    rmse_mean = c(0.22, 0.18, 0.28, 0.21),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_benchmark_rankings(rankings), "ggplot")

  purity <- data.frame(
    dataset_id = c("D1", "D2"),
    purity_margin_mean = c(0.6, 0.4),
    low_purity_fraction = c(0.05, 0.12),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_purity_qc(purity), "ggplot")

  enrichment <- data.frame(
    state = c("resident", "inflammatory"),
    pathway = c("ECM organization", "cytokine signaling"),
    q_value = c(0.001, 0.02),
    n_overlap = c(8, 5),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_pathway_enrichment(enrichment), "ggplot")

  marker_support <- data.frame(
    state = c("resident", "inflammatory"),
    log2_ratio_own_vs_other = c(0.4, 0.8),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_marker_support(marker_support), "ggplot")

  cluster_agreement <- data.frame(
    cluster = rep(c("0", "1"), each = 2),
    state = rep(c("resident", "inflammatory"), 2),
    fraction_within_cluster = c(0.8, 0.2, 0.3, 0.7),
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_cluster_state_agreement(cluster_agreement), "ggplot")
})

test_that("visualization functions validate required columns", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_state_composition(data.frame(dataset_id = "D1", state = "resident")),
    "missing required column"
  )
  expect_error(
    plot_cell_state_heatmap(data.frame(cell_id = "cell1", resident = 1)),
    "provide `state_cols`"
  )
})

test_that("fibroblast annotation plots work with Seurat objects", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("ggplot2")

  counts <- matrix(
    c(
      5, 0, 1, 0, 3, 0,
      0, 4, 0, 2, 0, 1,
      3, 1, 4, 1, 2, 1,
      1, 3, 1, 4, 1, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:4), paste0("cell", 1:6))
  )
  metadata <- data.frame(
    fibrodynmix_dominant_state = rep(c("resident", "inflammatory"), each = 3),
    sample = rep(c("s1", "s2"), 3),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    meta.data = metadata
  )
  embedding <- matrix(
    c(seq_len(6), rev(seq_len(6))),
    ncol = 2,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  expect_s3_class(plot_fibroblast_annotation(object), "ggplot")
  expect_s3_class(
    plot_fibroblast_annotation(object, annotation_col = "fibrodynmix_dominant_state", facet_col = "sample"),
    "ggplot"
  )
  expect_s3_class(
    plot_fibroblast_marker_dot(object, features = c("gene1", "gene2")),
    "ggplot"
  )
  expect_error(
    plot_fibroblast_marker_dot(object, features = "missing_gene"),
    "None of `features`"
  )
})
