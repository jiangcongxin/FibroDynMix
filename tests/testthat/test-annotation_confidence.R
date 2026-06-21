test_that("evaluate_fibrodynmix_annotation summarizes confidence and support", {
  weights <- data.frame(
    cell_id = paste0("cell", 1:6),
    resident = c(0.9, 0.8, 0.2, 0.1, 0.45, 0.3),
    inflammatory = c(0.1, 0.2, 0.8, 0.9, 0.55, 0.7),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    cell_id = weights$cell_id,
    cluster = c("0", "0", "1", "1", "1", "1"),
    condition = c("normal", "normal", "disease", "disease", "disease", "normal"),
    stringsAsFactors = FALSE
  )
  expression <- matrix(
    c(
      5, 4, 0, 0, 1, 1,
      3, 4, 0, 0, 1, 1,
      0, 0, 5, 4, 3, 4,
      0, 0, 3, 4, 2, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("DCN", "LUM", "IL6", "CCL2"), weights$cell_id)
  )
  markers <- list(resident = c("DCN", "LUM"), inflammatory = c("IL6", "CCL2"))

  evaluation <- evaluate_fibrodynmix_annotation(
    weights,
    metadata = metadata,
    metadata_cell_col = "cell_id",
    cluster_col = "cluster",
    condition_col = "condition",
    expression = expression,
    marker_index = markers,
    min_state_cells = 1,
    min_high_confidence_fraction = 0.1
  )

  expect_s3_class(evaluation, "FibroDynMixAnnotationEvaluation")
  expect_equal(nrow(evaluation$cell_summary), 6)
  expect_true(all(c("max_state_weight", "normalized_entropy", "confidence_class") %in% colnames(evaluation$cell_summary)))
  expect_equal(sum(evaluation$state_summary$n_cells), 6)
  expect_true(all(c("support_label", "formal_ready", "support_reason") %in% colnames(evaluation$state_summary)))
  expect_true(all(evaluation$state_summary$formal_ready))
  expect_true(all(c("cluster", "state", "fraction_within_cluster") %in% colnames(evaluation$cluster_agreement)))
  expect_true(all(c("condition", "state", "composition", "dominant_fraction") %in% colnames(evaluation$condition_composition)))
  expect_true(all(evaluation$marker_support$log2_ratio_own_vs_other > 0))
})

test_that("evaluate_fibrodynmix_annotation keeps sparse marker support sparse-friendly", {
  skip_if_not_installed("Matrix")
  weights <- data.frame(
    cell_id = paste0("cell", 1:4),
    resident = c(0.9, 0.8, 0.2, 0.1),
    inflammatory = c(0.1, 0.2, 0.8, 0.9),
    stringsAsFactors = FALSE
  )
  expression <- Matrix::Matrix(
    matrix(
      c(5, 4, 0, 0, 0, 1, 4, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("DCN", "IL6"), weights$cell_id)
    ),
    sparse = TRUE
  )

  evaluation <- evaluate_fibrodynmix_annotation(
    weights,
    expression = expression,
    marker_index = list(resident = "DCN", inflammatory = "IL6")
  )

  expect_true(all(is.finite(evaluation$marker_support$mean_score_own_state)))
})

test_that("evaluate_fibrodynmix_annotation flags weak state support", {
  weights <- data.frame(
    cell_id = paste0("cell", 1:5),
    resident = c(0.8, 0.75, 0.7, 0.6, 0.2),
    inflammatory = c(0.2, 0.25, 0.3, 0.4, 0.8),
    stringsAsFactors = FALSE
  )
  evaluation <- evaluate_fibrodynmix_annotation(
    weights,
    min_state_cells = 3,
    min_high_confidence_fraction = 0.5
  )

  resident <- evaluation$state_summary[evaluation$state_summary$state == "resident", ]
  inflammatory <- evaluation$state_summary[evaluation$state_summary$state == "inflammatory", ]
  expect_true(resident$formal_ready)
  expect_false(inflammatory$formal_ready)
  expect_equal(inflammatory$support_label, "exploratory")
  expect_match(inflammatory$support_reason, "low_cell_count")
})

test_that("evaluate_fibrodynmix_annotation validates marker inputs", {
  weights <- data.frame(cell_id = c("c1", "c2"), a = c(0.8, 0.2), b = c(0.2, 0.8))

  expect_error(
    evaluate_fibrodynmix_annotation(weights, expression = matrix(1, nrow = 1, ncol = 2)),
    "must be supplied together"
  )
  expect_error(
    evaluate_fibrodynmix_annotation(weights, high_threshold = 0.2, moderate_threshold = 0.4),
    "moderate_threshold"
  )
  expect_error(
    evaluate_fibrodynmix_annotation(weights, min_high_confidence_fraction = 1.1),
    "min_high_confidence_fraction"
  )
})
