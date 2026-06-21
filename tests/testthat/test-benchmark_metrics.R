test_that("evaluate_state_weights reports perfect recovery", {
  z <- matrix(
    c(0.8, 0.2, 0.1, 0.9, 0.4, 0.6),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("resident", "inflammatory"))
  )

  metrics <- evaluate_state_weights(z, z)

  expect_equal(metrics$rmse, 0)
  expect_equal(metrics$mean_absolute_error, 0)
  expect_equal(metrics$dominant_accuracy, 1)
  expect_equal(unname(metrics$per_state_rmse), c(0, 0))
  expect_equal(sum(metrics$dominant_confusion), nrow(z))
})

test_that("evaluate_state_weights rejects mismatched states", {
  z_true <- matrix(c(0.7, 0.3), nrow = 1, dimnames = list(NULL, c("a", "b")))
  z_pred <- matrix(c(0.7, 0.3), nrow = 1, dimnames = list(NULL, c("a", "c")))

  expect_error(evaluate_state_weights(z_true, z_pred), "column names")
})

test_that("evaluate_marker_recovery computes state-wise AUPRC", {
  truth <- matrix(
    c(1, 1, 0, 0, 0, 0, 1, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("state_a", "state_b"), paste0("g", 1:4))
  )
  scores <- matrix(
    c(0.9, 0.8, 0.2, 0.1, 0.1, 0.2, 0.8, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(truth)
  )

  metrics <- evaluate_marker_recovery(truth, scores)

  expect_equal(metrics$macro_auprc, 1)
  expect_equal(unname(metrics$per_state_auprc), c(1, 1))
})

test_that("evaluate_transition_detection handles rare positives", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  score <- c(0.9, 0.8, 0.2, 0.1, 0.05)

  metrics <- evaluate_transition_detection(truth, score)

  expect_equal(metrics$auroc, 1)
  expect_equal(metrics$auprc, 1)
  expect_equal(metrics$precision, 1)
  expect_equal(metrics$recall, 1)
  expect_equal(metrics$f1, 1)
})

test_that("evaluate_transition_detection rejects degenerate truth", {
  expect_error(evaluate_transition_detection(c(TRUE, TRUE), c(0.9, 0.8)), "both positive and negative")
})

test_that("evaluate_downstream_classification separates simple class centroids", {
  features <- matrix(
    c(
      0.9, 0.1,
      0.8, 0.2,
      0.85, 0.15,
      0.1, 0.9,
      0.2, 0.8,
      0.15, 0.85
    ),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("state_a", "state_b"))
  )
  labels <- rep(c("normal", "disease"), each = 3)

  metrics <- evaluate_downstream_classification(
    features = features,
    labels = labels,
    n_folds = 3,
    seed = 1
  )

  expect_equal(metrics$status, "ok")
  expect_equal(metrics$accuracy, 1)
  expect_equal(metrics$balanced_accuracy, 1)
  expect_equal(metrics$macro_f1, 1)
  expect_equal(metrics$n_folds, 3)
})

test_that("evaluate_downstream_classification respects grouped folds", {
  features <- matrix(
    c(
      0.9, 0.1,
      0.8, 0.2,
      0.1, 0.9,
      0.2, 0.8,
      0.88, 0.12,
      0.82, 0.18,
      0.12, 0.88,
      0.18, 0.82
    ),
    ncol = 2,
    byrow = TRUE
  )
  labels <- rep(c("normal", "disease", "normal", "disease"), each = 2)
  groups <- rep(c("donor1", "donor2", "donor3", "donor4"), each = 2)

  metrics <- evaluate_downstream_classification(
    features = features,
    labels = labels,
    groups = groups,
    n_folds = 2,
    seed = 2
  )

  expect_equal(metrics$status, "ok")
  expect_equal(metrics$n_folds, 2)
  expect_true(metrics$balanced_accuracy >= 0.5)
  expect_true(metrics$macro_f1 >= 0.5)
})
