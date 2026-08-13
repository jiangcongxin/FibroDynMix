test_that("fixed logistic-normal priors normalize named means and covariance", {
  eta_mean <- matrix(
    c(2, 1,
      4, 3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("cell_b", "cell_a"), c("eta_b", "eta_a"))
  )
  covariance <- matrix(
    c(4, 1,
      1, 9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("eta_b", "eta_a"), c("eta_b", "eta_a"))
  )

  prior <- FibroDynMix:::prepare_logistic_normal_prior(
    eta_mean = eta_mean,
    n_cells = 2,
    n_logits = 2,
    covariance = covariance,
    cell_names = c("cell_a", "cell_b"),
    logit_names = c("eta_a", "eta_b")
  )

  expect_equal(rownames(prior$mean), c("cell_a", "cell_b"))
  expect_equal(colnames(prior$mean), c("eta_a", "eta_b"))
  expect_equal(
    prior$mean,
    matrix(c(3, 1, 4, 2), 2, 2, dimnames = list(c("cell_a", "cell_b"), c("eta_a", "eta_b")))
  )
  expect_equal(
    prior$covariance,
    matrix(c(9, 1, 1, 4), 2, 2, dimnames = list(c("eta_a", "eta_b"), c("eta_a", "eta_b")))
  )
  expect_equal(unname(prior$precision %*% prior$covariance), diag(2), tolerance = 1e-12)
  expect_equal(prior$matrix_source, "covariance")
  expect_true("prepare_logistic_normal_prior" %in% getNamespaceExports("FibroDynMix"))
})

test_that("shared vector means and precision matrices are validated", {
  precision <- matrix(
    c(5, 1,
      1, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("eta_b", "eta_a"), c("eta_b", "eta_a"))
  )
  prior <- FibroDynMix:::prepare_logistic_normal_prior(
    eta_mean = c(eta_b = 0.4, eta_a = -0.2),
    n_cells = 3,
    n_logits = 2,
    precision = precision,
    cell_names = c("c1", "c2", "c3"),
    logit_names = c("eta_a", "eta_b")
  )

  expect_equal(unname(prior$mean[, "eta_a"]), rep(-0.2, 3))
  expect_equal(unname(prior$mean[, "eta_b"]), rep(0.4, 3))
  expect_equal(unname(prior$precision), matrix(c(2, 1, 1, 5), 2, 2))
  expect_equal(prior$matrix_source, "precision")
  expect_error(
    FibroDynMix:::prepare_logistic_normal_prior(
      eta_mean = c(0, 0),
      n_cells = 1,
      n_logits = 2,
      covariance = diag(2),
      precision = diag(2)
    ),
    "exactly one"
  )
  expect_error(
    FibroDynMix:::prepare_logistic_normal_prior(
      eta_mean = c(0, 0),
      n_cells = 1,
      n_logits = 2,
      covariance = matrix(c(1, 2, 2, 1), 2, 2)
    ),
    "positive definite"
  )
})

test_that("the logistic-normal quadratic uses per-cell means and aligned logits", {
  prior <- FibroDynMix:::prepare_logistic_normal_prior(
    eta_mean = matrix(
      c(0, 1,
        2, -1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("c1", "c2"), c("eta_a", "eta_b"))
    ),
    n_cells = 2,
    n_logits = 2,
    covariance = diag(2),
    cell_names = c("c1", "c2"),
    logit_names = c("eta_a", "eta_b")
  )
  eta <- matrix(
    c(0, 0,
      2, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c2", "c1"), c("eta_b", "eta_a"))
  )

  penalty <- FibroDynMix:::logistic_normal_prior_quadratic(eta, prior)
  vector_penalty <- FibroDynMix:::logistic_normal_prior_quadratic(
    c(eta_b = 1, eta_a = 2),
    prior,
    cell_index = "c2"
  )

  expect_equal(penalty, 3.5, tolerance = 1e-12)
  expect_equal(vector_penalty, 2, tolerance = 1e-12)
  expect_error(
    FibroDynMix:::logistic_normal_prior_quadratic(c(0, 1), prior),
    "cell_index"
  )
})

test_that("empirical group priors estimate group means and ridge-stabilized covariance", {
  z <- matrix(
    c(0.60, 0.30, 0.10,
      0.50, 0.40, 0.10,
      0.20, 0.30, 0.50,
      0.25, 0.25, 0.50),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("c1", "c2", "c3", "c4"), c("s1", "s2", "s3"))
  )
  group_id <- c(c4 = "B", c3 = "B", c2 = "A", c1 = "A")

  prior <- FibroDynMix:::estimate_group_logistic_normal_prior(
    z = z,
    group_id = group_id,
    covariance_ridge = 0.25,
    reference_state = "s3"
  )

  expected_eta <- log(z[, c("s1", "s2")] / z[, "s3"])
  expected_a <- colMeans(expected_eta[c("c1", "c2"), , drop = FALSE])
  expected_b <- colMeans(expected_eta[c("c3", "c4"), , drop = FALSE])
  expect_equal(unname(prior$mean[c("c1", "c2"), ]), matrix(rep(expected_a, each = 2), 2, 2))
  expect_equal(unname(prior$mean[c("c3", "c4"), ]), matrix(rep(expected_b, each = 2), 2, 2))
  expect_equal(prior$group_sizes, c(A = 2L, B = 2L))
  expect_equal(prior$reference_state, "s3")
  expect_true(all(eigen(prior$covariance, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(unname(prior$covariance - prior$empirical_covariance), diag(0.25, 2), tolerance = 1e-12)
  expect_equal(prior$residual_df, 2)
})

test_that("initializer logit anchors are centered at the supplied state coordinate", {
  z <- matrix(
    c(0.60, 0.30, 0.10,
      0.25, 0.25, 0.50),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("s1", "s2", "s3"))
  )
  prior <- FibroDynMix:::prepare_initializer_logit_anchor(z, logit_sd = 0.2)
  eta <- FibroDynMix:::simplex_matrix_to_reference_logits(z)

  expect_equal(prior$mean, eta, tolerance = 1e-12)
  expect_equal(unname(prior$covariance), diag(0.04, 2), tolerance = 1e-12)
  expect_equal(
    FibroDynMix:::logistic_normal_prior_quadratic(eta, prior),
    0,
    tolerance = 1e-12
  )
  expect_identical(prior$prior_origin, "cell_specific_marker_initializer_logit_anchor")
  expect_error(
    FibroDynMix:::prepare_initializer_logit_anchor(z, logit_sd = 0),
    "positive"
  )
})
