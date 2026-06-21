test_that("fibrodynmix_nb_loglik matches stats::dnbinom for a small model", {
  counts <- matrix(
    c(0, 2, 1, 3, 4, 0),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:2))
  )
  z <- matrix(
    c(0.7, 0.3, 0.2, 0.8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("s1", "s2"))
  )
  beta <- matrix(
    c(0.1, -0.2, 0.3, -0.1, 0.2, 0.05),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("s1", "s2"), rownames(counts))
  )
  alpha <- c(-4, -3.8, -4.2)
  phi <- c(5, 6, 7)
  library_size <- c(1000, 1200)

  eta <- z %*% beta
  eta <- sweep(eta, 2, alpha, "+")
  mu <- t(sweep(exp(eta), 1, library_size, "*"))
  expected <- sum(stats::dnbinom(
    as.vector(counts),
    size = rep(phi, times = ncol(counts)),
    mu = as.vector(mu),
    log = TRUE
  ))

  observed <- fibrodynmix_nb_loglik(counts, z, beta, alpha, phi, library_size)
  expect_equal(observed, expected)

  loglik_matrix <- fibrodynmix_nb_loglik(
    counts, z, beta, alpha, phi, library_size,
    return_matrix = TRUE
  )
  expect_equal(dim(loglik_matrix), dim(counts))
  expect_equal(sum(loglik_matrix), expected)
})

test_that("true simulator parameters improve NB log-likelihood over shuffled z", {
  sim <- simulate_fibrodynmix(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = 8,
    n_genes = 100,
    marker_genes_per_state = 5,
    seed = 12
  )
  true_loglik <- fibrodynmix_nb_loglik(
    counts = sim$counts,
    z = sim$z,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    study_effect = sim$parameters$study_effect,
    donor_effect = sim$parameters$donor_effect,
    study_id = sim$cell_metadata$study_id,
    donor_id = sim$cell_metadata$donor_id
  )

  shuffled_z <- sim$z[sample(seq_len(nrow(sim$z))), , drop = FALSE]
  rownames(shuffled_z) <- rownames(sim$z)
  shuffled_loglik <- fibrodynmix_nb_loglik(
    counts = sim$counts,
    z = shuffled_z,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    study_effect = sim$parameters$study_effect,
    donor_effect = sim$parameters$donor_effect,
    study_id = sim$cell_metadata$study_id,
    donor_id = sim$cell_metadata$donor_id
  )

  expect_gt(true_loglik, shuffled_loglik)
})

test_that("NB deviance and objective are finite and non-negative", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 5,
    n_genes = 90,
    marker_genes_per_state = 4,
    seed = 13
  )
  deviance <- fibrodynmix_nb_deviance(
    counts = sim$counts,
    z = sim$z,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    study_effect = sim$parameters$study_effect,
    donor_effect = sim$parameters$donor_effect,
    study_id = sim$cell_metadata$study_id,
    donor_id = sim$cell_metadata$donor_id
  )
  objective <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = sim$z,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    beta_l2 = 0.01,
    average = TRUE
  )

  expect_true(is.finite(deviance))
  expect_gte(deviance, 0)
  expect_true(is.finite(objective))
  expect_gt(objective, 0)
})

test_that("NB likelihood can score initializer output", {
  sim <- simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 2,
    cells_per_donor = 5,
    n_genes = 90,
    marker_genes_per_state = 4,
    seed = 14
  )
  fit <- fit_fibrodynmix_initializer(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    n_iter = 2
  )
  alpha_nb <- log((rowMeans(sim$counts) + 0.1) / mean(sim$cell_metadata$library_size))
  phi <- rep(10, nrow(sim$counts))
  names(phi) <- rownames(sim$counts)

  objective <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = fit$z_hat,
    beta = fit$beta_hat,
    alpha = alpha_nb,
    phi = phi,
    library_size = sim$cell_metadata$library_size,
    average = TRUE
  )

  expect_true(is.finite(objective))
  expect_gt(objective, 0)
})

test_that("NB likelihood validates simplex rows", {
  counts <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  z <- matrix(c(0.8, 0.8, 0.5, 0.5), nrow = 2, byrow = TRUE, dimnames = list(c("c1", "c2"), c("s1", "s2")))
  beta <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), c("g1", "g2")))

  expect_error(
    fibrodynmix_nb_loglik(
      counts = counts,
      z = z,
      beta = beta,
      alpha = c(0, 0),
      phi = c(1, 1),
      library_size = c(1, 1)
    ),
    "Rows of `z`"
  )
})

test_that("NB objective includes marker orientation penalty", {
  counts <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  z <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE, dimnames = list(c("c1", "c2"), c("s1", "s2")))
  beta <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), c("g1", "g2")))
  marker_target <- matrix(1, nrow = 2, ncol = 2, dimnames = dimnames(beta))

  unpenalized <- fibrodynmix_nb_objective(
    counts = counts,
    z = z,
    beta = beta,
    alpha = c(0, 0),
    phi = c(1, 1),
    library_size = c(1, 1)
  )
  penalized <- fibrodynmix_nb_objective(
    counts = counts,
    z = z,
    beta = beta,
    alpha = c(0, 0),
    phi = c(1, 1),
    library_size = c(1, 1),
    marker_target = marker_target,
    marker_l2 = 0.5
  )

  expect_gt(penalized, unpenalized)
})

test_that("NB objective supports donor-specific effect penalty", {
  counts <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  z <- matrix(
    c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("s1", "s2"))
  )
  beta <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), rownames(counts)))
  donor_effect <- matrix(
    0.2,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("d1", "d2"), rownames(counts))
  )
  donor_id <- c("d1", "d1", "d2")

  low_penalty <- fibrodynmix_nb_objective(
    counts = counts,
    z = z,
    beta = beta,
    alpha = c(0, 0),
    phi = c(1, 1),
    library_size = c(1, 1, 1),
    donor_effect = donor_effect,
    donor_id = donor_id,
    effect_l2 = 0,
    donor_effect_l2 = 0.01
  )
  high_penalty <- fibrodynmix_nb_objective(
    counts = counts,
    z = z,
    beta = beta,
    alpha = c(0, 0),
    phi = c(1, 1),
    library_size = c(1, 1, 1),
    donor_effect = donor_effect,
    donor_id = donor_id,
    effect_l2 = 0,
    donor_effect_l2 = 1
  )

  expect_gt(high_penalty, low_penalty)
})
