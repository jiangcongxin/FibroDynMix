nb_repair_api_ready <- function() {
  all(c("initial_state", "beta_constraint") %in% names(formals(fit_fibrodynmix_nb))) &&
    all(c("z_l2", "z_prior") %in% names(formals(fibrodynmix_nb_objective))) &&
    all(c("z_prior", "z_anchor", "z_anchor_sd") %in% names(formals(fit_fibrodynmix_nb)))
}

skip_if_nb_repair_api_is_unavailable <- function() {
  if (!nb_repair_api_ready()) {
    skip("Requires the repaired NB start-state, beta-contrast, and full-objective API.")
  }
}

make_nb_repair_fixture <- function(seed = 991) {
  simulate_fibrodynmix(
    n_studies = 1,
    donors_per_study = 1,
    cells_per_donor = 5,
    n_genes = 30,
    marker_genes_per_state = 3,
    seed = seed
  )
}

fit_nb_from_fixed_state <- function(sim, z_l2 = 0.37) {
  fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    initial_state = list(
      z_hat = sim$z,
      beta_hat = sim$parameters$beta_kg
    ),
    beta_constraint = "sum_to_zero",
    n_outer = 1,
    estimate_phi = FALSE,
    phi_init = sim$parameters$phi_g,
    beta_l2 = 0.11,
    marker_l2 = 0.23,
    z_l2 = z_l2,
    z_anchor = "none",
    maxit_beta = 3,
    maxit_z = 3,
    early_stopping = FALSE,
    rollback_to_best = TRUE
  )
}

test_that("NB objective includes the z-logit penalty used by the z update", {
  if (!"z_l2" %in% names(formals(fibrodynmix_nb_objective))) {
    skip("Requires `z_l2` in `fibrodynmix_nb_objective()`.")
  }

  counts <- matrix(
    c(1, 0, 2, 3, 1, 0),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
  )
  z <- matrix(
    c(0.2, 0.3, 0.5,
      0.6, 0.1, 0.3,
      0.4, 0.4, 0.2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("s1", "s2", "s3"))
  )
  beta <- matrix(
    c(0.2, -0.1,
      -0.3, 0.4,
      0.1, 0.2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(z), rownames(counts))
  )
  alpha <- c(g1 = -3.5, g2 = -3.2)
  phi <- c(g1 = 5, g2 = 6)
  library_size <- c(c1 = 100, c2 = 120, c3 = 110)
  z_l2 <- 0.7

  base <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    z_l2 = 0,
    average = FALSE
  )
  penalized <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    z_l2 = z_l2,
    average = FALSE
  )
  logits <- log(sweep(z[, -ncol(z), drop = FALSE], 1L, z[, ncol(z)], "/"))
  expected_penalty <- z_l2 * sum(logits^2)

  expect_equal(penalized - base, expected_penalty, tolerance = 1e-10)

  base_average <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    z_l2 = 0,
    average = TRUE
  )
  penalized_average <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    z_l2 = z_l2,
    average = TRUE
  )
  expect_equal(
    penalized_average - base_average,
    expected_penalty / length(counts),
    tolerance = 1e-10
  )
})

test_that("NB objective includes an aligned fixed logistic-normal prior", {
  counts <- matrix(
    c(1, 0, 2, 3, 1, 0),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
  )
  z <- matrix(
    c(0.2, 0.3, 0.5,
      0.6, 0.1, 0.3,
      0.4, 0.4, 0.2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("s1", "s2", "s3"))
  )
  beta <- matrix(
    c(0.2, -0.1,
      -0.3, 0.4,
      0.1, 0.2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(z), rownames(counts))
  )
  alpha <- c(g1 = -3.5, g2 = -3.2)
  phi <- c(g1 = 5, g2 = 6)
  library_size <- c(c1 = 100, c2 = 120, c3 = 110)
  prior <- prepare_logistic_normal_prior(
    eta_mean = matrix(
      0,
      nrow = nrow(z),
      ncol = ncol(z) - 1L,
      dimnames = list(rownames(z), colnames(z)[-ncol(z)])
    ),
    n_cells = nrow(z),
    n_logits = ncol(z) - 1L,
    covariance = diag(ncol(z) - 1L),
    cell_names = rownames(z),
    logit_names = colnames(z)[-ncol(z)]
  )

  unpenalized <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    average = FALSE
  )
  penalized <- fibrodynmix_nb_objective(
    counts, z, beta, alpha, phi, library_size,
    z_prior = prior,
    average = FALSE
  )
  expected <- logistic_normal_prior_quadratic(
    simplex_matrix_to_reference_logits(z),
    prior
  )

  expect_equal(penalized - unpenalized, expected, tolerance = 1e-10)
  expect_error(
    fibrodynmix_nb_objective(
      counts, z, beta, alpha, phi, library_size,
      z_prior = prepare_logistic_normal_prior(
        eta_mean = c(0, 0),
        n_cells = nrow(z),
        n_logits = ncol(z) - 1L,
        covariance = diag(ncol(z) - 1L),
        cell_names = rownames(z),
        logit_names = c("wrong_1", "wrong_2")
      )
    ),
    "incompatible reference-logit names"
  )
})

test_that("sum-to-zero fitting returns beta and marker targets in one contrast space", {
  skip_if_nb_repair_api_is_unavailable()
  sim <- make_nb_repair_fixture()
  fit <- fit_nb_from_fixed_state(sim)

  expect_equal(unname(colMeans(fit$beta_hat)), rep(0, ncol(fit$beta_hat)), tolerance = 1e-8)
  expect_equal(unname(colMeans(fit$marker_target)), rep(0, ncol(fit$marker_target)), tolerance = 1e-12)
})

test_that("a provided NB start state bypasses randomized initialization", {
  skip_if_nb_repair_api_is_unavailable()
  sim <- make_nb_repair_fixture()
  supplied_state <- list(z_hat = sim$z, beta_hat = sim$parameters$beta_kg)

  set.seed(1)
  fit_a <- fit_nb_from_fixed_state(sim)
  set.seed(999)
  fit_b <- fit_nb_from_fixed_state(sim)

  expect_equal(fit_a$initializer$z_hat, supplied_state$z_hat)
  expect_equal(fit_a$initializer$beta_hat, supplied_state$beta_hat)
  expect_equal(fit_a$nb_objective_trace, fit_b$nb_objective_trace, tolerance = 1e-10)
  expect_equal(fit_a$z_hat, fit_b$z_hat, tolerance = 1e-10)
  expect_equal(fit_a$beta_hat, fit_b$beta_hat, tolerance = 1e-10)

  expect_error(
    fit_fibrodynmix_nb(
      counts = sim$counts,
      marker_index = sim$parameters$marker_index,
      library_size = sim$cell_metadata$library_size,
      initial_state = list(z_hat = sim$z),
      n_outer = 1
    ),
    "beta_hat"
  )
})

test_that("the best reported NB objective equals the returned repaired fit", {
  skip_if_nb_repair_api_is_unavailable()
  sim <- make_nb_repair_fixture()
  fit <- fit_nb_from_fixed_state(sim, z_l2 = 0.43)

  recomputed <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = fit$z_hat,
    beta = fit$beta_hat,
    alpha = fit$alpha_hat,
    phi = fit$phi_hat,
    library_size = sim$cell_metadata$library_size,
    beta_l2 = 0.11,
    marker_target = fit$marker_target,
    marker_l2 = 0.23,
    z_l2 = 0.43,
    average = TRUE
  )

  expect_equal(fit$best_objective, recomputed, tolerance = 1e-8)
})

test_that("a fitted logistic-normal prior is used in the returned NB objective", {
  skip_if_nb_repair_api_is_unavailable()
  sim <- make_nb_repair_fixture()
  prior <- prepare_logistic_normal_prior(
    eta_mean = matrix(
      0,
      nrow = nrow(sim$z),
      ncol = ncol(sim$z) - 1L,
      dimnames = list(rownames(sim$z), colnames(sim$z)[-ncol(sim$z)])
    ),
    n_cells = nrow(sim$z),
    n_logits = ncol(sim$z) - 1L,
    covariance = diag(ncol(sim$z) - 1L),
    cell_names = rownames(sim$z),
    logit_names = colnames(sim$z)[-ncol(sim$z)]
  )
  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    initial_state = list(z_hat = sim$z, beta_hat = sim$parameters$beta_kg),
    beta_constraint = "sum_to_zero",
    n_outer = 1,
    estimate_phi = FALSE,
    phi_init = sim$parameters$phi_g,
    beta_l2 = 0.11,
    marker_l2 = 0.23,
    z_l2 = 0,
    z_prior = prior,
    maxit_beta = 3,
    maxit_z = 3,
    early_stopping = FALSE,
    rollback_to_best = TRUE
  )
  recomputed <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = fit$z_hat,
    beta = fit$beta_hat,
    alpha = fit$alpha_hat,
    phi = fit$phi_hat,
    library_size = sim$cell_metadata$library_size,
    beta_l2 = 0.11,
    marker_target = fit$marker_target,
    marker_l2 = 0.23,
    z_l2 = 0,
    z_prior = prior,
    average = TRUE
  )

  expect_s3_class(fit$z_prior, "fibrodynmix_logistic_normal_prior")
  expect_equal(fit$best_objective, recomputed, tolerance = 1e-8)
})

test_that("the default NB fit retains a marker-logit coordinate anchor", {
  skip_if_nb_repair_api_is_unavailable()
  sim <- make_nb_repair_fixture()
  fit <- fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    initial_state = list(z_hat = sim$z, beta_hat = sim$parameters$beta_kg),
    beta_constraint = "sum_to_zero",
    n_outer = 1,
    estimate_phi = FALSE,
    phi_init = sim$parameters$phi_g,
    beta_l2 = 0.11,
    marker_l2 = 0.23,
    z_l2 = 0,
    maxit_beta = 3,
    maxit_z = 3,
    early_stopping = FALSE,
    rollback_to_best = TRUE
  )
  initial_penalty <- logistic_normal_prior_quadratic(
    simplex_matrix_to_reference_logits(fit$initializer$z_hat),
    fit$z_prior
  )
  recomputed <- fibrodynmix_nb_objective(
    counts = sim$counts,
    z = fit$z_hat,
    beta = fit$beta_hat,
    alpha = fit$alpha_hat,
    phi = fit$phi_hat,
    library_size = sim$cell_metadata$library_size,
    beta_l2 = 0.11,
    marker_target = fit$marker_target,
    marker_l2 = 0.23,
    z_l2 = 0,
    z_prior = fit$z_prior,
    average = TRUE
  )

  expect_identical(fit$z_anchor, "initializer_logit")
  expect_equal(fit$z_anchor_sd, 0.1)
  expect_identical(fit$z_prior$prior_origin, "cell_specific_marker_initializer_logit_anchor")
  expect_equal(initial_penalty, 0, tolerance = 1e-12)
  expect_equal(fit$best_objective, recomputed, tolerance = 1e-8)
})
