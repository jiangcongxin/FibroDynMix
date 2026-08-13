#!/usr/bin/env Rscript

# Controlled oracle diagnostic for the FibroDynMix simulation model.
#
# This script deliberately lives outside the package API.  It compares state
# recovery when all generative nuisance parameters are known with recovery from
# the currently implemented joint NB fit.  It is a mechanism diagnostic, not
# a replacement for a real-data analysis or a production inference routine.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

has_flag <- function(name) {
  any(args %in% paste0("--", name))
}

parse_logical <- function(x, name) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(x)
  }
  value <- tolower(as.character(x))
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("`%s` must be true or false.", name), call. = FALSE)
}

parse_int_vector <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || anyNA(out) || any(out < 1L)) {
    stop("Integer vector arguments must contain positive comma-separated integers.", call. = FALSE)
  }
  sort(unique(out))
}

parse_char_vector <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  if (length(out) == 0L) {
    stop("Character vector arguments cannot be empty.", call. = FALSE)
  }
  unique(out)
}

print_help <- function() {
  cat(
    paste(
      "Controlled FibroDynMix oracle-identifiability diagnostic",
      "",
      "Usage:",
      "  Rscript scripts/run_oracle_identifiability_diagnostic.R [--name=value]",
      "",
      "Core options:",
      "  --out=analysis/oracle_identifiability_diagnostic",
      "  --scenarios=continuous,batch_confounding",
      "  --n-replicates=10 --seed=260803",
      "  --n-studies=2 --donors-per-study=2 --cells-per-donor=8",
      "  --n-genes=90 --marker-genes-per-state=4",
      "  --joint-n-outer=20 --fixed-beta-n-outer=20",
      "  --posterior-draws=512 --posterior-scenarios=continuous,batch_confounding",
      "",
      "Method switches:",
      "  --run-fixed-beta=true --run-joint=true --run-map=true --run-posterior=true",
      "",
      "The all-truth MAP/posterior methods use the exact conditional",
      "logistic-normal log-ratio prior used by the simulator.  By default they",
      "are limited to continuous and batch_confounding simulations, for which",
      "the latent simplex is not subsequently replaced by a discrete or rare",
      "transition construction.",
      sep = "\n"
    ),
    "\n"
  )
}

if (has_flag("help") || has_flag("h")) {
  print_help()
  quit(status = 0L)
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "scripts/run_oracle_identifiability_diagnostic.R"
}
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "oracle_identifiability_diagnostic"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT)

source_files <- c(
  "matrix_utils.R",
  "z_logistic_normal_prior.R",
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R"
)
source_paths <- file.path(ROOT, "R", source_files)
if (any(!file.exists(source_paths))) {
  stop("Could not find every required R source file.", call. = FALSE)
}
invisible(lapply(source_paths, source))

valid_scenarios <- c("continuous", "discrete", "batch_confounding", "rare_transition")
scenarios <- parse_char_vector(get_arg("scenarios"), c("continuous", "batch_confounding"))
unknown_scenarios <- setdiff(scenarios, valid_scenarios)
if (length(unknown_scenarios) > 0L) {
  stop(sprintf("Unsupported scenario(s): %s", paste(unknown_scenarios, collapse = ", ")), call. = FALSE)
}

map_scenarios <- parse_char_vector(get_arg("map-scenarios"), c("continuous", "batch_confounding"))
posterior_scenarios <- parse_char_vector(get_arg("posterior-scenarios"), c("continuous", "batch_confounding"))
if (length(setdiff(map_scenarios, valid_scenarios)) > 0L ||
    length(setdiff(posterior_scenarios, valid_scenarios)) > 0L) {
  stop("`map-scenarios` and `posterior-scenarios` must contain supported scenarios.", call. = FALSE)
}

n_replicates <- as.integer(get_arg("n-replicates", "10"))
base_seed <- as.integer(get_arg("seed", "260803"))
n_studies <- as.integer(get_arg("n-studies", "2"))
donors_per_study <- as.integer(get_arg("donors-per-study", "2"))
cells_per_donor <- as.integer(get_arg("cells-per-donor", "8"))
n_genes <- as.integer(get_arg("n-genes", "90"))
marker_genes_per_state <- as.integer(get_arg("marker-genes-per-state", "4"))
tau_high <- as.numeric(get_arg("tau-high", "1.0"))
tau_low <- as.numeric(get_arg("tau-low", "0.08"))
study_effect_sd <- as.numeric(get_arg("study-effect-sd", "0.15"))
donor_effect_sd <- as.numeric(get_arg("donor-effect-sd", "0.08"))
state_sd <- as.numeric(get_arg("state-sd", "0.35"))
library_meanlog <- as.numeric(get_arg("library-meanlog", as.character(log(5000))))
library_sdlog <- as.numeric(get_arg("library-sdlog", "0.35"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
oracle_maxit <- as.integer(get_arg("oracle-maxit", "120"))
joint_n_outer <- as.integer(get_arg("joint-n-outer", "20"))
fixed_beta_n_outer <- as.integer(get_arg("fixed-beta-n-outer", as.character(joint_n_outer)))
joint_maxit_beta <- as.integer(get_arg("joint-maxit-beta", "50"))
joint_maxit_z <- as.integer(get_arg("joint-maxit-z", "35"))
fixed_beta_maxit_alpha <- as.integer(get_arg("fixed-beta-maxit-alpha", as.character(joint_maxit_beta)))
fixed_beta_maxit_z <- as.integer(get_arg("fixed-beta-maxit-z", as.character(joint_maxit_z)))
beta_l2 <- as.numeric(get_arg("beta-l2", "0.01"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))
z_l2 <- as.numeric(get_arg("z-l2", "0.001"))
posterior_draws <- as.integer(get_arg("posterior-draws", "512"))
posterior_proposal_scale <- as.numeric(get_arg("posterior-proposal-scale", "1.25"))
run_fixed_beta <- parse_logical(get_arg("run-fixed-beta", "true"), "run-fixed-beta")
run_joint <- parse_logical(get_arg("run-joint", "true"), "run-joint")
run_map <- parse_logical(get_arg("run-map", "true"), "run-map")
run_posterior <- parse_logical(get_arg("run-posterior", "true"), "run-posterior")

positive_ints <- list(
  n_replicates = n_replicates,
  n_studies = n_studies,
  donors_per_study = donors_per_study,
  cells_per_donor = cells_per_donor,
  n_genes = n_genes,
  marker_genes_per_state = marker_genes_per_state,
  initializer_iter = initializer_iter,
  oracle_maxit = oracle_maxit,
  joint_n_outer = joint_n_outer,
  fixed_beta_n_outer = fixed_beta_n_outer,
  joint_maxit_beta = joint_maxit_beta,
  joint_maxit_z = joint_maxit_z,
  fixed_beta_maxit_alpha = fixed_beta_maxit_alpha,
  fixed_beta_maxit_z = fixed_beta_maxit_z,
  posterior_draws = posterior_draws
)
for (name in names(positive_ints)) {
  if (length(positive_ints[[name]]) != 1L || is.na(positive_ints[[name]]) || positive_ints[[name]] < 1L) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
}
if (marker_genes_per_state * 6L > n_genes) {
  stop("`n-genes` must accommodate all six state-specific marker sets.", call. = FALSE)
}
if (any(!is.finite(c(tau_high, tau_low, study_effect_sd, donor_effect_sd, state_sd, library_meanlog, library_sdlog))) ||
    any(c(tau_high, tau_low, study_effect_sd, donor_effect_sd, state_sd, library_sdlog) <= 0)) {
  stop("Simulation scale arguments must be finite and positive.", call. = FALSE)
}
if (any(!is.finite(c(beta_l2, marker_l2, z_l2))) || any(c(beta_l2, marker_l2, z_l2) < 0)) {
  stop("Penalty arguments must be finite and non-negative.", call. = FALSE)
}
if (!is.finite(posterior_proposal_scale) || posterior_proposal_scale <= 0) {
  stop("`posterior-proposal-scale` must be positive.", call. = FALSE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) NA_real_ else stats::sd(x)
}

log_sum_exp <- function(x) {
  if (length(x) == 0L || all(!is.finite(x))) {
    return(-Inf)
  }
  upper <- max(x[is.finite(x)])
  upper + log(sum(exp(x - upper)))
}

logits_matrix_to_simplex <- function(logits) {
  logits <- as.matrix(logits)
  with_reference <- cbind(logits, 0)
  with_reference <- with_reference - apply(with_reference, 1L, max)
  weights <- exp(with_reference)
  weights / rowSums(weights)
}

mvn_log_density <- function(x, mean, covariance, chol_covariance = NULL) {
  x <- as.numeric(x)
  mean <- as.numeric(mean)
  if (is.null(chol_covariance)) {
    chol_covariance <- chol(covariance)
  }
  centered <- x - mean
  solved <- backsolve(chol_covariance, centered, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(chol_covariance))) + sum(solved^2))
}

mvn_log_density_rows <- function(x, mean, covariance, chol_covariance = NULL) {
  x <- as.matrix(x)
  if (is.null(chol_covariance)) {
    chol_covariance <- chol(covariance)
  }
  centered <- sweep(x, 2L, mean, "-")
  solved <- backsolve(chol_covariance, t(centered), transpose = TRUE)
  -0.5 * (
    ncol(x) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      colSums(solved^2)
  )
}

rmvn_draws <- function(n, mean, covariance, chol_covariance = NULL) {
  if (is.null(chol_covariance)) {
    chol_covariance <- chol(covariance)
  }
  standard <- matrix(stats::rnorm(n * length(mean)), nrow = n, ncol = length(mean))
  sweep(standard %*% chol_covariance, 2L, mean, "+")
}

safe_inverse_hessian <- function(hessian, minimum_eigenvalue = 1e-6) {
  hessian <- (hessian + t(hessian)) / 2
  eigen_hessian <- eigen(hessian, symmetric = TRUE)
  eigenvalues <- pmax(eigen_hessian$values, minimum_eigenvalue)
  covariance <- eigen_hessian$vectors %*% diag(1 / eigenvalues, nrow = length(eigenvalues)) %*% t(eigen_hessian$vectors)
  covariance <- (covariance + t(covariance)) / 2
  list(
    covariance = covariance,
    raw_min_eigenvalue = min(eigen_hessian$values),
    raw_max_eigenvalue = max(eigen_hessian$values),
    condition_number = max(eigenvalues) / min(eigenvalues)
  )
}

oracle_logistic_normal_prior <- function(sim, state_sd) {
  state_names <- colnames(sim$z)
  n_states <- length(state_names)
  n_logit <- n_states - 1L
  donor_id <- as.character(sim$cell_metadata$donor_id)
  unique_donors <- unique(donor_id)
  donor_means <- sim$parameters$donor_state_logit_mean
  if (is.null(donor_means)) {
    donor_means <- vapply(
      unique_donors,
      function(id) {
        donor_state_mean(
          donor_row = sim$donor_metadata[id, , drop = FALSE],
          state_names = state_names,
          scenario = sim$parameters$scenario
        )
      },
      numeric(n_states)
    )
    donor_means <- t(donor_means)
    rownames(donor_means) <- unique_donors
  } else {
    donor_means <- as.matrix(donor_means[unique_donors, state_names, drop = FALSE])
  }
  mean_logits_by_cell <- donor_means[donor_id, seq_len(n_logit), drop = FALSE] - donor_means[donor_id, n_states]
  rownames(mean_logits_by_cell) <- rownames(sim$z)
  colnames(mean_logits_by_cell) <- paste0("eta_", seq_len(n_logit))
  covariance <- state_sd^2 * (diag(n_logit) + matrix(1, nrow = n_logit, ncol = n_logit))
  list(
    mean_logits_by_cell = mean_logits_by_cell,
    covariance = covariance,
    chol_covariance = chol(covariance),
    state_sd = state_sd,
    reference_state = state_names[n_states]
  )
}

new_oracle_context <- function(sim, prior = NULL) {
  list(
    counts = sim$counts,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    study_effect = sim$parameters$study_effect,
    donor_effect = sim$parameters$donor_effect,
    study_id = as.character(sim$cell_metadata$study_id),
    donor_id = as.character(sim$cell_metadata$donor_id),
    prior = prior
  )
}

cell_nb_loglik_from_logits <- function(logits, i, context) {
  z_i <- logits_to_simplex(logits)
  eta <- as.vector(z_i %*% context$beta) + context$alpha
  eta <- eta + context$study_effect[context$study_id[i], ]
  eta <- eta + context$donor_effect[context$donor_id[i], ]
  log_mu <- log(context$library_size[i]) + eta
  log_mu <- pmin(pmax(log_mu, -745), 700)
  sum(stats::dnbinom(
    x = as.numeric(context$counts[, i]),
    size = context$phi,
    mu = exp(log_mu),
    log = TRUE
  ))
}

cell_log_prior_from_logits <- function(logits, prior_mean, prior) {
  mvn_log_density(
    x = logits,
    mean = prior_mean,
    covariance = prior$covariance,
    chol_covariance = prior$chol_covariance
  )
}

fit_oracle_z <- function(context,
                         z_start,
                         prior = NULL,
                         maxit = 120L,
                         method_label = "oracle_z_mle") {
  n_cells <- ncol(context$counts)
  n_states <- nrow(context$beta)
  z_hat <- matrix(
    NA_real_, nrow = n_cells, ncol = n_states,
    dimnames = list(colnames(context$counts), rownames(context$beta))
  )
  diagnostics <- vector("list", n_cells)

  for (i in seq_len(n_cells)) {
    start <- simplex_to_logits(z_start[i, ])
    prior_mean <- if (!is.null(prior)) prior$mean_logits_by_cell[i, ] else NULL
    objective <- function(logits) {
      nll <- -cell_nb_loglik_from_logits(logits, i = i, context = context)
      if (!is.null(prior)) {
        nll <- nll - cell_log_prior_from_logits(logits, prior_mean = prior_mean, prior = prior)
      }
      nll
    }
    fit <- tryCatch(
      stats::optim(par = start, fn = objective, method = "BFGS", control = list(maxit = maxit)),
      error = function(error) list(
        par = start,
        value = objective(start),
        convergence = 999L,
        counts = c("function" = NA_integer_, "gradient" = NA_integer_),
        error = conditionMessage(error)
      )
    )
    z_hat[i, ] <- logits_to_simplex(fit$par)
    loglik <- cell_nb_loglik_from_logits(fit$par, i = i, context = context)
    logprior <- if (is.null(prior)) NA_real_ else cell_log_prior_from_logits(fit$par, prior_mean = prior_mean, prior = prior)
    counts_fit <- fit$counts
    diagnostics[[i]] <- data.frame(
      method = method_label,
      cell_id = colnames(context$counts)[i],
      converged = isTRUE(fit$convergence == 0L),
      convergence_code = as.integer(fit$convergence),
      function_evaluations = if (!is.null(counts_fit) && "function" %in% names(counts_fit)) as.integer(counts_fit[["function"]]) else NA_integer_,
      gradient_evaluations = if (!is.null(counts_fit) && "gradient" %in% names(counts_fit)) as.integer(counts_fit[["gradient"]]) else NA_integer_,
      log_likelihood = loglik,
      log_prior = logprior,
      objective = -loglik - ifelse(is.na(logprior), 0, logprior),
      error = if (!is.null(fit$error)) as.character(fit$error) else "",
      stringsAsFactors = FALSE
    )
  }
  list(
    z_hat = z_hat,
    cell_diagnostics = do.call(rbind, diagnostics),
    prior = prior
  )
}

matrix_nb_loglik_for_z <- function(z, context) {
  fibrodynmix_nb_loglik(
    counts = context$counts,
    z = z,
    beta = context$beta,
    alpha = context$alpha,
    phi = context$phi,
    library_size = context$library_size,
    study_effect = context$study_effect,
    donor_effect = context$donor_effect,
    study_id = context$study_id,
    donor_id = context$donor_id,
    return_matrix = FALSE
  )
}

matrix_log_prior_for_z <- function(z, prior) {
  logits <- t(apply(z, 1L, simplex_to_logits))
  values <- vapply(
    seq_len(nrow(logits)),
    function(i) cell_log_prior_from_logits(logits[i, ], prior$mean_logits_by_cell[i, ], prior),
    numeric(1)
  )
  sum(values)
}

posterior_mean_importance <- function(context,
                                      z_map,
                                      prior,
                                      n_draws,
                                      proposal_scale,
                                      seed_base,
                                      maxit) {
  n_cells <- ncol(context$counts)
  n_states <- nrow(context$beta)
  z_mean <- matrix(
    NA_real_, nrow = n_cells, ncol = n_states,
    dimnames = list(colnames(context$counts), rownames(context$beta))
  )
  diagnostics <- vector("list", n_cells)

  for (i in seq_len(n_cells)) {
    map_logits <- simplex_to_logits(z_map[i, ])
    prior_mean <- prior$mean_logits_by_cell[i, ]
    neg_log_posterior <- function(logits) {
      -cell_nb_loglik_from_logits(logits, i = i, context = context) -
        cell_log_prior_from_logits(logits, prior_mean = prior_mean, prior = prior)
    }
    hessian <- tryCatch(
      stats::optimHess(par = map_logits, fn = neg_log_posterior),
      error = function(error) diag(1 / prior$state_sd^2, length(map_logits))
    )
    inverse <- safe_inverse_hessian(hessian)
    proposal_covariance <- inverse$covariance * proposal_scale^2
    proposal_chol <- chol(proposal_covariance)
    posterior_seed <- as.integer(seed_base + i)
    set.seed(posterior_seed)
    draws <- rmvn_draws(
      n = n_draws,
      mean = map_logits,
      covariance = proposal_covariance,
      chol_covariance = proposal_chol
    )
    z_draws <- logits_matrix_to_simplex(draws)
    effect <- context$study_effect[context$study_id[i], ] + context$donor_effect[context$donor_id[i], ]
    eta <- z_draws %*% context$beta
    eta <- sweep(eta, 2L, context$alpha + effect + log(context$library_size[i]), "+")
    eta <- pmin(pmax(eta, -745), 700)
    mu <- exp(eta)
    y <- as.numeric(context$counts[, i])
    loglik <- rowSums(matrix(
      stats::dnbinom(
        x = rep(y, each = n_draws),
        size = rep(context$phi, each = n_draws),
        mu = as.vector(mu),
        log = TRUE
      ),
      nrow = n_draws,
      ncol = length(y)
    ))
    logprior <- mvn_log_density_rows(
      x = draws,
      mean = prior_mean,
      covariance = prior$covariance,
      chol_covariance = prior$chol_covariance
    )
    logproposal <- mvn_log_density_rows(
      x = draws,
      mean = map_logits,
      covariance = proposal_covariance,
      chol_covariance = proposal_chol
    )
    log_weights <- loglik + logprior - logproposal
    log_normalizer <- log_sum_exp(log_weights)
    weights <- exp(log_weights - log_normalizer)
    z_mean[i, ] <- colSums(z_draws * weights)
    ess <- 1 / sum(weights^2)
    diagnostics[[i]] <- data.frame(
      cell_id = colnames(context$counts)[i],
      posterior_seed = posterior_seed,
      n_draws = n_draws,
      importance_ess = ess,
      importance_ess_fraction = ess / n_draws,
      maximum_weight = max(weights),
      log_weight_range = max(log_weights) - min(log_weights),
      hessian_min_eigenvalue = inverse$raw_min_eigenvalue,
      hessian_max_eigenvalue = inverse$raw_max_eigenvalue,
      hessian_condition_number = inverse$condition_number,
      stringsAsFactors = FALSE
    )
  }
  list(z_mean = z_mean, diagnostics = do.call(rbind, diagnostics))
}

update_alpha_fixed_beta <- function(counts,
                                    z,
                                    beta,
                                    alpha,
                                    phi,
                                    library_size,
                                    maxit) {
  n_genes <- nrow(counts)
  base_eta <- z %*% beta
  alpha_new <- as.numeric(alpha)
  converged <- logical(n_genes)
  for (g in seq_len(n_genes)) {
    y <- as.numeric(counts[g, ])
    fit <- tryCatch(
      stats::optim(
        par = alpha_new[g],
        fn = function(value) {
          log_mu <- log(library_size) + base_eta[, g] + value
          log_mu <- pmin(pmax(log_mu, -745), 700)
          -sum(stats::dnbinom(y, size = phi[g], mu = exp(log_mu), log = TRUE))
        },
        method = "BFGS",
        control = list(maxit = maxit)
      ),
      error = function(error) list(par = alpha_new[g], convergence = 999L)
    )
    alpha_new[g] <- fit$par[1L]
    converged[g] <- isTRUE(fit$convergence == 0L)
  }
  names(alpha_new) <- rownames(counts)
  list(alpha = alpha_new, convergence_rate = mean(converged))
}

fixed_beta_mu <- function(z, beta, alpha, library_size) {
  eta <- z %*% beta
  eta <- sweep(eta, 2L, alpha, "+")
  mu <- exp(pmin(pmax(eta, -745), 700))
  mu <- sweep(mu, 1L, library_size, "*")
  t(mu)
}

fit_fixed_beta_current_nuisance <- function(counts,
                                             z_start,
                                             beta_truth,
                                             library_size,
                                             z_truth = NULL,
                                             n_outer,
                                             maxit_alpha,
                                             maxit_z,
                                             z_l2) {
  z_hat <- z_start
  alpha_hat <- initialize_nb_alpha(counts, library_size)
  phi_hat <- estimate_phi_moments(counts, library_size)
  trace <- vector("list", n_outer + 1L)
  record_trace <- function(iteration, alpha_rate = NA_real_, z_rate = NA_real_) {
    loglik <- fibrodynmix_nb_loglik(
      counts = counts,
      z = z_hat,
      beta = beta_truth,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size
    )
    logits <- t(apply(z_hat, 1L, simplex_to_logits))
    data.frame(
      iteration = iteration,
      fit_log_likelihood = loglik,
      fit_log_likelihood_per_entry = loglik / length(counts),
      z_logit_l2_penalty = z_l2 * sum(logits^2),
      pseudo_objective_per_entry = (-loglik + z_l2 * sum(logits^2)) / length(counts),
      alpha_convergence_rate = alpha_rate,
      z_convergence_rate = z_rate,
      rmse = if (is.null(z_truth)) NA_real_ else evaluate_state_weights(z_truth, z_hat)$rmse,
      stringsAsFactors = FALSE
    )
  }
  trace[[1L]] <- record_trace(0L)
  for (iteration in seq_len(n_outer)) {
    alpha_update <- update_alpha_fixed_beta(
      counts = counts,
      z = z_hat,
      beta = beta_truth,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size,
      maxit = maxit_alpha
    )
    alpha_hat <- alpha_update$alpha
    z_update <- update_z_nb(
      counts = counts,
      z = z_hat,
      beta = beta_truth,
      alpha = alpha_hat,
      phi = phi_hat,
      library_size = library_size,
      study_effect = NULL,
      study_id = NULL,
      donor_effect = NULL,
      donor_id = NULL,
      z_l2 = z_l2,
      optimizer = "BFGS",
      optimizer_control = list(),
      maxit = maxit_z
    )
    z_hat <- z_update$z
    phi_hat <- estimate_phi_moments(
      counts = counts,
      library_size = library_size,
      mu = fixed_beta_mu(z_hat, beta_truth, alpha_hat, library_size)
    )
    trace[[iteration + 1L]] <- record_trace(
      iteration = iteration,
      alpha_rate = alpha_update$convergence_rate,
      z_rate = z_update$convergence_rate
    )
  }
  list(
    z_hat = z_hat,
    alpha_hat = alpha_hat,
    phi_hat = phi_hat,
    trace = do.call(rbind, trace)
  )
}

truth_parameter_loglik <- function(sim, z) {
  fibrodynmix_nb_loglik(
    counts = sim$counts,
    z = z,
    beta = sim$parameters$beta_kg,
    alpha = sim$parameters$alpha_g,
    phi = sim$parameters$phi_g,
    library_size = sim$cell_metadata$library_size,
    study_effect = sim$parameters$study_effect,
    donor_effect = sim$parameters$donor_effect,
    study_id = sim$cell_metadata$study_id,
    donor_id = sim$cell_metadata$donor_id
  )
}

metric_row <- function(scenario,
                       replicate,
                       simulation_seed,
                       initializer_seed,
                       method,
                       parameter_setting,
                       z_hat,
                       sim,
                       prior = NULL,
                       fit_loglik = NA_real_,
                       status = "ok",
                       note = "") {
  recovered <- evaluate_state_weights(sim$z, z_hat)
  truth_ll <- truth_parameter_loglik(sim, z_hat)
  logit_error <- t(apply(z_hat, 1L, simplex_to_logits)) - t(apply(sim$z, 1L, simplex_to_logits))
  truth_prior <- if (is.null(prior)) NA_real_ else matrix_log_prior_for_z(z_hat, prior)
  data.frame(
    scenario = scenario,
    replicate = replicate,
    simulation_seed = simulation_seed,
    initializer_seed = initializer_seed,
    method = method,
    parameter_setting = parameter_setting,
    status = status,
    rmse = recovered$rmse,
    mean_absolute_error = recovered$mean_absolute_error,
    dominant_accuracy = recovered$dominant_accuracy,
    macro_dominant_f1 = recovered$macro_dominant_f1,
    mean_abs_logit_error = mean(abs(logit_error)),
    truth_parameter_log_likelihood = truth_ll,
    truth_parameter_log_likelihood_per_entry = truth_ll / length(sim$counts),
    truth_logistic_normal_log_prior = truth_prior,
    truth_log_posterior = if (is.na(truth_prior)) NA_real_ else truth_ll + truth_prior,
    method_fit_log_likelihood = fit_loglik,
    method_fit_log_likelihood_per_entry = if (is.na(fit_loglik)) NA_real_ else fit_loglik / length(sim$counts),
    n_cells = ncol(sim$counts),
    n_genes = nrow(sim$counts),
    note = note,
    stringsAsFactors = FALSE
  )
}

config <- data.frame(
  key = c(
    "analysis", "output_directory", "scenarios", "map_scenarios", "posterior_scenarios",
    "n_replicates", "base_seed", "n_studies", "donors_per_study", "cells_per_donor",
    "n_genes", "marker_genes_per_state", "tau_high", "tau_low", "study_effect_sd",
    "donor_effect_sd", "state_sd", "library_meanlog", "library_sdlog", "initializer_iter",
    "oracle_maxit", "joint_n_outer", "fixed_beta_n_outer", "joint_maxit_beta", "joint_maxit_z",
    "fixed_beta_maxit_alpha", "fixed_beta_maxit_z", "beta_l2", "marker_l2", "z_l2",
    "posterior_draws", "posterior_proposal_scale", "run_fixed_beta", "run_joint", "run_map", "run_posterior"
  ),
  value = as.character(c(
    "oracle_identifiability_diagnostic", OUT, paste(scenarios, collapse = ";"),
    paste(map_scenarios, collapse = ";"), paste(posterior_scenarios, collapse = ";"),
    n_replicates, base_seed, n_studies, donors_per_study, cells_per_donor,
    n_genes, marker_genes_per_state, tau_high, tau_low, study_effect_sd,
    donor_effect_sd, state_sd, library_meanlog, library_sdlog, initializer_iter,
    oracle_maxit, joint_n_outer, fixed_beta_n_outer, joint_maxit_beta, joint_maxit_z,
    fixed_beta_maxit_alpha, fixed_beta_maxit_z, beta_l2, marker_l2, z_l2,
    posterior_draws, posterior_proposal_scale, run_fixed_beta, run_joint, run_map, run_posterior
  )),
  stringsAsFactors = FALSE
)
write_tsv(config, file.path(OUT, "config.tsv"))
write_tsv(
  data.frame(
    source_file = source_files,
    md5 = unname(tools::md5sum(source_paths)),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "source_manifest.tsv")
)

metric_rows <- list()
oracle_cell_rows <- list()
posterior_rows <- list()
fixed_trace_rows <- list()
joint_trace_rows <- list()
seed_rows <- list()
metric_index <- 1L
oracle_index <- 1L
posterior_index <- 1L
fixed_trace_index <- 1L
joint_trace_index <- 1L
seed_index <- 1L

for (scenario in scenarios) {
  scenario_index <- match(scenario, valid_scenarios)
  for (replicate in seq_len(n_replicates)) {
    simulation_seed <- as.integer(base_seed + scenario_index * 100000L + replicate)
    initializer_seed <- as.integer(simulation_seed + 10000000L)
    posterior_seed_base <- as.integer(simulation_seed + 20000000L)
    sim <- simulate_fibrodynmix(
      n_studies = n_studies,
      donors_per_study = donors_per_study,
      cells_per_donor = cells_per_donor,
      n_genes = n_genes,
      marker_genes_per_state = marker_genes_per_state,
      scenario = scenario,
      tau_high = tau_high,
      tau_low = tau_low,
      study_effect_sd = study_effect_sd,
      donor_effect_sd = donor_effect_sd,
      state_sd = state_sd,
      library_meanlog = library_meanlog,
      library_sdlog = library_sdlog,
      seed = simulation_seed
    )
    seed_rows[[seed_index]] <- data.frame(
      scenario = scenario,
      replicate = replicate,
      simulation_seed = simulation_seed,
      initializer_seed = initializer_seed,
      posterior_seed_base = posterior_seed_base,
      stringsAsFactors = FALSE
    )
    seed_index <- seed_index + 1L

    # Resetting this seed before the joint fit makes its random initializer
    # exactly match the separately recorded initializer below.
    set.seed(initializer_seed)
    initializer <- fit_fibrodynmix_initializer(
      counts = sim$counts,
      marker_index = sim$parameters$marker_index,
      library_size = sim$cell_metadata$library_size,
      n_iter = initializer_iter
    )
    metric_rows[[metric_index]] <- metric_row(
      scenario = scenario,
      replicate = replicate,
      simulation_seed = simulation_seed,
      initializer_seed = initializer_seed,
      method = "initializer",
      parameter_setting = "marker_guided_log_normalized_initializer",
      z_hat = initializer$z_hat,
      sim = sim,
      note = "Shared deterministic start for oracle and current joint comparisons."
    )
    metric_index <- metric_index + 1L

    prior <- oracle_logistic_normal_prior(sim, state_sd = state_sd)
    oracle_context <- new_oracle_context(sim, prior = prior)

    oracle_mle <- fit_oracle_z(
      context = oracle_context,
      z_start = initializer$z_hat,
      prior = NULL,
      maxit = oracle_maxit,
      method_label = "oracle_z_mle_all_truth"
    )
    oracle_mle_ll <- matrix_nb_loglik_for_z(oracle_mle$z_hat, oracle_context)
    metric_rows[[metric_index]] <- metric_row(
      scenario = scenario,
      replicate = replicate,
      simulation_seed = simulation_seed,
      initializer_seed = initializer_seed,
      method = "oracle_z_mle_all_truth",
      parameter_setting = "true_beta_alpha_phi_study_and_donor_effects",
      z_hat = oracle_mle$z_hat,
      sim = sim,
      fit_loglik = oracle_mle_ll,
      note = "Conditional per-cell NB maximum likelihood; no z prior."
    )
    metric_index <- metric_index + 1L
    mle_cells <- oracle_mle$cell_diagnostics
    mle_cells$scenario <- scenario
    mle_cells$replicate <- replicate
    mle_cells$simulation_seed <- simulation_seed
    oracle_cell_rows[[oracle_index]] <- mle_cells
    oracle_index <- oracle_index + 1L

    oracle_map <- NULL
    if (isTRUE(run_map) && scenario %in% map_scenarios) {
      oracle_map <- fit_oracle_z(
        context = oracle_context,
        z_start = initializer$z_hat,
        prior = prior,
        maxit = oracle_maxit,
        method_label = "oracle_z_map_all_truth"
      )
      oracle_map_ll <- matrix_nb_loglik_for_z(oracle_map$z_hat, oracle_context)
      metric_rows[[metric_index]] <- metric_row(
        scenario = scenario,
        replicate = replicate,
        simulation_seed = simulation_seed,
        initializer_seed = initializer_seed,
        method = "oracle_z_map_all_truth",
        parameter_setting = "true_count_parameters_and_exact_conditional_logistic_normal_prior",
        z_hat = oracle_map$z_hat,
        sim = sim,
        prior = prior,
        fit_loglik = oracle_map_ll,
        note = "MAP under the simulator's donor-specific logistic-normal log-ratio prior."
      )
      metric_index <- metric_index + 1L
      map_cells <- oracle_map$cell_diagnostics
      map_cells$scenario <- scenario
      map_cells$replicate <- replicate
      map_cells$simulation_seed <- simulation_seed
      oracle_cell_rows[[oracle_index]] <- map_cells
      oracle_index <- oracle_index + 1L
    }

    if (isTRUE(run_posterior) && scenario %in% posterior_scenarios) {
      if (is.null(oracle_map)) {
        oracle_map <- fit_oracle_z(
          context = oracle_context,
          z_start = initializer$z_hat,
          prior = prior,
          maxit = oracle_maxit,
          method_label = "oracle_z_map_all_truth"
        )
      }
      posterior <- posterior_mean_importance(
        context = oracle_context,
        z_map = oracle_map$z_hat,
        prior = prior,
        n_draws = posterior_draws,
        proposal_scale = posterior_proposal_scale,
        seed_base = posterior_seed_base,
        maxit = oracle_maxit
      )
      posterior_ll <- matrix_nb_loglik_for_z(posterior$z_mean, oracle_context)
      metric_rows[[metric_index]] <- metric_row(
        scenario = scenario,
        replicate = replicate,
        simulation_seed = simulation_seed,
        initializer_seed = initializer_seed,
        method = "oracle_z_posterior_mean_all_truth",
        parameter_setting = "self_normalized_importance_posterior_mean_with_true_count_parameters_and_prior",
        z_hat = posterior$z_mean,
        sim = sim,
        prior = prior,
        fit_loglik = posterior_ll,
        note = "Self-normalized importance posterior mean; inspect ESS before interpreting."
      )
      metric_index <- metric_index + 1L
      posterior_diagnostics <- posterior$diagnostics
      posterior_diagnostics$scenario <- scenario
      posterior_diagnostics$replicate <- replicate
      posterior_diagnostics$simulation_seed <- simulation_seed
      posterior_rows[[posterior_index]] <- posterior_diagnostics
      posterior_index <- posterior_index + 1L
    }

    if (isTRUE(run_fixed_beta)) {
      fixed_beta_fit <- fit_fixed_beta_current_nuisance(
        counts = sim$counts,
        z_start = initializer$z_hat,
        beta_truth = sim$parameters$beta_kg,
        library_size = sim$cell_metadata$library_size,
        z_truth = sim$z,
        n_outer = fixed_beta_n_outer,
        maxit_alpha = fixed_beta_maxit_alpha,
        maxit_z = fixed_beta_maxit_z,
        z_l2 = z_l2
      )
      fixed_beta_ll <- fibrodynmix_nb_loglik(
        counts = sim$counts,
        z = fixed_beta_fit$z_hat,
        beta = sim$parameters$beta_kg,
        alpha = fixed_beta_fit$alpha_hat,
        phi = fixed_beta_fit$phi_hat,
        library_size = sim$cell_metadata$library_size
      )
      metric_rows[[metric_index]] <- metric_row(
        scenario = scenario,
        replicate = replicate,
        simulation_seed = simulation_seed,
        initializer_seed = initializer_seed,
        method = "fixed_beta_estimated_nuisance",
        parameter_setting = "true_beta_estimated_alpha_phi_no_study_or_donor_effect",
        z_hat = fixed_beta_fit$z_hat,
        sim = sim,
        fit_loglik = fixed_beta_ll,
        note = "True beta fixed; alpha and phi re-estimated; current no-effect nuisance specification retained."
      )
      metric_index <- metric_index + 1L
      fixed_trace <- fixed_beta_fit$trace
      fixed_trace$scenario <- scenario
      fixed_trace$replicate <- replicate
      fixed_trace$simulation_seed <- simulation_seed
      fixed_trace_rows[[fixed_trace_index]] <- fixed_trace
      fixed_trace_index <- fixed_trace_index + 1L
    }

    if (isTRUE(run_joint)) {
      set.seed(initializer_seed)
      joint_fit <- fit_fibrodynmix_nb(
        counts = sim$counts,
        marker_index = sim$parameters$marker_index,
        library_size = sim$cell_metadata$library_size,
        n_outer = joint_n_outer,
        initializer_args = list(n_iter = initializer_iter),
        beta_l2 = beta_l2,
        marker_l2 = marker_l2,
        z_l2 = z_l2,
        z_anchor = "none",
        maxit_beta = joint_maxit_beta,
        maxit_z = joint_maxit_z,
        rollback_to_best = TRUE
      )
      joint_ll <- fibrodynmix_nb_loglik(
        counts = sim$counts,
        z = joint_fit$z_hat,
        beta = joint_fit$beta_hat,
        alpha = joint_fit$alpha_hat,
        phi = joint_fit$phi_hat,
        library_size = sim$cell_metadata$library_size
      )
      metric_rows[[metric_index]] <- metric_row(
        scenario = scenario,
        replicate = replicate,
        simulation_seed = simulation_seed,
        initializer_seed = initializer_seed,
        method = "joint_current_nb",
        parameter_setting = "current_joint_nb_alpha_beta_phi_no_study_or_donor_effect",
        z_hat = joint_fit$z_hat,
        sim = sim,
        fit_loglik = joint_ll,
        note = "Legacy unanchored alternating NB implementation, deterministically initialized from the reported initializer."
      )
      metric_index <- metric_index + 1L
      joint_trace_rows[[joint_trace_index]] <- data.frame(
        scenario = scenario,
        replicate = replicate,
        simulation_seed = simulation_seed,
        iteration = seq_along(joint_fit$nb_objective_trace) - 1L,
        nb_objective = joint_fit$nb_objective_trace,
        best_iteration = joint_fit$best_iteration,
        executed_iterations = joint_fit$executed_iterations,
        stop_reason = joint_fit$stop_reason,
        convergence_beta = c(NA, joint_fit$convergence$beta),
        convergence_z = c(NA, joint_fit$convergence$z),
        rolled_back = c(NA, joint_fit$convergence$rolled_back),
        stringsAsFactors = FALSE
      )
      joint_trace_index <- joint_trace_index + 1L
    }
  }
}

metrics <- do.call(rbind, metric_rows)
write_tsv(metrics, file.path(OUT, "diagnostic_metrics.tsv"))
write_tsv(do.call(rbind, seed_rows), file.path(OUT, "seed_manifest.tsv"))

if (length(oracle_cell_rows) > 0L) {
  write_tsv(do.call(rbind, oracle_cell_rows), file.path(OUT, "oracle_cell_optimization.tsv"))
} else {
  write_tsv(data.frame(), file.path(OUT, "oracle_cell_optimization.tsv"))
}
if (length(posterior_rows) > 0L) {
  write_tsv(do.call(rbind, posterior_rows), file.path(OUT, "posterior_importance_diagnostics.tsv"))
} else {
  write_tsv(data.frame(), file.path(OUT, "posterior_importance_diagnostics.tsv"))
}
if (length(fixed_trace_rows) > 0L) {
  write_tsv(do.call(rbind, fixed_trace_rows), file.path(OUT, "fixed_beta_nuisance_trace.tsv"))
} else {
  write_tsv(data.frame(), file.path(OUT, "fixed_beta_nuisance_trace.tsv"))
}
if (length(joint_trace_rows) > 0L) {
  write_tsv(do.call(rbind, joint_trace_rows), file.path(OUT, "joint_current_nb_trace.tsv"))
} else {
  write_tsv(data.frame(), file.path(OUT, "joint_current_nb_trace.tsv"))
}

summary_keys <- unique(metrics[, c("scenario", "method", "parameter_setting"), drop = FALSE])
summary_rows <- lapply(seq_len(nrow(summary_keys)), function(i) {
  keep <- metrics$scenario == summary_keys$scenario[i] &
    metrics$method == summary_keys$method[i] &
    metrics$parameter_setting == summary_keys$parameter_setting[i] &
    metrics$status == "ok"
  current <- metrics[keep, , drop = FALSE]
  data.frame(
    scenario = summary_keys$scenario[i],
    method = summary_keys$method[i],
    parameter_setting = summary_keys$parameter_setting[i],
    n_replicates = nrow(current),
    rmse_mean = safe_mean(current$rmse),
    rmse_sd = safe_sd(current$rmse),
    mean_absolute_error_mean = safe_mean(current$mean_absolute_error),
    dominant_accuracy_mean = safe_mean(current$dominant_accuracy),
    macro_dominant_f1_mean = safe_mean(current$macro_dominant_f1),
    mean_abs_logit_error_mean = safe_mean(current$mean_abs_logit_error),
    truth_parameter_log_likelihood_per_entry_mean = safe_mean(current$truth_parameter_log_likelihood_per_entry),
    truth_log_posterior_mean = safe_mean(current$truth_log_posterior),
    method_fit_log_likelihood_per_entry_mean = safe_mean(current$method_fit_log_likelihood_per_entry),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write_tsv(summary, file.path(OUT, "diagnostic_summary.tsv"))

posterior_summary <- if (length(posterior_rows) > 0L) {
  posterior_all <- do.call(rbind, posterior_rows)
  stats::aggregate(
    cbind(importance_ess, importance_ess_fraction, maximum_weight, log_weight_range) ~ scenario,
    data = posterior_all,
    FUN = mean
  )
} else {
  data.frame()
}
write_tsv(posterior_summary, file.path(OUT, "posterior_importance_summary.tsv"))

readme <- c(
  "# FibroDynMix oracle-identifiability diagnostic",
  "",
  "This directory is produced by scripts/run_oracle_identifiability_diagnostic.R.",
  "It is a simulation mechanism diagnostic, not a real-data benchmark.",
  "",
  "Method definitions:",
  "- initializer: marker-guided log-normalized initializer.",
  "- oracle_z_mle_all_truth: only z is optimized while beta, alpha, phi, study effects, and donor effects are fixed at simulated truth; no prior is used.",
  "- oracle_z_map_all_truth: the same oracle count model with the simulator's donor-specific logistic-normal prior expressed in reference-logit coordinates.",
  "- oracle_z_posterior_mean_all_truth: self-normalized importance posterior mean around the oracle MAP; inspect ESS diagnostics before drawing conclusions.",
  "- fixed_beta_estimated_nuisance: beta is fixed at truth while alpha and phi are re-estimated with the current no-study/no-donor nuisance specification.",
  "- joint_current_nb: current alternating NB fit using its usual marker and ridge penalties, without study/donor effects unless the script is extended.",
  "",
  "Interpretation boundary:",
  "A recovery difference between the joint fit and the oracle methods localizes a problem to joint estimation, priors/regularization, nuisance handling, or optimization. It does not by itself identify an exact rotation or scaling mechanism. The discrete and rare-transition scenarios do not retain the exact continuous logistic-normal prior after simulation-time transformations; by default MAP/posterior methods are restricted to continuous and batch_confounding scenarios."
)
writeLines(readme, con = file.path(OUT, "README.md"))

run_manifest <- data.frame(
  analysis = "oracle_identifiability_diagnostic",
  status = "completed",
  scenarios = paste(scenarios, collapse = ";"),
  n_replicates = n_replicates,
  n_metric_rows = nrow(metrics),
  n_oracle_cell_rows = if (length(oracle_cell_rows) == 0L) 0L else nrow(do.call(rbind, oracle_cell_rows)),
  n_posterior_cell_rows = if (length(posterior_rows) == 0L) 0L else nrow(do.call(rbind, posterior_rows)),
  primary_comparison = "initializer vs conditional all-truth z-MLE vs all-truth logistic-normal MAP/posterior mean vs fixed-beta estimated-nuisance vs current joint NB",
  claim_boundary = "The script tests mechanisms under the project simulator. It cannot alone establish a general likelihood--recovery principle or a unique non-identifiability mechanism.",
  r_version = R.version.string,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(run_manifest, file.path(OUT, "run_manifest.tsv"))

message("Oracle identifiability diagnostic written to: ", OUT)
