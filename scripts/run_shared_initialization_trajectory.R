#!/usr/bin/env Rscript

# Shared-initialization trajectory audit for the raw-count NB optimizer.
#
# This script deliberately does not call fit_fibrodynmix_nb() repeatedly for
# different outer-iteration budgets. Instead, it generates one simulation and
# one marker-guided initializer per scenario/replicate, then advances each
# method along a single, saved alternating-optimization path. The resulting
# state snapshots make it possible to distinguish a true within-run trajectory
# from differences caused by independent initializer draws.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

has_flag <- function(name) {
  paste0("--", name) %in% args
}

parse_bool <- function(x, name) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("`%s` must be true or false.", name), call. = FALSE)
}

parse_positive_integer <- function(x, name) {
  out <- suppressWarnings(as.integer(x))
  if (length(out) != 1L || is.na(out) || out < 1L || out != as.numeric(x)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  out
}

parse_nonnegative_number <- function(x, name) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) != 1L || is.na(out) || !is.finite(out) || out < 0) {
    stop(sprintf("`%s` must be a non-negative finite number.", name), call. = FALSE)
  }
  out
}

parse_character_vector <- function(x, default, name) {
  if (is.null(x)) {
    return(default)
  }
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  if (length(out) == 0L) {
    stop(sprintf("`%s` cannot be empty.", name), call. = FALSE)
  }
  unique(out)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  stats::sd(x)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "scripts/run_shared_initialization_trajectory.R"
}
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (has_flag("help") || has_flag("h")) {
  cat(paste(
    "Usage:",
    "  Rscript scripts/run_shared_initialization_trajectory.R [options]",
    "",
    "Each scenario/replicate is simulated once and initialized once. Every method",
    "then advances from that exact saved initial NB state, writing state snapshots",
    "at iteration 0 and after each accepted or rolled-back outer iteration.",
    "",
    "Core options:",
    "  --out=PATH                         Output directory (must be empty).",
    "  --scenarios=a,b,c                  Defaults to all four simulation scenarios.",
    "  --methods=a,b                      Defaults to fibrodynmix_nb,fibrodynmix_nb_study.",
    "                                      Optional donor method: fibrodynmix_nb_study_donor.",
    "  --n-replicates=5                   Number of deterministic replicates per scenario.",
    "  --n-outer=20                       Number of successive outer updates per trajectory.",
    "  --seed=940                         Base simulation seed.",
    "  --early-stopping=false             Defaults to false so all requested states are saved.",
    "  --rollback-on-increase=true        Mirrors current NB objective rollback behavior.",
    "",
    "Simulation options:",
    "  --n-studies=2 --donors-per-study=2 --cells-per-donor=8",
    "  --n-genes=90 --marker-genes-per-state=4",
    "",
    "Optimizer options:",
    "  --initializer-iter=3 --maxit-beta=50 --maxit-z=35",
    "  --beta-l2=0.01 --marker-l2=0.05 --beta-constraint=sum_to_zero",
    "  --marker-l2-schedule=constant       constant or likelihood_normalized",
    "  --marker-l2-nll-fraction=0.01       Target marker penalty / NB NLL when normalized",
    "  --study-l2=5 --donor-l2=0.1 --z-l2=0.001",
    "  --z-prior=none                      none, empirical_donor, empirical_study, initializer_cell, or oracle_donor",
    "  --z-prior-covariance-ridge=0.1      Ridge for empirical pooled covariance",
    "  --z-prior-cell-sd=1                 Reference-logit SD for initializer_cell anchoring",
    "  --estimate-phi=true --min-delta=1e-5",
    sep = "\n"
  ), "\n")
  quit(status = 0L)
}

OUT <- get_arg("out", file.path(ROOT, "analysis", "shared_initialization_trajectory"))
OUT <- normalizePath(OUT, mustWork = FALSE)
if (dir.exists(OUT)) {
  existing <- list.files(OUT, all.files = TRUE, no.. = TRUE)
  if (length(existing) > 0L) {
    stop(
      sprintf("Output directory already contains files: %s. Choose a new `--out` directory.", OUT),
      call. = FALSE
    )
  }
} else {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
}

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
invisible(lapply(file.path(ROOT, "R", source_files), source))

valid_scenarios <- c("continuous", "discrete", "batch_confounding", "rare_transition")
scenarios <- parse_character_vector(
  get_arg("scenarios"),
  valid_scenarios,
  "scenarios"
)
unknown_scenarios <- setdiff(scenarios, valid_scenarios)
if (length(unknown_scenarios) > 0L) {
  stop(
    sprintf("Unsupported scenario(s): %s", paste(unknown_scenarios, collapse = ", ")),
    call. = FALSE
  )
}

method_specs <- list(
  fibrodynmix_nb = list(fit_study_effect = FALSE, fit_donor_effect = FALSE),
  fibrodynmix_nb_study = list(fit_study_effect = TRUE, fit_donor_effect = FALSE),
  fibrodynmix_nb_study_donor = list(fit_study_effect = TRUE, fit_donor_effect = TRUE)
)
methods <- parse_character_vector(
  get_arg("methods"),
  c("fibrodynmix_nb", "fibrodynmix_nb_study"),
  "methods"
)
unknown_methods <- setdiff(methods, names(method_specs))
if (length(unknown_methods) > 0L) {
  stop(
    sprintf("Unsupported method(s): %s", paste(unknown_methods, collapse = ", ")),
    call. = FALSE
  )
}

n_replicates <- parse_positive_integer(get_arg("n-replicates", "5"), "n-replicates")
n_outer <- parse_positive_integer(get_arg("n-outer", "20"), "n-outer")
base_seed <- parse_positive_integer(get_arg("seed", "940"), "seed")
n_studies <- parse_positive_integer(get_arg("n-studies", "2"), "n-studies")
donors_per_study <- parse_positive_integer(get_arg("donors-per-study", "2"), "donors-per-study")
cells_per_donor <- parse_positive_integer(get_arg("cells-per-donor", "8"), "cells-per-donor")
n_genes <- parse_positive_integer(get_arg("n-genes", "90"), "n-genes")
marker_genes_per_state <- parse_positive_integer(get_arg("marker-genes-per-state", "4"), "marker-genes-per-state")
initializer_iter <- parse_positive_integer(get_arg("initializer-iter", "3"), "initializer-iter")
maxit_beta <- parse_positive_integer(get_arg("maxit-beta", "50"), "maxit-beta")
maxit_z <- parse_positive_integer(get_arg("maxit-z", "35"), "maxit-z")
beta_l2 <- parse_nonnegative_number(get_arg("beta-l2", "0.01"), "beta-l2")
marker_l2 <- parse_nonnegative_number(get_arg("marker-l2", "0.05"), "marker-l2")
marker_l2_schedule <- match.arg(
  get_arg("marker-l2-schedule", "constant"),
  c("constant", "likelihood_normalized")
)
marker_l2_nll_fraction <- parse_nonnegative_number(
  get_arg("marker-l2-nll-fraction", "0.01"),
  "marker-l2-nll-fraction"
)
beta_constraint <- match.arg(
  get_arg("beta-constraint", "sum_to_zero"),
  c("sum_to_zero", "none")
)
study_l2 <- parse_nonnegative_number(get_arg("study-l2", "5"), "study-l2")
donor_l2 <- parse_nonnegative_number(get_arg("donor-l2", "0.1"), "donor-l2")
z_l2 <- parse_nonnegative_number(get_arg("z-l2", "0.001"), "z-l2")
z_prior_mode <- match.arg(
  get_arg("z-prior", "none"),
  c("none", "empirical_donor", "empirical_study", "initializer_cell", "oracle_donor")
)
z_prior_covariance_ridge <- parse_nonnegative_number(
  get_arg("z-prior-covariance-ridge", "0.1"),
  "z-prior-covariance-ridge"
)
z_prior_cell_sd <- parse_nonnegative_number(
  get_arg("z-prior-cell-sd", "1"),
  "z-prior-cell-sd"
)
if (identical(z_prior_mode, "initializer_cell") && z_prior_cell_sd <= 0) {
  stop("`z-prior-cell-sd` must be positive for `z-prior=initializer_cell`.", call. = FALSE)
}
estimate_phi <- parse_bool(get_arg("estimate-phi", "true"), "estimate-phi")
rollback_on_increase <- parse_bool(get_arg("rollback-on-increase", "true"), "rollback-on-increase")
early_stopping <- parse_bool(get_arg("early-stopping", "false"), "early-stopping")
patience <- parse_positive_integer(get_arg("patience", "2"), "patience")
stagnation_window <- parse_positive_integer(get_arg("stagnation-window", "5"), "stagnation-window")
min_delta <- parse_nonnegative_number(get_arg("min-delta", "1e-5"), "min-delta")
objective_rel_tol <- parse_nonnegative_number(get_arg("objective-rel-tol", "1e-6"), "objective-rel-tol")
objective_abs_tol <- parse_nonnegative_number(get_arg("objective-abs-tol", "1e-8"), "objective-abs-tol")

if (identical(marker_l2_schedule, "likelihood_normalized") && isTRUE(early_stopping)) {
  stop(
    "`early-stopping` must be false for `marker-l2-schedule=likelihood_normalized` because the reported objective changes its penalty weight between updates.",
    call. = FALSE
  )
}

if (marker_genes_per_state * 6L > n_genes) {
  stop("`n-genes` must be at least six times `marker-genes-per-state` for the default six-state simulator.", call. = FALSE)
}

relative_to_out <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(normalizePath(OUT, mustWork = TRUE), "/")
  sub(paste0("^", prefix), "", path)
}

write_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, version = 3L)
}

method_ids <- function(sim, spec) {
  list(
    study_id = if (isTRUE(spec$fit_study_effect)) as.character(sim$cell_metadata$study_id) else NULL,
    donor_id = if (isTRUE(spec$fit_donor_effect)) as.character(sim$cell_metadata$donor_id) else NULL
  )
}

build_trajectory_z_prior <- function(mode,
                                     initializer_z,
                                     sim,
                                     prior_scale) {
  if (identical(mode, "none")) {
    return(NULL)
  }
  if (identical(mode, "empirical_donor")) {
    prior <- estimate_group_logistic_normal_prior(
      z = initializer_z,
      group_id = as.character(sim$cell_metadata$donor_id),
      covariance_ridge = prior_scale
    )
    prior$prior_origin <- "empirical_marker_initializer_by_donor"
    return(prior)
  }
  if (identical(mode, "empirical_study")) {
    prior <- estimate_group_logistic_normal_prior(
      z = initializer_z,
      group_id = as.character(sim$cell_metadata$study_id),
      covariance_ridge = prior_scale
    )
    prior$prior_origin <- "empirical_marker_initializer_by_study"
    return(prior)
  }
  if (identical(mode, "initializer_cell")) {
    return(prepare_initializer_logit_anchor(
      z = initializer_z,
      logit_sd = prior_scale
    ))
  }
  if (!sim$parameters$scenario %in% c("continuous", "batch_confounding")) {
    stop(
      "`z-prior=oracle_donor` is valid only for continuous and batch_confounding simulations.",
      call. = FALSE
    )
  }
  state_names <- colnames(initializer_z)
  reference_index <- length(state_names)
  donor_means <- sim$parameters$donor_state_logit_mean[
    as.character(sim$cell_metadata$donor_id),
    state_names,
    drop = FALSE
  ]
  eta_mean <- donor_means[, seq_len(reference_index - 1L), drop = FALSE] -
    donor_means[, reference_index]
  rownames(eta_mean) <- rownames(initializer_z)
  colnames(eta_mean) <- state_names[seq_len(reference_index - 1L)]
  covariance <- sim$parameters$state_sd^2 * (
    diag(reference_index - 1L) +
      matrix(1, nrow = reference_index - 1L, ncol = reference_index - 1L)
  )
  prior <- prepare_logistic_normal_prior(
    eta_mean = eta_mean,
    n_cells = nrow(initializer_z),
    n_logits = reference_index - 1L,
    covariance = covariance,
    cell_names = rownames(initializer_z),
    logit_names = colnames(eta_mean)
  )
  prior$prior_origin <- "simulator_oracle_donor_logistic_normal"
  prior$reference_state <- state_names[reference_index]
  prior
}

state_objective_components <- function(state,
                                       sim,
                                       ids,
                                       marker_target,
                                       z_prior = NULL,
                                       marker_l2_value = marker_l2) {
  negative_loglik <- -fibrodynmix_nb_loglik(
    counts = sim$counts,
    z = state$z,
    beta = state$beta,
    alpha = state$alpha,
    phi = state$phi,
    library_size = sim$cell_metadata$library_size,
    study_effect = state$study_effect,
    donor_effect = state$donor_effect,
    study_id = ids$study_id,
    donor_id = ids$donor_id,
    return_matrix = FALSE
  )
  beta_penalty <- beta_l2 * sum(state$beta^2)
  marker_penalty <- marker_l2_value * sum((state$beta - marker_target)^2)
  study_penalty <- if (is.null(state$study_effect)) 0 else study_l2 * sum(state$study_effect^2)
  donor_penalty <- if (is.null(state$donor_effect)) 0 else donor_l2 * sum(state$donor_effect^2)
  z_penalty <- z_l2 * sum(vapply(
    seq_len(nrow(state$z)),
    function(i) sum(simplex_to_logits(state$z[i, ])^2),
    numeric(1)
  ))
  z_prior_penalty <- if (is.null(z_prior)) {
    0
  } else {
    logistic_normal_prior_quadratic(
      simplex_matrix_to_reference_logits(state$z),
      z_prior
    )
  }
  legacy_total_raw <- negative_loglik + beta_penalty + marker_penalty + study_penalty + donor_penalty
  complete_total_raw <- legacy_total_raw + z_penalty + z_prior_penalty
  entries <- matrix_n_entries(sim$counts)
  data.frame(
    negative_loglik_raw = negative_loglik,
    beta_l2_penalty_raw = beta_penalty,
    marker_l2_effective = marker_l2_value,
    marker_l2_penalty_raw = marker_penalty,
    study_l2_penalty_raw = study_penalty,
    donor_l2_penalty_raw = donor_penalty,
    z_l2_penalty_raw = z_penalty,
    z_prior_penalty_raw = z_prior_penalty,
    reported_penalty_raw = beta_penalty + marker_penalty + study_penalty + donor_penalty + z_penalty + z_prior_penalty,
    legacy_objective_raw = legacy_total_raw,
    legacy_objective_average = legacy_total_raw / entries,
    reported_objective_raw = complete_total_raw,
    reported_objective_average = complete_total_raw / entries,
    complete_objective_raw = complete_total_raw,
    complete_objective_average = complete_total_raw / entries,
    stringsAsFactors = FALSE
  )
}

resolve_marker_l2 <- function(state, sim, ids, marker_target) {
  if (identical(marker_l2_schedule, "constant")) {
    return(marker_l2)
  }
  marker_error <- sum((state$beta - marker_target)^2)
  if (!is.finite(marker_error) || marker_error <= .Machine$double.eps) {
    return(marker_l2)
  }
  negative_loglik <- -fibrodynmix_nb_loglik(
    counts = sim$counts,
    z = state$z,
    beta = state$beta,
    alpha = state$alpha,
    phi = state$phi,
    library_size = sim$cell_metadata$library_size,
    study_effect = state$study_effect,
    donor_effect = state$donor_effect,
    study_id = ids$study_id,
    donor_id = ids$donor_id,
    return_matrix = FALSE
  )
  if (!is.finite(negative_loglik) || negative_loglik <= 0) {
    return(marker_l2)
  }
  marker_l2_nll_fraction * negative_loglik / marker_error
}

state_metrics <- function(state, sim, initializer_state, beta_constraint) {
  z_metrics <- evaluate_state_weights(sim$z, state$z)
  truth_beta <- sim$parameters$beta_kg
  truth_alpha <- sim$parameters$alpha_g
  if (identical(beta_constraint, "sum_to_zero")) {
    truth_canonical <- canonicalize_nb_alpha_beta(alpha = truth_alpha, beta = truth_beta)
    truth_beta <- truth_canonical$beta
    truth_alpha <- truth_canonical$alpha
  }
  study_rmse <- if (is.null(state$study_effect)) {
    NA_real_
  } else {
    sqrt(mean((state$study_effect - sim$parameters$study_effect)^2))
  }
  donor_rmse <- if (is.null(state$donor_effect)) {
    NA_real_
  } else {
    sqrt(mean((state$donor_effect - sim$parameters$donor_effect)^2))
  }
  data.frame(
    z_rmse = z_metrics$rmse,
    z_mean_absolute_error = z_metrics$mean_absolute_error,
    z_dominant_accuracy = z_metrics$dominant_accuracy,
    z_mean_entropy = z_metrics$mean_entropy_pred,
    beta_rmse_direct = sqrt(mean((state$beta - truth_beta)^2)),
    alpha_rmse_direct = sqrt(mean((as.numeric(state$alpha) - as.numeric(truth_alpha))^2)),
    log_phi_rmse_direct = sqrt(mean((log(as.numeric(state$phi)) - log(as.numeric(sim$parameters$phi_g)))^2)),
    study_effect_rmse_direct = study_rmse,
    donor_effect_rmse_direct = donor_rmse,
    z_mean_abs_delta_from_initializer = mean(abs(state$z - initializer_state$z)),
    z_rmse_delta_from_initializer = sqrt(mean((state$z - initializer_state$z)^2)),
    beta_rmse_delta_from_initializer = sqrt(mean((state$beta - initializer_state$beta)^2)),
    stringsAsFactors = FALSE
  )
}

state_snapshot <- function(state,
                           scenario,
                           replicate,
                           sim_seed,
                           initializer_seed,
                           method,
                           iteration,
                           phase,
                           components,
                           metrics,
                           update,
                           simulation_path,
                           initializer_path) {
  list(
    format_version = "1.0",
    analysis = "shared_initialization_trajectory",
    scenario = scenario,
    replicate = replicate,
    simulation_seed = sim_seed,
    initializer_seed = initializer_seed,
    method = method,
    iteration = iteration,
    phase = phase,
    state = state,
    objective_components = components,
    recovery_metrics = metrics,
    update = update,
    simulation_path = simulation_path,
    shared_initializer_path = initializer_path
  )
}

row_id <- function(scenario, replicate, method, iteration) {
  sprintf("%s__replicate_%03d__%s__iteration_%03d", scenario, replicate, method, iteration)
}

metrics_rows <- list()
components_rows <- list()
state_index_rows <- list()
shared_audit_rows <- list()
metrics_index <- 1L
components_index <- 1L
state_index <- 1L
shared_audit_index <- 1L

for (scenario in scenarios) {
  for (replicate in seq_len(n_replicates)) {
    sim_seed <- base_seed + match(scenario, valid_scenarios) * 10000L + replicate
    initializer_seed <- sim_seed + 500000L
    sim <- simulate_fibrodynmix(
      n_studies = n_studies,
      donors_per_study = donors_per_study,
      cells_per_donor = cells_per_donor,
      n_genes = n_genes,
      marker_genes_per_state = marker_genes_per_state,
      scenario = scenario,
      seed = sim_seed
    )

    replicate_data_dir <- file.path(OUT, "data", scenario, sprintf("replicate_%03d", replicate))
    simulation_path <- file.path(replicate_data_dir, "simulation.rds")
    write_rds(sim, simulation_path)

    set.seed(initializer_seed)
    initializer <- fit_fibrodynmix_initializer(
      counts = sim$counts,
      marker_index = sim$parameters$marker_index,
      library_size = sim$cell_metadata$library_size,
      n_iter = initializer_iter
    )
    marker_target_raw <- build_marker_beta_target(
      marker_index = sim$parameters$marker_index,
      gene_names = rownames(sim$counts),
      state_names = colnames(initializer$z_hat)
    )
    initial_alpha <- initialize_nb_alpha(sim$counts, sim$cell_metadata$library_size)
    initial_beta <- initializer$beta_hat
    marker_target <- marker_target_raw
    if (identical(beta_constraint, "sum_to_zero")) {
      initial_canonical <- canonicalize_nb_alpha_beta(alpha = initial_alpha, beta = initial_beta)
      initial_alpha <- initial_canonical$alpha
      initial_beta <- initial_canonical$beta
      marker_target <- center_state_program_matrix(marker_target_raw)
    }
    shared_initial_state <- list(
      z = initializer$z_hat,
      beta = initial_beta,
      alpha = initial_alpha,
      phi = estimate_phi_moments(sim$counts, sim$cell_metadata$library_size),
      study_effect = NULL,
      donor_effect = NULL
    )
    shared_dir <- file.path(OUT, "states", scenario, sprintf("replicate_%03d", replicate))
    z_prior <- build_trajectory_z_prior(
      mode = z_prior_mode,
      initializer_z = shared_initial_state$z,
      sim = sim,
      prior_scale = if (identical(z_prior_mode, "initializer_cell")) z_prior_cell_sd else z_prior_covariance_ridge
    )
    z_prior_path <- if (is.null(z_prior)) {
      NA_character_
    } else {
      file.path(shared_dir, "z_prior.rds")
    }
    if (!is.null(z_prior)) {
      write_rds(z_prior, z_prior_path)
    }
    shared_initializer_path <- file.path(shared_dir, "shared_initializer.rds")
    write_rds(
      list(
        format_version = "1.0",
        analysis = "shared_initialization_trajectory",
        scenario = scenario,
        replicate = replicate,
        simulation_seed = sim_seed,
        initializer_seed = initializer_seed,
        initializer = initializer,
        nb_initial_state = shared_initial_state,
        marker_target = marker_target,
        marker_target_raw = marker_target_raw,
        beta_constraint = beta_constraint,
        z_prior_mode = z_prior_mode,
        z_prior_covariance_ridge = z_prior_covariance_ridge,
        z_prior_path = if (is.na(z_prior_path)) NA_character_ else relative_to_out(z_prior_path),
        simulation_path = relative_to_out(simulation_path)
      ),
      shared_initializer_path
    )

    for (method in methods) {
      spec <- method_specs[[method]]
      ids <- method_ids(sim, spec)
      state <- shared_initial_state
      if (isTRUE(spec$fit_study_effect)) {
        state$study_effect <- initialize_group_effect(ids$study_id, rownames(sim$counts))
      }
      if (isTRUE(spec$fit_donor_effect)) {
        state$donor_effect <- initialize_group_effect(ids$donor_id, rownames(sim$counts))
      }

      method_dir <- file.path(shared_dir, method)
      initial_marker_l2 <- resolve_marker_l2(state, sim, ids, marker_target)
      initial_components <- state_objective_components(
        state,
        sim,
        ids,
        marker_target,
        z_prior = z_prior,
        marker_l2_value = initial_marker_l2
      )
      initial_metrics <- state_metrics(state, sim, shared_initial_state, beta_constraint)
      initial_update <- list(
        accepted = TRUE,
        rollback = FALSE,
        rollback_reason = NA_character_,
        beta_converged = NA,
        study_converged = NA,
        donor_converged = NA,
        z_converged = NA,
        z_convergence_rate = NA_real_,
        z_nonconverged_cells = NA_integer_,
        marker_l2_effective = initial_marker_l2,
        proposed_reported_objective_average = initial_components$reported_objective_average,
        best_iteration = 0L,
        stop_reason = "not_stopped"
      )
      initial_state_path <- file.path(method_dir, "iteration_000.rds")
      write_rds(
        state_snapshot(
          state = state,
          scenario = scenario,
          replicate = replicate,
          sim_seed = sim_seed,
          initializer_seed = initializer_seed,
          method = method,
          iteration = 0L,
          phase = "shared_initializer",
          components = initial_components,
          metrics = initial_metrics,
          update = initial_update,
          simulation_path = relative_to_out(simulation_path),
          initializer_path = relative_to_out(shared_initializer_path)
        ),
        initial_state_path
      )

      current_id <- row_id(scenario, replicate, method, 0L)
      metrics_rows[[metrics_index]] <- cbind(
        data.frame(
          trajectory_id = current_id,
          scenario = scenario,
          replicate = replicate,
          simulation_seed = sim_seed,
          initializer_seed = initializer_seed,
          method = method,
          iteration = 0L,
          phase = "shared_initializer",
          accepted = TRUE,
          rollback = FALSE,
          rollback_reason = NA_character_,
          beta_converged = NA,
          study_converged = NA,
          donor_converged = NA,
          z_converged = NA,
          z_convergence_rate = NA_real_,
          z_nonconverged_cells = NA_integer_,
          best_iteration = 0L,
          stop_reason = "not_stopped",
          stringsAsFactors = FALSE
        ),
        initial_metrics
      )
      metrics_index <- metrics_index + 1L
      components_rows[[components_index]] <- cbind(
        data.frame(
          trajectory_id = current_id,
          scenario = scenario,
          replicate = replicate,
          method = method,
          iteration = 0L,
          phase = "shared_initializer",
          accepted = TRUE,
          proposed_reported_objective_average = initial_components$reported_objective_average,
          proposed_complete_objective_average = initial_components$complete_objective_average,
          previous_reported_objective_average = NA_real_,
          stringsAsFactors = FALSE
        ),
        initial_components
      )
      components_index <- components_index + 1L
      state_index_rows[[state_index]] <- data.frame(
        trajectory_id = current_id,
        scenario = scenario,
        replicate = replicate,
        method = method,
        iteration = 0L,
        phase = "shared_initializer",
        state_path = relative_to_out(initial_state_path),
        state_bytes = file.info(initial_state_path)$size,
        stringsAsFactors = FALSE
      )
      state_index <- state_index + 1L
      shared_audit_rows[[shared_audit_index]] <- data.frame(
        scenario = scenario,
        replicate = replicate,
        method = method,
        simulation_path = relative_to_out(simulation_path),
        shared_initializer_path = relative_to_out(shared_initializer_path),
        z_identical_to_shared_initializer = identical(state$z, shared_initial_state$z),
        beta_identical_to_shared_initializer = identical(state$beta, shared_initial_state$beta),
        alpha_identical_to_shared_initializer = identical(state$alpha, shared_initial_state$alpha),
        phi_identical_to_shared_initializer = identical(state$phi, shared_initial_state$phi),
        z_prior_mode = z_prior_mode,
        z_prior_path = if (is.na(z_prior_path)) NA_character_ else relative_to_out(z_prior_path),
        stringsAsFactors = FALSE
      )
      shared_audit_index <- shared_audit_index + 1L

      best_objective <- initial_components$reported_objective_average
      best_iteration <- 0L
      stale_iterations <- 0L
      objective_history <- initial_components$reported_objective_average
      stop_reason <- "not_stopped"

      for (iteration in seq_len(n_outer)) {
        previous_state <- state
        effective_marker_l2 <- resolve_marker_l2(state, sim, ids, marker_target)
        previous_components <- state_objective_components(
          state,
          sim,
          ids,
          marker_target,
          z_prior = z_prior,
          marker_l2_value = effective_marker_l2
        )

        beta_update <- update_alpha_beta_nb(
          counts = sim$counts,
          z = state$z,
          beta = state$beta,
          alpha = state$alpha,
          phi = state$phi,
          library_size = sim$cell_metadata$library_size,
          study_effect = state$study_effect,
          study_id = ids$study_id,
          donor_effect = state$donor_effect,
          donor_id = ids$donor_id,
          beta_l2 = beta_l2,
          marker_target = marker_target,
          marker_l2 = effective_marker_l2,
          beta_constraint = beta_constraint,
          optimizer = "BFGS",
          optimizer_control = list(),
          maxit = maxit_beta
        )
        state$alpha <- beta_update$alpha
        state$beta <- beta_update$beta

        study_update <- NULL
        if (isTRUE(spec$fit_study_effect)) {
          study_update <- update_group_effect_nb(
            counts = sim$counts,
            z = state$z,
            beta = state$beta,
            alpha = state$alpha,
            phi = state$phi,
            library_size = sim$cell_metadata$library_size,
            effect_id = ids$study_id,
            effect = state$study_effect,
            effect_l2 = study_l2,
            other_effect = state$donor_effect,
            other_id = ids$donor_id,
            optimizer = "BFGS",
            optimizer_control = list(),
            maxit = maxit_beta
          )
          state$study_effect <- study_update$effect
        }

        donor_update <- NULL
        if (isTRUE(spec$fit_donor_effect)) {
          donor_update <- update_group_effect_nb(
            counts = sim$counts,
            z = state$z,
            beta = state$beta,
            alpha = state$alpha,
            phi = state$phi,
            library_size = sim$cell_metadata$library_size,
            effect_id = ids$donor_id,
            effect = state$donor_effect,
            effect_l2 = donor_l2,
            other_effect = state$study_effect,
            other_id = ids$study_id,
            optimizer = "BFGS",
            optimizer_control = list(),
            maxit = maxit_beta
          )
          state$donor_effect <- donor_update$effect
        }

        z_update <- update_z_nb(
          counts = sim$counts,
          z = state$z,
          beta = state$beta,
          alpha = state$alpha,
          phi = state$phi,
          library_size = sim$cell_metadata$library_size,
          study_effect = state$study_effect,
          study_id = ids$study_id,
          donor_effect = state$donor_effect,
          donor_id = ids$donor_id,
          z_l2 = z_l2,
          z_prior = z_prior,
          optimizer = "BFGS",
          optimizer_control = list(),
          maxit = maxit_z
        )
        state$z <- z_update$z

        if (isTRUE(estimate_phi)) {
          mu <- fibrodynmix_nb_mu_public(
            counts = sim$counts,
            z = state$z,
            beta = state$beta,
            alpha = state$alpha,
            library_size = sim$cell_metadata$library_size,
            study_effect = state$study_effect,
            study_id = ids$study_id,
            donor_effect = state$donor_effect,
            donor_id = ids$donor_id
          )
          state$phi <- estimate_phi_moments(
            counts = sim$counts,
            library_size = sim$cell_metadata$library_size,
            mu = mu
          )
        }

        proposed_components <- state_objective_components(
          state,
          sim,
          ids,
          marker_target,
          z_prior = z_prior,
          marker_l2_value = effective_marker_l2
        )
        rollback <- FALSE
        rollback_reason <- NA_character_
        if (!is.finite(proposed_components$reported_objective_average)) {
          state <- previous_state
          rollback <- TRUE
          rollback_reason <- "non_finite_objective"
          stop_reason <- "non_finite_objective"
        } else if (isTRUE(rollback_on_increase) &&
                   proposed_components$reported_objective_average > previous_components$reported_objective_average + min_delta) {
          state <- previous_state
          rollback <- TRUE
          rollback_reason <- "objective_increase"
        }

        accepted_components <- state_objective_components(
          state,
          sim,
          ids,
          marker_target,
          z_prior = z_prior,
          marker_l2_value = effective_marker_l2
        )
        accepted_metrics <- state_metrics(state, sim, shared_initial_state, beta_constraint)
        if (accepted_components$reported_objective_average < best_objective - min_delta) {
          best_objective <- accepted_components$reported_objective_average
          best_iteration <- iteration
          stale_iterations <- 0L
        } else {
          stale_iterations <- stale_iterations + 1L
        }
        objective_history <- c(objective_history, accepted_components$reported_objective_average)

        if (isTRUE(early_stopping) &&
            !identical(stop_reason, "non_finite_objective") &&
            iteration >= stagnation_window) {
          window_start <- iteration + 1L - stagnation_window
          objective_window <- objective_history[window_start:length(objective_history)]
          window_gain <- objective_window[1L] - objective_window[length(objective_window)]
          window_relative_gain <- window_gain / max(abs(objective_window[1L]), .Machine$double.eps)
          if (window_gain <= objective_abs_tol || window_relative_gain <= objective_rel_tol) {
            stale_iterations <- stale_iterations + 1L
          } else {
            stale_iterations <- 0L
          }
          if (stale_iterations >= patience) {
            stop_reason <- "early_stopping"
          }
        }
        if (identical(stop_reason, "not_stopped") && iteration == n_outer) {
          stop_reason <- "max_iterations"
        }

        update_info <- list(
          accepted = !rollback,
          rollback = rollback,
          rollback_reason = rollback_reason,
          beta_converged = beta_update$converged,
          study_converged = if (is.null(study_update)) NA else study_update$converged,
          donor_converged = if (is.null(donor_update)) NA else donor_update$converged,
          z_converged = z_update$converged,
          z_convergence_rate = z_update$convergence_rate,
          z_nonconverged_cells = z_update$n_nonconverged,
          marker_l2_effective = effective_marker_l2,
          proposed_reported_objective_average = proposed_components$reported_objective_average,
          proposed_complete_objective_average = proposed_components$complete_objective_average,
          previous_reported_objective_average = previous_components$reported_objective_average,
          best_iteration = best_iteration,
          best_reported_objective_average = best_objective,
          stop_reason = stop_reason,
          z_cell_diagnostics = z_update$cell_diagnostics
        )
        state_path <- file.path(method_dir, sprintf("iteration_%03d.rds", iteration))
        write_rds(
          state_snapshot(
            state = state,
            scenario = scenario,
            replicate = replicate,
            sim_seed = sim_seed,
            initializer_seed = initializer_seed,
            method = method,
            iteration = iteration,
            phase = "outer_update",
            components = accepted_components,
            metrics = accepted_metrics,
            update = update_info,
            simulation_path = relative_to_out(simulation_path),
            initializer_path = relative_to_out(shared_initializer_path)
          ),
          state_path
        )

        current_id <- row_id(scenario, replicate, method, iteration)
        metrics_rows[[metrics_index]] <- cbind(
          data.frame(
            trajectory_id = current_id,
            scenario = scenario,
            replicate = replicate,
            simulation_seed = sim_seed,
            initializer_seed = initializer_seed,
            method = method,
            iteration = iteration,
            phase = "outer_update",
            accepted = !rollback,
            rollback = rollback,
            rollback_reason = rollback_reason,
            beta_converged = beta_update$converged,
            study_converged = if (is.null(study_update)) NA else study_update$converged,
            donor_converged = if (is.null(donor_update)) NA else donor_update$converged,
            z_converged = z_update$converged,
            z_convergence_rate = z_update$convergence_rate,
            z_nonconverged_cells = z_update$n_nonconverged,
            best_iteration = best_iteration,
            stop_reason = stop_reason,
            stringsAsFactors = FALSE
          ),
          accepted_metrics
        )
        metrics_index <- metrics_index + 1L
        components_rows[[components_index]] <- cbind(
          data.frame(
            trajectory_id = current_id,
            scenario = scenario,
            replicate = replicate,
            method = method,
            iteration = iteration,
            phase = "outer_update",
            accepted = !rollback,
            proposed_reported_objective_average = proposed_components$reported_objective_average,
            proposed_complete_objective_average = proposed_components$complete_objective_average,
            previous_reported_objective_average = previous_components$reported_objective_average,
            stringsAsFactors = FALSE
          ),
          accepted_components
        )
        components_index <- components_index + 1L
        state_index_rows[[state_index]] <- data.frame(
          trajectory_id = current_id,
          scenario = scenario,
          replicate = replicate,
          method = method,
          iteration = iteration,
          phase = "outer_update",
          state_path = relative_to_out(state_path),
          state_bytes = file.info(state_path)$size,
          stringsAsFactors = FALSE
        )
        state_index <- state_index + 1L

        if (identical(stop_reason, "non_finite_objective") || identical(stop_reason, "early_stopping")) {
          break
        }
      }
    }
  }
}

trajectory_metrics <- do.call(rbind, metrics_rows)
objective_components <- do.call(rbind, components_rows)
state_index_table <- do.call(rbind, state_index_rows)
shared_audit <- do.call(rbind, shared_audit_rows)

summary_keys <- unique(trajectory_metrics[, c("scenario", "method", "iteration"), drop = FALSE])
summary_rows <- lapply(seq_len(nrow(summary_keys)), function(i) {
  keep <- trajectory_metrics$scenario == summary_keys$scenario[i] &
    trajectory_metrics$method == summary_keys$method[i] &
    trajectory_metrics$iteration == summary_keys$iteration[i]
  metric_slice <- trajectory_metrics[keep, , drop = FALSE]
  component_slice <- objective_components[
    objective_components$scenario == summary_keys$scenario[i] &
      objective_components$method == summary_keys$method[i] &
      objective_components$iteration == summary_keys$iteration[i],
    ,
    drop = FALSE
  ]
  data.frame(
    scenario = summary_keys$scenario[i],
    method = summary_keys$method[i],
    iteration = summary_keys$iteration[i],
    n_replicates_observed = nrow(metric_slice),
    z_rmse_mean = safe_mean(metric_slice$z_rmse),
    z_rmse_sd = safe_sd(metric_slice$z_rmse),
    z_dominant_accuracy_mean = safe_mean(metric_slice$z_dominant_accuracy),
    beta_rmse_direct_mean = safe_mean(metric_slice$beta_rmse_direct),
    z_delta_from_initializer_mean = safe_mean(metric_slice$z_rmse_delta_from_initializer),
    reported_objective_average_mean = safe_mean(component_slice$reported_objective_average),
    reported_objective_average_sd = safe_sd(component_slice$reported_objective_average),
    complete_objective_average_mean = safe_mean(component_slice$complete_objective_average),
    rollback_rate = mean(metric_slice$rollback),
    z_convergence_rate_mean = safe_mean(metric_slice$z_convergence_rate),
    stringsAsFactors = FALSE
  )
})
trajectory_summary <- do.call(rbind, summary_rows)

write_tsv(trajectory_metrics, file.path(OUT, "trajectory_metrics.tsv"))
write_tsv(objective_components, file.path(OUT, "objective_components.tsv"))
write_tsv(state_index_table, file.path(OUT, "state_index.tsv"))
write_tsv(shared_audit, file.path(OUT, "shared_initialization_audit.tsv"))
write_tsv(trajectory_summary, file.path(OUT, "trajectory_summary.tsv"))

config <- data.frame(
  key = c(
    "analysis", "scenarios", "methods", "n_replicates", "n_outer", "base_seed",
    "n_studies", "donors_per_study", "cells_per_donor", "n_genes", "marker_genes_per_state",
    "initializer_iter", "maxit_beta", "maxit_z", "beta_l2", "marker_l2", "marker_l2_schedule", "marker_l2_nll_fraction", "beta_constraint", "study_l2",
    "donor_l2", "z_l2", "z_prior_mode", "z_prior_covariance_ridge", "z_prior_cell_sd", "estimate_phi", "rollback_on_increase", "early_stopping",
    "patience", "stagnation_window", "min_delta", "objective_rel_tol", "objective_abs_tol",
    "state_snapshot_policy", "simulation_snapshot_policy", "r_version", "generated_at"
  ),
  value = c(
    "shared_initialization_trajectory", paste(scenarios, collapse = ";"), paste(methods, collapse = ";"),
    n_replicates, n_outer, base_seed, n_studies, donors_per_study, cells_per_donor, n_genes,
    marker_genes_per_state, initializer_iter, maxit_beta, maxit_z, beta_l2, marker_l2, marker_l2_schedule, marker_l2_nll_fraction, beta_constraint, study_l2,
    donor_l2, z_l2, z_prior_mode, z_prior_covariance_ridge, z_prior_cell_sd, estimate_phi, rollback_on_increase, early_stopping, patience, stagnation_window,
    min_delta, objective_rel_tol, objective_abs_tol,
    "one RDS at iteration 0 and after each outer update; state contains z, beta, alpha, phi, and fitted effects",
    "one simulation.rds and one shared_initializer.rds per scenario/replicate",
    R.version.string, as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)
write_tsv(config, file.path(OUT, "run_config.tsv"))

manifest <- data.frame(
  analysis = "shared_initialization_trajectory",
  primary_claim = "Each method follows one same-data, same-initialization NB trajectory rather than a grid of independently initialized fits.",
  claim_boundary = "This audit records within-run optimization behavior. It does not by itself identify the cause of any recovery degradation.",
  n_scenarios = length(scenarios),
  n_replicates = n_replicates,
  n_methods = length(methods),
  n_outer_requested = n_outer,
  marker_l2_schedule = marker_l2_schedule,
  marker_l2_nll_fraction = marker_l2_nll_fraction,
  beta_constraint = beta_constraint,
  z_prior_mode = z_prior_mode,
  z_prior_covariance_ridge = z_prior_covariance_ridge,
  z_prior_cell_sd = z_prior_cell_sd,
  n_metric_rows = nrow(trajectory_metrics),
  n_objective_rows = nrow(objective_components),
  n_state_snapshots = nrow(state_index_table),
  n_shared_initializer_audits = nrow(shared_audit),
  all_initial_z_identical = all(shared_audit$z_identical_to_shared_initializer),
  all_initial_beta_identical = all(shared_audit$beta_identical_to_shared_initializer),
  all_initial_alpha_identical = all(shared_audit$alpha_identical_to_shared_initializer),
  all_initial_phi_identical = all(shared_audit$phi_identical_to_shared_initializer),
  total_rollbacks = sum(trajectory_metrics$rollback, na.rm = TRUE),
  seed = base_seed,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "shared_initialization_trajectory_manifest.tsv"))

message("Shared-initialization trajectory written to: ", normalizePath(OUT))
