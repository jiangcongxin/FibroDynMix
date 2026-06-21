#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_nb_convergence_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "nb_convergence_benchmark"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "matrix_utils.R",
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R"
)
invisible(lapply(file.path(ROOT, "R", source_files), source))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
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
  out <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- trimws(out[nzchar(trimws(out))])
  if (length(out) == 0L) {
    stop("Character vector arguments cannot be empty.", call. = FALSE)
  }
  out
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

convergence_methods <- parse_char_vector(get_arg("methods"), c("fibrodynmix_nb", "fibrodynmix_nb_study"))
unknown_methods <- setdiff(convergence_methods, c("fibrodynmix_nb", "fibrodynmix_nb_study"))
if (length(unknown_methods) > 0L) {
  stop(sprintf("Unsupported convergence method(s): %s", paste(unknown_methods, collapse = ", ")), call. = FALSE)
}

scenarios <- parse_char_vector(
  get_arg("scenarios"),
  c("continuous", "discrete", "batch_confounding", "rare_transition")
)
valid_scenarios <- c("continuous", "discrete", "batch_confounding", "rare_transition")
if (length(setdiff(scenarios, valid_scenarios)) > 0L) {
  stop(sprintf("Unsupported scenario(s): %s", paste(setdiff(scenarios, valid_scenarios), collapse = ", ")), call. = FALSE)
}

n_outer_grid <- parse_int_vector(get_arg("n-outer-grid"), c(2L, 5L, 10L, 20L))
reference_n_outer <- max(n_outer_grid)
n_replicates <- as.integer(get_arg("n-replicates", "2"))
base_seed <- as.integer(get_arg("seed", "940"))
n_studies <- as.integer(get_arg("n-studies", "2"))
donors_per_study <- as.integer(get_arg("donors-per-study", "2"))
cells_per_donor <- as.integer(get_arg("cells-per-donor", "8"))
n_genes <- as.integer(get_arg("n-genes", "90"))
marker_genes_per_state <- as.integer(get_arg("marker-genes-per-state", "4"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "16"))
maxit_z <- as.integer(get_arg("maxit-z", "16"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))
stagnation_window <- as.integer(get_arg("stagnation-window", "5"))
patience <- as.integer(get_arg("patience", "2"))
objective_rel_tol <- as.numeric(get_arg("objective-rel-tol", "1e-6"))
objective_abs_tol <- as.numeric(get_arg("objective-abs-tol", "1e-8"))

fit_one <- function(sim, method, n_outer) {
  fit_fibrodynmix_nb(
    counts = sim$counts,
    marker_index = sim$parameters$marker_index,
    library_size = sim$cell_metadata$library_size,
    study_id = if (identical(method, "fibrodynmix_nb_study")) sim$cell_metadata$study_id else NULL,
    fit_study_effect = identical(method, "fibrodynmix_nb_study"),
    n_outer = n_outer,
    initializer_args = list(n_iter = initializer_iter),
    maxit_beta = maxit_beta,
    maxit_z = maxit_z,
    study_l2 = study_l2,
    marker_l2 = marker_l2,
    stagnation_window = stagnation_window,
    patience = patience,
    objective_rel_tol = objective_rel_tol,
    objective_abs_tol = objective_abs_tol
  )
}

fit_row <- function(scenario, replicate, seed, method, n_outer, sim, fit) {
  state_metrics <- evaluate_state_weights(sim$z, fit$z_hat)
  downstream <- evaluate_downstream_classification(
    features = fit$z_hat,
    labels = sim$cell_metadata$disease,
    groups = sim$cell_metadata$donor_id,
    n_folds = 3,
    seed = seed
  )
  trace <- fit$nb_objective_trace
  last_gain <- if (length(trace) >= 2L) trace[length(trace) - 1L] - trace[length(trace)] else NA_real_
  rel_last_gain <- if (is.finite(last_gain) && length(trace) >= 2L) {
    last_gain / max(abs(trace[length(trace) - 1L]), .Machine$double.eps)
  } else {
    NA_real_
  }
  data.frame(
    scenario = scenario,
    replicate = replicate,
    seed = seed,
    method = method,
    n_outer_requested = n_outer,
    n_outer_executed = fit$executed_iterations,
    stop_reason = fit$stop_reason,
    converged_subproblems = fit$converged,
    best_iteration = fit$best_iteration,
    initial_objective = trace[1L],
    final_objective = trace[length(trace)],
    best_objective = fit$best_objective,
    objective_improvement = trace[1L] - fit$best_objective,
    last_outer_objective_gain = last_gain,
    last_outer_relative_gain = rel_last_gain,
    rollback_count = sum(fit$convergence$rolled_back),
    beta_convergence_rate = mean(fit$convergence$beta),
    z_convergence_rate = mean(fit$convergence$z),
    study_convergence_rate = if (length(fit$convergence$study) == 0L) NA_real_ else mean(fit$convergence$study),
    rmse = state_metrics$rmse,
    mean_absolute_error = state_metrics$mean_absolute_error,
    dominant_accuracy = state_metrics$dominant_accuracy,
    mean_entropy_pred = state_metrics$mean_entropy_pred,
    downstream_status = downstream$status,
    downstream_balanced_accuracy = downstream$balanced_accuracy,
    downstream_macro_f1 = downstream$macro_f1,
    downstream_macro_auroc = downstream$macro_auroc,
    stringsAsFactors = FALSE
  )
}

trace_rows <- function(scenario, replicate, seed, method, n_outer, fit) {
  data.frame(
    scenario = scenario,
    replicate = replicate,
    seed = seed,
    method = method,
    n_outer_requested = n_outer,
    iteration = seq_along(fit$nb_objective_trace) - 1L,
    objective = fit$nb_objective_trace,
    stringsAsFactors = FALSE
  )
}

z_delta_row <- function(scenario, replicate, seed, method, n_outer, fit, reference_fit) {
  z <- fit$z_hat
  z_ref <- reference_fit$z_hat
  state_names <- colnames(z_ref)
  pred_state <- state_names[max.col(z, ties.method = "first")]
  ref_state <- state_names[max.col(z_ref, ties.method = "first")]
  data.frame(
    scenario = scenario,
    replicate = replicate,
    seed = seed,
    method = method,
    n_outer_requested = n_outer,
    reference_n_outer = reference_n_outer,
    best_objective_gap_vs_reference = fit$best_objective - reference_fit$best_objective,
    mean_abs_z_delta_vs_reference = mean(abs(z - z_ref)),
    rmse_z_delta_vs_reference = sqrt(mean((z - z_ref)^2)),
    dominant_agreement_vs_reference = mean(pred_state == ref_state),
    stringsAsFactors = FALSE
  )
}

metric_rows <- list()
trace_list <- list()
z_delta_rows <- list()
row_idx <- 1L
trace_idx <- 1L
delta_idx <- 1L

for (scenario in scenarios) {
  for (replicate in seq_len(n_replicates)) {
    replicate_seed <- base_seed + match(scenario, valid_scenarios) * 10000L + replicate
    sim <- simulate_fibrodynmix(
      n_studies = n_studies,
      donors_per_study = donors_per_study,
      cells_per_donor = cells_per_donor,
      n_genes = n_genes,
      marker_genes_per_state = marker_genes_per_state,
      scenario = scenario,
      seed = replicate_seed
    )

    for (method in convergence_methods) {
      fits <- list()
      for (n_outer in n_outer_grid) {
        fit <- fit_one(sim, method, n_outer)
        fits[[as.character(n_outer)]] <- fit
        metric_rows[[row_idx]] <- fit_row(scenario, replicate, replicate_seed, method, n_outer, sim, fit)
        row_idx <- row_idx + 1L
        trace_list[[trace_idx]] <- trace_rows(scenario, replicate, replicate_seed, method, n_outer, fit)
        trace_idx <- trace_idx + 1L
      }

      reference_fit <- fits[[as.character(reference_n_outer)]]
      for (n_outer in n_outer_grid) {
        z_delta_rows[[delta_idx]] <- z_delta_row(
          scenario = scenario,
          replicate = replicate,
          seed = replicate_seed,
          method = method,
          n_outer = n_outer,
          fit = fits[[as.character(n_outer)]],
          reference_fit = reference_fit
        )
        delta_idx <- delta_idx + 1L
      }
    }
  }
}

metrics <- do.call(rbind, metric_rows)
traces <- do.call(rbind, trace_list)
z_deltas <- do.call(rbind, z_delta_rows)

write_tsv(metrics, file.path(OUT, "nb_convergence_metrics.tsv"))
write_tsv(traces, file.path(OUT, "nb_convergence_objective_traces.tsv"))
write_tsv(z_deltas, file.path(OUT, "nb_convergence_z_delta_vs_reference.tsv"))

summary_keys <- unique(metrics[, c("scenario", "method", "n_outer_requested"), drop = FALSE])
summary_rows <- lapply(seq_len(nrow(summary_keys)), function(i) {
  keep <- metrics$scenario == summary_keys$scenario[i] &
    metrics$method == summary_keys$method[i] &
    metrics$n_outer_requested == summary_keys$n_outer_requested[i]
  df <- metrics[keep, , drop = FALSE]
  delta <- z_deltas[keep, , drop = FALSE]
  data.frame(
    scenario = summary_keys$scenario[i],
    method = summary_keys$method[i],
    n_outer_requested = summary_keys$n_outer_requested[i],
    n_replicates = nrow(df),
    rmse_mean = mean(df$rmse),
    rmse_sd = safe_sd(df$rmse),
    dominant_accuracy_mean = mean(df$dominant_accuracy),
    downstream_balanced_accuracy_mean = mean(df$downstream_balanced_accuracy),
    downstream_macro_f1_mean = mean(df$downstream_macro_f1),
    best_objective_mean = mean(df$best_objective),
    objective_improvement_mean = mean(df$objective_improvement),
    last_outer_relative_gain_mean = safe_mean(df$last_outer_relative_gain),
    early_stopping_rate = mean(df$stop_reason == "early_stopping"),
    mean_abs_z_delta_vs_reference = mean(delta$mean_abs_z_delta_vs_reference),
    rmse_z_delta_vs_reference = mean(delta$rmse_z_delta_vs_reference),
    dominant_agreement_vs_reference = mean(delta$dominant_agreement_vs_reference),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write_tsv(summary, file.path(OUT, "nb_convergence_summary.tsv"))

delta_from_two <- merge(
  summary,
  summary[summary$n_outer_requested == min(n_outer_grid), c("scenario", "method", "rmse_mean", "dominant_accuracy_mean", "downstream_balanced_accuracy_mean", "best_objective_mean"), drop = FALSE],
  by = c("scenario", "method"),
  suffixes = c("", "_at_min_n_outer"),
  sort = FALSE
)
delta_from_two$rmse_delta_vs_min_n_outer <- delta_from_two$rmse_mean - delta_from_two$rmse_mean_at_min_n_outer
delta_from_two$dominant_accuracy_delta_vs_min_n_outer <- delta_from_two$dominant_accuracy_mean - delta_from_two$dominant_accuracy_mean_at_min_n_outer
delta_from_two$downstream_balanced_accuracy_delta_vs_min_n_outer <- delta_from_two$downstream_balanced_accuracy_mean - delta_from_two$downstream_balanced_accuracy_mean_at_min_n_outer
delta_from_two$best_objective_delta_vs_min_n_outer <- delta_from_two$best_objective_mean - delta_from_two$best_objective_mean_at_min_n_outer
write_tsv(delta_from_two, file.path(OUT, "nb_convergence_delta_vs_min_n_outer.tsv"))

non_reference <- z_deltas[z_deltas$n_outer_requested < reference_n_outer, , drop = FALSE]
min_grid <- min(n_outer_grid)
min_delta <- z_deltas[z_deltas$n_outer_requested == min_grid, , drop = FALSE]
manifest <- data.frame(
  analysis = "nb_convergence_benchmark",
  primary_claim = "FibroDynMix NB optimization was evaluated across an outer-iteration grid to test whether n_outer=2 is sufficient for benchmark conclusions.",
  claim_boundary = "Bounded simulation convergence sensitivity, not a full production-scale rerun. If metrics or z change materially at higher n_outer, manuscript NB conclusions should use converged settings or state the limitation.",
  scenarios = paste(scenarios, collapse = ";"),
  methods = paste(convergence_methods, collapse = ";"),
  n_outer_grid = paste(n_outer_grid, collapse = ";"),
  reference_n_outer = reference_n_outer,
  n_replicates = n_replicates,
  n_rows = nrow(metrics),
  n_trace_rows = nrow(traces),
  max_abs_best_objective_gap_vs_reference_before_reference = max(abs(non_reference$best_objective_gap_vs_reference), na.rm = TRUE),
  max_mean_abs_z_delta_at_min_n_outer = max(min_delta$mean_abs_z_delta_vs_reference, na.rm = TRUE),
  min_dominant_agreement_at_min_n_outer = min(min_delta$dominant_agreement_vs_reference, na.rm = TRUE),
  max_abs_rmse_delta_vs_min_n_outer = max(abs(delta_from_two$rmse_delta_vs_min_n_outer), na.rm = TRUE),
  max_abs_downstream_balanced_accuracy_delta_vs_min_n_outer = max(abs(delta_from_two$downstream_balanced_accuracy_delta_vs_min_n_outer), na.rm = TRUE),
  any_reference_early_stopping = any(metrics$n_outer_requested == reference_n_outer & metrics$stop_reason == "early_stopping"),
  seed = base_seed,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "nb_convergence_manifest.tsv"))

message("NB convergence benchmark written to: ", normalizePath(OUT))
