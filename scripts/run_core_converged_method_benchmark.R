#!/usr/bin/env Rscript

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

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_core_converged_method_benchmark.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- get_arg("out", file.path(ROOT, "analysis", "core_converged_method_benchmark"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  "matrix_utils.R",
  "simulate_fibrodynmix.R",
  "benchmark_metrics.R",
  "baseline_marker_scoring.R",
  "topic_nmf_baseline.R",
  "fibrodynmix_initializer.R",
  "nb_likelihood.R",
  "fit_nb_model.R",
  "bootstrap_uncertainty.R",
  "vi_posterior.R",
  "scvi_baseline.R",
  "simulation_benchmark.R"
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
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  if (length(out) == 0L) {
    stop("Character vector arguments cannot be empty.", call. = FALSE)
  }
  out
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
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

compute_paired_tests <- function(benchmark, metrics) {
  if (nrow(benchmark) == 0L) {
    return(data.frame())
  }
  required <- c("scenario", "replicate", "method", "nb_n_outer", metrics)
  missing <- setdiff(required, colnames(benchmark))
  if (length(missing) > 0L) {
    stop(
      sprintf("Cannot compute paired tests; missing columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  groups <- unique(benchmark[, c("nb_n_outer", "scenario"), drop = FALSE])
  rows <- list()
  row_idx <- 1L

  for (group_idx in seq_len(nrow(groups))) {
    subset <- benchmark[
      benchmark$nb_n_outer == groups$nb_n_outer[group_idx] &
        benchmark$scenario == groups$scenario[group_idx],
      ,
      drop = FALSE
    ]
    method_names <- sort(unique(subset$method))
    if (length(method_names) < 2L) {
      next
    }
    method_pairs <- utils::combn(method_names, 2L, simplify = FALSE)
    for (pair in method_pairs) {
      left <- subset[subset$method == pair[[1L]], , drop = FALSE]
      right <- subset[subset$method == pair[[2L]], , drop = FALSE]
      for (metric in metrics) {
        paired <- merge(
          left[, c("replicate", metric), drop = FALSE],
          right[, c("replicate", metric), drop = FALSE],
          by = "replicate",
          suffixes = c("_method_a", "_method_b")
        )
        metric_a <- paste0(metric, "_method_a")
        metric_b <- paste0(metric, "_method_b")
        keep <- is.finite(paired[[metric_a]]) & is.finite(paired[[metric_b]])
        paired <- paired[keep, , drop = FALSE]
        delta <- paired[[metric_b]] - paired[[metric_a]]
        wilcox_p <- NA_real_
        paired_t_p <- NA_real_
        if (nrow(paired) >= 2L && any(abs(delta) > .Machine$double.eps)) {
          wilcox_p <- tryCatch(
            stats::wilcox.test(
              paired[[metric_b]],
              paired[[metric_a]],
              paired = TRUE,
              exact = FALSE
            )$p.value,
            error = function(e) NA_real_
          )
          paired_t_p <- tryCatch(
            stats::t.test(
              paired[[metric_b]],
              paired[[metric_a]],
              paired = TRUE
            )$p.value,
            error = function(e) NA_real_
          )
        }
        rows[[row_idx]] <- data.frame(
          nb_n_outer = groups$nb_n_outer[group_idx],
          scenario = groups$scenario[group_idx],
          metric = metric,
          metric_direction = if (metric %in% c("rmse", "mean_absolute_error")) "lower_is_better" else "higher_is_better",
          method_a = pair[[1L]],
          method_b = pair[[2L]],
          n_pairs = nrow(paired),
          method_a_mean = safe_mean(paired[[metric_a]]),
          method_b_mean = safe_mean(paired[[metric_b]]),
          mean_delta_method_b_minus_a = safe_mean(delta),
          median_delta_method_b_minus_a = if (length(delta) == 0L) NA_real_ else stats::median(delta, na.rm = TRUE),
          wilcox_p = wilcox_p,
          paired_t_p = paired_t_p,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

n_outer_grid <- parse_int_vector(get_arg("n-outer-grid"), c(10L, 20L))
scenarios <- parse_char_vector(
  get_arg("scenarios"),
  c("continuous", "discrete", "batch_confounding", "rare_transition")
)
n_replicates <- as.integer(get_arg("n-replicates", "2"))
seed <- as.integer(get_arg("seed", "1060"))
n_studies <- as.integer(get_arg("n-studies", "2"))
donors_per_study <- as.integer(get_arg("donors-per-study", "2"))
cells_per_donor <- as.integer(get_arg("cells-per-donor", "8"))
n_genes <- as.integer(get_arg("n-genes", "90"))
marker_genes_per_state <- as.integer(get_arg("marker-genes-per-state", "4"))
initializer_iter <- as.integer(get_arg("initializer-iter", "3"))
maxit_beta <- as.integer(get_arg("maxit-beta", "12"))
maxit_z <- as.integer(get_arg("maxit-z", "12"))
study_l2 <- as.numeric(get_arg("study-l2", "5"))
marker_l2 <- as.numeric(get_arg("marker-l2", "0.05"))
stagnation_window <- as.integer(get_arg("stagnation-window", "5"))
patience <- as.integer(get_arg("patience", "2"))
objective_rel_tol <- as.numeric(get_arg("objective-rel-tol", "1e-6"))
objective_abs_tol <- as.numeric(get_arg("objective-abs-tol", "1e-8"))

methods <- parse_char_vector(
  get_arg("methods"),
  c("marker_scoring", "topic_nmf", "fibrodynmix_initializer", "fibrodynmix_nb", "fibrodynmix_nb_study")
)
if (has_flag("include-scvi") && !"scvi_latent" %in% methods) {
  methods <- c(methods, "scvi_latent")
}
if (has_flag("include-vi") && !"fibrodynmix_vi" %in% methods) {
  methods <- c(methods, "fibrodynmix_vi")
}
scvi_max_epochs <- as.integer(get_arg("scvi-max-epochs", "40"))

topic_backend <- if (requireNamespace("NMF", quietly = TRUE)) "nmf" else "multiplicative_update"
outer_sensitive_methods <- c(
  "fibrodynmix_nb",
  "fibrodynmix_nb_study",
  "fibrodynmix_nb_donor",
  "fibrodynmix_nb_study_donor",
  "fibrodynmix_vi"
)
outer_static_methods <- setdiff(methods, outer_sensitive_methods)
static_benchmark_cache <- NULL

benchmark_rows <- list()
optimizer_rows <- list()
summary_rows <- list()
row_idx <- 1L
opt_idx <- 1L
sum_idx <- 1L

for (outer_idx in seq_along(n_outer_grid)) {
  n_outer <- n_outer_grid[[outer_idx]]
  loop_methods <- if (outer_idx == 1L) {
    methods
  } else {
    intersect(methods, outer_sensitive_methods)
  }
  benchmark <- data.frame()
  if (length(loop_methods) > 0L) {
    benchmark <- run_simulation_benchmark(
      scenarios = scenarios,
      n_replicates = n_replicates,
      seed = seed,
      methods = loop_methods,
      simulation_args = list(
        n_studies = n_studies,
        donors_per_study = donors_per_study,
        cells_per_donor = cells_per_donor,
        n_genes = n_genes,
        marker_genes_per_state = marker_genes_per_state
      ),
      initializer_args = list(n_iter = initializer_iter),
      topic_nmf_args = list(
        backend = topic_backend,
        n_iter = 40
      ),
      scvi_args = list(max_epochs = scvi_max_epochs),
      nb_args = list(
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
      ),
      vi_args = list(
        n_draws = 18,
        n_elbo_draws = 4,
        n_vi_iter = 2,
        seed = seed + 1L,
        keep_draws = FALSE
      ),
      downstream_task = TRUE,
      downstream_label_col = "disease",
      downstream_group_col = "donor_id",
      downstream_folds = 3
    )
    benchmark$nb_n_outer <- n_outer
  }
  if (outer_idx == 1L && length(outer_static_methods) > 0L && nrow(benchmark) > 0L) {
    static_benchmark_cache <- benchmark[benchmark$method %in% outer_static_methods, , drop = FALSE]
  }
  if (outer_idx > 1L && length(outer_static_methods) > 0L && !is.null(static_benchmark_cache)) {
    static_rows <- static_benchmark_cache
    static_rows$nb_n_outer <- n_outer
    benchmark <- rbind(benchmark, static_rows)
  }
  benchmark_rows[[row_idx]] <- benchmark
  row_idx <- row_idx + 1L

  optimizer <- summarize_optimizer_diagnostics(benchmark)
  optimizer$nb_n_outer <- n_outer
  optimizer_rows[[opt_idx]] <- optimizer
  opt_idx <- opt_idx + 1L

  summary <- summarize_benchmark_results(
    benchmark,
    metrics = c(
      "rmse",
      "mean_absolute_error",
      "dominant_accuracy",
      "mean_entropy_pred",
      "downstream_balanced_accuracy",
      "downstream_macro_f1",
      "downstream_macro_auroc"
    )
  )
  summary$nb_n_outer <- n_outer
  summary_rows[[sum_idx]] <- summary
  sum_idx <- sum_idx + 1L
}

benchmark_all <- do.call(rbind, benchmark_rows)
summary_all <- do.call(rbind, summary_rows)
optimizer_all <- do.call(rbind, optimizer_rows)

write_tsv(benchmark_all, file.path(OUT, "core_converged_benchmark_metrics.tsv"))
write_tsv(summary_all, file.path(OUT, "core_converged_benchmark_summary.tsv"))
write_tsv(optimizer_all, file.path(OUT, "core_converged_optimizer_diagnostics.tsv"))

rankings <- do.call(rbind, lapply(split(summary_all, paste(summary_all$nb_n_outer, summary_all$scenario, sep = "|")), function(df) {
  df <- df[order(df$rmse_mean), , drop = FALSE]
  df$rmse_rank <- seq_len(nrow(df))
  df
}))
rankings <- rankings[, c("nb_n_outer", "scenario", "rmse_rank", "method", "rmse_mean", "dominant_accuracy_mean", "downstream_balanced_accuracy_mean", "downstream_macro_f1_mean"), drop = FALSE]
write_tsv(rankings, file.path(OUT, "core_converged_benchmark_rankings.tsv"))

nb_methods <- intersect(methods, c("fibrodynmix_nb", "fibrodynmix_nb_study", "fibrodynmix_vi"))
tradeoff_rows <- list()
tradeoff_idx <- 1L
if (length(n_outer_grid) >= 2L) {
  sorted_outer <- sort(n_outer_grid)
  for (method in nb_methods) {
    for (scenario in scenarios) {
      for (replicate in seq_len(n_replicates)) {
        df <- benchmark_all[
          benchmark_all$method == method &
            benchmark_all$scenario == scenario &
            benchmark_all$replicate == replicate,
          ,
          drop = FALSE
        ]
        df <- df[order(df$nb_n_outer), , drop = FALSE]
        if (nrow(df) < 2L) {
          next
        }
        for (i in seq_len(nrow(df) - 1L)) {
          current <- df[i, , drop = FALSE]
          next_row <- df[i + 1L, , drop = FALSE]
          tradeoff_rows[[tradeoff_idx]] <- data.frame(
            scenario = scenario,
            replicate = replicate,
            method = method,
            from_n_outer = current$nb_n_outer,
            to_n_outer = next_row$nb_n_outer,
            best_objective_delta = next_row$nb_best_objective - current$nb_best_objective,
            rmse_delta = next_row$rmse - current$rmse,
            dominant_accuracy_delta = next_row$dominant_accuracy - current$dominant_accuracy,
            downstream_balanced_accuracy_delta = next_row$downstream_balanced_accuracy - current$downstream_balanced_accuracy,
            lower_objective_worse_rmse = is.finite(next_row$nb_best_objective) &&
              next_row$nb_best_objective < current$nb_best_objective &&
              next_row$rmse > current$rmse,
            stringsAsFactors = FALSE
          )
          tradeoff_idx <- tradeoff_idx + 1L
        }
      }
    }
  }
}
tradeoff <- if (length(tradeoff_rows) == 0L) data.frame() else do.call(rbind, tradeoff_rows)
write_tsv(tradeoff, file.path(OUT, "core_converged_objective_rmse_tradeoff.tsv"))

paired_tests <- compute_paired_tests(
  benchmark_all,
  metrics = c("rmse", "dominant_accuracy", "downstream_balanced_accuracy", "downstream_macro_f1")
)
write_tsv(paired_tests, file.path(OUT, "core_converged_paired_tests.tsv"))

fibro_rows <- benchmark_all[benchmark_all$method %in% nb_methods, , drop = FALSE]
baseline_rows <- benchmark_all[benchmark_all$method %in% c("marker_scoring", "topic_nmf"), , drop = FALSE]
manifest <- data.frame(
  analysis = "core_converged_method_benchmark",
  primary_claim = "Core simulation methods were rerun under higher NB outer-iteration budgets to separate smoke settings from convergence-candidate settings.",
  claim_boundary = "Bounded simulation benchmark with small sampled datasets. n_outer=10/20 results test optimizer sensitivity but do not by themselves solve NB objective/model-selection mismatch.",
  scenarios = paste(scenarios, collapse = ";"),
  methods = paste(methods, collapse = ";"),
  n_outer_grid = paste(n_outer_grid, collapse = ";"),
  n_replicates = n_replicates,
  n_rows = nrow(benchmark_all),
  topic_backend = topic_backend,
  include_scvi = "scvi_latent" %in% methods,
  scvi_max_epochs = if ("scvi_latent" %in% methods) scvi_max_epochs else NA_integer_,
  include_vi = "fibrodynmix_vi" %in% methods,
  mean_fibrodynmix_rmse = mean(fibro_rows$rmse, na.rm = TRUE),
  mean_baseline_rmse = mean(baseline_rows$rmse, na.rm = TRUE),
  mean_fibrodynmix_downstream_balanced_accuracy = mean(fibro_rows$downstream_balanced_accuracy, na.rm = TRUE),
  mean_baseline_downstream_balanced_accuracy = mean(baseline_rows$downstream_balanced_accuracy, na.rm = TRUE),
  n_lower_objective_worse_rmse_pairs = if (nrow(tradeoff) == 0L) 0L else sum(tradeoff$lower_objective_worse_rmse, na.rm = TRUE),
  n_objective_rmse_tradeoff_pairs = nrow(tradeoff),
  n_paired_test_rows = nrow(paired_tests),
  seed = seed,
  generated_at = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUT, "core_converged_benchmark_manifest.tsv"))

message("Core converged method benchmark written to: ", normalizePath(OUT))
