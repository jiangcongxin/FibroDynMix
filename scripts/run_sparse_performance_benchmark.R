#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(FibroDynMix)
}

out_dir <- file.path("analysis", "performance_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260606)
n_genes <- as.integer(Sys.getenv("FIBRODYNMIX_BENCH_GENES", "800"))
n_cells <- as.integer(Sys.getenv("FIBRODYNMIX_BENCH_CELLS", "1200"))
n_states <- 4L
density <- as.numeric(Sys.getenv("FIBRODYNMIX_BENCH_DENSITY", "0.04"))

counts_sparse <- Matrix::rsparsematrix(
  nrow = n_genes,
  ncol = n_cells,
  density = density,
  rand.x = function(n) stats::rpois(n, lambda = 2) + 1
)
counts_sparse <- as(counts_sparse, "dgCMatrix")
rownames(counts_sparse) <- paste0("gene", seq_len(n_genes))
colnames(counts_sparse) <- paste0("cell", seq_len(n_cells))
counts_dense <- as.matrix(counts_sparse)
library_size <- as.numeric(Matrix::colSums(counts_sparse))

marker_index <- lapply(seq_len(n_states), function(i) {
  paste0("gene", seq.int((i - 1L) * 15L + 1L, i * 15L))
})
names(marker_index) <- paste0("state", seq_len(n_states))
cell_metadata <- data.frame(
  cell_id = colnames(counts_sparse),
  sample_id = rep(c("sample_a", "sample_b"), length.out = n_cells),
  row.names = colnames(counts_sparse),
  stringsAsFactors = FALSE
)

prepared <- prepare_fibrodynmix_data(
  counts = counts_sparse,
  cell_metadata = cell_metadata,
  marker_index = marker_index,
  cell_id_col = "cell_id",
  donor_col = "sample_id",
  min_cells_per_gene = 0,
  min_counts_per_gene = 0
)
fit <- fit_fibrodynmix_prepared(
  prepared,
  n_outer = 1,
  initializer_args = list(n_iter = 1),
  maxit_beta = 8,
  maxit_z = 10
)
warm_start <- score_marker_baseline(
  counts = prepared$counts,
  marker_index = prepared$marker_index,
  library_size = prepared$library_size
)$z_pred

bench_one <- function(label, expr) {
  gc()
  elapsed <- system.time(force(expr))[["elapsed"]]
  data.frame(label = label, elapsed_sec = elapsed, stringsAsFactors = FALSE)
}

results <- rbind(
  bench_one("marker_baseline_sparse", score_marker_baseline(counts_sparse, marker_index, library_size)),
  bench_one("marker_baseline_dense", score_marker_baseline(counts_dense, marker_index, library_size)),
  bench_one("nb_objective_sparse", fibrodynmix_nb_objective(counts_sparse, fit$z_hat, fit$beta_hat, fit$alpha_hat, fit$phi_hat, library_size = library_size)),
  bench_one("nb_objective_dense", fibrodynmix_nb_objective(counts_dense, fit$z_hat, fit$beta_hat, fit$alpha_hat, fit$phi_hat, library_size = library_size)),
  bench_one(
    "transfer_sparse_chunked",
    fit_fibrodynmix_transfer(
      counts_sparse,
      fit = fit,
      library_size = library_size,
      z_init = warm_start[, rownames(fit$beta_hat), drop = FALSE],
      chunk_size = 300,
      maxit_z = 20,
      return_cell_diagnostics = TRUE
    )
  )
)

results$n_genes <- n_genes
results$n_cells <- n_cells
results$density <- density
utils::write.table(
  results,
  file.path(out_dir, "sparse_performance_benchmark.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
print(results)
