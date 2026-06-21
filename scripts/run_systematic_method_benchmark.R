#!/usr/bin/env Rscript

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(FibroDynMix)
}

out_dir <- file.path("analysis", "systematic_method_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c("marker_scoring", "topic_nmf", "fibrodynmix_initializer", "fibrodynmix_nb")
if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)) {
  methods <- c("seurat_cluster", methods)
}

n_replicates <- as.integer(Sys.getenv("FIBRODYNMIX_BENCH_REPLICATES", "5"))
benchmark <- run_simulation_benchmark(
  scenarios = c("continuous", "discrete", "batch_confounding", "rare_transition"),
  n_replicates = n_replicates,
  seed = 20260606,
  methods = methods,
  simulation_args = list(
    n_studies = 2,
    donors_per_study = 2,
    cells_per_donor = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_CELLS_PER_DONOR", "35")),
    n_genes = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_N_GENES", "240")),
    marker_genes_per_state = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_MARKERS_PER_STATE", "8"))
  ),
  topic_nmf_args = list(n_iter = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_NMF_ITER", "50")), backend = "multiplicative_update"),
  initializer_args = list(n_iter = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_INIT_ITER", "5"))),
  nb_args = list(
    n_outer = 1,
    initializer_args = list(n_iter = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_NB_INIT_ITER", "4"))),
    maxit_beta = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_MAXIT_BETA", "10")),
    maxit_z = as.integer(Sys.getenv("FIBRODYNMIX_BENCH_MAXIT_Z", "10")),
    early_stopping = TRUE
  ),
  seurat_cluster_args = list(n_pcs = 10, resolution = 0.4)
)
summary <- summarize_benchmark_results(benchmark)
optimizer <- summarize_optimizer_diagnostics(benchmark)

utils::write.table(benchmark, file.path(out_dir, "method_benchmark_replicates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(summary, file.path(out_dir, "method_benchmark_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(optimizer, file.path(out_dir, "method_benchmark_optimizer_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

print(summary)
