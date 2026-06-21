\name{run_simulation_benchmark}
\alias{run_simulation_benchmark}
\title{Run simulation benchmarks across scenarios}
\usage{
run_simulation_benchmark(
  scenarios = c("continuous", "discrete", "batch_confounding", "rare_transition"),
  n_replicates = 3,
  seed = 1,
  methods = "marker_scoring",
  simulation_args = list(),
  temperature = 1,
  initializer_args = list(),
  topic_nmf_args = list(),
  scvi_args = list(),
  seurat_cluster_args = list(),
  nb_args = list(),
  vi_args = list(),
  downstream_task = TRUE,
  downstream_label_col = "disease",
  downstream_group_col = "donor_id",
  downstream_folds = 3,
  keep_results = FALSE
)
}
\arguments{
\item{scenarios}{Character vector of simulation scenarios.}
\item{n_replicates}{Number of replicates per scenario.}
\item{seed}{Optional base seed. Replicate seeds are deterministically derived from this value.}
\item{methods}{Character vector of methods to run. Implemented methods are \code{marker_scoring}, \code{seurat_cluster}, \code{topic_nmf}, \code{scvi_latent}, \code{fibrodynmix_initializer}, \code{fibrodynmix_nb}, \code{fibrodynmix_nb_study}, \code{fibrodynmix_nb_donor}, \code{fibrodynmix_nb_study_donor}, and \code{fibrodynmix_vi}.}
\item{simulation_args}{Named list of arguments passed to \code{simulate_fibrodynmix()}.}
\item{temperature}{Softmax temperature for the marker-scoring baseline.}
\item{initializer_args}{Named list of arguments passed to \code{fit_fibrodynmix_initializer()}.}
\item{topic_nmf_args}{Named list of arguments passed to \code{fit_topic_nmf_baseline()}.}
\item{scvi_args}{Named list of arguments passed to the optional \code{fit_scvi_latent_baseline()}.}
\item{seurat_cluster_args}{Named list of arguments passed to the Seurat cluster baseline.}
\item{nb_args}{Named list of arguments passed to \code{fit_fibrodynmix_nb()}.}
\item{vi_args}{Named list of arguments passed to \code{fit_fibrodynmix_vi()}.}
\item{downstream_task}{Whether to evaluate a downstream label-prediction benchmark from each inferred state-feature matrix.}
\item{downstream_label_col}{Column in simulated cell metadata used as the downstream prediction label.}
\item{downstream_group_col}{Optional cell metadata column used for grouped cross-validation to avoid donor leakage.}
\item{downstream_folds}{Requested number of downstream CV folds.}
\item{keep_results}{Whether to retain full simulation and baseline objects.}
}
\value{
A data frame with one row per scenario, replicate, and method. If
\code{keep_results = TRUE}, full run objects are attached as the \code{results}
attribute. The table includes NB optimizer diagnostics for NB methods and
\code{NA} diagnostics for non-NB methods.
}
\description{
Runs repeated FibroDynMix simulations and evaluates the marker-scoring baseline
across one or more scenarios.
}
