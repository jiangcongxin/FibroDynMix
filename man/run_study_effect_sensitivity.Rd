\name{run_study_effect_sensitivity}
\alias{run_study_effect_sensitivity}
\title{Run study-effect penalty sensitivity benchmarks}
\usage{
run_study_effect_sensitivity(
  study_l2_grid = c(0.01, 0.05, 0.1, 0.5, 1, 5),
  marker_l2_grid = c(0.01, 0.05, 0.1),
  n_replicates = 2,
  seed = 1,
  simulation_args = list(),
  nb_args = list(),
  include_no_study = TRUE
)
}
\arguments{
\item{study_l2_grid}{Numeric vector of study-effect ridge penalties.}
\item{marker_l2_grid}{Numeric vector of marker-orientation penalties.}
\item{n_replicates}{Number of simulation replicates per parameter setting.}
\item{seed}{Optional base seed.}
\item{simulation_args}{Named list of arguments passed to \code{simulate_fibrodynmix()}.}
\item{nb_args}{Named list of additional arguments passed to \code{fit_fibrodynmix_nb()}.}
\item{include_no_study}{Whether to include a non-study-effect NB baseline for each replicate.}
}
\value{
A data frame with one row per replicate and parameter setting.
}
\description{
Evaluates how ridge penalties on study effects and marker orientation affect
FibroDynMix NB performance in a batch-confounded simulation.
}
