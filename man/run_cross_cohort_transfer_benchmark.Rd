\name{run_cross_cohort_transfer_benchmark}
\alias{run_cross_cohort_transfer_benchmark}
\title{Run leave-study-out FibroDynMix transfer benchmark}
\usage{
run_cross_cohort_transfer_benchmark(
  n_replicates = 2,
  seed = 1,
  simulation_args = list(),
  train_fit_args = list(),
  transfer_args = list(),
  holdout_studies = NULL,
  keep_results = FALSE
)
}
\arguments{
\item{n_replicates}{Number of simulation replicates.}
\item{seed}{Optional base seed.}
\item{simulation_args}{Named list passed to \code{simulate_fibrodynmix()}.}
\item{train_fit_args}{Named list passed to \code{fit_fibrodynmix_nb()} for the training studies.}
\item{transfer_args}{Named list passed to \code{fit_fibrodynmix_transfer()}.}
\item{holdout_studies}{Optional study identifiers to hold out. If \code{NULL}, every simulated study is held out once.}
\item{keep_results}{Whether to retain full simulation, train fit, and transfer fit objects as an attribute.}
}
\value{
A data frame with transfer metrics and transfer convergence diagnostics per
replicate and held-out study.
}
\description{
Simulates a multi-study dataset, fits FibroDynMix on all but one study, and
transfers the fitted program to the held-out study by optimizing held-out cell
state mixtures under a fixed raw-count NB program.
}
