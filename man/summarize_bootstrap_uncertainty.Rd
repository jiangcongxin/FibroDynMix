\name{summarize_bootstrap_uncertainty}
\alias{summarize_bootstrap_uncertainty}
\title{Summarize FibroDynMix bootstrap uncertainty draws}
\usage{
summarize_bootstrap_uncertainty(
  cell_draws,
  sample_draws,
  marker_draws,
  probs = c(0.025, 0.975)
)
}
\arguments{
\item{cell_draws}{Cell-level bootstrap state/entropy draws.}
\item{sample_draws}{Sample-level bootstrap composition draws.}
\item{marker_draws}{State-gene program bootstrap draws.}
\item{probs}{Lower and upper quantiles used for intervals.}
}
\value{
A list with cell, sample, and marker summary tables.
}
\description{
Summarizes FibroDynMix bootstrap uncertainty draws into mean, standard
deviation, and quantile intervals.
}
