\name{compute_fpi}
\alias{compute_fpi}
\title{Compute Fibroblast Plasticity Index}
\usage{
compute_fpi(
  z,
  flow = NULL,
  transition_potential = NULL,
  lambda = 1,
  normalize_entropy = TRUE
)
}
\arguments{
\item{z}{Cell-by-state simplex matrix.}
\item{flow}{Optional state-by-state transition flow. If supplied, transition potential is computed from outgoing state flow.}
\item{transition_potential}{Optional named state-level transition potential vector. Used when \code{flow} is not supplied.}
\item{lambda}{Weight on transition potential.}
\item{normalize_entropy}{Whether to divide entropy by log(K).}
}
\value{
Data frame with cell-level entropy, transition potential, and FPI.
}
\description{
Computes a cell-level plasticity index from state entropy and transition
potential.
}
