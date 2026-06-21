\name{estimate_transition_flow}
\alias{estimate_transition_flow}
\title{Estimate condition-to-condition fibroblast transition flow}
\usage{
estimate_transition_flow(
  source,
  target,
  cost,
  lambda = 0.1,
  prior = NULL,
  max_iter = 1000,
  tol = 1e-08
)
}
\arguments{
\item{source}{Composition vector for the source condition.}
\item{target}{Composition vector for the target condition.}
\item{cost}{State-by-state transition cost matrix.}
\item{lambda}{Entropic regularization strength.}
\item{prior}{Optional prior transition matrix. If supplied, the kernel is multiplied by this prior before Sinkhorn scaling.}
\item{max_iter}{Maximum Sinkhorn iterations.}
\item{tol}{Convergence tolerance on row/column marginals.}
}
\value{
A list with transition matrix \code{flow}, normalized marginals, expected cost,
entropy, convergence flag, and iteration count.
}
\description{
Estimates an entropy-regularized optimal-transport flow from source state
composition to target state composition.
}
