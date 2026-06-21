\name{compute_state_cost}
\alias{compute_state_cost}
\title{Compute transcriptional cost between fibroblast states}
\usage{
compute_state_cost(beta, method = c("euclidean", "correlation"), scale_genes = TRUE)
}
\arguments{
\item{beta}{State-by-gene program matrix.}
\item{method}{Distance method. Currently \code{euclidean} or \code{correlation}.}
\item{scale_genes}{Whether to z-score genes before computing distances.}
}
\value{
State-by-state cost matrix with zero diagonal.
}
\description{
Builds a state-by-state transition cost matrix from state-gene programs.
}
