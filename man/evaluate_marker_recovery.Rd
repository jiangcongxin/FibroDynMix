\name{evaluate_marker_recovery}
\alias{evaluate_marker_recovery}
\title{Evaluate recovery of state-specific marker programs}
\usage{
evaluate_marker_recovery(marker_truth, marker_scores)
}
\arguments{
\item{marker_truth}{Logical or 0/1 matrix with states in rows and genes in columns. Non-zero entries indicate true state markers.}
\item{marker_scores}{Numeric matrix with the same dimensions as \code{marker_truth}; larger values indicate stronger predicted marker support.}
}
\value{
A list with macro AUPRC and per-state AUPRC.
}
\description{
Computes state-wise and macro-average AUPRC for predicted gene-state marker
scores against true marker assignments.
}
