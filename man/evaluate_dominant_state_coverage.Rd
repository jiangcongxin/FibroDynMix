\name{evaluate_dominant_state_coverage}
\alias{evaluate_dominant_state_coverage}
\title{Evaluate dominant-state classification among confident cells}
\usage{
evaluate_dominant_state_coverage(
  z_true,
  z_pred,
  confidence_threshold = 0.8
)
}
\arguments{
\item{z_true}{Numeric matrix of true cell-by-state weights.}

\item{z_pred}{Numeric matrix of predicted cell-by-state weights.}

\item{confidence_threshold}{Confidence cutoff on predicted dominant state
weight.}
}
\value{
A list with coverage and confident-cell dominant accuracy.
}
\description{
Computes coverage and dominant accuracy only on cells whose maximum predicted
state probability passes a confidence threshold.
}
