\name{evaluate_state_weights}
\alias{evaluate_state_weights}
\title{Evaluate inferred fibroblast state weights}
\usage{
evaluate_state_weights(z_true, z_pred)
}
\arguments{
\item{z_true}{Numeric matrix of true cell-by-state weights.}
\item{z_pred}{Numeric matrix of predicted cell-by-state weights.}
}
\value{
A list with scalar metrics, per-state RMSE, and a confusion matrix for
dominant states.
}
\description{
Computes state-weight recovery metrics against the true latent simplex.
}
