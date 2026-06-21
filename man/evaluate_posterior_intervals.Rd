\name{evaluate_posterior_intervals}
\alias{evaluate_posterior_intervals}
\title{Evaluate posterior interval calibration for state weights}
\usage{
evaluate_posterior_intervals(z_true, cell_summary_z)
}
\arguments{
\item{z_true}{Cell-by-state matrix of true simplex weights.}
\item{cell_summary_z}{Data frame containing \code{cell_id}, \code{state},
\code{lower}, \code{upper}, and optionally \code{mean} columns.}
}
\value{
A list with scalar interval coverage and width diagnostics.
}
\description{
Compares true simulated latent state weights with posterior credible intervals
returned by \code{fit_fibrodynmix_vi()}.
}
