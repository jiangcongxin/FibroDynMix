\name{calibrate_posterior_interval_scale}
\alias{calibrate_posterior_interval_scale}
\title{Calibrate posterior interval scale in simulation}
\usage{
calibrate_posterior_interval_scale(
  z_true,
  cell_summary_z,
  scale_grid = c(1, 1.5, 2, 3, 4, 6, 8, 12, 24, 100),
  target_coverage = 0.9
)
}
\arguments{
\item{z_true}{Cell-by-state matrix of true simplex weights.}
\item{cell_summary_z}{Data frame containing \code{cell_id}, \code{state},
\code{lower}, \code{upper}, and optionally \code{mean} columns.}
\item{scale_grid}{Positive interval-width multipliers to evaluate.}
\item{target_coverage}{Target empirical coverage.}
}
\value{
A list with the selected scale and calibrated interval diagnostics.
}
\description{
Uses known simulated latent state weights to choose a multiplicative expansion
of posterior intervals. This is intended for simulation calibration and method
diagnostics, not for real-data inference where latent truth is unknown.
}
