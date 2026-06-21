\name{run_marker_scoring_benchmark}
\alias{run_marker_scoring_benchmark}
\title{Run a marker-scoring simulation benchmark}
\usage{
run_marker_scoring_benchmark(..., temperature = 1, marker_weight = NULL)
}
\arguments{
\item{...}{Arguments passed to \code{simulate_fibrodynmix()}.}
\item{temperature}{Softmax temperature for \code{score_marker_baseline()}.}
\item{marker_weight}{Optional \code{score_marker_baseline()} marker weights.}
}
\value{
A list with \code{simulation}, \code{baseline}, and \code{metrics}.
}
\description{
Simulates one FibroDynMix scenario, applies the marker-based baseline, and
evaluates state-weight recovery.
}
