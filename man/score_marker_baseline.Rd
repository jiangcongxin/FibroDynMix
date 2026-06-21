\name{score_marker_baseline}
\alias{score_marker_baseline}
\title{Score cells with a marker-based baseline}
\usage{
score_marker_baseline(
  counts,
  marker_index,
  marker_weight = NULL,
  library_size = NULL,
  scale_factor = 10000,
  temperature = 1
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of marker genes for each state. Entries can be gene names or integer row indices.}
\item{marker_weight}{Optional per-state marker-gene reliability weights. Supported formats are a state-by-gene matrix, a named list keyed by state, or a numeric vector.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{scale_factor}{Library-size normalization scale factor.}
\item{temperature}{Softmax temperature applied to centered state scores.}
}
\value{
A list with \code{z_pred}, \code{scores}, and \code{normalized_expression}.
}
\description{
Computes a simple gene-set scoring baseline from raw counts. This function is
intended as a comparator for simulation benchmarks, not as the FibroDynMix
generative model.
}
