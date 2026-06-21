\name{bootstrap_fibrodynmix}
\alias{bootstrap_fibrodynmix}
\title{Bootstrap FibroDynMix state uncertainty}
\usage{
bootstrap_fibrodynmix(
  counts,
  marker_index,
  library_size = NULL,
  cell_metadata = NULL,
  sample_col = NULL,
  method = c("initializer", "nb", "nb_study"),
  n_boot = 20,
  seed = 1,
  fit_args = list(),
  keep_fits = FALSE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of weak prior marker genes for each state.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{cell_metadata}{Optional cell metadata. If supplied, \code{sample_col} can be used for sample-level composition intervals.}
\item{sample_col}{Column in \code{cell_metadata} used to aggregate sample-level composition. If \code{NULL}, all cells are treated as one sample.}
\item{method}{Fitting method: \code{initializer}, \code{nb}, or \code{nb_study}.}
\item{n_boot}{Number of bootstrap replicates.}
\item{seed}{Optional random seed.}
\item{fit_args}{Additional arguments passed to the selected fit function.}
\item{keep_fits}{Whether to retain full fit objects for each replicate.}
}
\value{
A list containing bootstrap draws and uncertainty summaries.
}
\description{
Runs a lightweight cell bootstrap to estimate uncertainty in cell-state
mixtures, entropy, sample-level composition, and marker-program stability.
}
