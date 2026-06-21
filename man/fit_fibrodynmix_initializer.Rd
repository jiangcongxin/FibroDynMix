\name{fit_fibrodynmix_initializer}
\alias{fit_fibrodynmix_initializer}
\title{Fit a count-based FibroDynMix initializer}
\usage{
fit_fibrodynmix_initializer(
  counts,
  marker_index,
  library_size = NULL,
  n_iter = 25,
  scale_factor = 10000,
  marker_boost = 1,
  ridge = 0.1,
  temperature = 1,
  prior_weight = 1,
  program_weight = 0.1,
  return_normalized_expression = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of weak prior marker genes for each state. Entries can be gene names or integer row indices.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{n_iter}{Number of alternating update iterations.}
\item{scale_factor}{Library-size normalization scale factor.}
\item{marker_boost}{Additive initialization boost for marker genes in their corresponding state program.}
\item{ridge}{Ridge penalty used in least-squares program updates.}
\item{temperature}{Softmax temperature for cell-state updates.}
\item{prior_weight}{Weight of the weak marker-prior logits retained during cell-state updates.}
\item{program_weight}{Weight of learned state-program reconstruction logits during cell-state updates.}
\item{return_normalized_expression}{Whether to retain the normalized gene-by-cell expression matrix in the returned object.}
\item{verbose}{Whether to print objective values.}
}
\value{
A list with \code{z_hat}, \code{beta_hat}, \code{alpha_hat},
\code{objective}, \code{marker_index}, and \code{prior_scores}. If
\code{return_normalized_expression = TRUE}, the returned object also contains
\code{normalized_expression}.
}
\description{
Learns an interpretable cell-state simplex and state-gene program from raw
counts using weak marker priors and alternating least-squares updates on
log-normalized expression. This is an initializer for later probabilistic
FibroDynMix inference, not the final Bayesian model.
}
