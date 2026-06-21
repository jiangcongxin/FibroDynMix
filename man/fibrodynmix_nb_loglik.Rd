\name{fibrodynmix_nb_loglik}
\alias{fibrodynmix_nb_loglik}
\title{FibroDynMix negative-binomial log-likelihood}
\usage{
fibrodynmix_nb_loglik(
  counts,
  z,
  beta,
  alpha,
  phi,
  library_size,
  study_effect = NULL,
  donor_effect = NULL,
  study_id = NULL,
  donor_id = NULL,
  return_matrix = FALSE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{z}{Cell-by-state simplex matrix.}
\item{beta}{State-by-gene expression program matrix.}
\item{alpha}{Gene-level baseline vector.}
\item{phi}{Gene-specific negative-binomial size/overdispersion vector.}
\item{library_size}{Positive vector with one library size per cell.}
\item{study_effect}{Optional study-by-gene effect matrix.}
\item{donor_effect}{Optional donor-by-gene effect matrix.}
\item{study_id}{Optional cell-level study identifiers matching rows of \code{study_effect}.}
\item{donor_id}{Optional cell-level donor identifiers matching rows of \code{donor_effect}.}
\item{return_matrix}{If \code{TRUE}, return a gene-by-cell log-likelihood matrix; otherwise return the scalar sum.}
}
\value{
Numeric scalar or gene-by-cell matrix.
}
\description{
Evaluates the raw-count negative-binomial likelihood under the FibroDynMix mean
model:
\code{log(mu_ig) = log(library_i) + alpha_g + sum_k z_ik beta_kg + effects_ig}.
}
