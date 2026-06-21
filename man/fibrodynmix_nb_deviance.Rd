\name{fibrodynmix_nb_deviance}
\alias{fibrodynmix_nb_deviance}
\title{FibroDynMix negative-binomial deviance}
\usage{
fibrodynmix_nb_deviance(
  counts,
  z,
  beta,
  alpha,
  phi,
  library_size,
  study_effect = NULL,
  donor_effect = NULL,
  study_id = NULL,
  donor_id = NULL
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
}
\value{
Numeric scalar deviance.
}
\description{
Computes twice the difference between the saturated negative-binomial
log-likelihood and the fitted FibroDynMix negative-binomial log-likelihood,
using the supplied gene-specific \code{phi}.
}
