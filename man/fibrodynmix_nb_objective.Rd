\name{fibrodynmix_nb_objective}
\alias{fibrodynmix_nb_objective}
\title{FibroDynMix negative-binomial optimization objective}
\usage{
fibrodynmix_nb_objective(
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
  beta_l2 = 0,
  marker_target = NULL,
  marker_l2 = 0,
  effect_l2 = 0,
  donor_effect_l2 = NULL,
  average = FALSE
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
\item{beta_l2}{Non-negative L2 penalty for \code{beta}.}
\item{marker_target}{Optional state-by-gene matrix used as a weak marker orientation target for \code{beta}.}
\item{marker_l2}{Non-negative L2 penalty around \code{marker_target}.}
\item{effect_l2}{Non-negative L2 penalty for study effects and, by default, donor effects.}
\item{donor_effect_l2}{Optional non-negative L2 penalty for donor effects. If \code{NULL}, \code{effect_l2} is used.}
\item{average}{If \code{TRUE}, divide the objective by the number of observed count entries.}
}
\value{
Numeric scalar objective to minimize.
}
\description{
Returns a minimization objective based on the negative raw-count
negative-binomial log-likelihood. Optional ridge penalties can be applied to
\code{beta}, \code{study_effect}, and \code{donor_effect}.
}
