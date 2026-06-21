\name{fit_fibrodynmix_vi}
\alias{fit_fibrodynmix_vi}
\title{Fit a lightweight logistic-normal variational posterior for FibroDynMix}
\usage{
fit_fibrodynmix_vi(
  counts,
  marker_index,
  library_size = NULL,
  nb_args = list(),
  n_draws = 50,
  n_elbo_draws = 10,
  n_vi_iter = 3,
  posterior_scale = 0.25,
  scale_grid = c(0.5, 1, 1.5),
  prior_sd = 2,
  probs = c(0.025, 0.975),
  seed = 1,
  cell_metadata = NULL,
  sample_col = NULL,
  keep_draws = TRUE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of weak prior marker genes for each state.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{nb_args}{Named list passed to \code{fit_fibrodynmix_nb()}.}
\item{n_draws}{Number of posterior draws used for summaries.}
\item{n_elbo_draws}{Number of Monte Carlo draws used per ELBO evaluation.}
\item{n_vi_iter}{Number of scale-refinement iterations for the diagonal logistic-normal posterior.}
\item{posterior_scale}{Initial logit posterior standard deviation.}
\item{scale_grid}{Multipliers evaluated around the current posterior scale.}
\item{prior_sd}{Logistic-normal prior standard deviation for cell logits.}
\item{probs}{Lower and upper probabilities for credible intervals.}
\item{seed}{Optional random seed.}
\item{cell_metadata}{Optional cell metadata for sample-level composition posterior summaries.}
\item{sample_col}{Optional column in \code{cell_metadata} used for sample-level composition summaries.}
\item{keep_draws}{Whether to retain full posterior cell-level draws.}
}
\value{
A \code{FibroDynMixVI} list containing the underlying NB fit, variational
parameters \code{eta_mean} and \code{eta_sd}, posterior mean state weights,
cell-level credible intervals, optional sample-level credible intervals,
posterior draws, and an ELBO-like trace.
}
\description{
Fits the current negative-binomial FibroDynMix optimizer, then constructs a
mean-field logistic-normal variational posterior over each cell's latent simplex
weights. \code{alpha}, \code{beta}, \code{phi}, study effects, and donor effects are
fixed at the fitted mode while \code{q(z_i)} is represented in cell-logit space.
}
\details{
This function is a lightweight posterior layer around the implemented raw-count
NB optimizer. It provides auditable posterior summaries and an ELBO-like
Monte Carlo evidence trace, but it is not yet an amortized neural inference
engine or a full hierarchical variational posterior over every model parameter.
}
