\name{fit_fibrodynmix_nb}
\alias{fit_fibrodynmix_nb}
\title{Fit the first FibroDynMix negative-binomial model}
\usage{
fit_fibrodynmix_nb(
  counts,
  marker_index,
  library_size = NULL,
  study_id = NULL,
  n_outer = 5,
  initializer_args = list(),
  estimate_phi = TRUE,
  phi_init = NULL,
  beta_l2 = 0.01,
  marker_l2 = 0.05,
  marker_weight = NULL,
  fit_study_effect = FALSE,
  study_l2 = 0.1,
  donor_id = NULL,
  fit_donor_effect = FALSE,
  donor_l2 = 0.1,
  z_l2 = 0.001,
  optimizer = c("BFGS", "L-BFGS-B"),
  optimizer_control = list(),
  maxit_beta = 50,
  maxit_z = 35,
  early_stopping = TRUE,
  patience = 2,
  min_delta = 1e-05,
  objective_rel_tol = 1e-06,
  stagnation_window = 5L,
  objective_abs_tol = 1e-08,
  rollback_to_best = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of weak prior marker genes for each state.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{study_id}{Optional cell-level study identifiers. When supplied with \code{fit_study_effect = TRUE}, a ridge-penalized study-by-gene effect is fit.}
\item{donor_id}{Optional cell-level donor identifiers. When supplied with \code{fit_donor_effect = TRUE}, a ridge-penalized donor-by-gene effect is fit.}
\item{n_outer}{Number of alternating NB optimization iterations.}
\item{initializer_args}{Named list of arguments passed to \code{fit_fibrodynmix_initializer()}.}
\item{estimate_phi}{Whether to update gene-specific \code{phi} by method-of-moments after each outer iteration.}
\item{phi_init}{Optional initial gene-specific NB size vector.}
\item{beta_l2}{Non-negative L2 penalty for state-gene program coefficients.}
\item{marker_l2}{Non-negative L2 penalty that keeps prior marker coefficients oriented toward positive state-specific programs.}
\item{marker_weight}{Optional marker reliability weights used by initialization and marker orientation.}
\item{fit_study_effect}{Whether to fit ridge-penalized study-by-gene effects.}
\item{study_l2}{Non-negative L2 penalty for study effects.}
\item{fit_donor_effect}{Whether to fit ridge-penalized donor-by-gene effects.}
\item{donor_l2}{Non-negative L2 penalty for donor effects.}
\item{z_l2}{Non-negative L2 penalty for cell logits in simplex updates.}
\item{optimizer}{Optimizer method passed to \code{optim()} for beta and z updates.}
\item{optimizer_control}{Additional control list passed to \code{optim()}.}
\item{maxit_beta}{Maximum \code{optim()} iterations for each gene-level \code{alpha/beta} update.}
\item{maxit_z}{Maximum \code{optim()} iterations for each cell-level \code{z} update.}
\item{early_stopping}{Whether to stop when the best objective has not improved for \code{patience} iterations.}
\item{patience}{Number of non-improving iterations tolerated when \code{early_stopping = TRUE}.}
\item{min_delta}{Minimum objective decrease required to count as an improvement.}
\item{objective_rel_tol}{Relative objective tolerance for convergence/stagnation checks.}
\item{stagnation_window}{Number of recent objective values used for stagnation checks.}
\item{objective_abs_tol}{Absolute objective tolerance for convergence/stagnation checks.}
\item{rollback_to_best}{Whether to return the best-seen parameters instead of the last parameters.}
\item{verbose}{Whether to print objective values.}
}
\value{
A list with fitted \code{z_hat}, \code{beta_hat}, \code{alpha_hat},
\code{phi_hat}, \code{nb_objective_trace}, \code{best_objective},
\code{best_iteration}, \code{stop_reason}, \code{converged}, and optimizer
metadata.
}
\description{
Fits a lightweight raw-count negative-binomial FibroDynMix model by alternating
optimization. By default this optimizer estimates
\code{log(mu_ig) = log(library_i) + alpha_g + sum_k z_ik beta_kg}. When
\code{fit_study_effect = TRUE} or \code{fit_donor_effect = TRUE}, it additionally
fits ridge-penalized study-by-gene and donor-by-gene effects. This is intended as
the bridge from the log-normalized initializer toward the full hierarchical
Bayesian/VI implementation.
}
