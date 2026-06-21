\name{select_fibrodynmix_nb_model}
\alias{select_fibrodynmix_nb_model}
\title{Select a FibroDynMix NB outer-iteration budget}
\usage{
select_fibrodynmix_nb_model(
  counts,
  marker_index,
  library_size = NULL,
  candidate_n_outer = c(2L, 5L, 10L, 20L),
  labels = NULL,
  groups = NULL,
  holdout_fraction = 0.25,
  seed = NULL,
  nb_args = list(),
  transfer_args = list(),
  selection_weights = c(heldout_nb_objective = 0.35, z_stability = 0.25,
    marker_gradient = 0.2, downstream_balanced_accuracy = 0.2),
  truth_z = NULL
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of marker genes for each state.}
\item{library_size}{Optional vector of per-cell library sizes. If omitted, column sums of \code{counts} are used.}
\item{candidate_n_outer}{Positive integer candidate outer-iteration budgets.}
\item{labels}{Optional cell labels used for downstream validation.}
\item{groups}{Optional grouping variable used for the holdout split and grouped downstream cross-validation.}
\item{holdout_fraction}{Fraction of cells or groups assigned to the held-out validation split.}
\item{seed}{Optional random seed for the validation split.}
\item{nb_args}{Named list of arguments passed to \code{fit_fibrodynmix_nb()}.}
\item{transfer_args}{Named list of arguments passed to \code{fit_fibrodynmix_transfer()}.}
\item{selection_weights}{Named numeric vector of weights for the validation metrics.}
\item{truth_z}{Optional simulated truth used only for audit metrics. It is not used in the selection score.}
}
\value{
A list with:
\itemize{
\item \code{selected_n_outer}: selected outer-iteration budget.
\item \code{candidate_scores}: validation metrics and selection scores for all candidates.
\item \code{selection_components}: per-state marker-gradient preservation components.
\item \code{split}: train/holdout split metadata.
\item \code{selection_weights}: weights used by the selection score.
\item \code{selection_rule}: short description of the rule.
}
}
\description{
Fits candidate FibroDynMix negative-binomial models on a training cell split,
transfers each fitted program to held-out cells, and ranks candidates by
validation metrics rather than by the training objective alone.
}
\details{
The default selection score combines held-out NB objective, z stability across
candidate budgets, marker-gradient preservation, and downstream validation
balanced accuracy when labels are supplied. The training objective is reported
for auditability but is not used as a standalone selector. When \code{truth_z}
is supplied for simulations, truth-based RMSE and dominant-state accuracy are
reported as audit columns only.
}
