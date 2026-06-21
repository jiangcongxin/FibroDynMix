\name{evaluate_downstream_classification}
\alias{evaluate_downstream_classification}
\title{Evaluate downstream classification from inferred state features}
\usage{
evaluate_downstream_classification(
  features,
  labels,
  groups = NULL,
  n_folds = 5,
  seed = NULL
)
}
\arguments{
\item{features}{Numeric observation-by-feature matrix.}

\item{labels}{Class labels, one per observation.}

\item{groups}{Optional grouping variable. When supplied, cross-validation folds
are assigned at the group level to avoid train/test leakage.}

\item{n_folds}{Requested number of cross-validation folds.}

\item{seed}{Optional seed for fold assignment.}
}
\value{
A list with cross-validated accuracy, balanced accuracy, macro-F1, macro
one-vs-rest AUROC, and fold metadata.
}
\description{
Runs a lightweight nearest-centroid classifier with stratified
cross-validation. The classifier is intentionally simple so benchmark
differences primarily reflect the supplied feature representation rather than
a heavily tuned downstream model.
}
