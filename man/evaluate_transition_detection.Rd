\name{evaluate_transition_detection}
\alias{evaluate_transition_detection}
\title{Evaluate rare transition detection}
\usage{
evaluate_transition_detection(transition_truth, transition_score, threshold = NULL)
}
\arguments{
\item{transition_truth}{Logical vector indicating true transition cells.}
\item{transition_score}{Numeric vector where larger values indicate higher transition evidence.}
\item{threshold}{Optional score threshold. If omitted, the top \code{sum(truth)} cells are called positive.}
}
\value{
A list with AUROC, AUPRC, threshold, precision, recall, and F1.
}
\description{
Computes threshold-free and thresholded detection metrics for rare
transition-like cells.
}
