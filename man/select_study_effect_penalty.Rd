\name{select_study_effect_penalty}
\alias{select_study_effect_penalty}
\title{Select study-effect and marker penalties from a sensitivity benchmark}
\usage{
select_study_effect_penalty(
  sensitivity_results,
  rmse_tolerance = 0.05,
  objective_tolerance = 0.03,
  weights = c(rmse = 0.45, objective = 0.25, study_effect = 0.2, rollback = 0.1)
)
}
\arguments{
\item{sensitivity_results}{Data frame returned by \code{run_study_effect_sensitivity()}.}
\item{rmse_tolerance}{Relative tolerance around the best RMSE. Candidates within \code{best_rmse * (1 + rmse_tolerance)} are preferred.}
\item{objective_tolerance}{Relative tolerance around the best NB objective.}
\item{weights}{Named numeric vector with weights for \code{rmse}, \code{objective}, \code{study_effect}, and \code{rollback}.}
}
\value{
A list with \code{recommended_study_l2}, \code{recommended_marker_l2},
\code{selection_reason}, and \code{tradeoff_table}.
}
\description{
Chooses a simulation-calibrated default penalty pair from
\code{run_study_effect_sensitivity()} output. The selection rule prioritizes low
RMSE while penalizing objective degradation, large study effects, and optimizer
rollbacks.
}
