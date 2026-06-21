\name{summarize_optimizer_diagnostics}
\alias{summarize_optimizer_diagnostics}
\title{Summarize FibroDynMix NB optimizer diagnostics}
\usage{
summarize_optimizer_diagnostics(benchmark_results)
}
\arguments{
\item{benchmark_results}{Data frame returned by \code{run_simulation_benchmark()}.}
}
\value{
A data frame summarizing objective improvement, rollback frequency, early
stopping frequency, and best iteration by scenario and method.
}
\description{
Aggregates NB optimizer diagnostics emitted by \code{run_simulation_benchmark()}.
Non-NB methods are retained with \code{NA} diagnostics so benchmark tables
remain rectangular.
}
