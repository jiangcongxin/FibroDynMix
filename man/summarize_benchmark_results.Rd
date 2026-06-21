\name{summarize_benchmark_results}
\alias{summarize_benchmark_results}
\title{Summarize benchmark results}
\usage{
summarize_benchmark_results(
  benchmark_results,
  metrics = c("rmse", "mean_absolute_error", "dominant_accuracy", "mean_entropy_pred")
)
}
\arguments{
\item{benchmark_results}{Data frame returned by \code{run_simulation_benchmark()}.}
\item{metrics}{Character vector of metric columns to summarize.}
}
\value{
A data frame with mean, standard deviation, and replicate count for each metric.
}
\description{
Aggregates replicate-level benchmark rows by scenario and method.
}
