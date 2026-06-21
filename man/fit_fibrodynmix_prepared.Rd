\name{fit_fibrodynmix_prepared}
\alias{fit_fibrodynmix_prepared}
\title{Fit FibroDynMix from a prepared real-data object}
\usage{
fit_fibrodynmix_prepared(data, fit_study_effect = NULL, fit_donor_effect = NULL, ...)
}
\arguments{
\item{data}{A \code{FibroDynMixData} object.}
\item{fit_study_effect}{Whether to fit study-by-gene effects. If \code{NULL}, inferred from \code{data$study_id}.}
\item{fit_donor_effect}{Whether to fit donor-by-gene effects. If \code{NULL}, inferred when donor identifiers are available and represent a finer level than study identifiers.}
\item{...}{Additional arguments passed to \code{fit_fibrodynmix_nb()}.}
}
\value{
A fitted FibroDynMix NB model.
}
\description{
Convenience wrapper around \code{fit_fibrodynmix_nb()} for objects returned by
\code{prepare_fibrodynmix_data()}.
}
