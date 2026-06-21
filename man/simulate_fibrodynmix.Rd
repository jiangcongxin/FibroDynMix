\name{simulate_fibrodynmix}
\alias{simulate_fibrodynmix}
\title{Simulate FibroDynMix count data}
\usage{
simulate_fibrodynmix(
  n_studies = 2,
  donors_per_study = 4,
  cells_per_donor = 100,
  n_genes = 1000,
  state_names = c("resident", "inflammatory", "myofibroblast", "ECM-remodeling",
    "antigen-presenting", "IFN-stress"),
  scenario = c("continuous", "discrete", "batch_confounding", "rare_transition"),
  marker_genes_per_state = 30,
  tau_high = 1,
  tau_low = 0.08,
  study_effect_sd = 0.15,
  donor_effect_sd = 0.08,
  state_sd = 0.35,
  library_meanlog = log(5000),
  library_sdlog = 0.35,
  seed = NULL
)
}
\arguments{
\item{n_studies}{Number of studies or cohorts.}
\item{donors_per_study}{Number of donors per study. Either a scalar or an integer vector of length \code{n_studies}.}
\item{cells_per_donor}{Number of cells per donor. Either a scalar or an integer vector with one value per donor.}
\item{n_genes}{Number of genes.}
\item{state_names}{Character vector of fibroblast state names.}
\item{scenario}{Simulation scenario: \code{continuous}, \code{discrete}, \code{batch_confounding}, or \code{rare_transition}.}
\item{marker_genes_per_state}{Number of prior marker genes assigned to each state in the data-generating beta matrix.}
\item{tau_high}{Standard deviation scale for marker state effects.}
\item{tau_low}{Standard deviation for non-marker state effects.}
\item{study_effect_sd}{Standard deviation of study-by-gene effects.}
\item{donor_effect_sd}{Standard deviation of donor-by-gene effects.}
\item{state_sd}{Cell-level logistic-normal standard deviation.}
\item{library_meanlog}{Meanlog for log-normal library sizes.}
\item{library_sdlog}{Sdlog for log-normal library sizes.}
\item{seed}{Optional random seed.}
}
\value{
A list with \code{counts}, \code{cell_metadata}, \code{donor_metadata},
\code{gene_metadata}, \code{z}, and \code{parameters}.
}
\description{
Generates raw UMI-like counts from the FibroDynMix negative-binomial
hierarchical mixture model. The simulator returns observed counts together
with the latent state simplex, donor/study metadata, and true generative
parameters used for benchmarking.
}
