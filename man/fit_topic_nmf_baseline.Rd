\name{fit_topic_nmf_baseline}
\alias{fit_topic_nmf_baseline}
\title{Fit a topic/NMF baseline from raw counts}
\usage{
fit_topic_nmf_baseline(
  counts,
  marker_index,
  n_topics = length(marker_index),
  n_iter = 100,
  seed = NULL,
  eps = 1e-08,
  backend = c("auto", "nmf", "multiplicative_update"),
  nmf_method = "brunet",
  nrun = 1
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}
\item{marker_index}{Named list of marker genes for each state. Entries can be gene names or integer row indices.}
\item{n_topics}{Number of NMF topics. Defaults to the number of states.}
\item{n_iter}{Number of multiplicative-update iterations.}
\item{seed}{Optional random seed for initialization.}
\item{eps}{Small positive value for numerical stability.}
\item{backend}{Factorization backend. \code{"auto"} uses the \code{NMF} package when installed and otherwise falls back to the internal KL multiplicative updates.}
\item{nmf_method}{Method passed to \code{NMF::nmf()} when \code{backend = "nmf"}.}
\item{nrun}{Number of NMF runs passed to \code{NMF::nmf()}.}
}
\value{
A list with topic gene loadings \code{w}, cell topic loadings \code{h},
state-aligned cell weights \code{z_pred}, topic-state scores, topic
assignment, and KL objective trace.
}
\description{
Fits a lightweight KL-divergence non-negative matrix factorization baseline to a
gene-by-cell raw count matrix. Topics are aligned to fibroblast states after
fitting by marker enrichment, so marker genes orient the output labels but do
not define the cell-topic factors.
}
