\name{fit_scvi_latent_baseline}
\alias{fit_scvi_latent_baseline}
\title{Fit an optional scVI latent baseline}
\usage{
fit_scvi_latent_baseline(
  counts,
  marker_index,
  library_size = NULL,
  n_latent = length(marker_index),
  max_epochs = 40,
  seed = NULL,
  projection_temperature = NULL
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix.}

\item{marker_index}{Named list of marker genes for each state.}

\item{library_size}{Optional vector of per-cell library sizes used only for
marker-score alignment.}

\item{n_latent}{scVI latent dimensionality.}

\item{max_epochs}{Training epochs passed to \code{scvi.model.SCVI.train()}.}

\item{seed}{Optional seed passed to Python and R clustering.}

\item{projection_temperature}{Optional distance softmax temperature. If
omitted, the median latent distance to cluster centers is used.}
}
\value{
A list with \code{latent}, \code{z_pred}, cluster assignments, cluster-state
alignment, and projection metadata.
}
\description{
Fits a small scVI model through \code{reticulate}, then projects the learned
latent representation to a state simplex by soft cluster membership. Clusters
are aligned to fibroblast states using marker-score enrichment, so marker genes
orient labels but do not define the deep latent representation.
}
\details{
This baseline is optional because it requires a Python environment with
\code{scvi-tools} and \code{anndata}. The simplex projection is a post hoc
readout from the scVI latent space, not a native scVI mixture-composition
estimate.
}
