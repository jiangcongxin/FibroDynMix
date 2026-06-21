\name{fit_fibrodynmix_transfer}
\alias{fit_fibrodynmix_transfer}
\title{Transfer a fitted FibroDynMix program to held-out cells}
\usage{
fit_fibrodynmix_transfer(
  counts,
  fit,
  library_size = NULL,
  z_l2 = 0.001,
  maxit_z = 35,
  optimizer = c("BFGS", "L-BFGS-B"),
  optimizer_control = list(),
  z_init = NULL,
  chunk_size = NULL,
  parallel = FALSE,
  n_workers = 2,
  return_cell_diagnostics = TRUE
)
}
\arguments{
\item{counts}{Gene-by-cell raw count matrix for held-out cells.}
\item{fit}{A fitted object returned by \code{fit_fibrodynmix_nb()}.}
\item{library_size}{Optional held-out library sizes. If omitted, column sums of \code{counts} are used.}
\item{z_l2}{Non-negative L2 penalty for held-out cell logits.}
\item{maxit_z}{Maximum \code{optim()} iterations per held-out cell.}

\item{optimizer}{Optimizer method passed to \code{optim()} for held-out cell z updates.}

\item{optimizer_control}{Additional control list passed to \code{optim()}.}

\item{z_init}{Optional warm-start state-weight matrix for held-out cells. Rows
must match cells and columns must match fitted states.}

\item{chunk_size}{Optional number of cells to optimize per chunk. Chunking
bounds memory use and makes long transfers easier to monitor.}

\item{parallel}{Whether to optimize cells within each chunk in parallel.}

\item{n_workers}{Number of parallel workers when \code{parallel = TRUE}.}

\item{return_cell_diagnostics}{Whether to include a per-cell diagnostics table
with convergence code, iteration counts, and objective values.}
}
\value{
A list with \code{z_hat}, \code{nb_objective}, \code{nb_loglik}, shared-gene
metadata, and transfer diagnostics including \code{z_convergence_rate},
\code{n_nonconverged_cells}, and \code{cell_converged}.
}
\description{
Estimates held-out cell state mixtures under a fixed fitted FibroDynMix program.
The transfer step keeps \code{alpha}, \code{beta}, and \code{phi} fixed and optimizes only
held-out cell simplex weights under the raw-count negative-binomial likelihood.
}
