\name{prepare_fibrodynmix_data}
\alias{prepare_fibrodynmix_data}
\title{Prepare real single-cell count data for FibroDynMix}
\usage{
prepare_fibrodynmix_data(
  counts,
  cell_metadata,
  marker_index,
  cell_id_col = NULL,
  study_col = NULL,
  donor_col = NULL,
  library_size_col = NULL,
  min_cells_per_gene = 1,
  min_counts_per_gene = 1,
  drop_zero_library_cells = TRUE,
  require_all_states = TRUE
)
}
\arguments{
\item{counts}{Gene-by-cell raw UMI count matrix.}
\item{cell_metadata}{Data frame with one row per cell.}
\item{marker_index}{Named list of weak prior marker genes for each state.}
\item{cell_id_col}{Optional metadata column containing cell identifiers.}
\item{study_col}{Optional metadata column containing study/batch identifiers.}
\item{donor_col}{Optional metadata column containing donor/sample identifiers.}
\item{library_size_col}{Optional metadata column containing precomputed library sizes.}
\item{min_cells_per_gene}{Keep genes detected in at least this many cells.}
\item{min_counts_per_gene}{Keep genes with at least this total raw count.}
\item{drop_zero_library_cells}{Whether to remove cells with zero library size.}
\item{require_all_states}{Whether every state must retain at least one marker after filtering.}
}
\value{
A \code{FibroDynMixData} list with filtered counts, aligned metadata, marker priors,
library sizes, optional study/donor identifiers, and audit summaries.
}
\description{
Validates and aligns a raw gene-by-cell count matrix, cell metadata, and weak
marker priors before fitting FibroDynMix.
}
