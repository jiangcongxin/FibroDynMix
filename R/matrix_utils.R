is_matrix_like <- function(x) {
  length(dim(x)) == 2L && (is.numeric(x) || inherits(x, "Matrix"))
}

matrix_col_sums <- function(x) {
  if (inherits(x, "sparseMatrix") || inherits(x, "Matrix")) {
    return(as.numeric(Matrix::colSums(x)))
  }
  colSums(x)
}

matrix_row_sums <- function(x) {
  if (inherits(x, "sparseMatrix") || inherits(x, "Matrix")) {
    return(as.numeric(Matrix::rowSums(x)))
  }
  rowSums(x)
}

matrix_col_means <- function(x) {
  if (inherits(x, "sparseMatrix") || inherits(x, "Matrix")) {
    return(as.numeric(Matrix::colMeans(x)))
  }
  colMeans(x)
}

matrix_row_means_core <- function(x) {
  if (inherits(x, "sparseMatrix") || inherits(x, "Matrix")) {
    return(as.numeric(Matrix::rowMeans(x)))
  }
  rowMeans(x)
}

matrix_n_entries <- function(x) {
  nrow(x) * ncol(x)
}

matrix_gene_vector <- function(x, row_index) {
  as.numeric(x[row_index, ])
}

matrix_cell_vector <- function(x, col_index) {
  as.numeric(x[, col_index])
}

matrix_is_nonnegative_integerish <- function(x) {
  if (anyNA(x)) {
    return(FALSE)
  }
  if (inherits(x, "sparseMatrix")) {
    values <- x@x
    return(all(is.finite(values)) && all(values >= 0) && all(values == round(values)))
  }
  all(is.finite(x)) && all(x >= 0) && all(x == round(x))
}

log_normalize_counts <- function(counts, library_size, scale_factor = 10000) {
  scale <- scale_factor / library_size
  if (inherits(counts, "sparseMatrix") || inherits(counts, "Matrix")) {
    return(log1p(counts %*% Matrix::Diagonal(x = scale)))
  }
  log1p(sweep(counts, 2L, library_size, "/") * scale_factor)
}
