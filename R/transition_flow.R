#' Compute transcriptional cost between fibroblast states
#'
#' Builds a state-by-state transition cost matrix from state-gene programs.
#'
#' @param beta State-by-gene program matrix.
#' @param method Distance method. Currently `euclidean` or `correlation`.
#' @param scale_genes Whether to z-score genes before computing distances.
#'
#' @return State-by-state cost matrix with zero diagonal.
#' @export
compute_state_cost <- function(beta,
                               method = c("euclidean", "correlation"),
                               scale_genes = TRUE) {
  method <- match.arg(method)
  beta <- as.matrix(beta)
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be a numeric state-by-gene matrix without NA values.", call. = FALSE)
  }
  if (is.null(rownames(beta))) {
    rownames(beta) <- sprintf("state_%d", seq_len(nrow(beta)))
  }

  x <- beta
  if (isTRUE(scale_genes)) {
    x <- t(scale(t(x)))
    x[is.na(x)] <- 0
  }

  if (method == "euclidean") {
    cost <- as.matrix(stats::dist(x, method = "euclidean"))
  } else {
    corr <- stats::cor(t(x), use = "pairwise.complete.obs")
    corr[is.na(corr)] <- 0
    cost <- 1 - corr
  }
  diag(cost) <- 0
  dimnames(cost) <- list(rownames(beta), rownames(beta))
  cost
}

#' Estimate condition-to-condition fibroblast transition flow
#'
#' Estimates an entropy-regularized optimal-transport flow from source state
#' composition to target state composition.
#'
#' @param source Composition vector for the source condition.
#' @param target Composition vector for the target condition.
#' @param cost State-by-state transition cost matrix.
#' @param lambda Entropic regularization strength.
#' @param prior Optional prior transition matrix. If supplied, the kernel is
#'   multiplied by this prior before Sinkhorn scaling.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Convergence tolerance on row/column marginals.
#'
#' @return A list with transition matrix `flow`, normalized marginals, expected
#'   cost, entropy, convergence flag, and iteration count.
#' @export
estimate_transition_flow <- function(source,
                                     target,
                                     cost,
                                     lambda = 0.1,
                                     prior = NULL,
                                     max_iter = 1000,
                                     tol = 1e-8) {
  source <- normalize_composition(source, "source")
  target <- normalize_composition(target, "target")
  cost <- as.matrix(cost)
  if (!identical(names(source), rownames(cost)) || !identical(names(target), colnames(cost))) {
    stop("Names of `source` and `target` must match row/column names of `cost`.", call. = FALSE)
  }
  if (length(lambda) != 1L || is.na(lambda) || lambda <= 0) {
    stop("`lambda` must be a positive numeric scalar.", call. = FALSE)
  }
  assert_positive_integer(max_iter, "max_iter")
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a positive numeric scalar.", call. = FALSE)
  }

  kernel <- exp(-cost / lambda)
  if (!is.null(prior)) {
    prior <- as.matrix(prior)
    if (!identical(dim(prior), dim(cost))) {
      stop("`prior` must have the same dimensions as `cost`.", call. = FALSE)
    }
    kernel <- kernel * pmax(prior, .Machine$double.eps)
  }
  kernel <- pmax(kernel, .Machine$double.eps)

  u <- rep(1, length(source))
  v <- rep(1, length(target))
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    u <- source / as.vector(kernel %*% v)
    v <- target / as.vector(t(kernel) %*% u)
    flow <- diag(u, nrow = length(u)) %*% kernel %*% diag(v, nrow = length(v))
    row_error <- max(abs(rowSums(flow) - source))
    col_error <- max(abs(colSums(flow) - target))
    if (max(row_error, col_error) < tol) {
      converged <- TRUE
      break
    }
  }
  dimnames(flow) <- dimnames(cost)
  flow <- pmax(flow, 0)
  flow <- flow / sum(flow)

  entropy <- -sum(flow * log(pmax(flow, .Machine$double.eps)))
  list(
    flow = flow,
    source = source,
    target = target,
    cost = cost,
    expected_cost = sum(flow * cost),
    entropy = entropy,
    converged = converged,
    iterations = iter
  )
}

#' Compute Fibroblast Plasticity Index
#'
#' Computes a cell-level plasticity index from state entropy and transition
#' potential.
#'
#' @param z Cell-by-state simplex matrix.
#' @param flow Optional state-by-state transition flow. If supplied,
#'   transition potential is computed from outgoing state flow.
#' @param transition_potential Optional named state-level transition potential
#'   vector. Used when `flow` is not supplied.
#' @param lambda Weight on transition potential.
#' @param normalize_entropy Whether to divide entropy by log(K).
#'
#' @return Data frame with cell-level entropy, transition potential, and FPI.
#' @export
compute_fpi <- function(z,
                        flow = NULL,
                        transition_potential = NULL,
                        lambda = 1,
                        normalize_entropy = TRUE) {
  z <- as.matrix(z)
  if (!is.numeric(z) || anyNA(z) || any(z < 0)) {
    stop("`z` must be a non-negative numeric matrix.", call. = FALSE)
  }
  if (any(abs(rowSums(z) - 1) > 1e-4)) {
    stop("Rows of `z` must sum to one.", call. = FALSE)
  }
  if (is.null(colnames(z))) {
    colnames(z) <- sprintf("state_%d", seq_len(ncol(z)))
  }
  if (is.null(rownames(z))) {
    rownames(z) <- sprintf("cell_%d", seq_len(nrow(z)))
  }
  if (length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("`lambda` must be a non-negative numeric scalar.", call. = FALSE)
  }

  entropy <- -rowSums(pmax(z, .Machine$double.eps) * log(pmax(z, .Machine$double.eps)))
  if (isTRUE(normalize_entropy)) {
    entropy <- entropy / log(ncol(z))
  }

  if (!is.null(flow)) {
    flow <- as.matrix(flow)
    if (!identical(colnames(z), rownames(flow))) {
      stop("`colnames(z)` must match `rownames(flow)`.", call. = FALSE)
    }
    transition_potential <- rowSums(flow * (1 - diag(nrow(flow))))
    names(transition_potential) <- rownames(flow)
  }
  if (is.null(transition_potential)) {
    transition_potential <- rep(0, ncol(z))
    names(transition_potential) <- colnames(z)
  }
  transition_potential <- as.numeric(transition_potential[colnames(z)])
  if (anyNA(transition_potential)) {
    stop("`transition_potential` must be named for every state in `z`.", call. = FALSE)
  }
  cell_transition <- as.vector(z %*% transition_potential)
  fpi <- entropy + lambda * cell_transition
  data.frame(
    cell_id = rownames(z),
    entropy = entropy,
    transition_potential = cell_transition,
    fpi = fpi,
    stringsAsFactors = FALSE
  )
}

normalize_composition <- function(x, name) {
  x_names <- names(x)
  x <- as.numeric(x)
  names(x) <- x_names
  if (is.null(names(x))) {
    stop(sprintf("`%s` must be a named composition vector.", name), call. = FALSE)
  }
  if (anyNA(x) || any(x < 0) || sum(x) <= 0) {
    stop(sprintf("`%s` must contain non-negative values with positive sum.", name), call. = FALSE)
  }
  x / sum(x)
}
