#' Fit a topic/NMF baseline from raw counts
#'
#' Fits a lightweight KL-divergence non-negative matrix factorization baseline
#' to a gene-by-cell raw count matrix. Topics are aligned to fibroblast states
#' after fitting by marker enrichment, so marker genes orient the output labels
#' but do not define the cell-topic factors.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param marker_index Named list of marker genes for each state. Entries can be
#'   gene names or integer row indices.
#' @param n_topics Number of NMF topics. Defaults to the number of states.
#' @param n_iter Number of multiplicative-update iterations.
#' @param seed Optional random seed for initialization.
#' @param eps Small positive value for numerical stability.
#' @param backend Factorization backend. `"auto"` uses the `NMF` package when
#'   installed and otherwise falls back to the internal KL multiplicative
#'   updates.
#' @param nmf_method Method passed to `NMF::nmf()` when `backend = "nmf"`.
#' @param nrun Number of NMF runs passed to `NMF::nmf()`.
#'
#' @return A list with topic gene loadings `w`, cell topic loadings `h`,
#'   state-aligned cell weights `z_pred`, topic-state scores, topic assignment,
#'   and KL objective trace.
#' @export
fit_topic_nmf_baseline <- function(counts,
                                   marker_index,
                                   n_topics = length(marker_index),
                                   n_iter = 100,
                                   seed = NULL,
                                   eps = 1e-8,
                                   backend = c("auto", "nmf", "multiplicative_update"),
                                   nmf_method = "brunet",
                                   nrun = 1) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || length(dim(counts)) != 2L) {
    stop("`counts` must be a numeric matrix.", call. = FALSE)
  }
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("`counts` must contain finite non-negative values.", call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- sprintf("gene_%d", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("cell_%d", seq_len(ncol(counts)))
  }
  if (!is.list(marker_index) || is.null(names(marker_index))) {
    stop("`marker_index` must be a named list.", call. = FALSE)
  }
  assert_positive_integer(n_topics, "n_topics")
  assert_positive_integer(n_iter, "n_iter")
  if (length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("`eps` must be a positive numeric scalar.", call. = FALSE)
  }
  backend <- match.arg(backend)
  if (backend == "auto") {
    backend <- if (requireNamespace("NMF", quietly = TRUE)) "nmf" else "multiplicative_update"
  }
  assert_positive_integer(nrun, "nrun")

  if (backend == "nmf" && !requireNamespace("NMF", quietly = TRUE)) {
    stop("`backend = \"nmf\"` requires the NMF package.", call. = FALSE)
  }

  if (!is.null(seed) && backend == "multiplicative_update") {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  v <- counts + eps
  n_genes <- nrow(v)
  n_cells <- ncol(v)
  if (backend == "nmf") {
    nmf_seed <- if (is.null(seed)) "random" else seed
    nmf_fit <- suppressPackageStartupMessages(suppressMessages(
      NMF::nmf(v, rank = n_topics, method = nmf_method, nrun = nrun, seed = nmf_seed)
    ))
    w <- as.matrix(NMF::basis(nmf_fit))
    h <- as.matrix(NMF::coef(nmf_fit))
    rownames(w) <- rownames(counts)
    colnames(w) <- sprintf("topic_%d", seq_len(ncol(w)))
    rownames(h) <- colnames(w)
    colnames(h) <- colnames(counts)
    objective_trace <- kl_nmf_objective(v, w %*% h + eps)
    executed_iterations <- suppressWarnings(as.integer(NMF::niter(nmf_fit)))
  } else {
    w <- matrix(stats::rgamma(n_genes * n_topics, shape = 1, rate = 1), nrow = n_genes, ncol = n_topics)
    h <- matrix(stats::rgamma(n_topics * n_cells, shape = 1, rate = 1), nrow = n_topics, ncol = n_cells)
    rownames(w) <- rownames(counts)
    colnames(w) <- sprintf("topic_%d", seq_len(n_topics))
    rownames(h) <- colnames(w)
    colnames(h) <- colnames(counts)

    objective_trace <- numeric(n_iter)
    for (iter in seq_len(n_iter)) {
      wh <- w %*% h + eps
      h <- h * ((t(w) %*% (v / wh)) / pmax(matrix(colSums(w), nrow = n_topics, ncol = n_cells), eps))
      h <- pmax(h, eps)

      wh <- w %*% h + eps
      w <- w * (((v / wh) %*% t(h)) / pmax(matrix(rowSums(h), nrow = n_genes, ncol = n_topics, byrow = TRUE), eps))
      w <- pmax(w, eps)

      scale <- pmax(colSums(w), eps)
      w <- sweep(w, 2, scale, "/")
      h <- sweep(h, 1, scale, "*")
      objective_trace[iter] <- kl_nmf_objective(v, w %*% h + eps)
    }
    executed_iterations <- n_iter
  }

  alignment <- align_nmf_topics_to_states(w, marker_index)
  state_names <- names(marker_index)
  z_pred <- matrix(eps, nrow = n_cells, ncol = length(state_names), dimnames = list(colnames(counts), state_names))
  for (state in state_names) {
    topic <- alignment$state_to_topic[[state]]
    if (!is.na(topic)) {
      z_pred[, state] <- h[topic, ]
    }
  }
  z_pred <- sweep(z_pred, 1, pmax(rowSums(z_pred), eps), "/")

  list(
    w = w,
    h = h,
    z_pred = z_pred,
    topic_state_scores = alignment$topic_state_scores,
    topic_to_state = alignment$topic_to_state,
    state_to_topic = alignment$state_to_topic,
    objective_trace = objective_trace,
    final_objective = utils::tail(objective_trace, 1),
    backend = backend,
    nmf_method = if (backend == "nmf") nmf_method else NA_character_,
    n_iter = executed_iterations
  )
}

align_nmf_topics_to_states <- function(w, marker_index) {
  state_names <- names(marker_index)
  topic_names <- colnames(w)
  scores <- matrix(0, nrow = ncol(w), ncol = length(marker_index), dimnames = list(topic_names, state_names))
  for (state in state_names) {
    marker_rows <- resolve_marker_rows(marker_index[[state]], rownames(w))
    if (length(marker_rows) > 0L) {
      scores[, state] <- colMeans(w[marker_rows, , drop = FALSE])
    }
  }

  assignments <- data.frame(
    topic = character(),
    state = character(),
    score = numeric(),
    stringsAsFactors = FALSE
  )
  available_topics <- topic_names
  available_states <- state_names
  while (length(available_topics) > 0L && length(available_states) > 0L) {
    sub_scores <- scores[available_topics, available_states, drop = FALSE]
    best <- which(sub_scores == max(sub_scores), arr.ind = TRUE)[1, , drop = FALSE]
    topic <- rownames(sub_scores)[best[1, "row"]]
    state <- colnames(sub_scores)[best[1, "col"]]
    assignments <- rbind(
      assignments,
      data.frame(topic = topic, state = state, score = scores[topic, state], stringsAsFactors = FALSE)
    )
    available_topics <- setdiff(available_topics, topic)
    available_states <- setdiff(available_states, state)
  }

  topic_to_state <- stats::setNames(rep(NA_character_, length(topic_names)), topic_names)
  state_to_topic <- stats::setNames(rep(NA_character_, length(state_names)), state_names)
  for (i in seq_len(nrow(assignments))) {
    topic_to_state[[assignments$topic[i]]] <- assignments$state[i]
    state_to_topic[[assignments$state[i]]] <- assignments$topic[i]
  }

  list(
    topic_state_scores = scores,
    topic_to_state = topic_to_state,
    state_to_topic = state_to_topic
  )
}

kl_nmf_objective <- function(v, wh) {
  v_safe <- pmax(v, .Machine$double.eps)
  wh_safe <- pmax(wh, .Machine$double.eps)
  sum(v_safe * log(v_safe / wh_safe) - v_safe + wh_safe) / length(v_safe)
}
