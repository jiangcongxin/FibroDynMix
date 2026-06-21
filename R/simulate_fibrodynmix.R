#' Simulate FibroDynMix count data
#'
#' Generates raw UMI-like counts from the FibroDynMix negative-binomial
#' hierarchical mixture model. The simulator returns observed counts together
#' with the latent state simplex, donor/study metadata, and true generative
#' parameters used for benchmarking.
#'
#' @param n_studies Number of studies or cohorts.
#' @param donors_per_study Number of donors per study. Either a scalar or an
#'   integer vector of length `n_studies`.
#' @param cells_per_donor Number of cells per donor. Either a scalar or an
#'   integer vector with one value per donor.
#' @param n_genes Number of genes.
#' @param state_names Character vector of fibroblast state names.
#' @param scenario Simulation scenario: `continuous`, `discrete`,
#'   `batch_confounding`, or `rare_transition`.
#' @param marker_genes_per_state Number of prior marker genes assigned to each
#'   state in the data-generating beta matrix.
#' @param tau_high Standard deviation for true marker state effects.
#' @param tau_low Standard deviation for true non-marker state effects.
#' @param study_effect_sd Standard deviation of study-by-gene effects.
#' @param donor_effect_sd Standard deviation of donor-by-gene effects.
#' @param state_sd Cell-level logistic-normal standard deviation.
#' @param library_meanlog Meanlog for log-normal library sizes.
#' @param library_sdlog Sdlog for log-normal library sizes.
#' @param seed Optional random seed.
#'
#' @return A list with `counts`, `cell_metadata`, `gene_metadata`, `z`, and
#'   `parameters`.
#' @export
simulate_fibrodynmix <- function(n_studies = 2,
                                 donors_per_study = 4,
                                 cells_per_donor = 100,
                                 n_genes = 1000,
                                 state_names = c(
                                   "resident",
                                   "inflammatory",
                                   "myofibroblast",
                                   "ECM-remodeling",
                                   "antigen-presenting",
                                   "IFN-stress"
                                 ),
                                 scenario = c(
                                   "continuous",
                                   "discrete",
                                   "batch_confounding",
                                   "rare_transition"
                                 ),
                                 marker_genes_per_state = 30,
                                 tau_high = 1.0,
                                 tau_low = 0.08,
                                 study_effect_sd = 0.15,
                                 donor_effect_sd = 0.08,
                                 state_sd = 0.35,
                                 library_meanlog = log(5000),
                                 library_sdlog = 0.35,
                                 seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  assert_positive_integer(n_studies, "n_studies")
  assert_positive_integer(n_genes, "n_genes")
  assert_positive_integer(marker_genes_per_state, "marker_genes_per_state")

  state_names <- as.character(state_names)
  n_states <- length(state_names)
  if (n_states < 2) {
    stop("`state_names` must contain at least two states.", call. = FALSE)
  }
  if (marker_genes_per_state * n_states > n_genes) {
    stop(
      "`marker_genes_per_state * length(state_names)` must be <= `n_genes`.",
      call. = FALSE
    )
  }

  donors_per_study <- recycle_integer(donors_per_study, n_studies, "donors_per_study")
  n_donors <- sum(donors_per_study)
  cells_per_donor <- recycle_integer(cells_per_donor, n_donors, "cells_per_donor")
  n_cells <- sum(cells_per_donor)

  gene_ids <- sprintf("Gene%04d", seq_len(n_genes))
  study_ids <- sprintf("Study%02d", seq_len(n_studies))
  donor_ids <- sprintf("Donor%03d", seq_len(n_donors))

  donor_metadata <- build_donor_metadata(study_ids, donors_per_study, scenario)
  rownames(donor_metadata) <- donor_ids

  marker_state <- rep(NA_character_, n_genes)
  marker_index <- vector("list", n_states)
  names(marker_index) <- state_names
  cursor <- 1L
  for (k in seq_len(n_states)) {
    idx <- cursor:(cursor + marker_genes_per_state - 1L)
    marker_index[[k]] <- idx
    marker_state[idx] <- state_names[k]
    cursor <- cursor + marker_genes_per_state
  }

  alpha_g <- rnorm(n_genes, mean = -8.6, sd = 0.7)
  beta_kg <- matrix(
    rnorm(n_states * n_genes, mean = 0, sd = tau_low),
    nrow = n_states,
    dimnames = list(state_names, gene_ids)
  )
  for (k in seq_len(n_states)) {
    beta_kg[k, marker_index[[k]]] <- rnorm(marker_genes_per_state, 1.2, tau_high / 4)
  }

  phi_g <- rgamma(n_genes, shape = 10, rate = 2)
  study_effect <- matrix(
    rnorm(n_studies * n_genes, 0, study_effect_sd),
    nrow = n_studies,
    dimnames = list(study_ids, gene_ids)
  )
  donor_effect <- matrix(
    rnorm(n_donors * n_genes, 0, donor_effect_sd),
    nrow = n_donors,
    dimnames = list(donor_ids, gene_ids)
  )

  if (scenario == "batch_confounding") {
    disease_genes <- unlist(marker_index[c("inflammatory", "myofibroblast")], use.names = FALSE)
    disease_genes <- disease_genes[!is.na(disease_genes)]
    study_effect[n_studies, disease_genes] <- study_effect[n_studies, disease_genes] + 0.45
  }

  counts <- matrix(0L, nrow = n_genes, ncol = n_cells, dimnames = list(gene_ids, NULL))
  z <- matrix(NA_real_, nrow = n_cells, ncol = n_states, dimnames = list(NULL, state_names))
  cell_metadata <- data.frame(
    cell_id = sprintf("Cell%06d", seq_len(n_cells)),
    study_id = character(n_cells),
    donor_id = character(n_cells),
    disease = character(n_cells),
    tissue = character(n_cells),
    library_size = numeric(n_cells),
    dominant_state = character(n_cells),
    is_transition = logical(n_cells),
    stringsAsFactors = FALSE
  )

  cell_cursor <- 1L
  for (d in seq_len(n_donors)) {
    donor_id <- donor_ids[d]
    study_id <- donor_metadata$study_id[d]
    study_idx <- match(study_id, study_ids)
    donor_cells <- cells_per_donor[d]
    cell_idx <- cell_cursor:(cell_cursor + donor_cells - 1L)

    donor_mean <- donor_state_mean(donor_metadata[d, ], state_names, scenario)
    eta <- matrix(
      rnorm(donor_cells * n_states, 0, state_sd),
      nrow = donor_cells,
      ncol = n_states
    )
    eta <- sweep(eta, 2, donor_mean, "+")
    z_i <- softmax_rows(eta)

    if (scenario == "discrete") {
      dominant <- max.col(z_i, ties.method = "first")
      z_i[,] <- 0
      z_i[cbind(seq_len(donor_cells), dominant)] <- 1
    }

    transition_cells <- rep(FALSE, donor_cells)
    if (scenario == "rare_transition") {
      rare_n <- max(1L, round(0.04 * donor_cells))
      rare_idx <- sample.int(donor_cells, rare_n)
      transition_cells[rare_idx] <- TRUE
      transition_profile <- rep(0.02, n_states)
      names(transition_profile) <- state_names
      transition_profile[intersect(c("resident", "inflammatory", "myofibroblast"), state_names)] <-
        c(0.35, 0.35, 0.22)[seq_along(intersect(c("resident", "inflammatory", "myofibroblast"), state_names))]
      transition_profile <- transition_profile / sum(transition_profile)
      z_i[rare_idx, ] <- matrix(
        transition_profile,
        nrow = rare_n,
        ncol = n_states,
        byrow = TRUE
      )
    }

    library_size <- round(rlnorm(donor_cells, library_meanlog, library_sdlog))
    linear_predictor <- z_i %*% beta_kg
    linear_predictor <- sweep(linear_predictor, 2, alpha_g, "+")
    linear_predictor <- sweep(linear_predictor, 2, study_effect[study_idx, ], "+")
    linear_predictor <- sweep(linear_predictor, 2, donor_effect[d, ], "+")
    mu <- exp(linear_predictor)
    mu <- sweep(mu, 1, library_size, "*")

    counts[, cell_idx] <- vapply(
      seq_len(donor_cells),
      function(i) stats::rnbinom(n_genes, size = phi_g, mu = mu[i, ]),
      numeric(n_genes)
    )

    z[cell_idx, ] <- z_i
    cell_metadata[cell_idx, "study_id"] <- study_id
    cell_metadata[cell_idx, "donor_id"] <- donor_id
    cell_metadata[cell_idx, "disease"] <- donor_metadata$disease[d]
    cell_metadata[cell_idx, "tissue"] <- donor_metadata$tissue[d]
    cell_metadata[cell_idx, "library_size"] <- library_size
    cell_metadata[cell_idx, "dominant_state"] <- state_names[max.col(z_i, ties.method = "first")]
    cell_metadata[cell_idx, "is_transition"] <- transition_cells

    cell_cursor <- cell_cursor + donor_cells
  }

  colnames(counts) <- cell_metadata$cell_id
  rownames(z) <- cell_metadata$cell_id

  gene_metadata <- data.frame(
    gene_id = gene_ids,
    marker_state = marker_state,
    alpha = alpha_g,
    phi = phi_g,
    stringsAsFactors = FALSE
  )

  list(
    counts = counts,
    cell_metadata = cell_metadata,
    donor_metadata = donor_metadata,
    gene_metadata = gene_metadata,
    z = z,
    parameters = list(
      state_names = state_names,
      alpha_g = alpha_g,
      beta_kg = beta_kg,
      phi_g = phi_g,
      study_effect = study_effect,
      donor_effect = donor_effect,
      marker_index = marker_index,
      scenario = scenario
    )
  )
}

assert_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive integer scalar.", name), call. = FALSE)
  }
}

recycle_integer <- function(x, n, name) {
  if (length(x) == 1L) {
    x <- rep(x, n)
  }
  if (length(x) != n || any(is.na(x)) || any(x < 1) || any(x != as.integer(x))) {
    stop(sprintf("`%s` must be a positive integer scalar or length-%d vector.", name, n), call. = FALSE)
  }
  as.integer(x)
}

softmax_rows <- function(eta) {
  eta <- eta - apply(eta, 1, max)
  exp_eta <- exp(eta)
  exp_eta / rowSums(exp_eta)
}

build_donor_metadata <- function(study_ids, donors_per_study, scenario) {
  study_id <- rep(study_ids, donors_per_study)
  n_donors <- length(study_id)
  disease <- rep(c("normal", "disease"), length.out = n_donors)
  if (scenario == "batch_confounding" && length(study_ids) > 1L) {
    disease <- ifelse(study_id == tail(study_ids, 1), "disease", "normal")
  }
  tissue <- rep(c("skin", "lung", "colon"), length.out = n_donors)
  data.frame(
    study_id = study_id,
    disease = disease,
    tissue = tissue,
    stringsAsFactors = FALSE
  )
}

donor_state_mean <- function(donor_row, state_names, scenario) {
  state_mean <- rep(-1.5, length(state_names))
  names(state_mean) <- state_names
  state_mean["resident"] <- 1.2

  if (donor_row$disease == "disease") {
    state_mean[intersect(c("inflammatory", "myofibroblast", "ECM-remodeling"), state_names)] <-
      c(0.7, 0.9, 0.6)[seq_along(intersect(c("inflammatory", "myofibroblast", "ECM-remodeling"), state_names))]
    state_mean["resident"] <- 0.1
  }

  if (scenario == "continuous") {
    state_mean <- state_mean * 0.75
  }
  if (scenario == "rare_transition") {
    state_mean[intersect(c("inflammatory", "myofibroblast"), state_names)] <-
      state_mean[intersect(c("inflammatory", "myofibroblast"), state_names)] + 0.25
  }

  state_mean[is.na(state_mean)] <- -1.5
  state_mean
}
