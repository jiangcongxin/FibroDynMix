#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "scripts/run_transition_flow.R"
ROOT <- normalizePath(file.path(dirname(script_path), ".."))
OUT <- file.path(ROOT, "analysis", "transition_flow")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "R", "simulate_fibrodynmix.R"))
source(file.path(ROOT, "R", "transition_flow.R"))

sim <- simulate_fibrodynmix(
  n_studies = 2,
  donors_per_study = 3,
  cells_per_donor = 30,
  n_genes = 300,
  marker_genes_per_state = 12,
  scenario = "rare_transition",
  seed = 610
)

normal <- colMeans(sim$z[sim$cell_metadata$disease == "normal", , drop = FALSE])
disease <- colMeans(sim$z[sim$cell_metadata$disease == "disease", , drop = FALSE])
cost <- compute_state_cost(sim$parameters$beta_kg)
flow <- estimate_transition_flow(normal, disease, cost, lambda = 0.5)
fpi <- compute_fpi(sim$z, flow = flow$flow)
fpi <- cbind(sim$cell_metadata[, c("cell_id", "disease", "study_id", "donor_id", "is_transition")], fpi[, -1, drop = FALSE])

write.table(as.data.frame(cost), file.path(OUT, "state_transition_cost.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(flow$flow), file.path(OUT, "state_transition_flow.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(
  data.frame(state = names(normal), normal = normal, disease = disease),
  file.path(OUT, "condition_state_composition.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(fpi, file.path(OUT, "cell_fpi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(
    expected_cost = flow$expected_cost,
    entropy = flow$entropy,
    converged = flow$converged,
    iterations = flow$iterations
  ),
  file.path(OUT, "transition_flow_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Transition flow written to: ", OUT)
